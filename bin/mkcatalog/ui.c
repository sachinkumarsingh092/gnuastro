/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"
#include "mkcatalog.h"

#include "ui.h"
#include "columns.h"
#include "authors-cite.h"





/**************************************************************/
/*********      Argp necessary global entities     ************/
/**************************************************************/
/* Definition parameters for the Argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "ASTRdata";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will create a catalog from "
  "an input, labeled, and noise images.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;

/* Option groups particular to this program. */
enum program_args_groups
{
  ARGS_GROUP_UPPERLIMIT = GAL_OPTIONS_GROUP_AFTER_COMMON,
  ARGS_GROUP_COLUMNS,
};

















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct mkcatalogparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->program_struct     = p;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->poptions           = program_options;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;

  /* Specific to this program. */
  p->nsigmag             = NAN;
  p->upnsigma            = NAN;
  p->zeropoint           = NAN;
  p->threshold           = NAN;
  p->upsigmaclip[0]      = NAN;
  p->upsigmaclip[1]      = NAN;


  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_IGNORECASE:
          cp->coptions[i].flags=OPTION_HIDDEN;
          cp->coptions[i].mandatory=GAL_OPTIONS_NOT_MANDATORY;
          break;

        case GAL_OPTIONS_KEY_TABLEFORMAT:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;
        }

      /* Select by group. */
      switch(cp->coptions[i].group)
        {
        case GAL_OPTIONS_GROUP_TESSELLATION:
          cp->coptions[i].doc=NULL; /* Necessary to remove title. */
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }
    }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct mkcatalogparams *p = state->input;

  /* Pass `gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then `arg`
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign (`=`). For short "
               "options, `=` should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {

    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      if(p->inputname)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->inputname=arg;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}





/* Read the user's desired columns. Because the types of these options are
   `GAL_TYPE_INVALID', this function will not be called when printing the
   full list of parameters and their values. */
void *
ui_column_codes_ll(struct argp_option *option, char *arg,
                   char *filename, size_t lineno, void *params)
{
  struct mkcatalogparams *p=(struct mkcatalogparams *)params;

  /* These options don't take arguments on the command-line but in the
     configuration files they can take values of 0 or 1. In the latter
     case, the column shouldn't be added if the value is 0. */
  if(arg)
    {
      if( arg[0]=='0' && arg[1]=='\0' ) return NULL;
      else if ( !(arg[0]=='1' && arg[1]=='\0')  )
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "`%s' is not a "
                      "valid value to the `%s' option: (\"%s\").\n\n`%s' is "
                      "an on/off option specifying if you want this column "
                      "in the output catalog, it doesn't take any "
                      "arguments. In a configuration file, it can only take "
                      "a value of `0' (to be ignored) or `1'", arg,
                      option->name, option->doc, option->name);
    }


  /* The user wants this column, so add it to the list. Note that the `ids'
     column means three columns. */
  if(option->key==UI_KEY_IDS)
    {
      gal_list_i32_add(&p->columnids, UI_KEY_OBJID);
      gal_list_i32_add(&p->columnids, UI_KEY_HOSTOBJID);
      gal_list_i32_add(&p->columnids, UI_KEY_IDINHOSTOBJ);
    }
  else
    gal_list_i32_add(&p->columnids, option->key);

  /* Return NULL */
  return NULL;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in `ui_check_options_and_arguments'.
static void
ui_read_check_only_options(struct mkcatalogparams *p)
{
}
*/




static void
ui_check_options_and_arguments(struct mkcatalogparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->inputname)
    {
      if( gal_fits_name_is_fits(p->inputname) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified. When the input is a FITS "
              "file, a HDU must also be specified, you can use the `--hdu' "
              "(`-h') option and give it the HDU number (starting from "
              "zero), extension name, or anything acceptable by CFITSIO");

    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_preparations_read_inputs(struct mkcatalogparams *p)
{
  size_t one;
  gal_list_i32_t *tmp;
  int readclumps=0, readwcs=0;
  char *namestypes, **strarr=NULL;
  gal_data_t *zero, *key=gal_data_array_calloc(1);
  char *skyfile=p->skyfile ? p->skyfile : p->inputname;
  char *stdfile=p->stdfile ? p->stdfile : p->inputname;
  char *clumpsfile=p->clumpsfile ? p->clumpsfile : p->inputname;
  char *objectsfile=p->objectsfile ? p->objectsfile : p->inputname;


  /* Read the input image. */
  p->input=gal_fits_img_read_to_type(p->inputname, p->cp.hdu,
                                     GAL_TYPE_FLOAT32, p->cp.minmapsize);


  /* Currently MakeCatalog is only implemented for 2D images. */
  if(p->input->ndim!=2)
    error(EXIT_FAILURE, 0, "%s (%s) has %zu dimensions, MakeCatalog "
          "currently only supports 2D inputs", p->inputname, p->cp.hdu,
          p->input->ndim);


  /* See if the input has blank pixels and set the flags appropriately. */
  p->hasblank = gal_blank_present(p->input, 1);


  /* Read the object label image and check its size. */
  p->objects = gal_fits_img_read(objectsfile, p->objectshdu,
                                 p->cp.minmapsize);
  if( gal_data_dsize_is_different(p->input, p->objects) )
    error(EXIT_FAILURE, 0, "`%s' (hdu: %s) and `%s' (hdu: %s) have a"
          "different dimension/size", objectsfile, p->objectshdu,
          p->inputname, p->cp.hdu);


  /* Read the Sky image and check its size. */
  p->sky=gal_fits_img_read_to_type(skyfile, p->skyhdu, GAL_TYPE_FLOAT32,
                                   p->cp.minmapsize);
  if( gal_data_dsize_is_different(p->input, p->sky) )
    error(EXIT_FAILURE, 0, "`%s' (hdu: %s) and `%s' (hdu: %s) have a"
          "different dimension/size", skyfile, p->skyhdu, p->inputname,
          p->cp.hdu);


  /* Read the Sky standard deviation image and check its size. */
  p->std=gal_fits_img_read_to_type(stdfile, p->stdhdu, GAL_TYPE_FLOAT32,
                                   p->cp.minmapsize);
  if( gal_data_dsize_is_different(p->input, p->std) )
    error(EXIT_FAILURE, 0, "`%s' (hdu: %s) and `%s' (hdu: %s) have a"
          "different dimension/size", stdfile, p->stdhdu, p->inputname,
          p->cp.hdu);


  /* If an upper-limit mask image was given, read it. */
  if(p->upmaskfile)
    {
      /* Read the mask image. */
      p->upmask = gal_fits_img_read(p->upmaskfile, p->upmaskhdu,
                                    p->cp.minmapsize);
      if( gal_data_dsize_is_different(p->input, p->upmask) )
        error(EXIT_FAILURE, 0, "`%s' (hdu: %s) and `%s' (hdu: %s) have a"
              "different dimension/size", p->upmaskfile, p->upmaskhdu,
          p->inputname, p->cp.hdu);

      /* If it isn't an integer type, report an error, otherwise, convert
         it to a uint8_t: with a 1 for all non-zero pixels and 0 for zero
         pixels. */
      zero=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &one, NULL, 1, -1,
                          NULL, NULL, NULL);
      p->upmask=gal_arithmetic(GAL_ARITHMETIC_OP_NE,
                               ( GAL_ARITHMETIC_INPLACE | GAL_ARITHMETIC_FREE
                                 | GAL_ARITHMETIC_NUMOK ), p->upmask, zero);
    }


  /* Check to see if the objects extension has a `WCLUMPS' keyword and that
     its value is 'yes', `y', or `1'. */
  key->name="WCLUMPS";
  key->type=GAL_TYPE_STRING;
  gal_fits_key_read(objectsfile, p->objectshdu, key, 0, 0);
  if(key->status==KEY_NO_EXIST) readclumps=0;
  else
    {
      if(key->status)
        gal_fits_io_error(key->status, "CFITSIO error while reading "
                          "WCLUMPS keyword");
      else
        {
          strarr=key->array;
          if( !strcasecmp(strarr[0], "yes") || !strcasecmp(strarr[0], "y")
              || !strcmp(strarr[0], "1") ) readclumps=1;
        }
    }


  /* Read the clumps array if necessary. */
  if(readclumps)
    {
      /* Make sure the user did indeed give the clumps HDU. */
      if(p->clumpshdu==NULL)
        error(EXIT_FAILURE, 0, "no `--clumpshdu' given! The WCLUMPS keyword "
              "in %s (hdu: %s) has a value of `%s', so MakeCatalog expects "
              "a clump image.\n\nYou can use the optional `--clumpsfile' "
              "option to give the filename and the mandatory `--clumpshdu' "
              "to specify the extension. If `--clumpsfile' is not given, "
              "MakeCatalog will look into the input file for the given "
              "extension. Alternatively, you can modify/remove this "
              "keyword using Gnuastro's Fits program, please run `$ info "
              "astfits' for more information (press `SPACE' to go down and "
              "`q' to return to the command-line).", objectsfile,
              p->objectshdu, strarr[0]);

      /* Read the clumps image and check its size. */
      p->clumps = gal_fits_img_read(clumpsfile, p->clumpshdu,
                                    p->cp.minmapsize);
      if( gal_data_dsize_is_different(p->input, p->std) )
        error(EXIT_FAILURE, 0, "`%s' (hdu: %s) and `%s' (hdu: %s) have a"
              "different dimension/size", clumpsfile, p->clumpshdu,
              p->inputname, p->cp.hdu);
    }


  /* Make sure the object and clump label images have an integer type, then
     to be safe that they have the correct integer type, run
     `gal_data_copy_to_new_type' on them. */
  if( p->objects->type==GAL_TYPE_FLOAT32
      || p->objects->type==GAL_TYPE_FLOAT64
      || (p->clumps && ( p->clumps->type==GAL_TYPE_FLOAT32
                         || p->clumps->type==GAL_TYPE_FLOAT64 ) ) )
    {
      if(p->clumps)
        asprintf(&namestypes, "However, `%s' (hdu: %s) and `%s' (hdu: %s) "
                 "have types of `%s' and `%s' respectively", objectsfile,
                 p->objectshdu, clumpsfile, p->clumpshdu,
                 gal_type_name(p->objects->type, 1),
                 gal_type_name(p->clumps->type, 1) );
      else
        asprintf(&namestypes, "However, %s (hdu: %s) has a type of %s",
                 objectsfile, p->objectshdu,
                 gal_type_name(p->objects->type, 1));
      error(EXIT_FAILURE, 0, "labeled images (for objects or clumps) must "
            "have an integer datatype. %s.\n\n"
            "If you are sure the images contain only integer values but "
            "are just stored in a floating point container, you can "
            "put them in an integer container with Gnuastro's Arithmetic "
            "program using a command like below:\n\n"
            "    $ astarithmetic img.fits int32", namestypes);
    }
  p->objects=gal_data_copy_to_new_type_free(p->objects, GAL_TYPE_INT32);
  if(p->clumps)
    p->clumps=gal_data_copy_to_new_type_free(p->clumps, GAL_TYPE_INT32);


  /* If any WCS related parameter is requested then read the input's WCS
     structure. */
  for(tmp=p->columnids; tmp!=NULL; tmp=tmp->next)
    if(readwcs==0)
      switch(tmp->v)
        {
        case UI_KEY_RA:
        case UI_KEY_DEC:
        case UI_KEY_GEORA:
        case UI_KEY_GEODEC:
        case UI_KEY_CLUMPSRA:
        case UI_KEY_CLUMPSDEC:
        case UI_KEY_CLUMPSGEORA:
        case UI_KEY_CLUMPSGEODEC:
          readwcs=1;
          break;
        }
  if(readwcs)
    p->input->wcs=gal_wcs_read(p->inputname, p->cp.hdu, 0, 0, &p->input->nwcs);


  /* Clean up. */
  key->name=NULL;
  gal_data_array_free(key, 1, 1);
}





/* The input images can have extensions to speed up the processing. */
static void
ui_preparations_read_keywords(struct mkcatalogparams *p)
{
  gal_data_t *tmp;
  gal_data_t *keys=gal_data_array_calloc(2);
  char *stdfile=p->stdfile ? p->stdfile : p->inputname;
  char *clumpsfile=p->clumpsfile ? p->clumpsfile : p->inputname;
  char *objectsfile=p->objectsfile ? p->objectsfile : p->inputname;

  /* Read the keywords from the standard deviation image. */
  keys[0].next=&keys[1];
  keys[0].name="MINSTD";                    keys[1].name="MEDSTD";
  keys[0].type=GAL_TYPE_FLOAT32;            keys[1].type=GAL_TYPE_FLOAT32;
  keys[0].array=&p->minstd;                 keys[1].array=&p->medstd;
  gal_fits_key_read(stdfile, p->stdhdu, keys, 0, 0);

  /* If the two keywords couldn't be read, calculate them. */
  if(keys[0].status)
    {
      /* Calculate the minimum STD. */
      tmp=gal_statistics_minimum(p->std);
      p->minstd=*((float *)(tmp->array));
      gal_data_free(tmp);
    }
  if(keys[1].status)
    {
      /* Calculate the median STD. */
      tmp=gal_statistics_median(p->std, 0);
      p->medstd=*((float *)(tmp->array));
      gal_data_free(tmp);

      /* Alert the user if it wasn't calculated from a header keyword. */
      fprintf(stderr, "---------------\n"
              "Warning: Could not find the `MEDSTD' keyword in `%s' "
              "(hdu: %s). The median standard deviation is thus found on "
              "the (interpolated) standard deviation image. NoiseChisel "
              "finds the median before interpolation which is more "
              "accurate. Ho the reported value in the final catalog will "
              "not be accurate: it will depend on how many tiles were "
              "blank and their spatial position relative to the non-blank "
              "ones.\n"
              "---------------\n", stdfile, p->stdhdu);
    }
  p->cpscorr = p->minstd>1 ? 1.0f : p->minstd;


  /* Read the keywords from the objects image. */
  keys[0].name="DETSN";                     keys[1].name="NUMLABS";
  keys[0].type=GAL_TYPE_FLOAT32;            keys[1].type=GAL_TYPE_SIZE_T;
  keys[0].array=&p->detsn;                  keys[1].array=&p->numobjects;
  gal_fits_key_read(objectsfile, p->objectshdu, keys, 0, 0);
  if(keys[0].status) p->detsn=NAN;
  if(keys[1].status)
    {
      tmp=gal_statistics_maximum(p->objects);
      p->numobjects=*((int32_t *)(tmp->array)); /* objects is in int32_t. */
      gal_data_free(tmp);
    }


  /* If there were no objects in the input, then inform the user with an
     error (no catalog was built). */
  if(p->numobjects==0)
    error(EXIT_FAILURE, 0, "no object labels (non-zero pixels) in "
          "%s (hdu %s). To make a catalog, labeled regions must be defined",
          objectsfile, p->objectshdu);


  /* Read the keywords from the clumps image if necessary. */
  if(p->clumps)
    {
      keys[0].name="CLUMPSN";               keys[1].name="NUMLABS";
      keys[0].type=GAL_TYPE_FLOAT32;        keys[1].type=GAL_TYPE_SIZE_T;
      keys[0].array=&p->clumpsn;            keys[1].array=&p->numclumps;
      gal_fits_key_read(clumpsfile, p->clumpshdu, keys, 0, 0);
      if(keys[0].status) p->clumpsn=NAN;
      if(keys[1].status)
        error(EXIT_FAILURE, 0, "couldn't find NCLUMPS in the header of "
              "%s (hdu: %s).", p->clumpsfile, p->clumpshdu);
    }


  /* If there were no clumps, then free the clumps array and set it to
     NULL, so for the rest of the processing, MakeCatalog things that no
     clumps image was given. */
  if(p->numclumps==0)
    {
      gal_data_free(p->clumps);
      p->clumps=NULL;
    }


  /* Clean up. */
  keys[0].name=keys[1].name=NULL;
  keys[0].array=keys[1].array=NULL;
  gal_data_array_free(keys, 2, 1);
}





/* To make the catalog processing more scalable (and later allow for
   over-lappping regions), we will define a tile for each object. */
void
ui_one_tile_per_object(struct mkcatalogparams *p)
{
  size_t ndim=p->input->ndim;

  int32_t *l, *lf, *start;
  size_t i, d, *min, *max, width=2*ndim;
  size_t *minmax=gal_data_malloc_array(GAL_TYPE_SIZE_T,
                                       width*p->numobjects, __func__,
                                       "minmax");
  size_t *coord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                      "coord");


  /* Initialize the minimum and maximum position for each tile/object. So,
     we'll initialize the minimum coordinates to the maximum possible
     `size_t' value (in `GAL_BLANK_SIZE_T') and the maximums to zero. */
  for(i=0;i<p->numobjects;++i)
    for(d=0;d<ndim;++d)
      {
        minmax[ i * width +        d ] = GAL_BLANK_SIZE_T; /* Minimum. */
        minmax[ i * width + ndim + d ] = 0;                /* Maximum. */
      }

  /* Go over the objects label image and correct the minimum and maximum
     coordinates. */
  start=p->objects->array;
  lf=(l=p->objects->array)+p->objects->size;
  do
    if(*l>0)
      {
        /* Get the coordinates of this pixel. */
        gal_dimension_index_to_coord(l-start, ndim, p->objects->dsize, coord);

        /* Check to see this coordinate is the smallest/largest found so
           far for this label. Note that labels start from 1, while indexs
           here start from zero. */
        min = &minmax[ (*l-1) * width        ];
        max = &minmax[ (*l-1) * width + ndim ];
        for(d=0;d<ndim;++d)
          {
            if( coord[d] < min[d] ) min[d] = coord[d];
            if( coord[d] > max[d] ) max[d] = coord[d];
          }
      }
  while(++l<lf);

  /* For a check.
  for(i=0;i<p->numobjects;++i)
    printf("%zu: (%zu, %zu) --> (%zu, %zu)\n", i+1, minmax[i*width],
           minmax[i*width+1], minmax[i*width+2], minmax[i*width+3]);
  */

  /* Make the tiles. */
  p->tiles=gal_tile_series_from_minmax(p->input, minmax, p->numobjects);

  /* Clean up. */
  free(coord);
  free(minmax);
}





/* When both catalogs need to be made, we need a separator, the output
   names will either be built based on the input name or output name (if
   given). In both cases, the operations are the same, just the base name
   differs. So to keep things clean, we have defined this function. */
static void
ui_preparations_both_names(struct mkcatalogparams *p, char *basename)
{
  char *end, suffix[50];

  /* Set the objects name */
  end="_o";
  sprintf(suffix, "%s.%s", end,
          p->cp.tableformat==GAL_TABLE_FORMAT_TXT ? "txt" : "fits");
  p->objectsout=gal_checkset_automatic_output(&p->cp, p->inputname,
                                              suffix);

  /* Set the clumps name */
  end="_c";
  sprintf(suffix, "%s.%s", end,
          p->cp.tableformat==GAL_TABLE_FORMAT_TXT ? "txt" : "fits");
  p->clumpsout=gal_checkset_automatic_output(&p->cp, p->inputname,
                                             suffix);
}





/* Set the output name. */
static void
ui_preparations_outnames(struct mkcatalogparams *p)
{
  char *suffix;

  /* Set the output filename */
  if(p->cp.output)
    {
      /* Make sure the given output filename corresponds to the
         `tableformat' option. */
      gal_tableintern_check_fits_format(p->cp.output, p->cp.tableformat);

      /* If a clumps image has been read (and there are any clumps in the
         image), then we have two outputs. */
      if(p->clumps)
        ui_preparations_both_names(p, p->cp.output);
      else
        p->objectsout=p->cp.output;
    }
  else
    {
      /* Both clumps and object catalogs are necessary. */
      if(p->clumps)
        ui_preparations_both_names(p, p->inputname);

      /* We only need one objects catalog. */
      else
        {
          suffix = ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                     ? ".txt" : ".fits" );
          p->objectsout=gal_checkset_automatic_output(&p->cp, p->inputname,
                                                      suffix);
        }
    }
}





/* Sanity checks and preparations for the upper-limit magnitude. */
static void
ui_preparations_upperlimit(struct mkcatalogparams *p)
{
  /* Check the number of random samples. */
  if( p->upnum < MKCATALOG_UPPERLIMIT_MINIMUM_NUM )
    error(EXIT_FAILURE, 0, "%zu not acceptable as `--upnum'. The minimum "
          "acceptable number of random samples for the upper limit "
          "magnitude is %d", p->upnum, MKCATALOG_UPPERLIMIT_MINIMUM_NUM);

  /* Check if sigma-clipping parameters have been given. */
  if( isnan(p->upsigmaclip[0]) )
    error(EXIT_FAILURE, 0, "`--upsigmaclip' is mandatory for measuring "
          "the upper-limit magnitude. It takes two numbers separated by "
          "a comma. The first is the multiple of sigma and the second is "
          "the aborting criteria: <1: tolerance level, >1: number of "
          "clips");

  /* Check if the sigma multiple is given. */
  if( isnan(p->upnsigma) )
    error(EXIT_FAILURE, 0, "`--upnsigma' is mandatory for measuring the "
          "upperlimit magnitude. Its value is the multiple of final sigma "
          "that is reported as the upper-limit");

  /* Set the random number generator. */
  gsl_rng_env_setup();
  p->rng=gsl_rng_alloc(gsl_rng_default);
  p->seed = ( p->envseed
              ? gsl_rng_default_seed
              : gal_timing_time_based_rng_seed() );
  if(p->envseed) gsl_rng_set(p->rng, p->seed);
  p->rngname=gsl_rng_name(p->rng);
}









void
ui_preparations(struct mkcatalogparams *p)
{
  /* If no columns are requested, then inform the user. */
  if(p->columnids==NULL)
    error(EXIT_FAILURE, 0, "no columns requested, please run again with "
          "`--help' for the full list of columns you can ask for");


  /* Read the inputs. */
  ui_preparations_read_inputs(p);


  /* Read the helper keywords from the inputs and if they aren't present
     then calculate the necessary parameters. */
  ui_preparations_read_keywords(p);


  /* Set the output filename(s). */
  ui_preparations_outnames(p);


  /* Allocate the output columns to fill up with the program. */
  columns_define_alloc(p);


  /* Make the tiles that cover each object. */
  ui_one_tile_per_object(p);


  /* Allocate the reference random number generator and seed values. It
     will be cloned once for every thread. If the user hasn't called
     `envseed', then we want it to be different for every run, so we need
     to re-set the seed. */
  if(p->upperlimit) ui_preparations_upperlimit(p);

  if( p->hasmag && isnan(p->zeropoint) )
    error(EXIT_FAILURE, 0, "no zeropoint specified");
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkcatalogparams *p)
{
  struct gal_options_common_params *cp=&p->cp;
  char *skyfile, *stdfile, *clumpsfile, *objectsfile;


  /* Include the parameters necessary for argp from this program (`args.h')
     and for the common options to all Gnuastro (`commonopts.h'). We want
     to directly put the pointers to the fields in `p' and `cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <gnuastro-internal/commonopts.h>
#include "args.h"


  /* Initialize the options and necessary information.  */
  ui_initialize_options(p, program_options, gal_commonopts_options);


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Read the configuration files and set the common values. */
  gal_options_read_config_set(&p->cp);


  /* Read the options into the program's structure, and check them and
     their relations prior to printing.
  ui_read_check_only_options(p);
  */

  /* Print the option values if asked. Note that this needs to be done
     after the option checks so un-sane values are not printed in the
     output state. */
  gal_options_print_state(&p->cp);


  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);


  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);


  /* Inform the user. */
  if(!p->cp.quiet)
    {
      /* Set the names for easy reading. */
      skyfile     = p->skyfile     ? p->skyfile     : p->inputname;
      stdfile     = p->stdfile     ? p->stdfile     : p->inputname;
      clumpsfile  = p->clumpsfile  ? p->clumpsfile  : p->inputname;
      objectsfile = p->objectsfile ? p->objectsfile : p->inputname;

      /* Write the information. */
      printf(PROGRAM_NAME" started on %s", ctime(&p->rawtime));
      printf("  - Using %zu CPU thread%s\n", p->cp.numthreads,
             p->cp.numthreads==1 ? "." : "s.");
      printf("  - Input:   %s (hdu: %s)\n", p->inputname, p->cp.hdu);
      printf("  - Objects: %s (hdu: %s)\n", objectsfile,  p->objectshdu);
      if(p->clumps)
        printf("  - Clumps:  %s (hdu: %s)\n", clumpsfile, p->clumpshdu);
      printf("  - Sky:     %s (hdu: %s)\n", skyfile, p->skyhdu);
      printf("  - Sky STD: %s (hdu: %s)\n", stdfile, p->stdhdu);
      if(p->upmaskfile)
        printf("  - Upper limit magnitude mask: %s (hdu: %s)\n",
               p->upmaskfile, p->cp.hdu);
      if(p->upperlimit)
        {
          printf("  - Random number generator name: %s\n", p->rngname);
          printf("  - Random number generator seed: %"PRIu64"\n", p->seed);
        }
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct mkcatalogparams *p, struct timeval *t1)
{
  size_t d;

  /* The temporary arrays for RA and Dec. */
  if(p->rd_vo || p->rd_vc || p->rd_go || p->rd_gc || p->rd_vcc || p->rd_gcc)
    {
      /* First free the separate dimensions in each value. */
      for(d=0;d<p->input->ndim;++d)
        {
          if(p->rd_vo)  free(p->rd_vo[d]);
          if(p->rd_vc)  free(p->rd_vc[d]);
          if(p->rd_go)  free(p->rd_go[d]);
          if(p->rd_gc)  free(p->rd_gc[d]);
          if(p->rd_vcc) free(p->rd_vcc[d]);
          if(p->rd_gcc) free(p->rd_gcc[d]);
        }

      /* Then free the container arrays. */
      if(p->rd_vo)  free(p->rd_vo);
      if(p->rd_vc)  free(p->rd_vc);
      if(p->rd_go)  free(p->rd_go);
      if(p->rd_gc)  free(p->rd_gc);
      if(p->rd_vcc) free(p->rd_vcc);
      if(p->rd_gcc) free(p->rd_gcc);
    }

  /* If a random number generator was allocated, free it. */
  if(p->rng) gsl_rng_free(p->rng);

  /* Free the allocated arrays: */
  free(p->skyhdu);
  free(p->stdhdu);
  free(p->cp.hdu);
  free(p->oiflag);
  free(p->ciflag);
  free(p->skyfile);
  free(p->stdfile);
  free(p->cp.output);
  free(p->clumpshdu);
  free(p->objectshdu);
  free(p->clumpsfile);
  free(p->objectsfile);
  gal_data_free(p->sky);
  gal_data_free(p->std);
  gal_data_free(p->input);
  gal_data_free(p->upmask);
  gal_data_free(p->clumps);
  gal_data_free(p->objects);
  gal_data_array_free(p->tiles, p->numobjects, 0);

  /* Print the final message. */
  if(!p->cp.quiet)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
