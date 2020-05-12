/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <gnuastro/array.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
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
  p->medstd         = NAN;
  p->sfmagnsigma    = NAN;
  p->sfmagarea      = NAN;
  p->upnsigma       = NAN;
  p->zeropoint      = NAN;
  p->sigmaclip[0]   = NAN;
  p->sigmaclip[1]   = NAN;
  p->upsigmaclip[0] = NAN;
  p->upsigmaclip[1] = NAN;
  p->checkuplim[0]  = GAL_BLANK_INT32;
  p->checkuplim[1]  = GAL_BLANK_INT32;

  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_LOG:
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_IGNORECASE:
        case GAL_OPTIONS_KEY_WORKOVERCH:
        case GAL_OPTIONS_KEY_STDINTIMEOUT:
        case GAL_OPTIONS_KEY_INTERPNUMNGB:
        case GAL_OPTIONS_KEY_INTERPONLYBLANK:
          cp->coptions[i].flags=OPTION_HIDDEN;
          cp->coptions[i].mandatory=GAL_OPTIONS_NOT_MANDATORY;
          break;

        case GAL_OPTIONS_KEY_TABLEFORMAT:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;
        }
    }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct mkcatalogparams *p = state->input;

  /* Pass 'gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then 'arg'
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign ('='). For short "
               "options, '=' should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {

    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      if(p->objectsfile)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->objectsfile=arg;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}





/* Read the user's desired columns. Because the types of these options are
   'GAL_TYPE_INVALID', this function will not be called when printing the
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
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "'%s' is not a "
                      "valid value to the '%s' option: ('%s').\n\n'%s' is "
                      "an on/off option specifying if you want this column "
                      "in the output catalog, it doesn't take any "
                      "arguments. In a configuration file, it can only take "
                      "a value of '0' (to be ignored) or '1'", arg,
                      option->name, option->doc, option->name);
    }


  /* The user wants this column, so add it to the list. Note that the 'ids'
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





/* Prepare the upper-limit distribution parameters. */
void *
ui_check_upperlimit(struct argp_option *option, char *arg,
                    char *filename, size_t lineno, void *params)
{
  size_t i;
  char *str;
  double *d;
  gal_data_t *raw;
  struct mkcatalogparams *p=(struct mkcatalogparams *)params;

  /* Write. */
  if(lineno==-1)
    {
      if(p->checkuplim[1]==GAL_BLANK_INT32)
        {
          if( asprintf(&str, "%d", p->checkuplim[0])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        if( asprintf(&str, "%d,%d", p->checkuplim[0], p->checkuplim[1])<0 )
          error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      return str;
    }

  /* Read */
  else
    {
      /* If the option is already set, just return. */
      if(option->set) return NULL;

      /* Read the list of numbers as an array. */
      raw=gal_options_parse_list_of_numbers(arg, filename, lineno);

      /* Make sure there is at most only two numbers given. */
      if(raw->size>2)
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "'%s' (value to "
                      "'--%s') contains %zu numbers, but only one or two "
                      "are acceptable.\n\n"
                      "With this option MakeCatalog will write all the "
                      "positions and values of the random distribution for "
                      "one particular labeled region into a table. The "
                      "given value(s) is(are) the label identifier.\n\n"
                      "With one value the distribution for an object will "
                      "be printed: the givne number will be interpretted as "
                      "the requested object's label. With two values, the "
                      "distribution for a specific clump will be written. "
                      "The first will be interpretted as the clump's host "
                      "object label and the second as the clump's label "
                      "within the object", arg, option->name, raw->size);

      /* Make sure the given values are integers and that they are larger
         than zero. */
      d=raw->array;
      for(i=0;i<raw->size;++i)
        {
          if( ceil(d[i]) != d[i])
            error_at_line(EXIT_FAILURE, 0, filename, lineno, "%g (value "
                          "number %zu given to '--%s') is not an "
                          "integer. The value(s) to this option are "
                          "object/clump labels/identifiers, so they must be "
                          "integers", d[i], i+1, option->name);
          if( d[i]<=0 )
            error_at_line(EXIT_FAILURE, 0, filename, lineno, "%g (value "
                          "number %zu given to '--%s') is not positive. "
                          "The value(s) to this option are object/clump "
                          "labels/identifiers, so they must be positive "
                          "integers", d[i], i+1, option->name);
        }

      /* Write the values in. */
      p->checkuplim[0] = d[0];
      p->checkuplim[1] = raw->size==2 ? d[1] : GAL_BLANK_INT32;

      /* For no un-used variable warning. This function doesn't need the
         pointer.*/
      return NULL;
    }
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct mkcatalogparams *p)
{
  float tmp;
  size_t one=1;
  char *tailptr;

  /* If an upper-limit check table is requested with a specific clump, but
     no clump catalog has been requested, then abort and inform the
     user. */
  if( p->checkuplim[1]!=GAL_BLANK_INT32 && p->clumpscat==0 )
    error(EXIT_FAILURE, 0, "no clumps catalog is requested, hence "
          "'--checkuplim' is only available for objects (one value "
          "must be given to it).\n\n"
          "To ask for a clumps catalog, please append '--clumpscat' to the "
          "command calling MakeCatalog.\n\n"
          "If you want the upperlimit check table for an object, only give "
          "one value (the object's label) to '--checkuplim'.");

  /* See if '--skyin' is a filename or a value. When the string is ONLY a
     number (and nothing else), 'tailptr' will point to the end of the
     string ('\0'). */
  if(p->skyfile)
    {
      tmp=strtod(p->skyfile, &tailptr);
      if(*tailptr=='\0')
        {
          /* Allocate the data structure. */
          p->sky=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0, -1,
                                1, NULL, NULL, NULL);

          /* Write the value inside it. */
          *((float *)(p->sky->array))=tmp;
        }
    }

  /* Similar to the case for Sky above. */
  if(p->stdfile)
    {
      tmp=strtod(p->stdfile, &tailptr);
      if(*tailptr=='\0')
        {
          /* Allocate the data structure. */
          p->std=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0, -1,
                                1, NULL, NULL, NULL);

          /* Write the value inside it. */
          *((float *)(p->std->array))=tmp;
        }
    }
}




static void
ui_check_options_and_arguments(struct mkcatalogparams *p)
{
  /* Make sure the main input file name (for the object labels) was given
     and if it was a FITS file, that a HDU is also given. */
  if(p->objectsfile)
    {
      if( gal_fits_name_is_fits(p->objectsfile) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified. When the input is a FITS "
              "file, a HDU must also be specified, you can use the '--hdu' "
              "('-h') option and give it the HDU number (starting from "
              "zero), extension name, or anything acceptable by CFITSIO");

    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* If the user hasn't explicitly specified a filename for input,
   MakeCatalog will use other given file names. */
static void
ui_set_filenames(struct mkcatalogparams *p)
{
  p->usedclumpsfile = p->clumpsfile ? p->clumpsfile : p->objectsfile;

  p->usedvaluesfile = p->valuesfile ? p->valuesfile : p->objectsfile;

  p->usedskyfile    = ( p->skyfile
                       ? p->skyfile
                       : ( p->valuesfile ? p->valuesfile : p->objectsfile ) );

  p->usedstdfile    = ( p->stdfile
                       ? p->stdfile
                       : ( p->valuesfile ? p->valuesfile : p->objectsfile ) );
}





/* The clumps and objects images must be integer type, so we'll use this
   function to avoid having to write the error message two times. */
static void
ui_check_type_int(char *filename, char *hdu, uint8_t type)
{
  if( type==GAL_TYPE_FLOAT32 || type==GAL_TYPE_FLOAT64 )
    error(EXIT_FAILURE, 0, "%s (hdu: %s): type %s not acceptable as "
          "labels input. labeled images must have an integer datatype.\n\n"
          "If you are sure the extension contains only integer values but "
          "it is just stored in a floating point container, you can "
          "put it in an integer container with Gnuastro's Arithmetic "
          "program using this command:\n\n"
          "    $ astarithmetic %s int32 -h%s", filename, hdu,
          gal_type_name(type, 1), filename, hdu);
}





/* If a WCS structure is present, then read its basic information to use in
   the table meta-data. */
static void
ui_wcs_info(struct mkcatalogparams *p)
{
  char *c;
  size_t i;

  /* Read the WCS meta-data. */
  p->objects->wcs=gal_wcs_read(p->objectsfile, p->cp.hdu, 0, 0,
                               &p->objects->nwcs);

  /* Read the basic WCS information. */
  if(p->objects->wcs)
    {
      /* Allocate space for the array of strings. */
      errno=0;
      p->ctype=malloc(p->objects->ndim * sizeof *p->ctype);
      if(p->ctype==NULL)
        error(EXIT_FAILURE, 0, "%s: %zu bytes for 'p->ctype'", __func__,
              p->objects->ndim * sizeof *p->ctype);

      /* Fill in the values. */
      for(i=0;i<p->objects->ndim;++i)
        {
          /* CTYPE might contain '-' characters, we just want the first
             non-dash characters. The loop will either stop either at the end
             or where there is a dash. So we can just replace it with an
             end-of-string character. */
          gal_checkset_allocate_copy(p->objects->wcs->ctype[i], &p->ctype[i]);
          c=p->ctype[i]; while(*c!='\0' && *c!='-') ++c;
          *c='\0';
        }
    }
}





static size_t
ui_num_clumps(struct mkcatalogparams *p)
{
  char *basename;
  int keepinputdir;
  size_t i, counter, numclumps=0;
  gal_list_i32_t *tmp, **labsinobj;
  int32_t *o=p->objects->array, *of=o+p->objects->size, *c=p->clumps->array;

  /* Allocate array of lists to keep the unique labels within each object. */
  errno=0;
  labsinobj=calloc(p->numobjects+1, sizeof *labsinobj);
  if(labsinobj==NULL)
    error(EXIT_FAILURE, 0, "%s: couldn't allocate %zu bytes for "
          "'labsinobj'", __func__, p->numobjects * sizeof *labsinobj);

  /* Go over each pixel and find the unique clump labels within each
     object. */
  do
    {
      /* Do the steps if we are on a clump. */
      if(*o>0 && *c>0)
        {
          /* See if the label has already been found. */
          for(tmp=labsinobj[*o];tmp!=NULL;tmp=tmp->next) if(tmp->v==*c) break;

          /* When it wasn't found, 'tmp==NULL'. */
          if(tmp==NULL)
            {
              ++numclumps;
              gal_list_i32_add(&labsinobj[*o], *c);
            }
        }

      /* Increment the clumps pointer.*/
      ++c;
    }
  while(++o<of);

  /* Re-write the clump values so their numbering is contiguous, since this
     is assumed during the later steps. */
  o=p->objects->array;
  c=p->clumps->array;
  do
    {
      /* Do the steps if we are on a clump. */
      if(*o>0 && *c>0)
        {
          counter=0;
          for(tmp=labsinobj[*o];tmp!=NULL;tmp=tmp->next)
            { counter++; if(tmp->v==*c) {*c=counter; break;} }
        }

      /* Increment the clumps pointer.*/
      ++c;
    }
  while(++o<of);

  /* Write the created file into a file for the user to inspect. */
  keepinputdir=p->cp.keepinputdir;
  p->cp.keepinputdir = p->cp.output ? 1 : 0;
  basename = p->cp.output ? p->cp.output : p->objectsfile;
  p->relabclumps=gal_checkset_automatic_output(&p->cp, basename,
                                               "-clumps-relab.fits");
  gal_fits_img_write(p->clumps, p->relabclumps, NULL, PROGRAM_STRING);
  p->cp.keepinputdir=keepinputdir;

  /* Clean up. */
  for(i=0;i<p->numobjects;++i)
    if(labsinobj[i]) gal_list_i32_free(labsinobj[i]);
  free(labsinobj);

  /* Return the number of clumps. */
  return numclumps;
}





/* To make the catalog processing more scalable (and later allow for
   over-lappping regions), we will define a tile for each object. */
static void
ui_one_tile_per_object_correct_numobjects(struct mkcatalogparams *p)
{
  size_t ndim=p->objects->ndim;

  uint8_t *rarray=NULL;
  int32_t *l, *lf, *start;
  gal_data_t *rowsremove=NULL;
  size_t i, j, d, no, *min, *max, exists, width=2*ndim;
  size_t *minmax=gal_pointer_allocate(GAL_TYPE_SIZE_T,
                                      width*p->numobjects, 0, __func__,
                                      "minmax");
  size_t *coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                      "coord");

  /* Initialize the minimum and maximum position for each tile/object. So,
     we'll initialize the minimum coordinates to the maximum possible
     'size_t' value (in 'GAL_BLANK_SIZE_T') and the maximums to zero. */
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
        gal_dimension_index_to_coord(l-start, ndim, p->objects->dsize,
                                     coord);

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

  /* If a label doesn't exist in the image, then write over it and define
     the unique labels to use for the next steps. To over-write, we have
     two counters: 'i' (for the position in the input array) and 'no' (or
     'num-objects' for the counter in the output table). In the end, 'no'
     counts the total number of unique labels in the input. */
  no=0;
  for(i=0;i<p->numobjects;++i)
    {
      /* Make sure a pixel with this label exists in all dimensions. */
      exists=0;
      for(d=0;d<ndim;++d)
        if ( minmax[ i * width + d ] == GAL_BLANK_SIZE_T
             && minmax[ i * width + ndim + d ] == 0 )
          {
            /* When the object doesn't exist, but the user wants a row
               anyway, make all the minimums and maximums of all
               coordinates 0, note that the maximum is already zero. */
            if(p->inbetweenints)
              minmax[ i * width + d ] = 0;
          }
        else
          {
            /* Write over the blank elements when necessary
               (i!=j). When i==j, then these statements are
               redundant. */
            minmax[no*width+d]=minmax[i*width+d];
            minmax[no*width+ndim+d]=minmax[i*width+ndim+d];

            /* Set the checked flag. */
            exists=1;
          }

      /* If it does (or if the user wants to keep all integers), then
         increment the output counter.*/
      if(p->inbetweenints || exists) ++no;
      else
        {
          /* If 'rarray' isn't defined yet, define it. */
          if(rarray==NULL)
            {
              /* Note that by initializing with zeros, all (the possibly
                 existing) previous rows that shouldn't be removed are
                 flagged as zero in this array. */
              rowsremove=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1,
                                        &p->numobjects, NULL, 1,
                                        p->cp.minmapsize, p->cp.quietmmap,
                                        NULL, NULL, NULL);
              rarray=rowsremove->array;
            }
          rarray[i]=1;
        }
    }

  /* If 'rarray!=NULL', then there are elements to remove and we need to
     make some modifications/corrections. */
  if(rarray)
    {
      /* Build an array to keep the real ID of each tile. */
      j=0;
      p->outlabs=gal_pointer_allocate(GAL_TYPE_INT32, no, 0, __func__,
                                      "p->outlabs");
      for(i=0;i<p->numobjects;++i) if(rarray[i]==0) p->outlabs[j++]=i+1;

      /* Correct numobjects and clean up. */
      p->numobjects=no;
      gal_data_free(rowsremove);
    }

  /* For a check.
  for(i=0;i<p->numobjects;++i)
    printf("%zu: (%zu, %zu) --> (%zu, %zu)\n",
           p->outlabs ? p->outlabs[i] : i+1, minmax[i*width],
           minmax[i*width+1], minmax[i*width+2], minmax[i*width+3]);
  */

  /* Make the tiles. */
  p->tiles=gal_tile_series_from_minmax(p->objects, minmax, p->numobjects);

  /* Clean up. */
  free(coord);
  free(minmax);
}





/* The only mandatory input is the objects image, so first read that and
   make sure its type is correct. */
static void
ui_read_labels(struct mkcatalogparams *p)
{
  gal_data_t *tmp, *keys=gal_data_array_calloc(2);

  /* Read it into memory. */
  p->objects = gal_array_read_one_ch(p->objectsfile, p->cp.hdu, NULL,
                                     p->cp.minmapsize, p->cp.quietmmap);
  p->objects->ndim=gal_dimension_remove_extra(p->objects->ndim,
                                              p->objects->dsize, NULL);


  /* Make sure it has an integer type. */
  ui_check_type_int(p->objectsfile, p->cp.hdu, p->objects->type);


  /* Convert it to 'int32' type (if it already isn't). */
  p->objects=gal_data_copy_to_new_type_free(p->objects, GAL_TYPE_INT32);


  /* Currently MakeCatalog is only implemented for 2D images or 3D cubes. */
  if(p->objects->ndim!=2 && p->objects->ndim!=3)
    error(EXIT_FAILURE, 0, "%s (hdu %s) has %zu dimensions, MakeCatalog "
          "currently only supports 2D or 3D datasets", p->objectsfile,
          p->cp.hdu, p->objects->ndim);

  /* Make sure the '--spectrum' option is not given on a 2D image.  */
  if(p->spectrum && p->objects->ndim!=3)
    error(EXIT_FAILURE, 0, "%s (hdu %s) has %zu dimensions, but '--spectrum' "
          "is currently only defined on 3D datasets", p->objectsfile,
          p->cp.hdu, p->objects->ndim);

  /* See if the total number of objects is given in the header keywords. */
  keys[0].name="NUMLABS";
  keys[0].type=GAL_TYPE_SIZE_T;
  keys[0].array=&p->numobjects;
  gal_fits_key_read(p->objectsfile, p->cp.hdu, keys, 0, 0);
  if(keys[0].status) /* status!=0: the key couldn't be read by CFITSIO. */
    {
      tmp=gal_statistics_maximum(p->objects);
      p->numobjects=*((int32_t *)(tmp->array)); /*numobjects is in int32_t.*/
      gal_data_free(tmp);
    }

  /* If there were no objects in the input, then inform the user with an
     error (it is pointless to build a catalog). */
  if(p->numobjects==0)
    error(EXIT_FAILURE, 0, "no object labels (non-zero pixels) in "
          "%s (hdu %s). To make a catalog, labeled regions must be defined",
          p->objectsfile, p->cp.hdu);


  /* See if the labels image has blank pixels and set the flags
     appropriately. */
  p->hasblank = gal_blank_present(p->objects, 1);


  /* Prepare WCS information for final table meta-data. */
  ui_wcs_info(p);


  /* Make the tiles that cover each object and also correct the total
     number of objects based on the parsing of the image. */
  ui_one_tile_per_object_correct_numobjects(p);


  /* Read the clumps array if necessary. */
  if(p->clumpscat)
    {
      /* Make sure the HDU is also given. */
      if(p->clumpshdu==NULL)
        error(EXIT_FAILURE, 0, "%s: no HDU/extension provided for the "
              "CLUMPS dataset. Please use the '--clumpshdu' option to "
              "give a specific HDU using its number (counting from zero) "
              "or name. If the dataset is in another file, please use "
              "'--clumpsfile' to give the filename. If you don't want any "
              "clumps catalog output, remove the '--clumpscat' option from "
              "the command-line or give it a value of zero in a "
              "configuration file", p->usedclumpsfile);

      /* Read the clumps image. */
      p->clumps = gal_array_read_one_ch(p->usedclumpsfile, p->clumpshdu,
                                        NULL, p->cp.minmapsize,
                                        p->cp.quietmmap);
      p->clumps->ndim=gal_dimension_remove_extra(p->clumps->ndim,
                                                 p->clumps->dsize, NULL);

      /* Check its size. */
      if( gal_dimension_is_different(p->objects, p->clumps) )
        error(EXIT_FAILURE, 0, "'%s' (hdu: %s) and '%s' (hdu: %s) have a"
              "different dimension/size", p->usedclumpsfile, p->clumpshdu,
              p->objectsfile, p->cp.hdu);

      /* Check its type. */
      ui_check_type_int(p->usedclumpsfile, p->clumpshdu, p->clumps->type);
      p->clumps=gal_data_copy_to_new_type_free(p->clumps, GAL_TYPE_INT32);

      /* See if there are keywords to help in finding the number. */
      keys[0].next=&keys[1];
      keys[0].status=keys[1].status=0;
      keys[0].name="CLUMPSN";               keys[1].name="NUMLABS";
      keys[0].type=GAL_TYPE_FLOAT32;        keys[1].type=GAL_TYPE_SIZE_T;
      keys[0].array=&p->clumpsn;            keys[1].array=&p->numclumps;
      gal_fits_key_read(p->usedclumpsfile, p->clumpshdu, keys, 0, 0);
      if(keys[0].status) p->clumpsn=NAN;
      if(keys[1].status) p->numclumps=ui_num_clumps(p);

      /* If there were no clumps, then free the clumps array and set it to
         NULL, so for the rest of the processing, MakeCatalog things that
         no clumps image was given. */
      if(p->numclumps==0)
        {
          /* Just as a sanity check, see if there are any clumps (positive
             valued pixels) in the array. If there are, then 'NUMCLUMPS'
             wasn't set properly and we should abort with an error. */
          tmp=gal_statistics_maximum(p->clumps);
          if( *((int32_t *)(p->clumps->array))>0 )
            error(EXIT_FAILURE, 0, "%s (hdu: %s): the 'NUMCLUMPS' header "
                  "keyword has a value of zero, but there are positive "
                  "pixels in the array, showing that there are clumps in "
                  "image. This is a wrong usage of the 'NUMCLUMPS' keyword."
                  "It must contain the total number of clumps (irrespective "
                  "of how many objects there are). Please correct this issue "
                  "and run MakeCatalog again", p->usedclumpsfile,
                  p->clumpshdu);

          /* Since there are no clumps, we won't bother creating a clumps
             catalog and from this step onward, we'll act as if no clumps
             catalog was requested. In order to not confuse the user in the
             end, we'll print a warning first. */
          fprintf(stderr, "WARNING: %s (hdu %s): there are no clumps "
                  "in the image, therefore no clumps catalog will be "
                  "created.\n", p->usedclumpsfile, p->clumpshdu);
          gal_data_free(p->clumps);
          p->clumps=NULL;
        }
    }

  /* Clean up. */
  keys[0].name=keys[1].name=NULL;
  keys[0].array=keys[1].array=NULL;
  gal_data_array_free(keys, 2, 1);
}





/* See which inputs are necessary. Ultimate, there are only three extra
   inputs: a values image, a sky image and a sky standard deviation
   image. However, there are many raw column measurements. So to keep
   things clean, we'll just put a value of '1' in the three 'values', 'sky'
   and 'std' pointers everytime a necessary input is found. */
static void
ui_necessary_inputs(struct mkcatalogparams *p, int *values, int *sky,
                    int *std)
{
  size_t i;

  /* Set necessary inputs based on options. */
  if(p->forcereadstd) *std=1;
  if(p->upperlimit) *values=1;
  if(p->spectrum) *values=*std=1;

  /* Go over all the object columns. Note that the objects and clumps (if
     the '--clumpcat' option is given) inputs are mandatory and it is not
     necessary to specify it here. */
  for(i=0; i<OCOL_NUMCOLS; ++i)
    if(p->oiflag[i])
      switch(i)
        {
        case OCOL_NUMALL:             /* Only object labels. */    break;
        case OCOL_NUMALLXY:           /* Only object labels. */    break;
        case OCOL_NUM:                *values        = 1;          break;
        case OCOL_NUMXY:              *values        = 1;          break;
        case OCOL_SUM:                *values        = 1;          break;
        case OCOL_SUM_VAR:            *values = *std = 1;          break;
        case OCOL_MEDIAN:             *values        = 1;          break;
        case OCOL_SIGCLIPNUM:         *values        = 1;          break;
        case OCOL_SIGCLIPMEDIAN:      *values        = 1;          break;
        case OCOL_SIGCLIPMEAN:        *values        = 1;          break;
        case OCOL_SIGCLIPSTD:         *values        = 1;          break;
        case OCOL_VX:                 *values        = 1;          break;
        case OCOL_VY:                 *values        = 1;          break;
        case OCOL_VZ:                 *values        = 1;          break;
        case OCOL_VXX:                *values        = 1;          break;
        case OCOL_VYY:                *values        = 1;          break;
        case OCOL_VXY:                *values        = 1;          break;
        case OCOL_SUMSKY:             *sky           = 1;          break;
        case OCOL_SUMVAR:             *std           = 1;          break;
        case OCOL_SUMWHT:             *values        = 1;          break;
        case OCOL_NUMWHT:             *values        = 1;          break;
        case OCOL_GX:                 /* Only object labels. */    break;
        case OCOL_GY:                 /* Only object labels. */    break;
        case OCOL_GZ:                 /* Only object labels. */    break;
        case OCOL_GXX:                /* Only object labels. */    break;
        case OCOL_GYY:                /* Only object labels. */    break;
        case OCOL_GXY:                /* Only object labels. */    break;
        case OCOL_UPPERLIMIT_B:       *values        = 1;          break;
        case OCOL_UPPERLIMIT_S:       *values        = 1;          break;
        case OCOL_UPPERLIMIT_Q:       *values        = 1;          break;
        case OCOL_UPPERLIMIT_SKEW:    *values        = 1;          break;
        case OCOL_C_NUMALL:           /* Only clump labels.  */    break;
        case OCOL_C_NUM:              *values        = 1;          break;
        case OCOL_C_SUM:              *values        = 1;          break;
        case OCOL_C_VX:               *values        = 1;          break;
        case OCOL_C_VY:               *values        = 1;          break;
        case OCOL_C_VZ:               *values        = 1;          break;
        case OCOL_C_GX:               /* Only clump labels. */     break;
        case OCOL_C_GY:               /* Only clump labels. */     break;
        case OCOL_C_GZ:               /* Only clump labels. */     break;
        case OCOL_C_SUMWHT:           *values        = 1;          break;
        case OCOL_C_NUMWHT:           *values        = 1;          break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                "fix the problem. The code %zu is not a recognized "
                "intermediate OBJECT columns", __func__, PACKAGE_BUGREPORT,
                i);
        }

  /* Check the clump elements also. */
  if(p->clumps)
    for(i=0; i<CCOL_NUMCOLS; ++i)
      if(p->ciflag[i])
        switch(i)
          {
          case CCOL_NUMALL:           /* Only clump labels. */     break;
          case CCOL_NUMALLXY:         /* Only clump labels. */     break;
          case CCOL_NUM:              *values        = 1;          break;
          case CCOL_NUMXY:            *values        = 1;          break;
          case CCOL_SUM:              *values        = 1;          break;
          case CCOL_SUM_VAR:          *values = *std = 1;          break;
          case CCOL_MEDIAN:           *values        = 1;          break;
          case CCOL_SIGCLIPNUM:    *values        = 1;          break;
          case CCOL_SIGCLIPMEDIAN:    *values        = 1;          break;
          case CCOL_SIGCLIPMEAN:      *values        = 1;          break;
          case CCOL_SIGCLIPSTD:       *values        = 1;          break;
          case CCOL_RIV_NUM:          /* Only clump labels. */     break;
          case CCOL_RIV_SUM:          *values        = 1;          break;
          case CCOL_RIV_SUM_VAR:      *values = *std = 1;          break;
          case CCOL_VX:               *values        = 1;          break;
          case CCOL_VY:               *values        = 1;          break;
          case CCOL_VZ:               *values        = 1;          break;
          case CCOL_VXX:              *values        = 1;          break;
          case CCOL_VYY:              *values        = 1;          break;
          case CCOL_VXY:              *values        = 1;          break;
          case CCOL_SUMSKY:           *sky           = 1;          break;
          case CCOL_SUMVAR:           *std           = 1;          break;
          case CCOL_SUMWHT:           *values        = 1;          break;
          case CCOL_NUMWHT:           *values        = 1;          break;
          case CCOL_GX:               /* Only clump labels. */     break;
          case CCOL_GY:               /* Only clump labels. */     break;
          case CCOL_GZ:               /* Only clump labels. */     break;
          case CCOL_GXX:              /* Only clump labels. */     break;
          case CCOL_GYY:              /* Only clump labels. */     break;
          case CCOL_GXY:              /* Only clump labels. */     break;
          case CCOL_MINX:             /* Only clump labels. */     break;
          case CCOL_MAXX:             /* Only clump labels. */     break;
          case CCOL_MINY:             /* Only clump labels. */     break;
          case CCOL_MAXY:             /* Only clump labels. */     break;
          case CCOL_MINZ:             /* Only clump labels. */     break;
          case CCOL_MAXZ:             /* Only clump labels. */     break;
          case CCOL_UPPERLIMIT_B:     *values        = 1;          break;
          case CCOL_UPPERLIMIT_S:     *values        = 1;          break;
          case CCOL_UPPERLIMIT_Q:     *values        = 1;          break;
          case CCOL_UPPERLIMIT_SKEW:  *values        = 1;          break;
          default:
            error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                  "fix the problem. The code %zu is not a recognized "
                  "intermediate CLUMP column", __func__, PACKAGE_BUGREPORT,
                  i);
          }
}





/* When the Sky and its standard deviation are given as tiles, we need to
   define a tile structure. */
static void
ui_preparation_check_size_read_tiles(struct mkcatalogparams *p,
                                     gal_data_t *in, char *filename,
                                     char *hdu)
{
  struct gal_tile_two_layer_params *tl=&p->cp.tl;

  /* See if we should treat this dataset as tile values or not. */
  if( gal_dimension_is_different(p->objects, in) )
    {
      /* The 'tl' structure is initialized here. But this function may be
         called multiple times. So, first check if the 'tl' structure has
         already been initialized and if so, don't repeat it. */
      if(tl->ndim==0)
        {
          gal_tile_full_sanity_check(p->objectsfile, p->cp.hdu, p->objects,
                                     tl);
          gal_tile_full_two_layers(p->objects, tl);
          gal_tile_full_permutation(tl);
        }

      /* See if the size of the 'in' dataset corresponds to the
         tessellation. */
      if(in->size!=tl->tottiles)
        error(EXIT_FAILURE, 0, "%s (hdu: %s): doesn't have the right "
              "size (%zu elements or pixels).\n\n"
              "It must either be the same size as '%s' (hdu: '%s'), or "
              "it must have the same number of elements as the total "
              "number of tiles in the tessellation (%zu). In the latter "
              "case, each pixel is assumed to be a fixed value for a "
              "complete tile.\n\n"
              "Run with '-P' to see the (tessellation) options/settings "
              "and their values). For more information on tessellation in "
              "Gnuastro, please run the following command (use the arrow "
              "keys for up and down and press 'q' to return to the "
              "command-line):\n\n"
              "    $ info gnuastro tessellation",
              filename, hdu, in->size, p->objectsfile, p->cp.hdu,
              tl->tottiles);
    }
}





/* Subtract 'sky' from the input dataset depending on its size (it may be
   the whole array or a tile-values array).. */
static void
ui_subtract_sky(struct mkcatalogparams *p)
{
  size_t tid;
  gal_data_t *tile;
  float *s, *f, *ff, *skyarr=p->sky->array;
  struct gal_tile_two_layer_params *tl=&p->cp.tl;

  /* It is the same size as the input or a single value. */
  if( gal_dimension_is_different(p->values, p->sky)==0 || p->sky->size==1)
    {
      s=p->sky->array;
      ff = (f=p->values->array) + p->values->size;
      if(p->sky->size==1) { if(*s!=0.0) do *f-=*s;   while(++f<ff); }
      else                              do *f-=*s++; while(++f<ff);
    }

  /* It is the same size as the number of tiles. */
  else if( tl->tottiles==p->sky->size )
    {
      /* Go over all the tiles. */
      for(tid=0; tid<tl->tottiles; ++tid)
        {
          /* For easy reading. */
          tile=&tl->tiles[tid];

          /* Subtract the Sky value from the input image. */
          GAL_TILE_PARSE_OPERATE(tile, NULL, 0, 0, {*i-=skyarr[tid];});
        }
    }

  /* The size must have been checked before, so if control reaches here, we
     have a bug! */
  else
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. For some reason, the size doesn't match", __func__,
          PACKAGE_BUGREPORT);

  /* Inform the user that this operation is done (if necessary). */
  if(!p->cp.quiet)
    printf("  - Sky subtracted from input values.\n");
}





static void
ui_preparations_read_inputs(struct mkcatalogparams *p)
{
  size_t one=1;
  gal_data_t *zero;
  gal_data_t *column;
  int need_values=0, need_sky=0, need_std=0;

  /* See which inputs are necessary. */
  ui_necessary_inputs(p, &need_values, &need_sky, &need_std);


  /* If the values dataset is necessary, read it in and set the units of
     the columns from it (if it has any). */
  if(need_values)
    {
      /* Make sure the HDU is also given. */
      if(p->valueshdu==NULL)
        error(EXIT_FAILURE, 0, "%s: no HDU/extension provided for the "
              "VALUES dataset. Atleast one column needs this dataset. "
              "Please use the '--valueshdu' option to give a specific HDU "
              "using its number (counting from zero) or name. If the "
              "dataset is in another file, please use '--valuesfile' to "
              "give the filename", p->usedvaluesfile);

      /* Read the values dataset. */
      p->values=gal_array_read_one_ch_to_type(p->usedvaluesfile, p->valueshdu,
                                              NULL, GAL_TYPE_FLOAT32,
                                              p->cp.minmapsize,
                                              p->cp.quietmmap);
      p->values->ndim=gal_dimension_remove_extra(p->values->ndim,
                                                 p->values->dsize, NULL);

      /* Make sure it has the correct size. */
      if( gal_dimension_is_different(p->objects, p->values) )
        error(EXIT_FAILURE, 0, "'%s' (hdu: %s) and '%s' (hdu: %s) have a"
              "different dimension/size", p->usedvaluesfile, p->valueshdu,
              p->objectsfile, p->cp.hdu);

      /* Initially, 'p->hasblank' was set based on the objects image, but
         it may happen that the objects image only has zero values for
         blank pixels, so we'll also do a check on the input image. */
      p->hasblank = gal_blank_present(p->values, 1);

      /* Reset the units of the value-based columns if the input dataset
         has defined units. */
      if(p->values->unit)
        {
          for(column=p->objectcols; column!=NULL; column=column->next)
            if( !strcmp(column->unit, MKCATALOG_NO_UNIT) )
              { free(column->unit);
                gal_checkset_allocate_copy(p->values->unit, &column->unit); }
          for(column=p->clumpcols; column!=NULL; column=column->next)
            if( !strcmp(column->unit, MKCATALOG_NO_UNIT) )
              { free(column->unit);
                gal_checkset_allocate_copy(p->values->unit, &column->unit); }
        }
    }



  /* Read the Sky image and check its size. */
  if(p->subtractsky || need_sky)
    {
      /* If it wasn't a number, read the dataset into memory. */
      if(p->sky==NULL)
        {
          /* Make sure the HDU is also given. */
          if(p->skyhdu==NULL)
            error(EXIT_FAILURE, 0, "%s: no HDU/extension provided for the "
                  "SKY dataset. Atleast one column needs this dataset, or "
                  "you have asked to subtract the Sky from the values.\n\n"
                  "Please use the '--skyhdu' option to give a specific HDU "
                  "using its number (counting from zero) or name. If the "
                  "dataset is in another file, please use '--skyin' to "
                  "give the filename", p->usedskyfile);

          /* Read the Sky dataset. */
          p->sky=gal_array_read_one_ch_to_type(p->usedskyfile, p->skyhdu,
                                               NULL, GAL_TYPE_FLOAT32,
                                               p->cp.minmapsize,
                                               p->cp.quietmmap);
          p->sky->ndim=gal_dimension_remove_extra(p->sky->ndim,
                                                  p->sky->dsize, NULL);

          /* Check its size and prepare tile structure. */
          ui_preparation_check_size_read_tiles(p, p->sky, p->usedskyfile,
                                               p->skyhdu);
        }

      /* Subtract the Sky value. */
      if(p->subtractsky) ui_subtract_sky(p);
    }


  /* Read the Sky standard deviation dataset (if it wasn't already given as
     a number) and check its size. */
  if(need_std && p->std==NULL)
    {
      /* Make sure the HDU is also given. */
      if(p->stdhdu==NULL)
        error(EXIT_FAILURE, 0, "%s: no HDU/extension provided for the "
              "SKY STANDARD DEVIATION dataset.\n\n"
              "Atleast one column needs this dataset. Please use the "
              "'--stdhdu' option to give a specific HDU using its number "
              "(counting from zero) or name. If the dataset is in another "
              "file, please use '--stdin' to give the filename",
              p->usedstdfile);

      /* Read the Sky standard deviation image into memory. */
      p->std=gal_array_read_one_ch_to_type(p->usedstdfile, p->stdhdu,
                                           NULL, GAL_TYPE_FLOAT32,
                                           p->cp.minmapsize, p->cp.quietmmap);
      p->std->ndim=gal_dimension_remove_extra(p->std->ndim,
                                              p->std->dsize, NULL);

      /* Check its size and prepare tile structure. */
      ui_preparation_check_size_read_tiles(p, p->std, p->usedstdfile,
                                           p->stdhdu);
    }



  /* Sanity checks on upper-limit measurements. */
  if(p->upperlimit)
    {
      /* If an upperlimit check was requested, make sure the object number
         is not larger than the maximum number of labels. */
      if(p->checkuplim[0] != GAL_BLANK_INT32
         && p->checkuplim[0] > p->numobjects)
        error(EXIT_FAILURE, 0, "%d (object identifier for the "
              "'--checkuplim' option) is larger than the number of "
              "objects in the input labels (%zu)", p->checkuplim[0],
              p->numobjects);

      /* Read the mask file if it was given. */
      if(p->upmaskfile)
        {
          /* Make sure the HDU for the mask image is given. */
          if(p->upmaskhdu==NULL)
            error(EXIT_FAILURE, 0, "%s: no HDU/extension provided, please "
                  "use the '--upmaskhdu' option to specify a specific HDU "
                  "using its number (counting from zero) or name",
                  p->upmaskfile);

          /* Read the mask image. */
          p->upmask = gal_array_read_one_ch(p->upmaskfile, p->upmaskhdu,
                                            NULL, p->cp.minmapsize,
                                            p->cp.quietmmap);
          p->upmask->ndim=gal_dimension_remove_extra(p->upmask->ndim,
                                                     p->upmask->dsize,
                                                     NULL);

          /* Check its size. */
          if( gal_dimension_is_different(p->objects, p->upmask) )
            error(EXIT_FAILURE, 0, "'%s' (hdu: %s) and '%s' (hdu: %s) have a"
                  "different dimension/size", p->upmaskfile, p->upmaskhdu,
                  p->objectsfile, p->cp.hdu);

          /* If it isn't an integer type, report an error. */
          if( p->upmask->type==GAL_TYPE_FLOAT32
              || p->upmask->type==GAL_TYPE_FLOAT64 )
            error(EXIT_FAILURE, 0, "%s (hdu: %s) has a %s numerical data "
                  "type. Only integer type inputs are acceptable as a mask."
                  "If the values are indeed integers, only placed in a "
                  "floating point container, you can use Gnuastro's "
                  "Arithmetic program to conver the numeric data type",
                  p->upmaskfile, p->upmaskhdu,
                  gal_type_name(p->upmask->type, 1));

          /* Convert the mask to a uint8_t: with a 1 for all non-zero
             pixels and 0 for zero pixels. */
          zero=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &one, NULL, 1, -1,
                              1, NULL, NULL, NULL);
          p->upmask=gal_arithmetic(GAL_ARITHMETIC_OP_NE, 1,
                                   ( GAL_ARITHMETIC_INPLACE
                                     | GAL_ARITHMETIC_FREE
                                     | GAL_ARITHMETIC_NUMOK ),
                                   p->upmask, zero);
        }
    }
}





/* The necessary keywords from the objects or clumps image were read when
   we were reading them. They were necessary during the
   pre-processing. Here, we'll read the image from  */
static void
ui_preparations_read_keywords(struct mkcatalogparams *p)
{
  gal_data_t *tmp;
  float std, minstd;
  gal_data_t *keys=NULL;

  /* Set the counts-per-second correction. */
  if(p->std)
    {
      if(p->std->size>1)
        {
          /* Read the keywords from the standard deviation image. */
          keys=gal_data_array_calloc(2);
          keys[0].next=&keys[1];
          keys[0].name="MINSTD";              keys[1].name="MEDSTD";
          keys[0].type=GAL_TYPE_FLOAT32;      keys[1].type=GAL_TYPE_FLOAT32;
          keys[0].array=&minstd;              keys[1].array=&p->medstd;
          gal_fits_key_read(p->usedstdfile, p->stdhdu, keys, 0, 0);

          /* If the two keywords couldn't be read. We don't want to slow down
             the user for the median (which needs sorting). So we'll just
             calculate the minimum which is necessary for the 'p->cpscorr'. */
          if(keys[1].status) p->medstd=NAN;
          if(keys[0].status)
            {
              /* Calculate the minimum STD. */
              tmp=gal_statistics_minimum(p->std);
              minstd=*((float *)(tmp->array));
              gal_data_free(tmp);

              /* If the units are in variance, then take the square root. */
              if(p->variance) minstd=sqrt(minstd);
            }
          p->cpscorr = minstd>1 ? 1.0f : minstd;

          /* Clean up. */
          keys[0].name=keys[1].name=NULL;
          keys[0].array=keys[1].array=NULL;
          gal_data_array_free(keys, 2, 1);
        }
      else
        {
          std=((float *)(p->std->array))[0];
          p->cpscorr=std>1 ? 1.0f : std;
        }
    }
}





/* When both catalogs need to be made, we need a separator, the output
   names will either be built based on the input name or output name (if
   given). In both cases, the operations are the same, just the base name
   differs. So to keep things clean, we have defined this function. */
static void
ui_preparations_both_names(struct mkcatalogparams *p)
{
  char *basename, *suffix=".fits";
  uint8_t keepinputdir=p->cp.keepinputdir;  /* See below. */

  /* Set the type ending. */
  if(p->cp.output)
    {
      /* When the user has specified a name, any possible directories in
         that name must be respected. So we have kept the actual
         'keepinputdir' value in a temporary variable above and set it to 1
         only for this operation. Later we set it back to what it was. */
      p->cp.keepinputdir=1;

      /* Set the base name (if necessary). */
      basename = p->cp.output;

      /* FITS speicifc preparations. */
      if( gal_fits_name_is_fits(p->cp.output) )
        {
          /* The output file name that the user has given supersedes the
             'tableformat' argument. In this case, the filename is a FITS
             file, so if 'tableformat' is a text file, we will change it to
             a default binary FITS table. */
          if( p->cp.tableformat==GAL_TABLE_FORMAT_TXT )
            p->cp.tableformat=GAL_TABLE_FORMAT_BFITS;
        }
    }
  else
    {
      /* Note that the suffix is not used in the text table outputs, so it
         doesn't matter if the output table is not FITS. */
      suffix="_cat.fits";
      basename = p->objectsfile;
    }


  /* Set the final filename. If the output is a text file, we need two
     files. But when its a FITS file we want to make a multi-extension FITS
     file. */
  if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
    {
      p->objectsout=gal_checkset_automatic_output(&p->cp, basename, "_o.txt");
      p->clumpsout=gal_checkset_automatic_output(&p->cp, basename, "_c.txt");
    }
  else
    {
      p->objectsout=gal_checkset_automatic_output(&p->cp, basename, suffix);
      p->clumpsout=p->objectsout;
    }

  /* Revert 'keepinputdir' to what it was. */
  p->cp.keepinputdir=keepinputdir;
}





/* Set the output name. */
static void
ui_preparations_outnames(struct mkcatalogparams *p)
{
  char *suffix;
  uint8_t keepinputdir=p->cp.keepinputdir;

  /* The process differs if an output filename has been given. */
  if(p->cp.output)
    {
      /* If the output name is a FITS file, then
         'gal_tableintern_check_fits_format' will see if the tableformat
         corresponds to a FITS table or not. If the output name isn't a
         FITS file then the current value of 'p->cp.tableformat' is
         irrelevant and it must be set to text. We use this value in the
         end to determine specific features. */
      if( gal_fits_name_is_fits(p->cp.output) )
        gal_tableintern_check_fits_format(p->cp.output, p->cp.tableformat);
      else
        p->cp.tableformat=GAL_TABLE_FORMAT_TXT;

      /* If a clumps image is present, then we have two outputs. */
      if(p->clumps) ui_preparations_both_names(p);
      else
        {
          gal_checkset_writable_remove(p->cp.output, 0, p->cp.dontdelete);
          gal_checkset_allocate_copy(p->cp.output, &p->objectsout);
        }
    }
  else
    {
      /* Both clumps and object catalogs are necessary. */
      if(p->clumps) ui_preparations_both_names(p);

      /* We only need one objects catalog. */
      else
        {
          suffix = ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                     ? "_cat.txt" : "_cat.fits" );
          p->objectsout=gal_checkset_automatic_output(&p->cp, p->objectsfile,
                                                      suffix);
        }
    }

  /* If an upperlimit check image is requsted, then set its filename. */
  if(p->checkuplim[0]!=GAL_BLANK_INT32)
    {
      /* See if the directory should be respected. */
      p->cp.keepinputdir = p->cp.output ? 1 : p->cp.keepinputdir;

      /* Set the suffix. */
      suffix = ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                 ? "_upcheck.txt" : "_upcheck.fits" );

      /* Set the file name. */
      p->upcheckout=gal_checkset_automatic_output(&p->cp,
                                                  ( p->cp.output
                                                    ? p->cp.output
                                                    : p->objectsfile),
                                                  suffix);

      /* Set 'keepinputdir' to what it was before. */
      p->cp.keepinputdir=keepinputdir;
    }

  /* Just to avoid bugs ('p->cp.output' must no longer be used), we'll free
     it and set it to NULL.*/
  free(p->cp.output);
  p->cp.output=NULL;
}





/* When a spectrum is requested, the slice information (slice number and
   slice WCS) is common to all different spectra. So instead of calculating
   it every time, we'll just make it once here, then copy it for every
   object.

   The Slice information is going to be written in every spectrum. So we
   don't want it to take too much space. Therefore, only when the number of
   slices is less than 65535 (2^16-1), will we actually use a 32-bit
   integer type for the slice number column.
*/
static void
ui_preparations_spectrum_wcs(struct mkcatalogparams *p)
{
  double *xarr, *yarr, *zarr;
  gal_data_t *x, *y, *z, *coords;
  size_t i, numslices=p->objects->dsize[0];
  size_t slicenumtype=numslices>=65535 ? GAL_TYPE_UINT32 : GAL_TYPE_UINT16;

  /* A small sanity check. */
  if(p->objects->ndim!=3)
    error(EXIT_FAILURE, 0, "%s (hdu %s) is a %zuD dataset, but '--spectrum' "
          "is currently only defined on 3D datasets", p->objectsfile,
          p->cp.hdu, p->objects->ndim);

  /* Allocate space for the slice number as well as the X and Y positions
     for WCS conversion. Note that the 'z' axis is going to be converted to
     WCS later, so we'll just give it the basic information now.*/
  x=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 0,
                   p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);
  y=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 0,
                   p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);
  z=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numslices, NULL, 0,
                   p->cp.minmapsize, p->cp.quietmmap, p->ctype[2],
                   p->objects->wcs->cunit[2], "Slice WCS coordinates.");

  /* Write values into the 3 coordinates. */
  xarr=x->array; yarr=y->array; zarr=z->array;
  for(i=0;i<numslices;++i) { zarr[i]=i+1; xarr[i]=yarr[i]=1; }


  /* Convert the coordinates to WCS. We are doing this inplace to avoid too
     much memory/speed consumption. */
  coords=x;
  coords->next=y;
  coords->next->next=z;
  gal_wcs_img_to_world(coords, p->objects->wcs, 1);

  /* For a check.
  for(i=0;i<numslices;++i)
    printf("%g, %g, %g\n", xarr[i], yarr[i], zarr[i]);
  exit(0);
  */

  /* Allocate the slice counter array (we are doing it again because we
     want it to be in integer type now). */
  p->specsliceinfo=gal_data_alloc(NULL, slicenumtype, 1, &numslices, NULL, 0,
                                  p->cp.minmapsize, p->cp.quietmmap, "SLICE",
                                  "counter",
                                  "Slice number in cube (counting from 1).");
  if(p->specsliceinfo->type==GAL_TYPE_UINT16)
    for(i=0;i<numslices;++i) ((uint16_t *)(p->specsliceinfo->array))[i]=i+1;
  else
    for(i=0;i<numslices;++i) ((uint32_t *)(p->specsliceinfo->array))[i]=i+1;

  /* Set the slice WCS column information. Note that 'z' is now the WCS
     coordinate value of the third dimension, and to avoid wasting extra
     space (this column is repeated one very object's spectrum), we'll
     convert it to a 32-bit floating point dataset. */
  p->specsliceinfo->next=gal_data_copy_to_new_type(z, GAL_TYPE_FLOAT32);

  /* For a final check.
  gal_table_write(p->specsliceinfo, NULL, GAL_TABLE_FORMAT_BFITS,
                  "specsliceinfo.fits", "test-debug",0);
  */

  /* Clean up. */
  gal_data_free(x);
  gal_data_free(y);
  gal_data_free(z);
}





/* Sanity checks and preparations for the upper-limit magnitude. */
static void
ui_preparations_upperlimit(struct mkcatalogparams *p)
{
  size_t i, c=0;

  /* Check if the given range has the same number of elements as dimensions
     in the input. */
  if(p->uprange)
    {
      for(i=0;p->uprange[i]!=-1;++i) ++c;
      if(c!=p->objects->ndim)
        error(EXIT_FAILURE, 0, "%zu values given to '--uprange', but input "
              "has %zu dimensions", c, p->objects->ndim);
    }

  /* Check the number of random samples. */
  if( p->upnum < MKCATALOG_UPPERLIMIT_MINIMUM_NUM )
    error(EXIT_FAILURE, 0, "%zu not acceptable as '--upnum'. The minimum "
          "acceptable number of random samples for the upper limit "
          "magnitude is %d", p->upnum, MKCATALOG_UPPERLIMIT_MINIMUM_NUM);

  /* Check if sigma-clipping parameters have been given. */
  if( isnan(p->upsigmaclip[0]) )
    error(EXIT_FAILURE, 0, "'--upsigmaclip' is mandatory for measuring "
          "the upper-limit magnitude. It takes two numbers separated by "
          "a comma. The first is the multiple of sigma and the second is "
          "the aborting criteria: <1: tolerance level, >1: number of "
          "clips");

  /* Check if the sigma multiple is given. */
  if( isnan(p->upnsigma) )
    error(EXIT_FAILURE, 0, "'--upnsigma' is mandatory for measuring the "
          "upperlimit magnitude. Its value is the multiple of final sigma "
          "that is reported as the upper-limit");

  /* Set the random number generator. */
  p->rng=gal_checkset_gsl_rng(p->envseed, &p->rng_name, &p->rng_seed);

  /* Keep the minimum and maximum values of the random number generator. */
  p->rngmin=gsl_rng_min(p->rng);
  p->rngdiff=gsl_rng_max(p->rng)-p->rngmin;
}









void
ui_preparations(struct mkcatalogparams *p)
{
  /* If no columns are requested, then inform the user. */
  if(p->columnids==NULL && p->spectrum==0)
    error(EXIT_FAILURE, 0, "no measurements requested! Please run again "
          "with '--help' for the possible list of measurements");


  /* Set the actual filenames to use. */
  ui_set_filenames(p);


  /* Read the main input (the objects image). */
  ui_read_labels(p);


  /* Prepare the output columns. */
  columns_define_alloc(p);


  /* Read the inputs. */
  ui_preparations_read_inputs(p);


  /* Read the helper keywords from the inputs and if they aren't present
     then calculate the necessary parameters. */
  ui_preparations_read_keywords(p);


  /* Set the output filename(s). */
  ui_preparations_outnames(p);


  /* If a spectrum is requested, generate the two WCS columns. */
  if(p->spectrum)
    {
      ui_preparations_spectrum_wcs(p);
      p->spectra=gal_data_array_calloc(p->numobjects);
    }


  /* Allocate the reference random number generator and seed values. It
     will be cloned once for every thread. If the user hasn't called
     'envseed', then we want it to be different for every run, so we need
     to re-set the seed. */
  if(p->upperlimit) ui_preparations_upperlimit(p);

  if( p->hasmag && isnan(p->zeropoint) )
    error(EXIT_FAILURE, 0, "no zeropoint specified");


  /* Prepare the two internal arrays necessary to sort the clumps catalog
     by object and clump IDs. We are allocating and filling these in
     separately (and not using the actual output columns that have the same
     values), because playing with the output columns can cause bad
     bugs. If the user wants performance, they are encouraged to run
     MakeCatalog with '--noclumpsort' and avoid the whole process all
     together. */
  if(p->clumps && !p->noclumpsort && p->cp.numthreads>1)
    {
      p->hostobjid_c=gal_pointer_allocate(GAL_TYPE_SIZE_T,
                                          p->clumpcols->size, 0, __func__,
                                          "p->hostobjid_c");
      p->numclumps_c=gal_pointer_allocate(GAL_TYPE_SIZE_T,
                                          p->objectcols->size, 0, __func__,
                                          "p->numclumps_c");
    }
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkcatalogparams *p)
{
  char *tmp;
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program ('args.h')
     and for the common options to all Gnuastro ('commonopts.h'). We want
     to directly put the pointers to the fields in 'p' and 'cp', so we are
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
     their relations prior to printing. */
  ui_read_check_only_options(p);


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


  /* If the output is a FITS table, prepare all the options as FITS
     keywords to write in output later. */
  if(gal_fits_name_is_fits(p->objectsout))
      gal_options_as_fits_keywords(&p->cp);


  /* Inform the user. */
  if(!p->cp.quiet)
    {
      /* Write the information. */
      printf(PROGRAM_NAME" "PACKAGE_VERSION" started on %s",
             ctime(&p->rawtime));
      printf("  - Using %zu CPU thread%s\n", p->cp.numthreads,
             p->cp.numthreads==1 ? "." : "s.");
      printf("  - Objects: %s (hdu: %s)\n", p->objectsfile, p->cp.hdu);
      if(p->clumps)
        printf("  - Clumps:  %s (hdu: %s)\n", p->usedclumpsfile,
               p->clumpshdu);
      if(p->relabclumps)
        printf("  - RELABELED CLUMPS (no NUMLABS in original): %s\n",
               p->relabclumps);
      if(p->values)
        printf("  - Values:  %s (hdu: %s)\n", p->usedvaluesfile,
               p->valueshdu);

      if(p->subtractsky || p->sky)
        {
          if(p->sky->size==1)
            printf("  - Sky: %g\n", *((float *)(p->sky->array)) );
          else
            printf("  - Sky: %s (hdu: %s)\n", p->usedskyfile, p->skyhdu);
        }

      if(p->std)
        {
          tmp = p->variance ? "VAR" : "STD";
          if(p->std->size==1)
            printf("  - Sky %s: %g\n", tmp, *((float *)(p->std->array)) );
          else
            printf("  - Sky %s: %s (hdu: %s)\n", tmp, p->usedstdfile,
                   p->stdhdu);
        }

      if(p->upmaskfile)
        printf("  - Upper limit magnitude mask: %s (hdu: %s)\n",
               p->upmaskfile, p->cp.hdu);
      if(p->upperlimit)
        {
          printf("  - Random number generator name: %s\n", p->rng_name);
          printf("  - Random number generator seed: %lu\n", p->rng_seed);
        }
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct mkcatalogparams *p, struct timeval *t1)
{
  size_t d, i;

  /* The temporary arrays for WCS coordinates. */
  if(p->wcs_vo ) gal_list_data_free(p->wcs_vo);
  if(p->wcs_vc ) gal_list_data_free(p->wcs_vc);
  if(p->wcs_go ) gal_list_data_free(p->wcs_go);
  if(p->wcs_gc ) gal_list_data_free(p->wcs_gc);
  if(p->wcs_vcc) gal_list_data_free(p->wcs_vcc);
  if(p->wcs_gcc) gal_list_data_free(p->wcs_gcc);

  /* Free the types of the WCS coordinates (for catalog meta-data). */
  if(p->ctype)
    {
      for(d=0;d<p->objects->ndim;++d)
        free(p->ctype[d]);
      free(p->ctype);
    }

  /* If a random number generator was allocated, free it. */
  if(p->rng) gsl_rng_free(p->rng);

  /* Free output names. */
  if(p->clumpsout && p->clumpsout!=p->objectsout)
    free(p->clumpsout);
  free(p->objectsout);

  /* Free the allocated arrays: */
  free(p->skyhdu);
  free(p->stdhdu);
  free(p->cp.hdu);
  free(p->oiflag);
  free(p->ciflag);
  free(p->skyfile);
  free(p->stdfile);
  free(p->clumpshdu);
  free(p->valueshdu);
  free(p->clumpsfile);
  free(p->valuesfile);
  free(p->hostobjid_c);
  free(p->numclumps_c);
  gal_data_free(p->sky);
  gal_data_free(p->std);
  gal_data_free(p->values);
  gal_data_free(p->upmask);
  gal_data_free(p->clumps);
  gal_data_free(p->objects);
  if(p->outlabs) free(p->outlabs);
  gal_list_data_free(p->clumpcols);
  gal_list_data_free(p->objectcols);
  gal_list_data_free(p->specsliceinfo);
  if(p->upcheckout) free(p->upcheckout);
  gal_data_array_free(p->tiles, p->numobjects, 0);

  /* Clean up the spectra. */
  if(p->spectra)
    {
      /* Note that each element of the array is the first node in a list of
         datasets. So we can't free the first one with
         'gal_list_data_free', we'll delete all the nodes after it in the
         loop. */
      for(i=0;i<p->numobjects;++i)
        {
          gal_list_data_free( p->spectra[i].next );
          p->spectra[i].next=NULL;
          gal_data_free_contents(&p->spectra[i]);
        }
      gal_data_array_free(p->spectra, p->numobjects, 0);
    }

  /* If the Sky or its STD image were given in tiles, then we defined a
     tile structure to deal with them. The initialization of the tile
     structure is checked with its 'ndim' element. */
  if(p->cp.tl.ndim) gal_tile_full_free_contents(&p->cp.tl);

  /* If an upper limit range warning is necessary, print it here. */
  if(p->uprangewarning)
    fprintf(stderr, "\nMore on the WARNING-UPPERLIMIT(s) above: "
            "In order to obtain a good/robust random distribution (and "
            "thus a reliable upper-limit measurement), it is necessary "
            "to have a sufficienty wide enough range (in each dimension). "
            "As mentioned in the warning(s) above, the available "
            "range for random sampling of some of the labels in this "
            "input is less than double their length. If the input is taken "
            "from a larger dataset, this issue can be solved by using a "
            "larger part of it. You can also run MakeCatalog with "
            "'--checkuplim' to see the distribution for a special "
            "object or clump as a table and personally inspect its "
            "reliability. \n\n");

  /* Print the final message. */
  if(!p->cp.quiet)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
