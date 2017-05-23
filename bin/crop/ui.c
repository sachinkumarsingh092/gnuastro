/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro/wcs.h>
#include <gnuastro/list.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "onecrop.h"
#include "wcsmode.h"
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
args_doc[] = "[ASCIIcatalog] ASTRdata ...";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will create cutouts, "
  "thumbnails, postage stamps or crops of region(s) from input image(s) "
  "using image or celestial coordinates. If muliple crops are desired, a "
  "catalog must be provided. When in WCS mode, if the cut out covers more "
  "than one input image, all overlapping input images will be stitched in "
  "the output.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Option groups particular to this program. */
enum program_args_groups
{
  ARGS_GROUP_CENTER_GENERAL = GAL_OPTIONS_GROUP_AFTER_COMMON,
  ARGS_GROUP_CENTER_SINGLE,
  ARGS_GROUP_CENTER_CATALOG,
  ARGS_GROUP_REGION,
};




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct cropparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->poptions           = program_options;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;


  /* Initalize necessary parameters. */
  p->xc=p->yc=p->ra=p->dec=NAN;
  p->mode         = IMGCROP_MODE_INVALID;
  cp->searchin    = GAL_TABLE_SEARCH_INVALID;


  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_HDU:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          cp->coptions[i].doc="Extension name or number of (all) input(s).";
          break;

        case GAL_OPTIONS_KEY_MINMAPSIZE:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;

        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_IGNORECASE:
          cp->coptions[i].group=ARGS_GROUP_CENTER_CATALOG;
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
  struct cropparams *p = state->input;

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
      gal_list_str_add(&p->inputs, arg, 0);
      ++p->numin;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct cropparams *p)
{
  int checksum;

  /* Read the mode from the string the user specified. */
  if(p->modestr)
    {
      if      (!strcmp(p->modestr, "img"))   p->mode=IMGCROP_MODE_IMG;
      else if (!strcmp(p->modestr, "wcs"))   p->mode=IMGCROP_MODE_WCS;
      else
        error(EXIT_FAILURE, 0, "`%s' (value to `--mode') not recognized as "
              "an input mode. Recognized values are `img' and `wcs'. This "
              "option is necessary when the inputs are not sufficient to "
              "identify the nature of the coordinates.\n\n"
              "Please run the following command for more information "
              "(press the `SPACE' key to go down and `q' to return to the "
              "command-line):\n\n"
              "    $ info gnuastro \"Crop modes\"\n", p->modestr);
    }




  /* Check that if one of the coordinates is given, the other is also
     given. Note that if both are not given, their sum will be 0, if both
     are given, then the sum will be 2. If only one is given, then the sum
     will be 1. */
  if( (isnan(p->ra) + isnan(p->dec)) == 1 )
    error(EXIT_FAILURE, 0, "no `--%s' given, it is mandatory with the `--%s' "
          "option", isnan(p->ra)?"ra":"dec", isnan(p->ra)?"dec":"ra");

  if( (isnan(p->xc) + isnan(p->yc)) == 1 )
    error(EXIT_FAILURE, 0, "no `--%s' given, it is mandatory with the `--%s' "
          "option", isnan(p->xc)?"xc":"yc", isnan(p->xc)?"yc":"xc");

  if( ((p->racol!=NULL) + (p->deccol!=NULL)) == 1 )
    error(EXIT_FAILURE, 0, "no `--%scol' given, it is mandatory with the "
          "`--%scol' option", p->racol?"dec":"ra",
          p->racol?"ra":"dec");

  if( ((p->xcol!=NULL) + (p->ycol!=NULL)) == 1 )
    error(EXIT_FAILURE, 0, "no `--%scol' given, it is mandatory with the "
          "`--%scol' option", p->xcol?"y":"x", p->xcol?"x":"y");



  /* Make sure that the single crop modes are not called together. */
  checksum = ( (!isnan(p->xc)) + (!isnan(p->ra)) + (p->catname!=NULL)
               + (p->section!=NULL) + (p->polygon!=NULL) );
  switch(checksum)
    {
    case 0:
      error(EXIT_FAILURE, 0, "no crop definition. You can use any of the "
            "following options to define the crop(s): (`--xc', `--yc'), "
            "(`--ra', `--dec'), `--catalog', `--section', `--polygon'. "
            "Please run this command for more information:\n\n"
            "    $ info gnuastro \"Crop modes\"\n");
    case 1:
      /* Everything is ok, just ignore the switch structure. */
      break;
    default:
      error(EXIT_FAILURE, 0, "more than one crop type specified. In each "
            "run, only one crop definition is acceptable on the "
            "command-line, or in configuration files. You have called "
            "%s%s%s%s%s\b\b.",
            !isnan(p->xc) ? "(`--xc', `--yc'), " : "",
            !isnan(p->ra) ? "(`--ra', `--dec'), " : "",
            p->catname!=NULL ? "`--catalog', " : "",
            p->section!=NULL ? "`--section', " : "",
            p->polygon!=NULL ? "`--polygon', " : "");
    }



  /* Sanity checks and mode setting based on the desired crop. */
  if(p->catname)
    {
      /* If the searchin option has been given. */
      if(p->cp.searchin==GAL_TABLE_SEARCH_INVALID)
        error(EXIT_FAILURE, 0, "%s: no field specified to search for "
              "columns. Please use the `--searchin' option to specify "
              "which column meta-data you would like to search in: `name', "
              "`unit' and `comment'. You may also select columns by their "
              "number, which won't use this option, but for complentess its "
              "best for this option to have a value", p->catname);

      /* If it is a FITS file, we need the HDU. */
      if( gal_fits_name_is_fits(p->catname) && p->cathdu==NULL )
        error(EXIT_FAILURE, 0, "%s: no hdu given. Please use the `--cathdu' "
              "option to specify which extension contains the table",
              p->catname);

      /* Atleast one of the (X,Y), and (RA,Dec) set of columns are
         necessary. Note that we have checked that they are together if
         given, so we only need to check one of the two in each couple. */
      if(p->xcol==NULL && p->racol==NULL)
        error(EXIT_FAILURE, 0, "no crop center position column given "
              "to read from the input catalog (`%s'). Please use either of "
              "(`--xcol', `--ycol') or (`--racol', `--deccol'). For more "
              "information on how to select columns in Gnuastro, please run "
              "the following command:\n\n"
              "    $ info gnuastro \"Selecting table columns\"", p->catname);

      /* If only image columns are specified, then we have Image mode, if
         only WCS columns are specified, we have WSC mode. If both are
         specified, the mode is mandatory */
      if(p->xcol && p->racol)
        {
          if(p->mode==IMGCROP_MODE_INVALID)
            error(EXIT_FAILURE, 0, "both image and WCS coordinate columns "
                  "are specified to read the center of the crops in the "
                  "input catalog (%s). You can use the `--mode=img', or "
                  "`--mode=wcs' options to specify which set of columns "
                  "should be used", p->catname);
        }
      else
        p->mode = p->xcol ? IMGCROP_MODE_IMG : IMGCROP_MODE_WCS;
    }
  else if(p->polygon)
    {
      if(p->mode==IMGCROP_MODE_INVALID)
        error(EXIT_FAILURE, 0, "the polygon option works in both image and "
              "wcs mode, please specify how the vertices should be "
              "interpreted with `--mode=img', or `--mode=wcs' options to "
              "specify which set of columns should be used");
    }
  else if(!isnan(p->xc))
    p->mode = IMGCROP_MODE_IMG;
  else if(!isnan(p->ra))
    p->mode = IMGCROP_MODE_WCS;
  else if( p->section)
    p->mode = IMGCROP_MODE_IMG;


  /* Parse the polygon vertices if they are given to make sure that it is
     in the proper format. */
  if(p->polygon)
    {
      crop_polygonparser(p);
      if(p->nvertices<3)
        error(EXIT_FAILURE, 0, "a polygon has to have 3 or more vertices, "
              "you have only given %zu (%s)", p->nvertices, p->polygon);
      if(p->outpolygon && p->numin>1)
        error(EXIT_FAILURE, 0, "currently in WCS mode, outpolygon can only "
              "be set to zero when there is one image, you have given %zu "
              "images. For multiple images the region will be very large. "
              "It is best if you first crop out the larger region you want "
              "into one image, then mask the polygon", p->numin);
    }
  else
    p->wpolygon=p->ipolygon=NULL;


  /* If we are in WCS mode, noblanks must be off */
  if(p->mode==IMGCROP_MODE_WCS && p->noblank)
    error(EXIT_FAILURE, 0, "`--noblanks` (`-b`) is only for image mode. "
          "You have called it with WCS mode");
}





static void
ui_check_options_and_arguments(struct cropparams *p)
{
  /* Make sure we do actually have inputs. */
  if(p->inputs==NULL)
    error(EXIT_FAILURE, 0, "no input file given");

  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->cp.hdu==NULL )
    error(EXIT_FAILURE, 0, "no HDU specified. When the input is a FITS "
          "file, a HDU must also be specified, you can use the `--hdu' "
          "(`-h') option and give it the HDU number (starting from "
          "zero), extension name, or anything acceptable by CFITSIO");

  /* If in image mode, there should only be one input image. */
  if(p->mode==IMGCROP_MODE_IMG && p->numin>1)
    error(EXIT_FAILURE, 0, "in image mode, only one input image may be "
          "specified");

  /* If no output name is given, set it to the current directory. */
  if(p->cp.output==NULL)
    gal_checkset_allocate_copy("./", &p->cp.output);

  /* Only catalog mode needs multiple threads and a directory for the
     output. */
  if(p->catname)
    {
      /* When multiple threads need to access a file, CFITSIO needs to be
         configured with the `--enable-reentrant` option. */
      if(p->cp.numthreads>1 && fits_is_reentrant()==0)
        error(EXIT_FAILURE, 0, "CFITSIO was not configured with the "
              "`--enable-reentrant` option but you have asked to crop "
              "on %zu threads.\n\nPlease configure, make and install CFITSIO "
              "again with this flag. Alternatively, to avoid this error "
              "you can set the number of threads to 1 by adding the "
              "`--numthreads=1` or `-N1` options. Please run the following "
              "command to learn more about configuring CFITSIO:\n\n"
              "    $ info gnuastro CFITSIO", p->cp.numthreads);

      /* Make sure the given output is a directory. */
      gal_checkset_check_dir_write_add_slash(&p->cp.output);
    }
  else
    {
      p->cp.numthreads=1;
      p->outnameisfile=gal_checkset_dir_0_file_1(p->cp.output,
                                                 p->cp.dontdelete);
    }
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_read_cols(struct cropparams *p)
{
  char *colname;
  size_t counter=0;
  int toomanycols=0;
  gal_list_str_t *colstrs=NULL;
  gal_data_t *cols, *tmp, *corrtype=NULL;
  char *ax1col = p->mode==IMGCROP_MODE_IMG ? p->xcol : p->racol;
  char *ax2col = p->mode==IMGCROP_MODE_IMG ? p->ycol : p->deccol;

  /* Specify the order of columns. */
  if(p->namecol)
    gal_list_str_add(&colstrs, p->namecol, 0);
  gal_list_str_add(&colstrs, ax1col, 0);
  gal_list_str_add(&colstrs, ax2col, 0);

  /* Read the desired columns from the file. */
  cols=gal_table_read(p->catname, p->cathdu, colstrs, p->cp.searchin,
                      p->cp.ignorecase, p->cp.minmapsize);

  /* Set the number of objects. */
  p->numout=cols->size;

  /* For a sanity check, make sure that the total number of columns read is
     the same as those that were wanted (it might be more). */
  while(cols!=NULL)
    {
      /* Pop out the top node. */
      tmp=gal_list_data_pop(&cols);

      /* Note that the input was a linked list, so the output order is the
         inverse of the input order. For the position, we will store the
         values into the `x' and `y' arrays even if they are RA/Dec. */
      switch(++counter)
        {
        case 3:
          colname="crop name prefix";
          if(p->namecol)
            {
              if(tmp->type==GAL_TYPE_STRING)
                {
                  p->name=tmp->array;
                  tmp->array=NULL;
                  gal_data_free(tmp);
                }
              else
                {
                  corrtype=gal_data_copy_to_new_type_free(tmp,
                                            GAL_TYPE_STRING);
                  p->name=corrtype->array;
                }
            }
          else
            toomanycols=1;
          break;

        case 2:
          colname="first axis position";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
          p->c1=corrtype->array;
          break;

        case 1:
          colname="second axis position";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
          p->c2=corrtype->array;
          break;

        /* If the index isn't recognized, then it is larger, showing that
           there was more than one match for the given criteria */
        default:
          toomanycols=1;
        }

      /* Print an error if there were too many columns. */
      if(toomanycols)
        gal_tableintern_error_col_selection(p->catname, p->cathdu, "too many "
                                            "columns were selected by the "
                                            "given values to the options "
                                            "ending in `col'.");

      /* Sanity check and clean up.  Note that it might happen that the
         input structure is already freed. In that case, `corrtype' will be
         NULL. */
      if(corrtype)
        {
          /* Make sure there are no blank values in this column. */
          if( gal_blank_present(corrtype, 1) )
            error(EXIT_FAILURE, 0, "%s column has blank values. "
                  "Input columns cannot contain blank values", colname);

          /* Free the unnecessary sturcture information. The correct-type
             (`corrtype') data structure's array is necessary for later
             steps, so its pointer has been copied in the main program's
             structure. Hence, we should set the structure's pointer to
             NULL so the important data isn't freed.*/
          corrtype->array=NULL;
          gal_data_free(corrtype);
          corrtype=NULL;
        }
    }
}





/* Add all the columns of the log file. Just note that since this is a
   linked list, we have to add them in the opposite order. */
static void
ui_make_log(struct cropparams *p)
{
  char *comment;

  /* Return if no long file is to be created. */
  if(p->cp.log==0) return;

  /* If central pixels are filled. */
  asprintf(&comment, "Are the central pixels filled? (1: yes, 0: no, "
           "%u: not checked)", GAL_BLANK_UINT8);
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT8, 1, &p->numout,
                          NULL, 1, p->cp.minmapsize, "CENTER_FILLED",
                          "bool", comment);
  free(comment);

  /* Number of images used. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT16, 1, &p->numout,
                          NULL, 1, p->cp.minmapsize, "NUM_INPUTS", "count",
                          "Number of input images used to make this crop.");

  /* Row number in input catalog. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_STRING, 1, &p->numout,
                          NULL, 1, p->cp.minmapsize, "CROP_NAME", "name",
                          "File name of crop.");
}





void
ui_preparations(struct cropparams *p)
{
  char *msg;
  fitsfile *tmpfits;
  struct timeval t1;
  size_t input_counter;
  struct inputimgs *img;
  int status, firsttype=0;


  /* Set the initial iwidth. */
  p->iwidth[0]=p->iwidth[1]=p->iwidthin;


  /* For polygon and section, there should be no center checking. */
  if(p->polygon || p->section)
    p->checkcenter=0;


  /* Read the columns if there is a catalog, otherwise, set the number
     of images (crops) to 1.*/
  if(p->catname)
    ui_read_cols(p);
  else
    p->numout=1;


  /* Everything is ready, notify the user of the program starting. */
  if(!p->cp.quiet)
    {
      gettimeofday(&t1, NULL);
      printf(PROGRAM_NAME" started on %s", ctime(&p->rawtime));
    }


  /* Allocate space for all the input images. This is done here because
     WCSLIB is unfortunately not thread-safe when reading the WCS
     information from the FITS files. In cases where the number of cropped
     images are more than the input images, this can also be a preformance
     boost because each image information is only read once.

     The images are filled in opposite order because we used a linked list
     to read them in, which is a first in first out structure.*/
  errno=0;
  p->imgs=malloc(p->numin*sizeof *p->imgs);
  if(p->imgs==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %zu bytes for p->imgs",
          p->numin*sizeof *p->imgs);


  /* Fill in the WCS information of each image. */
  input_counter=p->numin;
  while(p->inputs)
    {
      /* Pop from the list of input images and get the info. */
      status=0;
      img=&p->imgs[--input_counter];
      img->name=gal_list_str_pop(&p->inputs);
      tmpfits=gal_fits_hdu_open_format(img->name, p->cp.hdu, 0);
      gal_fits_img_info(tmpfits, &p->type, &img->ndim, &img->dsize);
      img->wcs=gal_wcs_read_fitsptr(tmpfits, p->hstartwcs, p->hendwcs,
                                    &img->nwcs);
      if(img->wcs)
        {
          gal_wcs_decompose_pc_cdelt(img->wcs);
          status=wcshdo(0, img->wcs, &img->nwcskeys, &img->wcstxt);
          if(status)
            error(EXIT_FAILURE, 0, "wcshdo ERROR %d: %s", status,
                  wcs_errmsg[status]);
        }
      else
        if(p->mode==IMGCROP_MODE_WCS)
          error(EXIT_FAILURE, 0, "the WCS structure of %s (hdu: %s) "
                "image is not recognized. So RA and Dec cannot be used "
                "as input. You can try with pixel coordinates in the "
                "Image Mode (note that the crops will lack WCS "
                "header information)", img->name, p->cp.hdu);
      fits_close_file(tmpfits, &status);
      gal_fits_io_error(status, NULL);

      /* Make sure all the images have the same type. */
      if(firsttype==0)
        {
          firsttype=p->type;
          p->bitnul = gal_blank_alloc_write(p->type);
        }
      else
        {
          if(firsttype!=p->type)
            error(EXIT_FAILURE, 0, "%s: type is `%s' while revious image(s) "
                  "were `%s' type. All inputs must have the same pixel data "
                  "type.\n\nYou can use Gnuastro's Arithmetic program to "
                  "convert `%s' to `%s', please run this command for more "
                  "information (press `SPACE' for going down and `q' to "
                  "return to the command-line):\n\n"
                  "    $ info Arithmetic\n",
                  img->name, gal_type_name(p->type, 1),
                  gal_type_name(firsttype, 1), img->name,
                  gal_type_name(p->type, 1));
        }

      /* In WCS mode, Check resolution and get the first pixel
         positions. */
      if(p->mode==IMGCROP_MODE_WCS) wcs_check_prepare(p, img);
    }


  /* Report timing: */
  if(!p->cp.quiet)
    {
      asprintf(&msg, "Read metadata of %zu image%s.", p->numin,
              p->numin>1 ? "s" : "");
      gal_timing_report(&t1, msg, 1);
    }


  /* Prepare the log file if the user has asked for it. */
  ui_make_log(p);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct cropparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


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
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct cropparams *p, struct timeval *t1)
{
  size_t i;

  /* Free the simple arrays (if they were set). */
  if(p->c1) free(p->c1);
  if(p->c2) free(p->c2);
  if(p->cp.hdu) free(p->cp.hdu);
  if(p->cathdu) free(p->cathdu);
  if(p->catname) free(p->catname);

  /* The arguments (note that the values were not allocated). */
  gal_list_str_free(p->inputs, 0);

  /* Free the name/ array.  */
  if(p->name)
    {
      for(i=0;i<p->numout;++i)
        free(p->name[i]);
      free(p->name);
    }

  /* Free the log information. */
  if(p->cp.log) gal_list_data_free(p->log);

  /* Print the final message. */
  if(!p->cp.quiet)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
