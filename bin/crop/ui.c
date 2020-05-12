/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro/wcs.h>
#include <gnuastro/list.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>

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
args_doc[] = "[Crop-Identifier] ASTRdata ...";

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
  cp->program_struct     = p;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->poptions           = program_options;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;


  /* Initalize necessary parameters. */
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
          cp->coptions[i].group=UI_GROUP_CENTER_CATALOG;
          break;

        case GAL_OPTIONS_KEY_STDINTIMEOUT:
          cp->coptions[i].group=OPTION_HIDDEN;
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
      gal_list_str_add(&p->inputs, arg, 0);
      ++p->numin;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}




/* Parse the mode to interpret the given coordinates. */
void *
ui_parse_coordinate_mode(struct argp_option *option, char *arg,
                         char *filename, size_t lineno, void *junk)
{
  char *outstr;

  /* We want to print the stored values. */
  if(lineno==-1)
    {
      gal_checkset_allocate_copy( *(int *)(option->value)==IMGCROP_MODE_IMG
                                  ? "img" : "wcs", &outstr );
      return outstr;
    }
  else
    {
      if      (!strcmp(arg, "img")) *(int *)(option->value)=IMGCROP_MODE_IMG;
      else if (!strcmp(arg, "wcs")) *(int *)(option->value)=IMGCROP_MODE_WCS;
      else
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "'%s' (value to "
                      "'--mode') not recognized as an input mode. "
                      "Recognized values are 'img' and 'wcs'. This option "
                      "is necessary to identify the nature of your input "
                      "coordinates.\n\n"
                      "Please run the following command for more "
                      "information (press the 'SPACE' key to go down and "
                      "'q' to return to the command-line):\n\n"
                      "    $ info gnuastro \"Crop modes\"\n", arg);
      return NULL;
    }
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct cropparams *p)
{
  double *darray;
  int i, checksum;

  /* Make sure that only one of the crop definitions is given. */
  checksum = ( (p->center!=NULL)
               + (p->catname!=NULL)
               + (p->section!=NULL)
               + (p->polygon!=NULL) );
  switch(checksum)
    {
    case 0:
      error(EXIT_FAILURE, 0, "no crop definition. You can use any of the "
            "following options to define the crop(s): '--center', "
            "'--catalog', '--section', or '--polygon'. Please run this "
            "command for more information:\n\n"
            "    $ info gnuastro \"Crop modes\"\n");
    case 1:
      /* Everything is ok, just ignore the switch structure. */
      break;
    default:
      error(EXIT_FAILURE, 0, "more than one crop type specified. In each "
            "run, only one crop definition is acceptable on the "
            "command-line or in configuration files. You have called: "
            "%s%s%s%s\b\b.",
            p->center!=NULL  ? "'--center', " : "",
            p->catname!=NULL ? "'--catalog', " : "",
            p->section!=NULL ? "'--section', " : "",
            p->polygon!=NULL ? "'--polygon', " : "");
    }


  /* The width values must not be negative. */
  if(p->width)
    {
      darray=p->width->array;
      for(i=0;i<p->width->size;++i)
        if(darray[i]<=0.0f)
          error(EXIT_FAILURE, 0, "%g is <=0. The values to the '--width' "
                "option must be larger than zero. %g is input number %d to "
                "this option", darray[i], darray[i], i+1);
    }


  /* Checkcenter sanity check. */
  if(p->incheckcenter)
    {
      /* We only want a single number. */
      if(p->incheckcenter->size>1)
        error(EXIT_FAILURE, 0, "%zu values given to '--checkcenter'. This "
              "option only takes one value currently",
              p->incheckcenter->size);

      darray=p->incheckcenter->array;
      if(*darray<0.0f)
        error(EXIT_FAILURE, 0, "negative value (%f) given to "
              "'--checkcenter'. This option only takes positive values",
              *darray);
    }


  /* Section is currently only defined in Image mode. */
  if(p->section && p->mode!=IMGCROP_MODE_IMG)
    error(EXIT_FAILURE, 0, "The '--section' option is only available in "
          "image coordinate mode, currently it doesn't work with WCS mode. "
          "Please run with '--mode=img' and if necessary, change the "
          "values accordingly");


  /* Sanity checks and mode setting based on the desired crop. */
  if(p->catname)
    {
      /* If the searchin option has been given. */
      if(p->cp.searchin==GAL_TABLE_SEARCH_INVALID)
        error(EXIT_FAILURE, 0, "%s: no field specified to search for "
              "columns. Please use the '--searchin' option to specify "
              "which column meta-data you would like to search in: 'name', "
              "'unit' and 'comment'. You may also select columns by their "
              "number, which won't use this option, but for complentess its "
              "best for this option to have a value", p->catname);

      /* If it is a FITS file, we need the HDU. */
      if( gal_fits_name_is_fits(p->catname) && p->cathdu==NULL )
        error(EXIT_FAILURE, 0, "%s: no hdu given. Please use the '--cathdu' "
              "option to specify which extension contains the table",
              p->catname);

      /* Atleast one of the (X,Y), and (RA,Dec) set of columns are
         necessary. Note that we have checked that they are together if
         given, so we only need to check one of the two in each couple. */
      if(p->coordcol==NULL)
        error(EXIT_FAILURE, 0, "no crop center columns given to read from "
              "the input catalog ('%s'). Please use '--coordcol' several "
              "times (depending on dimensionality) to specify the column "
              "keeping the center position the respective dimension.\n\n"
              "For more information on how to select columns in Gnuastro, "
              "please run the following command:\n\n"
              "    $ info gnuastro \"Selecting table columns\"", p->catname);
    }


  /* Parse the polygon vertices if they are given to make sure that it is
     in the proper format. */
  if(p->polygon)
    {
      /* The number of vertices is half the total number of given values
         (currently only 2D spaces are considered so each vertice has 2
         coordinates. */
      p->nvertices=p->polygon->size/2;

      /* Basic sanity checks. */
      if(p->nvertices<3)
        error(EXIT_FAILURE, 0, "a polygon has to have 3 or more vertices, "
              "you have only given %zu", p->nvertices);
      if(p->polygonout && p->numin>1)
        error(EXIT_FAILURE, 0, "currently in WCS mode, '--polygonout' can "
              "only be set to zero when there is one image, you have given "
              "%zu images. For multiple images the region will be very "
              "large. It is best if you first crop out the larger region "
              "you want into one image, then mask the polygon", p->numin);

      /* Put the coordinates into an array while reversing their order so
         they correspond to the user's order, then put it in the right
         place.*/
      darray=p->polygon->array;
      if(p->mode==IMGCROP_MODE_IMG) {p->ipolygon=darray; p->wpolygon=NULL;  }
      else                          {p->ipolygon=NULL;   p->wpolygon=darray;}

      /* We know that the cropped region is not defined by its center. So
         it makes no sense to check if the center is blank. */
      p->checkcenter=0;
    }
  else
    p->wpolygon=p->ipolygon=NULL;


  /* If we are in WCS mode, noblanks must be off */
  if(p->mode==IMGCROP_MODE_WCS && p->noblank)
    error(EXIT_FAILURE, 0, "'--noblanks' ('-b') is only for image mode. "
          "You have called it with WCS mode");
}





static void
ui_check_options_and_arguments(struct cropparams *p)
{
  /* Make sure we actually have inputs. */
  if(p->inputs==NULL)
    error(EXIT_FAILURE, 0, "no input file given");

  /* Make sure that a HDU is also given. */
  if(p->cp.hdu==NULL )
    error(EXIT_FAILURE, 0, "no HDU specified. When the input is a FITS "
          "file, a HDU must also be specified, you can use the '--hdu' "
          "('-h') option and give it the HDU number (starting from "
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
         configured with the '--enable-reentrant' option and we can only
         know from the 'fits_is_reentrant' function that came from CFITSIO
         version 3.30. */
#if GAL_CONFIG_HAVE_FITS_IS_REENTRANT == 1
      if(p->cp.numthreads>1 && fits_is_reentrant()==0)
        {
          fprintf(stderr, "WARNING: CFITSIO was not configured with the "
                  "'--enable-reentrant' option but you have asked to crop "
                  "on %zu threads. Therefore only one thread will be used.\n\n"
                  "Please run the following command to learn more about "
                  "configuring CFITSIO:\n\n"
                  "    $ info gnuastro CFITSIO", p->cp.numthreads);
          p->cp.numthreads=1;
        }
#else
      if(p->cp.numthreads>1)
        {
          fprintf(stderr, "WARNING: the installed CFITSIO version doesn't "
                  "have 'fits_is_reentrant' function (it is older than "
                  "version 3.30). But you have asked to crop on %zu threads."
                  "Therefore only one thread will be used.\n\n"
                  "To avoid this warning, you can set the number of threads "
                  "to one with '-N1' or update your installation of CFITSIO.",
                  p->cp.numthreads);
          p->cp.numthreads=1;
        }
#endif

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
/* When the crop is defined by its center, the final width that we need
   must be in actual number of pixels (an integer). But the user's values
   can be in WCS mode or even in image mode, they may be non-integers. */
static void
ui_set_img_sizes(struct cropparams *p)
{
  gal_data_t *newwidth;
  size_t i, ndim=p->imgs->ndim;
  double pwidth, pcheckcenter, *warray;

  /* Make sure a width value is actually given. */
  if(p->width==NULL)
    error(EXIT_FAILURE, 0, "no crop width specified. When crops are "
          "defined by their center (with '--center' or '--catalog') a "
          "width is necessary (using the '--width' option)");

  /* Make sure that the width array only has one element or the same number
     of elements as the input's dimensions. */
  if(p->width->size!=ndim && p->width->size!=1)
    error(EXIT_FAILURE, 0, "%zu values give to '--width', but input is %zu "
          "dimensions. It can only take either one value (same width in all "
          "dimensions), or the same number as the input's dimensions",
          p->width->size, ndim);

  /* If the width array has only one value, that single value should be
     used for all dimensions. */
  if(p->width->size==1)
    {
      /* Allocate the new width dataset. */
      newwidth=gal_data_alloc(NULL, p->width->type, 1, &ndim, NULL, 0, -1,
                              1, NULL, NULL, NULL);

      /* Fill the new width. */
      warray=newwidth->array;
      for(i=0;i<ndim;++i) warray[i] = *(double *)(p->width->array);

      /* Free the old (single element) width dataset and put the new one in
         its place. */
      gal_data_free(p->width);
      p->width=newwidth;
    }
  else warray=p->width->array;

  /* WCS mode. */
  if(p->mode==IMGCROP_MODE_WCS)
    {
      /* Fill in the widths depending on the mode. */
      for(i=0;i<ndim;++i)
        {
          /* Convert the width in units of the input's WCS into pixels. */
          pwidth = warray[i]/p->pixscale[i];
          if(pwidth<3 || pwidth>50000)
            error(EXIT_FAILURE, 0, "value %g (requested width along "
                  "dimension %zu) translates to %.0f pixels on this "
                  "dataset. This is probably not what you wanted. Note "
                  "that the dataset's resolution in this dimension is "
                  "%g.\n\n"
                  "You can do the conversion to the dataset's WCS units "
                  "prior to calling Crop. Alternatively, you can specify "
                  "all the coordinates/sizes in image (not WCS) units and "
                  "use the '--mode=img' option", warray[i], i+1, pwidth,
                  p->pixscale[i]);

          /* Write the single valued width in WCS and the image width for
             this dimension. */
          p->iwidth[i]=GAL_DIMENSION_FLT_TO_INT(pwidth);
          if(p->iwidth[i]%2==0)
            {
              p->iwidth[i] += 1;
              warray[i]    += p->pixscale[i];
            }
        }

      /* Checkcenter: */
      if(p->incheckcenter)
        pcheckcenter=((double *)(p->incheckcenter->array))[0]/p->pixscale[0];
    }
  /* Image mode. */
  else
    {
      /* The width (along each dimension). */
      for(i=0;i<ndim;++i)
        {
          p->iwidth[i]=GAL_DIMENSION_FLT_TO_INT(warray[i]);
          if(p->iwidth[i]%2==0) p->iwidth[i] += 1;
        }

      /* Checkcenter: */
      if(p->incheckcenter)
        {
          /* Write the double value into the temporary variable */
          pcheckcenter=((double *)(p->incheckcenter->array))[0];

          /* In image-mode it has to be an integer. */
          if( ceilf(pcheckcenter)!=pcheckcenter )
            error(EXIT_FAILURE, 0, "%g is not an integer. When cropping in "
                  "image-mode, the number of pixels to check in the "
                  "center must be an integer", pcheckcenter);
        }
    }

  /* Finalize the number of central pixels to check. */
  if(p->incheckcenter)
    {
      /* Convert the floating point value to an integer. */
      p->checkcenter=GAL_DIMENSION_FLT_TO_INT(pcheckcenter);

      /* If 'checkcenter' isn't zero, but is even, convert it to an odd
         number (so the actual center can be checked). */
      if(p->checkcenter && p->checkcenter%2==0) p->checkcenter += 1;
    }

  /* For a check:
  printf("Width: "); for(i=0;i<ndim;++i) printf("\t%ld\n\t", p->iwidth[i]);
  exit(0);
  */
}





static void
ui_read_cols(struct cropparams *p)
{
  char colname[100];
  gal_data_t *cols, *tmp, *corrtype=NULL;
  gal_list_str_t *colstrs=NULL, *extracolstr, *lastcolstr;
  size_t i, ncoordcols, counter=0, dcounter=0, ndim=p->imgs->ndim;

  /* See if the number of columns given for coordinates corresponds to the
     number of dimensions of the input dataset. */
  if(p->coordcol)
    {
      /* Check if the number of columns given for coordinates is the same
         as the number of dimensions in the input dataset(s). */
      ncoordcols=gal_list_str_number(p->coordcol);
      if( ncoordcols < ndim)
        error(EXIT_FAILURE, 0, "'--coordcol' was called %zu times, but the "
              "input dataset%s %zu dimensions. Recall that through "
              "'--coordcol' you are specifying the columns containing the "
              "coordinates of the center of the crop in a catalog",
              ncoordcols, (p->numin==1?" has":"s have"), ndim);

      /* If the number of given columns is more than the input's
         dimensions, then we'll just delete all the unnecessary columns. */
      else if( ncoordcols > ndim )
        {
          /* Go over the columns find the last, but first initialize the
             two ('lastcolstr' to avoid compiler warnings). */
          lastcolstr=extracolstr=p->coordcol;
          for(i=0;i<ndim;++i)
            {
              /* Keep the last node if on the last (useful) column. */
              if(i==ndim-1) lastcolstr=extracolstr;

              /* Go onto the next one. */
              extracolstr=extracolstr->next;
            }

          /* Set the 'next' element of the last node to NULL and free the
             extra ones. */
          lastcolstr->next=NULL;
          gal_list_str_free(extracolstr, 1);
        }
    }
  else
    error(EXIT_FAILURE, 0, "no coordinate columns specified. When a catalog"
          "is given, it is necessary to identify which columns identify "
          "the coordinate values in which dimension.\n\n"
          "You can do this by calling '--coordcol' multiple times, the "
          "order must be in the same order as the input's dimensions. "
          "For more information on how to select columns in Gnuastro, "
          "please run the following command:\n\n"
          "    $ info gnuastro \"Selecting table columns\"");


  /* If a name column was also given, the put that as the first column to
     read, otherwise just use the given set of columns (in the same
     order). */
  if(p->namecol)
    {
      gal_list_str_add(&colstrs, p->namecol, 0);
      colstrs->next=p->coordcol;
    }
  else colstrs=p->coordcol;


  /* Read the desired columns from the file. */
  cols=gal_table_read(p->catname, p->cathdu, NULL, colstrs, p->cp.searchin,
                      p->cp.ignorecase, p->cp.minmapsize, p->cp.quietmmap,
                      NULL);
  if(cols==NULL)
    error(EXIT_FAILURE, 0, "%s: is empty! No usable information "
          "(un-commented lines) could be read from this file",
          gal_fits_name_save_as_string(p->catname, p->cathdu));


  /* Set the number of objects (rows in each column). */
  p->numout=cols->size;


  /* Make sure more columns were not read (the name matchings might result
     in more than one column being read from the inputs). */
  if( gal_list_data_number(cols) != ndim + (p->namecol!=NULL) )
    gal_tableintern_error_col_selection(p->catname, p->cathdu, "too many "
                                        "columns were selected by the given "
                                        "values to the options ending in "
                                        "'col'.");


  /* Put the information in each column in the proper place. */
  while(cols!=NULL)
    {
      /* Pop out the top node. */
      tmp=gal_list_data_pop(&cols);

      /* See which column it is. */
      switch(++counter)
        {
        case 1:
          if(p->namecol)
            {
              sprintf(colname, "crop name prefix");
              corrtype = (tmp->type==GAL_TYPE_STRING ? tmp
                          : gal_data_copy_to_new_type_free(tmp,
                                                            GAL_TYPE_STRING));
              p->name=corrtype->array;
            }
          else
            {
              sprintf(colname, "position in dimension %zu", dcounter+1);
              corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
              p->centercoords[ dcounter++ ]=corrtype->array;
            }
          break;

        default:
          sprintf(colname, "position in dimension %zu", dcounter+1);
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
          p->centercoords[ dcounter++ ] = corrtype->array;
        }

      /* Sanity check and clean up.  Note that it might happen that the
         input structure is already freed. In that case, 'corrtype' will be
         NULL. */
      if(corrtype)
        {
          /* Make sure there are no blank values in this column. */
          if( gal_blank_present(corrtype, 1) )
            error(EXIT_FAILURE, 0, "%s: column with %s has blank values. "
                  "Input columns must not contain blank values",
                  gal_fits_name_save_as_string(p->catname, p->cathdu),
                  colname);

          /* Free the unnecessary sturcture information. The correct-type
             ('corrtype') data structure's array is necessary for later
             steps, so its pointer has been copied in the main program's
             structure. Hence, we should set the structure's pointer to
             NULL so the important data isn't freed.*/
          corrtype->array=NULL;
          gal_data_free(corrtype);
          corrtype=NULL;
        }
    }
}





static void
ui_prepare_center(struct cropparams *p)
{
  double *carray;
  size_t i, ndim=p->imgs->ndim;

  /* Allocate space to keep the central positions. */
  errno=0;
  p->centercoords=malloc(ndim*sizeof *p->centercoords);
  if( p->centercoords==NULL )
    error(EXIT_FAILURE, 0, "%s: %zu bytes for 'p->centercoords'",
          __func__, ndim*sizeof *p->centercoords);


  /* Set the integer widths of the crop(s) when defined by center. */
  ui_set_img_sizes(p);

  /* For a catalog, we have a separate function, but for a single center
     value, put the center values into an array. This will essentially
     simulate a catalog with one row. So from this point on, there is no
     difference between a catalog input and a central position input. */
  if(p->catname)
    ui_read_cols(p);
  else
    {
      carray=p->center->array;
      for(i=0;i<ndim;++i)
        {
          p->centercoords[i]=gal_pointer_allocate(GAL_TYPE_FLOAT64,
                                                  1, 0, __func__,
                                                  "p->centercoords[i]");
          p->centercoords[i][0]=carray[i];
        }
    }
}





/* Add all the columns of the log file. Just note that since this is a
   linked list, we have to add them in the opposite order. */
static void
ui_make_log(struct cropparams *p)
{
  char *comment;

  /* Return if no long file was requested. */
  if(p->cp.log==0) return;

  /* Column to specify if the central pixels are filled. */
  if( asprintf(&comment, "Are the central pixels filled? (1: yes, 0: no, "
               "%u: not checked)", GAL_BLANK_UINT8)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT8, 1, &p->numout,
                          NULL, 1, p->cp.minmapsize, p->cp.quietmmap,
                          "CENTER_FILLED", "bool", comment);
  free(comment);

  /* Column for number of datasets used in this crop. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT16, 1, &p->numout,
                          NULL, 1, p->cp.minmapsize, p->cp.quietmmap,
                          "NUM_INPUTS", "count",
                          "Number of input datasets used to make this crop.");

  /* Filename of crop. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_STRING, 1, &p->numout,
                          NULL, 1, p->cp.minmapsize, p->cp.quietmmap,
                          "CROP_NAME", "name", "File name of crop.");
}





void
ui_preparations(struct cropparams *p)
{
  fitsfile *tmpfits;
  struct inputimgs *img;
  int status, firsttype=0;
  size_t input_counter, firstndim=0;


  /* For polygon and section, there should be no center checking. */
  if(p->polygon || p->section)
    p->checkcenter=0;


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
    error(EXIT_FAILURE, errno, "%s: %zu bytes for p->imgs",
          __func__, p->numin*sizeof *p->imgs);


  /* Fill in the WCS information of each image. */
  input_counter=p->numin;
  while(p->inputs)
    {
      /* Pop from the list of input images and get the info. */
      status=0;
      img=&p->imgs[--input_counter];
      img->name=gal_list_str_pop(&p->inputs);
      tmpfits=gal_fits_hdu_open_format(img->name, p->cp.hdu, 0);
      gal_fits_img_info(tmpfits, &p->type, &img->ndim, &img->dsize,
                        NULL, NULL);
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
                "image is not recognized. So WCS mode cannot be used "
                "as input coordinates. You can try with pixel coordinates "
                "with '--mode=img'", img->name, p->cp.hdu);
      fits_close_file(tmpfits, &status);
      gal_fits_io_error(status, NULL);

      /* Make sure all the images have the same type and dimensions. */
      if(firsttype==0)
        {
          /* Set the basic information. */
          firsttype = p->type;
          firstndim = img->ndim;
          p->blankptrread  = gal_blank_alloc_write(p->type);
          p->blankptrwrite = gal_fits_key_img_blank(p->type);

          /* Make sure the number of dimensions is supported. */
          if(firstndim>MAXDIM)
            error(EXIT_FAILURE, 0, "%s: is as %zu dimensional dataset, Crop "
                  "currently only supports a maximum of %d dimensions",
                  img->name, firstndim, MAXDIM);

          /* Make sure the number of coordinates given for center
             correspond to the dimensionality of the data. */
          if(p->center && p->center->size!=firstndim)
            error(EXIT_FAILURE, 0, "%s (hdu %s) has %zu dimensions, but "
                  "%zu coordinates were given to '--center'", img->name,
                  p->cp.hdu, firstndim, p->center->size);
        }
      else
        {
          if(firsttype!=p->type)
            error(EXIT_FAILURE, 0, "%s: type is '%s' while previous input(s) "
                  "were '%s' type. All inputs must have the same pixel data "
                  "type.\n\nYou can use Gnuastro's Arithmetic program to "
                  "convert '%s' to '%s', please run this command for more "
                  "information (press 'SPACE' for going down and 'q' to "
                  "return to the command-line):\n\n"
                  "    $ info Arithmetic\n",
                  img->name, gal_type_name(p->type, 1),
                  gal_type_name(firsttype, 1), img->name,
                  gal_type_name(p->type, 1));
          if(firstndim!=img->ndim)
            error(EXIT_FAILURE, 0, "%s: type has %zu dimensions, while "
                  "previous input(s) had %zu dimensions. All inputs must "
                  "have the same number of dimensions", img->name, img->ndim,
                  firstndim);
        }

      /* In WCS mode, we need some additional preparations. */
      if(p->mode==IMGCROP_MODE_WCS) wcsmode_check_prepare(p, img);
    }


  /* Polygon cropping is currently only supported on 2D */
  if(p->imgs->ndim!=2 && p->polygon)
    error(EXIT_FAILURE, 0, "%s: polygon cropping is currently only "
          "supported on 2D datasets (images), not %zuD datasets",
          p->imgs->name, p->imgs->ndim);


  /* Unify central crop methods into 'p->centercoords'. */
  if(p->catname || p->center)
    ui_prepare_center(p);


  /* 'ui_read_cols' set the number of output crops when a catalog was
     given, in all other cases, we only have one output. */
  if(p->catname==NULL) p->numout=1;


  /* Prepare the log file if the user has asked for it. */
  ui_make_log(p);
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct cropparams *p)
{
  char *msg;
  struct timeval t1;
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


  /* To see how long it takes to read meta-data. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);


  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);


  /* Report timing: */
  if(!p->cp.quiet)
    {
      printf(PROGRAM_NAME" "PACKAGE_VERSION" started on %s",
             ctime(&p->rawtime));
      if(p->cp.numthreads>1)
        printf("  - Using %zu CPU thread%s\n", p->cp.numthreads,
               p->cp.numthreads==1 ? "." : "s.");
      if(p->checkcenter)
        printf("  - Number of central pixels to check for blank: %zu\n",
               p->checkcenter);
      if( asprintf(&msg, "Read metadata of %zu dataset%s.", p->numin,
                   p->numin>1 ? "s" : "")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(&t1, msg, 1);
      if(p->numout>1)
        {
          if( asprintf(&msg, "Will try making %zu crops (from catalog).",
                       p->numout)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(NULL, msg, 1);
        }
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct cropparams *p, struct timeval *t1)
{
  size_t i;

  /* Free the simple arrays (if they were set). */
  free(p->blankptrread);
  free(p->blankptrwrite);
  gal_data_free(p->center);
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
