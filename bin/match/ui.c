/*********************************************************************
Match - A program to match catalogs and WCS warps
Match is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2018, Free Software Foundation, Inc.

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

#include <gnuastro/fits.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"

#include "ui.h"
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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" matches catalogs of objects "
  "and (by default) will return the re-arranged matching inputs. The "
  "optional log file will return low-level information about the match "
  "(indexs and distances).\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct matchparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;


  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_HDU:
          cp->coptions[i].doc="Extension name or number of first input.";
          break;
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_NUMTHREADS:
          cp->coptions[i].flags=OPTION_HIDDEN;
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
  struct matchparams *p = state->input;

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
      if(p->input1name)
        {
          if(p->input2name)
            argp_error(state, "only two arguments (input files) should be "
                       "given");
          else
            p->input2name=arg;
        }
      else
        p->input1name=arg;
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
ui_read_check_only_options(struct matchparams *p)
{
  if(p->outcols && p->notmatched)
    error(EXIT_FAILURE, 0, "`--outcols' and `--notmatched' cannot be called "
          "at the same time. The former is only for cases when the matches "
          "are required");
}






static void
ui_check_options_and_arguments(struct matchparams *p)
{
  /* Make sure two input file names were given and if they a FITS file,
     that a HDU is also given for each. */
  if(p->input1name)
    {
      if( gal_fits_name_is_fits(p->input1name) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU for first input. When the input is "
              "a FITS file, a HDU must also be specified, you can use the "
              "`--hdu' (`-h') option and give it the HDU number (starting "
              "from zero), extension name, or anything acceptable by "
              "CFITSIO");
    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified: two inputs are "
          "necessary");

  if(p->input2name)
    {
      if( gal_fits_name_is_fits(p->input2name) && p->hdu2==NULL )
        error(EXIT_FAILURE, 0, "no HDU for second input. Please use the "
              "`--hdu2' (`-H') option and give it the HDU number (starting "
              "from zero), extension name, or anything acceptable by "
              "CFITSIO");
    }
  else
    error(EXIT_FAILURE, 0, "second input file not specified: two inputs are "
          "necessary");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_set_mode(struct matchparams *p)
{
  /* Check if we are in image or catalog mode. We will base the mode on the
     first input, then check with the second. */
  if( gal_fits_name_is_fits(p->input1name) )
    p->mode = ( (gal_fits_hdu_format(p->input1name, p->cp.hdu) == IMAGE_HDU)
                ? MATCH_MODE_WCS
                : MATCH_MODE_CATALOG );
  else
    p->mode=MATCH_MODE_CATALOG;

  /* Now that the mode is set, check the second input's type. */
  if( gal_fits_name_is_fits(p->input2name) )
    {
      if(gal_fits_hdu_format(p->input2name, p->hdu2) == IMAGE_HDU)
        {
          if( p->mode==MATCH_MODE_CATALOG)
            error(EXIT_FAILURE, 0, "%s is a catalog, while %s is an image. "
                  "Both inputs have to be images or catalogs",
                  gal_checkset_dataset_name(p->input1name, p->cp.hdu),
                  gal_checkset_dataset_name(p->input2name, p->hdu2) );
        }
      else
        {
          if( p->mode==MATCH_MODE_WCS)
            error(EXIT_FAILURE, 0, "%s is an image, while %s is a catalog. "
                  "Both inputs have to be images or catalogs",
                  gal_checkset_dataset_name(p->input1name, p->cp.hdu),
                  gal_checkset_dataset_name(p->input2name, p->hdu2));
        }
    }
  else
    if(p->mode==MATCH_MODE_WCS)
      error(EXIT_FAILURE, 0, "%s is an image, while %s is a catalog! Both "
            "inputs have to be images or catalogs",
            gal_checkset_dataset_name(p->input1name, p->cp.hdu),
            gal_checkset_dataset_name(p->input2name, p->hdu2));
}





/* The final aperture must have the following values:

       p->aperture[0]: Major axis length.
       p->aperture[1]: Axis ratio.
       p->aperture[2]: Position angle (relative to first dim).     */
static void
ui_read_columns_aperture_2d(struct matchparams *p)
{
  size_t apersize=3;
  gal_data_t *newaper=NULL;
  double *naper, *oaper=p->aperture->array;

  /* A general sanity check: the first two elements of aperture cannot be
     zero or negative. */
  if( oaper[0]<=0 )
    error(EXIT_FAILURE, 0, "the first value of `--aperture' cannot be "
          "zero or negative");
  if( p->aperture->size>1 && oaper[1]<=0 )
    error(EXIT_FAILURE, 0, "the second value of `--aperture' cannot be "
          "zero or negative");

  /* Will be needed in more than one case. */
  if(p->aperture->size!=3)
    {
      newaper=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &apersize, NULL,
                             0, -1, NULL, NULL, NULL);
      naper=newaper->array;
    }

  /* Different based on  */
  switch(p->aperture->size)
    {
    case 1:
      naper[0]=oaper[0];
      naper[1]=1;
      naper[2]=0;
      break;

    case 2:
      naper[0] = oaper[0]>oaper[1] ? oaper[0]          : oaper[1];
      naper[1] = oaper[0]>oaper[1] ? oaper[1]/oaper[0] : oaper[0]/oaper[1];
      naper[2] = oaper[0]>oaper[1] ? 0                 : 90;
      break;

    case 3:
      if(oaper[1]>1)
        error(EXIT_FAILURE, 0, "second value to `--aperture' is larger "
              "than one. When three numbers are given to this option, the "
              "second is the axis ratio (which must always be less than 1).");
      break;

    default:
      error(EXIT_FAILURE, 0, "%zu values given to `--aperture'. In 2D, this "
            "option can only take 1, 2, or 3 values", p->aperture->size);
    }

  /* If a new aperture was defined, then replace it with the exitsting
     one. */
  if(newaper)
    {
      gal_data_free(p->aperture);
      p->aperture=newaper;
    }
}





/* We want to keep the columns as double type. So what-ever their original
   type is, convert it. */
static gal_data_t *
ui_read_columns_to_double(struct matchparams *p, char *filename, char *hdu,
                          gal_list_str_t *cols, size_t numcols)
{
  gal_data_t *tmp, *ttmp, *tout, *out=NULL;
  struct gal_options_common_params *cp=&p->cp;
  char *diff_cols_error="%s: the number of columns matched (%zu) "
    "differs from the number of usable calls to `--ccol1' (%zu). "
    "Please give more specific values to `--ccol1' (column "
    "numberes are the only identifiers guaranteed to be unique).";

  /* Read the columns. */
  tout=gal_table_read(filename, hdu, cols, cp->searchin, cp->ignorecase,
                      cp->minmapsize, NULL);

  /* A small sanity check. */
  if(gal_list_data_number(tout)!=numcols)
    error(EXIT_FAILURE, 0, diff_cols_error,
          gal_checkset_dataset_name(filename, hdu),
          gal_list_data_number(tout), numcols);

  /* Go over the columns and see if they are double or not. To keep things
     simple, we'll keep a new list even if all the types are float64.*/
  tmp=tout;
  while(tmp!=NULL)
    {
      /* We need ot set the `next' pointer  */
      ttmp=tmp->next;
      tmp->next=NULL;

      /* Correct the type if necessary. */
      if(tmp->type==GAL_TYPE_FLOAT64)
        gal_list_data_add(&out, tmp);
      else
        gal_list_data_add(&out,
                          gal_data_copy_to_new_type_free(tmp,
                                                         GAL_TYPE_FLOAT64) );

      /* Set `tmp' to the initial `next pointer. */
      tmp=ttmp;
    }

  /* The `out' above is in reverse, so correct it and return */
  gal_list_data_reverse(&out);
  return out;
}





/* Read catalog columns */
static void
ui_read_columns(struct matchparams *p)
{
  size_t i;
  size_t ccol1n=p->ccol1->size;
  size_t ccol2n=p->ccol2->size;
  gal_list_str_t *cols1=NULL, *cols2=NULL;
  char **strarr1=p->ccol1->array, **strarr2=p->ccol2->array;

  /* Make sure the same number of columns is given to both. */
  if(ccol1n!=ccol2n)
    error(EXIT_FAILURE, 0, "the number of values given to `--ccol1' and "
          "`--ccol2' (%zu and %zu) are not equal", ccol1n, ccol2n);


  /* Read/check the aperture values. */
  if(p->aperture)
    switch(ccol1n)
      {
      case 1:
        if(p->aperture->size>1)
          error(EXIT_FAILURE, 0, "%zu values given to `--aperture'. In a 1D "
                "match, this option can only take one value",
                p->aperture->size);
        break;

      case 2:
        ui_read_columns_aperture_2d(p);
        break;

      default:

        error(EXIT_FAILURE, 0, "%zu dimensional matches are not currently "
              "supported (maximum is 2 dimensions). The number of "
              "dimensions is deduced from the number of values given to "
              "`--ccol1' and `--ccol2'", ccol1n);
      }
  else
    error(EXIT_FAILURE, 0, "no matching aperture specified. Please use "
          "the `--aperture' option to define the acceptable aperture for "
          "matching the coordinates (in the same units as each "
          "dimension). Please run the following command for more "
          "information.\n\n    $ info %s\n", PROGRAM_EXEC);


  /* Convert the array of strings to a list of strings for the column
     names. */
  for(i=0;i<ccol1n;++i)
    {
      gal_list_str_add(&cols1, strarr1[i], 0);
      gal_list_str_add(&cols2, strarr2[i], 0);
      strarr1[i]=strarr2[i]=NULL;  /* So they are not freed later. */
    }
  gal_list_str_reverse(&cols1);
  gal_list_str_reverse(&cols2);


  /* Read the columns. */
  if(p->cp.searchin)
    {
      /* Read the first dataset. */
      p->cols1=ui_read_columns_to_double(p, p->input1name, p->cp.hdu,
                                         cols1, ccol1n);
      p->cols2=ui_read_columns_to_double(p, p->input2name, p->hdu2,
                                         cols2, ccol2n);
    }
  else
    error(EXIT_FAILURE, 0, "no `--searchin' option specified. Please run "
          "the following command for more information:\n\n"
          "    $ info gnuastro \"selecting table columns\"\n");

  /* Free the extra spaces. */
  gal_list_str_free(cols1, 1);
  gal_list_str_free(cols2, 1);
  gal_data_free(p->ccol1);
  gal_data_free(p->ccol2);
  p->ccol1=p->ccol2=NULL;
}





static void
ui_preparations_out_cols(struct matchparams *p)
{
  size_t i;
  char **strarr=p->outcols->array;

  /* Go over all the values and put the respective column identifier in the
     proper list. */
  for(i=0;i<p->outcols->size;++i)
    switch(strarr[i][0])
      {
      case 'a': gal_list_str_add(&p->acols, strarr[i]+1, 0); break;
      case 'b': gal_list_str_add(&p->bcols, strarr[i]+1, 0); break;
      default:
        error(EXIT_FAILURE, 0, "`%s' is not a valid value for `--outcols'. "
              "The first character of each value to this option must be "
              "either `a' or `b'. The former specifies a column from the "
              "first input and the latter a column from the second. The "
              "characters after them can be any column identifier (number, "
              "name, or regular expression). For more on column selection, "
              "please run this command:\n\n"
              "    $ info gnuastro \"Selecting table columns\"\n",
              strarr[i]);
      }

  /* Revere the lists so they correspond to the input order. */
  gal_list_str_reverse(&p->acols);
  gal_list_str_reverse(&p->bcols);
}





static void
ui_preparations_out_name(struct matchparams *p)
{
  if(p->logasoutput)
    {
      /* Set the logname (as output). */
      if(p->cp.output)
        gal_checkset_allocate_copy(p->cp.output, &p->logname);
      else
        {
          if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
            p->logname=gal_checkset_automatic_output(&p->cp, p->input1name,
                                                     "_matched.txt");
          else
            p->logname=gal_checkset_automatic_output(&p->cp, p->input1name,
                                                     "_matched.fits");
        }

      /* Make sure a file with this name doesn't exist. */
      gal_checkset_writable_remove(p->out1name, 0, p->cp.dontdelete);
    }
  else
    {
      if(p->outcols)
        {
          if(p->cp.output==NULL)
            p->cp.output = gal_checkset_automatic_output(&p->cp,
                 p->input1name, ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                                  ? "_matched.txt" : "_matched.fits") );
          gal_checkset_writable_remove(p->cp.output, 0, p->cp.dontdelete);
        }
      else
        {
          /* Set `p->out1name' and `p->out2name'. */
          if(p->cp.output)
            {
              if( gal_fits_name_is_fits(p->cp.output) )
                {
                  gal_checkset_allocate_copy(p->cp.output, &p->out1name);
                  gal_checkset_allocate_copy(p->cp.output, &p->out2name);
                }
              else
                {
                  p->out1name=gal_checkset_automatic_output(&p->cp,
                                                            p->cp.output,
                                                            "_matched_1.txt");
                  p->out2name=gal_checkset_automatic_output(&p->cp,
                                                            p->cp.output,
                                                            "_matched_2.txt");
                }
            }
          else
            {
              if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
                {
                  p->out1name=gal_checkset_automatic_output(&p->cp,
                                                            p->input1name,
                                                            "_matched_1.txt");
                  p->out2name=gal_checkset_automatic_output(&p->cp,
                                                            p->input2name,
                                                            "_matched_2.txt");
                }
              else
                {
                  p->out1name=gal_checkset_automatic_output(&p->cp,
                                                            p->input1name,
                                                            "_matched.fits");
                  gal_checkset_allocate_copy(p->out1name, &p->out2name);
                }
            }

          /* Make sure no file with these names exists. */
          gal_checkset_writable_remove(p->out1name, 0, p->cp.dontdelete);
          gal_checkset_writable_remove(p->out2name, 0, p->cp.dontdelete);
        }

      /* If a log file is necessary, set its name here. */
      if(p->cp.log)
        {
          p->logname = ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                         ? PROGRAM_EXEC".txt"
                         : PROGRAM_EXEC".fits" );
          gal_checkset_writable_remove(p->logname, 0, p->cp.dontdelete);
        }
    }
}





static void
ui_preparations(struct matchparams *p)
{
  /* Set the mode of the program. */
  ui_set_mode(p);

  /* Currently Match only works on catalogs. */
  if(p->mode==MATCH_MODE_WCS)
    error(EXIT_FAILURE, 0, "currently Match only works on catalogs, we will "
          "implement the WCS matching routines later");
  else
    {
      ui_read_columns(p);
      if(p->outcols) ui_preparations_out_cols(p);
    }

  /* Set the output filename. */
  ui_preparations_out_name(p);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct matchparams *p)
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
ui_free_report(struct matchparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->out1name);
  free(p->out2name);
  free(p->cp.output);

  /* Print the final message.
  if(!p->cp.quiet)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
  */
}
