/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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

#include <gnuastro/txt.h>
#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/tile.h>
#include <gnuastro/array.h>
#include <gnuastro/qsort.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/threads.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>
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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will do statistical "
  "analysis on the input dataset (table column or image). All blank "
  "pixels or pixels outside of the given range are ignored. You can "
  "either directly ask for certain statistics in one line/row as shown "
  "below with the same order as requested, or get tables of different "
  "statistical measures like the histogram, cumulative frequency style "
  "and etc. If no particular statistic is requested, some basic "
  "information about the dataset is printed on the command-line.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;






















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct statisticsparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;

  /* Set the necessary common parameters structure. */
  cp->program_struct     = p;
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;
  cp->numthreads         = gal_threads_number();
  cp->tl.remainderfrac   = NAN;

  /* Program-specific initializers */
  p->lessthan            = NAN;
  p->onebinstart         = NAN;
  p->greaterequal        = NAN;
  p->quantmin            = NAN;
  p->quantmax            = NAN;
  p->mirror              = NAN;
  p->mirrordist          = NAN;
  p->meanmedqdiff        = NAN;
  p->sclipparams[0]      = NAN;
  p->sclipparams[1]      = NAN;

  /* Set the mandatory common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    switch(cp->coptions[i].key)
      {
      case GAL_OPTIONS_KEY_LOG:
      case GAL_OPTIONS_KEY_TYPE:
        cp->coptions[i].flags=OPTION_HIDDEN;
        break;

      case GAL_OPTIONS_KEY_SEARCHIN:
      case GAL_OPTIONS_KEY_MINMAPSIZE:
      case GAL_OPTIONS_KEY_TABLEFORMAT:
        cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
        break;
      }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct statisticsparams *p = state->input;

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





static void *
ui_add_to_single_value(struct argp_option *option, char *arg,
                      char *filename, size_t lineno, void *params)
{
  size_t i;
  double *d;
  gal_data_t *inputs=NULL;
  struct statisticsparams *p=(struct statisticsparams *)params;

  /* In case of printing the option values. */
  if(lineno==-1)
    error(EXIT_FAILURE, 0, "currently the options to be printed in one row "
          "(like '--number', '--mean', and etc) do not support printing "
          "with the '--printparams' ('-P'), or writing into configuration "
          "files due to lack of time when implementing these features. "
          "You can put them into configuration files manually. Please get "
          "in touch with us at '%s', so we can implement it",
          PACKAGE_BUGREPORT);

  /* Some of these options take values and some don't. */
  if(option->type==GAL_OPTIONS_NO_ARG_TYPE)
    {
      /* If this option is given in a configuration file, then 'arg' will not
         be NULL and we don't want to do anything if it is '0'. */
      if(arg)
        {
          /* Make sure the value is only '0' or '1'. */
          if( arg[1]!='\0' && *arg!='0' && *arg!='1' )
            error_at_line(EXIT_FAILURE, 0, filename, lineno, "the '--%s' "
                          "option takes no arguments. In a configuration "
                          "file it can only have the values '1' or '0', "
                          "indicating if it should be used or not",
                          option->name);

          /* Only proceed if the (possibly given) argument is 1. */
          if(arg[0]=='0' && arg[1]=='\0') return NULL;
        }

      /* Add this option to the print list. */
      gal_list_i32_add(&p->singlevalue, option->key);
    }
  else
    {
      /* Read the string of numbers. */
      inputs=gal_options_parse_list_of_numbers(arg, filename, lineno);
      if(inputs->size==0)
        error(EXIT_FAILURE, 0, "'--%s' needs a value", option->name);

      /* Do the appropriate operations with the  */
      d=inputs->array;
      switch(option->key)
        {
        case UI_KEY_QUANTILE:
        case UI_KEY_QUANTFUNC:
          /* For the quantile and the quantile function, its possible to
             give any number of arguments, so add the operation index and
             the argument once for each given number. */
          for(i=0;i<inputs->size;++i)
            {
              if(option->key==UI_KEY_QUANTILE && (d[i]<0 || d[i]>1) )
                error_at_line(EXIT_FAILURE, 0, filename, lineno, "values "
                              "to '--quantile' ('-u') must be between 0 "
                              "and 1, you had asked for %g (read from '%s')",
                              d[i], arg);
              gal_list_f64_add(&p->tp_args, d[i]);
              gal_list_i32_add(&p->singlevalue, option->key);
            }
          break;

        default:
          error_at_line(EXIT_FAILURE, 0, filename, lineno, "a bug! please "
                        "contact us at %s so we can address the problem. "
                        "the option given to 'ui_add_to_print_in_row' is "
                        "marked as requiring a value, but is not recognized",
                        PACKAGE_BUGREPORT);
        }
    }

  return NULL;
}





static void *
ui_read_quantile_range(struct argp_option *option, char *arg,
                       char *filename, size_t lineno, void *params)
{
  char *str;
  gal_data_t *in;
  struct statisticsparams *p=(struct statisticsparams *)params;

  /* For the '--printparams' ('-P') option:*/
  if(lineno==-1)
    {
      if( isnan(p->quantmax) )
        {
          if( asprintf(&str, "%g", p->quantmin)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        {
          if( asprintf(&str, "%g,%g", p->quantmin, p->quantmax)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      return str;
    }

  /* Parse the inputs. */
  in=gal_options_parse_list_of_numbers(arg, filename, lineno);

  /* Check if there was only two numbers. */
  if(in->size!=1 && in->size!=2)
    error_at_line(EXIT_FAILURE, 0, filename, lineno, "the '--%s' "
                  "option takes one or two values values (separated by "
                  "a comma) to define the range of used values with "
                  "quantiles. However, %zu numbers were read in the "
                  "string '%s' (value to this option).\n\n"
                  "If there is only one number as input, it will be "
                  "interpretted as the lower quantile (Q) range. The "
                  "higher range will be set to the quantile (1-Q). "
                  "When two numbers are given, they will be used as the "
                  "lower and higher quantile range respectively",
                  option->name, in->size, arg);

  /* Read the values in. */
  p->quantmin = ((double *)(in->array))[0];
  if(in->size==2) p->quantmax = ((double *)(in->array))[1];

  /* Make sure the values are between 0 and 1. */
  if( (p->quantmin<0 || p->quantmin>1)
      || ( !isnan(p->quantmax) && (p->quantmax<0 || p->quantmax>1) ) )
    error_at_line(EXIT_FAILURE, 0, filename, lineno, "values to the "
                  "'--quantrange' option must be between 0 and 1 "
                  "(inclusive). Your input was: '%s'", arg);

  /* When only one value is given, make sure it is less than 0.5. */
  if( !isnan(p->quantmax) && p->quantmin>0.5 )
    error(EXIT_FAILURE, 0, "%g>=0.5! When only one value is given to the "
          "'--%s' option, the range is defined as Q and 1-Q. Thus, the "
          "value must be less than 0.5", p->quantmin, option->name);

  /* Clean up and return. */
  gal_data_free(in);
  return NULL;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct statisticsparams *p)
{
  gal_list_i32_t *tmp;
  struct gal_tile_two_layer_params *tl=&p->cp.tl;

  /* Check if the format of the output table is valid, given the type of
     the output. */
  gal_tableintern_check_fits_format(p->cp.output, p->cp.tableformat);


  /* If in tile-mode, we must have at least one single valued option. */
  if(p->ontile && p->singlevalue==NULL)
    error(EXIT_FAILURE, 0, "at least one of the single-value measurements "
          "(for example '--median') must be requested with the '--ontile' "
          "option: there is no value to put in each tile");

  /* Tessellation related options. */
  if( p->ontile || p->sky )
    {
      /* The tile or sky mode cannot be called with any other modes. */
      if(p->asciihist || p->asciicfp || p->histogram || p->cumulative
         || p->sigmaclip || !isnan(p->mirror) )
        error(EXIT_FAILURE, 0, "'--ontile' or '--sky' cannot be called with "
              "any of the 'particular' calculation options, for example "
              "'--histogram'. This is because the latter work over the whole "
              "dataset and element positions are changed, but in the former "
              "positions are significant");

      /* Make sure the tessellation defining options are given. */
      if( tl->tilesize==NULL || tl->numchannels==NULL
          || isnan(tl->remainderfrac) )
         error(EXIT_FAILURE, 0, "'--tilesize', '--numchannels', and "
               "'--remainderfrac' are mandatory options when dealing with "
               "a tessellation (in '--ontile' or '--sky' mode). Atleast "
               "one of these options wasn't given a value.");
    }


  /* In Sky mode, several options are mandatory. */
  if( p->sky )
    {
      /* Mandatory options. */
      if( isnan(p->meanmedqdiff) || isnan(p->sclipparams[0])
          || p->cp.interpmetric==0 || p->cp.interpnumngb==0 )
        error(EXIT_FAILURE, 0, "'--meanmedqdiff', '--sclipparams', "
              "'--interpmetric' and '--interpnumngb' are mandatory when "
              "requesting Sky measurement ('--sky')");

      /* If mode and median distance is a reasonable value. */
      if(p->meanmedqdiff>0.5)
        error(EXIT_FAILURE, 0, "%f not acceptable for '--meanmedqdiff'. It "
              "cannot take values larger than 0.5 (quantile of median)",
              p->meanmedqdiff);

      /* If a kernel name has been given, we need the HDU. */
      if(p->kernelname && gal_fits_name_is_fits(p->kernelname)
         && p->khdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified for the kernel image. When "
              "A HDU is necessary for FITS files. You can use the '--khdu' "
              "('-u') option and give it the HDU number (starting from "
              "zero), extension name, or anything acceptable by CFITSIO");
    }


  /* Sigma-clipping needs 'sclipparams'. */
  if(p->sigmaclip && isnan(p->sclipparams[0]))
    error(EXIT_FAILURE, 0, "'--sclipparams' is necessary with '--sigmaclip'. "
          "'--sclipparams' takes two values (separated by a comma) for "
          "defining the sigma-clip: the multiple of sigma, and tolerance "
          "(<1) or number of clips (>1).");


  /* If any of the mode measurements are requested, then 'mirrordist' is
     mandatory. */
  for(tmp=p->singlevalue; tmp!=NULL; tmp=tmp->next)
    switch(tmp->v)
      {
      case UI_KEY_MODE:
      case UI_KEY_MODESYM:
      case UI_KEY_MODEQUANT:
      case UI_KEY_MODESYMVALUE:
        if( isnan(p->mirrordist) )
          error(EXIT_FAILURE, 0, "'--mirrordist' is required for the "
                "mode-related single measurements ('--mode', '--modequant', "
                "'--modesym', and '--modesymvalue')");
        break;
      case UI_KEY_SIGCLIPSTD:
      case UI_KEY_SIGCLIPMEAN:
      case UI_KEY_SIGCLIPNUMBER:
      case UI_KEY_SIGCLIPMEDIAN:
        if( isnan(p->sclipparams[0]) )
          error(EXIT_FAILURE, 0, "'--sclipparams' is necessary with "
                "sigma-clipping measurements.\n\n"
                "'--sclipparams' takes two values (separated by a comma) for "
                "defining the sigma-clip: the multiple of sigma, and tolerance "
                "(<1) or number of clips (>1).");
        break;
      }


  /* If less than and greater than are both given, make sure that the value
     to greater than is smaller than the value to less-than. */
  if( !isnan(p->lessthan) && !isnan(p->greaterequal)
      && p->lessthan < p->greaterequal )
    error(EXIT_FAILURE, 0, "the value to '--lessthan' (%g) must be larger "
          "than the value to '--greaterequal' (%g)", p->lessthan,
          p->greaterequal);


  /* Less-than and greater-equal cannot be called together with
     quantrange. */
  if( ( !isnan(p->lessthan) || !isnan(p->greaterequal) )
      && !isnan(p->quantmin) )
    error(EXIT_FAILURE, 0, "'--lessthan' and/or '--greaterequal' cannot "
          "be called together with '--quantrange'");


  /* When binned outputs are requested, make sure that 'numbins' is set. */
  if( (p->histogram || p->cumulative || !isnan(p->mirror)) && p->numbins==0)
    error(EXIT_FAILURE, 0, "'--numbins' isn't set. When the histogram or "
          "cumulative frequency plots are requested, the number of bins "
          "('--numbins') is necessary");


  /* If an ascii plot is requested, check if the ascii number of bins and
     height are given. */
  if( (p->asciihist || p->asciicfp)
      && (p->numasciibins==0 || p->asciiheight==0) )
    error(EXIT_FAILURE, 0, "when an ascii plot is requested, "
          "'--numasciibins' and '--asciiheight' are mandatory, but atleast "
          "one of these has not been given");


  /* Reverse the list of statistics to print in one row and also the
     arguments, so it has the same order the user wanted. */
  gal_list_f64_reverse(&p->tp_args);
  gal_list_i32_reverse(&p->singlevalue);
}





static void
ui_check_options_and_arguments(struct statisticsparams *p)
{
  if(p->inputname)
    {
      /* If input is FITS. */
      if( (p->isfits=gal_fits_name_is_fits(p->inputname)) )
        {
          /* Check if a HDU is given. */
          if( p->cp.hdu==NULL )
            error(EXIT_FAILURE, 0, "no HDU specified. When the input is a "
                  "FITS file, a HDU must also be specified, you can use "
                  "the '--hdu' ('-h') option and give it the HDU number "
                  "(starting from zero), extension name, or anything "
                  "acceptable by CFITSIO");

          /* If its an image, make sure column isn't given (in case the
             user confuses an image with a table). */
          p->hdu_type=gal_fits_hdu_format(p->inputname, p->cp.hdu);
          if(p->hdu_type==IMAGE_HDU && p->column)
            error(EXIT_FAILURE, 0, "%s (hdu: %s): is a FITS image "
                  "extension. The '--column' option is only applicable "
                  "to tables.", p->inputname, p->cp.hdu);
        }
    }
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_out_of_range_to_blank(struct statisticsparams *p)
{
  size_t one=1;
  unsigned char flags=GAL_ARITHMETIC_NUMOK;
  unsigned char flagsor = ( GAL_ARITHMETIC_FREE
                            | GAL_ARITHMETIC_INPLACE
                            | GAL_ARITHMETIC_NUMOK );
  gal_data_t *tmp, *cond_g=NULL, *cond_l=NULL, *cond, *blank, *ref;


  /* Set the dataset that should be used for the condition. */
  ref = p->reference ? p->reference : p->input;


  /* If the user has given a quantile range, then set the 'greaterequal'
     and 'lessthan' values. */
  if( !isnan(p->quantmin) )
    {
      /* If only one value was given, set the maximum quantile range. */
      if( isnan(p->quantmax) ) p->quantmax = 1 - p->quantmin;

      /* Set the greater-equal value. */
      tmp=gal_statistics_quantile(ref, p->quantmin, 1);
      tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
      p->greaterequal=*((float *)(tmp->array));

      /* Set the lower-than value. */
      tmp=gal_statistics_quantile(ref, p->quantmax, 1);
      tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
      p->lessthan=*((float *)(tmp->array));
    }


  /* Set the condition. Note that the 'greaterequal' name is for the data
     we want. So we will set the condition based on those that are
     less-than  */
  if(!isnan(p->greaterequal))
    {
      tmp=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0, -1, 1,
                         NULL, NULL, NULL);
      *((float *)(tmp->array)) = p->greaterequal;
      cond_g=gal_arithmetic(GAL_ARITHMETIC_OP_LT, 1, flags, ref, tmp);
      gal_data_free(tmp);
    }


  /* Same reasoning as above for 'p->greaterthan'. */
  if(!isnan(p->lessthan))
    {
      tmp=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0, -1, 1,
                         NULL, NULL, NULL);
      *((float *)(tmp->array)) = p->lessthan;
      cond_l=gal_arithmetic(GAL_ARITHMETIC_OP_GE, 1, flags, ref, tmp);
      gal_data_free(tmp);
    }


  /* Now, set the final condition. If both values were specified, then use
     the GAL_ARITHMETIC_OP_OR to merge them into one. */
  switch( !isnan(p->greaterequal) + !isnan(p->lessthan) )
    {
    case 0: return;             /* No condition was specified, return.  */
    case 1:                     /* Only one condition was specified.    */
      cond = isnan(p->greaterequal) ? cond_l : cond_g;
      break;
    case 2:
      cond = gal_arithmetic(GAL_ARITHMETIC_OP_OR, 1, flagsor, cond_l, cond_g);
      break;
    }


  /* Allocate a blank value to mask all pixels that don't satisfy the
     condition. */
  blank=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL,
                       0, -1, 1, NULL, NULL, NULL);
  *((float *)(blank->array)) = NAN;


  /* Set all the pixels that satisfy the condition to blank. Note that a
     blank value will be used in the proper type of the input in the
     'where' operator.*/
  gal_arithmetic(GAL_ARITHMETIC_OP_WHERE, 1, flagsor, p->input, cond, blank);


  /* Reset the blank flags so they are checked again if necessary. */
  p->input->flag &= ~GAL_DATA_FLAG_BLANK_CH;
  p->input->flag &= ~GAL_DATA_FLAG_HASBLANK;
}





/* Check if a sorted array is necessary and if so, then make a sorted
   array. */
static void
ui_make_sorted_if_necessary(struct statisticsparams *p)
{
  int is_necessary=0;
  gal_list_i32_t *tmp;

  /* Check in the one-row outputs. */
  for(tmp=p->singlevalue; tmp!=NULL; tmp=tmp->next)
    switch(tmp->v)
      {
      case UI_KEY_MODE:
      case UI_KEY_MEDIAN:
      case UI_KEY_QUANTILE:
      case UI_KEY_QUANTFUNC:
      case UI_KEY_SIGCLIPSTD:
      case UI_KEY_SIGCLIPMEAN:
      case UI_KEY_SIGCLIPNUMBER:
      case UI_KEY_SIGCLIPMEDIAN:
        is_necessary=1;
        break;
      }

  /* Check in the rest of the outputs. */
  if( is_necessary==0 && ( p->sigmaclip || !isnan(p->mirror) ) )
    is_necessary=1;

  /* Do the sorting, we will keep the sorted array in a separate space,
     since the unsorted nature of the original dataset will help decrease
     floating point errors. If the input is already sorted, we'll just
     point it to the input.*/
  if(is_necessary)
    {
      if( gal_statistics_is_sorted(p->input, 1) )
        p->sorted=p->input;
      else
        {
          p->sorted=gal_data_copy(p->input);
          gal_statistics_sort_increasing(p->sorted);
        }
    }
}





void
ui_read_columns(struct statisticsparams *p)
{
  int toomanycols=0, tformat;
  gal_list_str_t *column=NULL;
  gal_data_t *cols, *tmp, *cinfo;
  size_t size, ncols, nrows, counter=0;
  gal_list_str_t *lines=gal_options_check_stdin(p->inputname,
                                                p->cp.stdintimeout, "input");

  /* If a reference column is also given, add it to the list of columns to
     read. */
  if(p->refcol)
    gal_list_str_add(&column, p->refcol, 0);

  /* If no column is specified, Statistics will abort and an error will be
     printed when the table has more than one column. If there is only one
     column, there is no need to specify any, so Statistics will use it. */
  if(p->column==NULL)
    {
      /* Get the basic table information. */
      cinfo=gal_table_info(p->inputname, p->cp.hdu, lines, &ncols, &nrows,
                           &tformat);
      gal_data_array_free(cinfo, ncols, 1);

      /* See how many columns it has and take the proper action. */
      switch(ncols)
        {
        case 0:
          error(EXIT_FAILURE, 0, "%s contains no usable information",
                ( p->inputname
                  ? gal_checkset_dataset_name(p->inputname, p->cp.hdu)
                  : "Standard input" ));
        case 1:
          gal_checkset_allocate_copy("1", &p->column);
          break;
        default:
          error(EXIT_FAILURE, 0, "%s is a table containing more than one "
                "column. However, the specific column to work on isn't "
                "specified.\n\n"
                "Please use the '--column' ('-c') option to specify a "
                "column. You can either give it the column number "
                "(couting from 1), or a match/search in its meta-data (e.g., "
                "column names).\n\n"
                "For more information, please run the following command "
                "(press the 'SPACE' key to go down and 'q' to return to the "
                "command-line):\n\n"
                "    $ info gnuastro \"Selecting table columns\"\n",
                ( p->inputname
                  ? gal_checkset_dataset_name(p->inputname, p->cp.hdu)
                  : "Standard input" ));
        }

    }
  gal_list_str_add(&column, p->column, 0);

  /* Read the desired column(s). */
  cols=gal_table_read(p->inputname, p->cp.hdu, lines, column, p->cp.searchin,
                      p->cp.ignorecase, p->cp.minmapsize, p->cp.quietmmap,
                      NULL);
  gal_list_str_free(lines, 1);

  /* If the input was from standard input, we can actually write this into
     it (for future reporting). */
  if(p->inputname==NULL)
    gal_checkset_allocate_copy("statistics", &p->inputname);

  /* Put the columns into the proper gal_data_t. */
  size=cols->size;
  while(cols!=NULL)
    {
      /* Pop out the top column. */
      tmp=gal_list_data_pop(&cols);

      /* Make sure it has the proper size. */
      if(tmp->size!=size)
        error(EXIT_FAILURE, 0, " read column number %zu has a %zu elements, "
              "while previous column(s) had %zu", counter, tmp->size, size);

      /* Make sure it is a usable datatype. */
      switch(tmp->type)
        {
        case GAL_TYPE_BIT:
        case GAL_TYPE_STRLL:
        case GAL_TYPE_STRING:
        case GAL_TYPE_COMPLEX32:
        case GAL_TYPE_COMPLEX64:
          error(EXIT_FAILURE, 0, " read column number %zu has a %s type, "
                "which is not currently supported by %s", counter,
                gal_type_name(tmp->type, 1), PROGRAM_NAME);
        }

      /* Put the column into the proper pointer. */
      switch(++counter)
        {
        case 1: p->input=tmp;                                         break;
        case 2: if(p->refcol) p->reference=tmp; else toomanycols=1;   break;
        default: toomanycols=1;
        }

      /* Print an error if there are too many columns: */
      if(toomanycols)
        gal_tableintern_error_col_selection(p->inputname, p->cp.hdu, "too "
                                            "many columns were selected by "
                                            "the given values to the "
                                            "'--column' and/or '--refcol' "
                                            "options. Only one is "
                                            "acceptable for each.");
    }

  /* Clean up. */
  gal_list_str_free(column, 0);
}





void
ui_preparations(struct statisticsparams *p)
{
  gal_data_t *check;
  int keepinputdir=p->cp.keepinputdir;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;
  char *checkbasename = p->cp.output ? p->cp.output : p->inputname;

  /* Change 'keepinputdir' based on if an output name was given. */
  p->cp.keepinputdir = p->cp.output ? 1 : 0;

  /* Read the input. */
  if(p->isfits && p->hdu_type==IMAGE_HDU)
    {
      p->inputformat=INPUT_FORMAT_IMAGE;
      p->input=gal_array_read_one_ch(p->inputname, cp->hdu, NULL,
                                     cp->minmapsize, p->cp.quietmmap);
      p->input->wcs=gal_wcs_read(p->inputname, cp->hdu, 0, 0,
                                 &p->input->nwcs);
      p->input->ndim=gal_dimension_remove_extra(p->input->ndim,
                                                p->input->dsize,
                                                p->input->wcs);
    }
  else
    {
      ui_read_columns(p);
      p->inputformat=INPUT_FORMAT_TABLE;
    }

  /* Read the convolution kernel if necessary. */
  if(p->sky && p->kernelname)
    {
      p->kernel=gal_fits_img_read_kernel(p->kernelname, p->khdu,
                                         cp->minmapsize, p->cp.quietmmap);
      p->kernel->ndim=gal_dimension_remove_extra(p->kernel->ndim,
                                                 p->kernel->dsize, NULL);
    }

  /* Tile and channel sanity checks and preparations. */
  if(p->ontile || p->sky)
    {
      /* Check the tiles and make the tile structure. */
      gal_tile_full_sanity_check(p->inputname, p->cp.hdu, p->input, tl);
      gal_tile_full_two_layers(p->input, tl);
      gal_tile_full_permutation(tl);

      /* Make the tile check image if requested. */
      if(tl->checktiles)
        {
          tl->tilecheckname=gal_checkset_automatic_output(cp, checkbasename,
                                                          "_tiled.fits");
          check=gal_tile_block_check_tiles(tl->tiles);
          if(p->inputformat==INPUT_FORMAT_IMAGE)
            gal_fits_img_write(check, tl->tilecheckname, NULL, PROGRAM_NAME);
          else
            {
              gal_checkset_writable_remove(tl->tilecheckname, 0,
                                           cp->dontdelete);
              gal_table_write(check, NULL, cp->tableformat, tl->tilecheckname,
                              "TABLE", 0);
            }
          gal_data_free(check);
        }

      /* Set the steps image name. */
      if(p->sky && p->checksky)
        p->checkskyname=gal_checkset_automatic_output(cp, checkbasename,
                                                      "_sky_steps.fits");
    }

  /* Set the out-of-range values in the input to blank. */
  ui_out_of_range_to_blank(p);

  /* If we are not to work on tiles, then re-order and change the input. */
  if(p->ontile==0 && p->sky==0 && p->contour==NULL)
    {
      /* Only keep the elements we want. */
      gal_blank_remove(p->input);

      /* Make sure there actually are any (non-blank) elements left. */
      if(p->input->size==0)
        error(EXIT_FAILURE, 0, "%s: all elements are blank",
              gal_fits_name_save_as_string(p->inputname, cp->hdu));

      p->input->flag &= ~GAL_DATA_FLAG_HASBLANK ;
      p->input->flag |= GAL_DATA_FLAG_BLANK_CH ;

      /* Make sure there is data remaining: */
      if(p->input->size==0)
        error(EXIT_FAILURE, 0, "%s: no data, maybe the '--greaterequal' or "
              "'--lessthan' options need to be adjusted",
              gal_fits_name_save_as_string(p->inputname, cp->hdu) );

      /* Make the sorted array if necessary. */
      ui_make_sorted_if_necessary(p);

      /* Set the number of output files. */
      if( !isnan(p->mirror) )             ++p->numoutfiles;
      if( p->histogram || p->cumulative ) ++p->numoutfiles;
    }

  /* Reset 'keepinputdir' to what it originally was. */
  p->cp.keepinputdir=keepinputdir;
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct statisticsparams *p)
{
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


  /* Prepare all the options as FITS keywords to write in output
     later. Note that in some modes, there is no output file, and
     'ui_add_to_single_value' isn't yet prepared. */
  if( (p->singlevalue && p->ontile) || p->sky || p->histogram \
      || p->cumulative)
    gal_options_as_fits_keywords(&p->cp);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct statisticsparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
  gal_data_free(p->input);
  gal_data_free(p->reference);
  gal_list_f64_free(p->tp_args);
  gal_list_i32_free(p->singlevalue);
  gal_tile_full_free_contents(&p->cp.tl);
  if(p->sorted!=p->input) gal_data_free(p->sorted);
}
