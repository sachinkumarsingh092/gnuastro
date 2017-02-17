/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

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
#include <gnuastro/qsort.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/linkedlist.h>
#include <gnuastro/statistics.h>

#include <timing.h>
#include <options.h>
#include <checkset.h>
#include <fixedstringmacros.h>

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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will print the basic "
  "statistics of the input image pixel flux distribution. All blank pixels "
  "or pixels specified by a mask image will be ignored.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Option groups particular to this program. */
enum program_args_groups
{
  ARGS_GROUP_IN_ONE_ROW = GAL_OPTIONS_GROUP_AFTER_COMMON,
  ARGS_GROUP_PARTICULAR_STAT,
  ARGS_GROUP_HIST_CFP,
};




/* Available letters for short options:

   a b c d e f i j k p r u v w x y z
   B E F G J L R W X Y Z
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_COLUMN       = 'c',
  ARGS_OPTION_KEY_GREATEREQUAL = 'g',
  ARGS_OPTION_KEY_LESSTHAN     = 'l',
  ARGS_OPTION_KEY_QUANTRANGE   = 'Q',
  ARGS_OPTION_KEY_MEAN         = 'm',
  ARGS_OPTION_KEY_STD          = 't',
  ARGS_OPTION_KEY_MEDIAN       = 'M',
  ARGS_OPTION_KEY_MODE         = 'O',
  ARGS_OPTION_KEY_ASCIIHIST    = 'A',
  ARGS_OPTION_KEY_HISTOGRAM    = 'H',
  ARGS_OPTION_KEY_CUMULATIVE   = 'C',
  ARGS_OPTION_KEY_SIGMACLIP    = 's',
  ARGS_OPTION_KEY_NORMALIZE    = 'n',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  ARGS_OPTION_KEY_NUMBER       = 1000,
  ARGS_OPTION_KEY_MINIMUM,
  ARGS_OPTION_KEY_MAXIMUM,
  ARGS_OPTION_KEY_SUM,
  ARGS_OPTION_KEY_ASCIICFP,
  ARGS_OPTION_KEY_MIRROR,
  ARGS_OPTION_KEY_NUMBINS,
  ARGS_OPTION_KEY_LOWERBIN,
  ARGS_OPTION_KEY_ONEBINSTART,
  ARGS_OPTION_KEY_MAXBINONE,
};



















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

  /* Program-specific initializers */
  p->lessthan            = NAN;
  p->onebinstart         = NAN;
  p->mirrorquant         = NAN;
  p->greaterequal        = NAN;
  p->quantilerange       = NAN;

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





static void *
ui_add_to_print_in_row(struct argp_option *option, char *arg,
                       char *filename, size_t lineno, void *params)
{
  struct statisticsparams *p=(struct statisticsparams *)params;

  /* If this option is given in a configuration file, then `arg' will not
     be NULL and we don't want to do anything if it is `0'. */
  if( arg && arg[1]!='\0' && *arg!='0' && *arg!='1' )
    error_at_line(EXIT_FAILURE, 0, filename, lineno, "the `--%s' "
                  "option takes no arguments. In a configuration file "
                  "it can only have the values `1' or `0', indicating "
                  "if it should be used or not", option->name);

  /* Only proceed if the (possibly given) argument is 1. */
  if(arg && *arg=='0') return NULL;

  /* Add this option to the print list. */
  gal_linkedlist_add_to_ill(&p->toprint, option->key);

  return NULL;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct statisticsparams *p)
{
  /* Check if the format of the output table is valid, given the type of
     the output. */
  gal_table_check_fits_format(p->cp.output, p->cp.tableformat);


  /* If less than and greater than are both given, make sure that the value
     to greater than is smaller than the value to less-than. */
  if( !isnan(p->lessthan) && !isnan(p->greaterequal)
      && p->lessthan < p->greaterequal )
    error(EXIT_FAILURE, 0, "the value to `--lessthan' (%g) must be larger "
          "than the value to `--greaterequal' (%g)", p->lessthan,
          p->greaterequal);


  /* Less-than and greater-equal cannot be called together with
     quantrange. */
  if( ( !isnan(p->lessthan) || !isnan(p->greaterequal) )
      && !isnan(p->quantilerange) )
    error(EXIT_FAILURE, 0, "`--lessthan' and/or `--greaterequal' cannot "
          "be called together with `--quantrange'");


  /* The quantile range only makes sense with value less than 0.5. */
  if( !isnan(p->quantilerange) && p->quantilerange>=0.5)
    error(EXIT_FAILURE, 0, "%g>=0.5! The quantile range must be less than "
          "0.5. Recall the the quantile (Q) range is defined to be: Q to "
          "1-Q", p->quantilerange);


  /* When binned outputs are requested, make sure that `numbins' is set. */
  if(p->histogram || p->cumulative)
    error(EXIT_FAILURE, 0, "`--numbins' isn't set. When the histogram or "
          "cumulative frequency plots are requested, the number of bins "
          "(`--numbins') is necessary");
}





static void
ui_check_options_and_arguments(struct statisticsparams *p)
{
  char *name=NULL;

  if(p->inputname)
    {
      /* If input is FITS. */
      if( (p->isfits=gal_fits_name_is_fits(p->inputname)) )
        {
          /* Check if a HDU is given. */
          if( p->cp.hdu==NULL )
            error(EXIT_FAILURE, 0, "no HDU specified. When the input is a "
                  "FITS file, a HDU must also be specified, you can use "
                  "the `--hdu' (`-h') option and give it the HDU number "
                  "(starting from zero), extension name, or anything "
                  "acceptable by CFITSIO");

          /* If its a table, make sure a column is also specified. */
          p->hdu_type=gal_fits_hdu_type(p->inputname, p->cp.hdu);
          if(p->hdu_type==IMAGE_HDU)
            {
              if(p->column)
                error(EXIT_FAILURE, 0, "%s (hdu: %s): is a FITS image "
                      "extension. The `--column' option is only applicable "
                      "to tables.", p->inputname, p->cp.hdu);
            }
          else if(p->column==NULL)
            asprintf(&name, "%s (hdu: %s)", p->inputname, p->cp.hdu);
        }

      /* If its not FITS, it must be a table. */
      else
        {
          if(p->column==NULL) name=p->inputname;
        }

      /* If a column was necessary, but not given, print an error. */
      if(name)
        error(EXIT_FAILURE, 0, "%s is a table but no column is "
              "specified. Please use the `--column' (`-c') option to "
              "specify a column.\n\nYou can either give it the column number "
              "(couting from 1), or a match/search in its meta-data (e.g., "
              "column names). For more information, please run the "
              "following command (press the `SPACE' key to go down and "
              "`q' to return to the command-line):\n\n"
              "    $ info gnuastro \"Selecting table columns\"\n", name);
    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_print_one_number(gal_data_t *out)
{
  char *toprint=gal_data_write_to_string(out->array, out->type, 0);
  printf("%s ", toprint);
  gal_data_free(out);
  free(toprint);
}





void
ui_print_one_row(struct statisticsparams *p)
{
  struct gal_linkedlist_ill *tmp;

  /* Reverse the list to have the same order the user wanted. */
  gal_linkedlist_reverse_ill(&p->toprint);

  /* Print the numbers. */
  for(tmp=p->toprint; tmp!=NULL; tmp=tmp->next)
    switch(tmp->v)
      {
      case ARGS_OPTION_KEY_NUMBER:
        ui_print_one_number( gal_statistics_number(p->input) );      break;
      case ARGS_OPTION_KEY_MINIMUM:
        ui_print_one_number( gal_statistics_minimum(p->input) );     break;
      case ARGS_OPTION_KEY_MAXIMUM:
        ui_print_one_number( gal_statistics_maximum(p->input) );     break;
      case ARGS_OPTION_KEY_SUM:
        ui_print_one_number( gal_statistics_sum(p->input) );         break;
      case ARGS_OPTION_KEY_MEAN:
        ui_print_one_number( gal_statistics_mean(p->input) );        break;
      case ARGS_OPTION_KEY_STD:
        ui_print_one_number( gal_statistics_std(p->input) );         break;
      case ARGS_OPTION_KEY_MEDIAN:
        ui_print_one_number( gal_statistics_median(p->input) );      break;
      case ARGS_OPTION_KEY_MODE:
        error(EXIT_FAILURE, 0, "mode isn't implemented yet!");
        break;
      default:
        error(EXIT_FAILURE, 0, "A bug! Operation code %d not recognized in "
              "`ui_print_row'. Please contact us at %s so we can address "
              "the problem", tmp->v, PACKAGE_BUGREPORT);
      }

  /* Print a new line. */
  printf("\n");
}





static void
ui_out_of_range_to_blank(struct statisticsparams *p)
{
  size_t one=1;
  unsigned char flags=GAL_ARITHMETIC_NUMOK;
  gal_data_t *tmp, *cond_g=NULL, *cond_l=NULL, *cond;
  unsigned char flagsor = ( GAL_ARITHMETIC_FREE
                            | GAL_ARITHMETIC_INPLACE
                            | GAL_ARITHMETIC_NUMOK );

  /* Set the condition. Note that the `greaterequal' name is for the data
     we want. So we will set the condition based on those that are
     less-than  */
  if(!isnan(p->greaterequal))
    {
      tmp=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT32, 1, &one, NULL, 0, -1,
                        NULL, NULL, NULL);
      *((float *)(tmp->array)) = p->greaterequal;
      cond_g=gal_arithmetic(GAL_ARITHMETIC_OP_LT, flags, p->input, tmp);
      gal_data_free(tmp);
    }

  /* Same reasoning as above for `p->greaterthan'. */
  if(!isnan(p->lessthan))
    {
      tmp=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT32, 1, &one, NULL, 0, -1,
                        NULL, NULL, NULL);
      *((float *)(tmp->array)) = p->lessthan;
      cond_l=gal_arithmetic(GAL_ARITHMETIC_OP_GE, flags, p->input, tmp);
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
      cond = gal_arithmetic(GAL_ARITHMETIC_OP_OR, flagsor, cond_l, cond_g);
      break;
    }

  /* Allocate a blank value to mask all input pixels. */
  tmp=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT32, 1, &one, NULL,
                     0, -1, NULL, NULL, NULL);
  *((float *)(tmp->array)) = NAN;

  /* Set all the pixels that satisfy the condition to blank. */
  gal_arithmetic(GAL_ARITHMETIC_OP_WHERE, flagsor, p->input, cond, tmp);
}





void
ui_preparations(struct statisticsparams *p)
{
  gal_data_t *tmp;
  char *errorstring;
  size_t numcolmatches=0;
  struct gal_linkedlist_stll *column=NULL;

  /* Read the input: no matter if its an image or a table column. */
  if(p->isfits && p->hdu_type==IMAGE_HDU)
    p->input=gal_fits_img_read(p->inputname, p->cp.hdu, p->cp.minmapsize);
  else
    {
      /* Read the input column. */
      gal_linkedlist_add_to_stll(&column, p->column, 0);
      p->input=gal_table_read(p->inputname, p->cp.hdu, column, p->cp.searchin,
                              p->cp.ignorecase, p->cp.minmapsize);
      if(p->input->next)
        {
          for(tmp=p->input;tmp!=NULL;tmp=tmp->next) ++numcolmatches;
          asprintf(&errorstring, "%zu columns were selected with `%s' "
                   "(value to `--column' option). In this context, "
                   "Statistics can only work on one data-set (column in a "
                   "table).", numcolmatches, p->column);
          gal_table_error_col_selection(p->inputname, p->cp.hdu, errorstring);
        }

      /* Clean up. */
      gal_linkedlist_free_stll(column, 0);
    }

  /* Set the out-of-range values in the input to blank. */
  ui_out_of_range_to_blank(p);

  /* Only keep the numbers we want. */
  gal_blank_remove(p->input);

  /* Print the one-row numbers if the user asked for them. We want this to
     be done before converting the array to a float, since these operations
     can be done on any type and if the user has just asked for these
     operations, we don't want to waste their time for nothing. */
  if(p->toprint) ui_print_one_row(p);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct statisticsparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program (`args.h')
     and for the common options to all Gnuastro (`commonopts.h'). We want
     to directly put the pointers to the fields in `p' and `cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <commonopts.h>
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
ui_free_report(struct statisticsparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
}
