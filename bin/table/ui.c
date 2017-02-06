/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/linkedlist.h>

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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" can be used to view the "
  "information, select columns, or convert tables. The inputs and outputs "
  "can be plain text (with whitespace or comma as delimiters), FITS ascii, "
  "or FITS binary tables. The output columns can either be selected by "
  "number (counting from 1), name or using regular expressions. For regular "
  "expressions, enclose the value to the `--column' (`-c') option in "
  "slashes (`\\', as in `-c\\^mag\\'). To print the selected columns on the "
  "command-line, don't specify an output file.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b d e f g j k l m n p r s t u v w x y z
   A B C E F G H J L M O Q R T U W X Y Z  */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_COLUMN      = 'c',
  ARGS_OPTION_KEY_INFORMATION = 'i',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};



















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct tableparams *p,
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


  /* Set the mandatory common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    switch(cp->coptions[i].key)
      {
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
  struct tableparams *p = state->input;

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
      if(p->filename)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->filename=arg;
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
ui_read_check_only_options(struct tableparams *p)
{

  /* Check if the format of the output table is valid, given the type of
     the output. */
  gal_table_check_fits_format(p->cp.output, p->cp.tableformat);

}





static void
ui_check_options_and_arguments(struct tableparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->filename)
    {
      if( gal_fits_name_is_fits(p->filename) && p->cp.hdu==NULL )
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
void
ui_preparations(struct tableparams *p)
{
  char *numstr;
  int tableformat;
  gal_data_t *allcols;
  size_t i, numcols, numrows;
  struct gal_options_common_params *cp=&p->cp;

  /* If there were no columns specified, we want the full set of
     columns. */
  if(p->columns==NULL)
    {
      /* Read the table information for the number of columns and rows. */
      allcols=gal_table_info(p->filename, cp->hdu, &numcols,
                             &numrows, &tableformat);

      /* If there was no actual data in the file, then inform the user */
      if(allcols==NULL)
        error(EXIT_FAILURE, 0, "%s: no usable data rows", p->filename);


      /* If the user just wanted information, then print it. */
      if(p->information)
        {
          /* Print the file information. */
          printf("--------\n");
          printf("%s", p->filename);
          if(gal_fits_name_is_fits(p->filename))
            printf(" (hdu: %s)\n", cp->hdu);
          else
            printf("\n");

          /* Print each column's information. */
          gal_table_print_info(allcols, numcols, numrows);
        }


      /* Free the information from all the columns. */
      for(i=0;i<numcols;++i)
        gal_data_free_contents(&allcols[i]);
      free(allcols);


      /* If the user just wanted information, then free the allocated
         spaces and exit. Otherwise, add the number of columns to the list
         if the user wanted to print the columns (didn't just want their
         information. */
      if(p->information)
        {
          ui_free_report(p);
          exit(EXIT_SUCCESS);
        }
      else
        for(i=1;i<=numcols;++i)
          {
            asprintf(&numstr, "%zu", i);
            gal_linkedlist_add_to_stll(&p->columns, numstr, 0);
          }
    }

  /* Reverse the list of column search criteria that we are looking for
     (since this is a last-in-first-out linked list, the order that
     elements were added to the list is the reverse of the order that they
     will be popped). */
  gal_linkedlist_reverse_stll(&p->columns);
  p->table=gal_table_read(p->filename, cp->hdu, p->columns, cp->searchin,
                          cp->ignorecase, cp->minmapsize);

  /* If there was no actual data in the file, then inform the user and
     abort. */
  if(p->table==NULL)
    error(EXIT_FAILURE, 0, "%s: no usable data rows (non-commented and "
          "non-blank lines)", p->filename);

  /* Now that the data columns are ready, we can free the string linked
     list. */
  gal_linkedlist_free_stll(p->columns, 1);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct tableparams *p)
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
ui_free_report(struct tableparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
  free(p->cp.searchinstr);
  free(p->cp.tableformatstr);
  gal_data_free_ll(p->table);
}
