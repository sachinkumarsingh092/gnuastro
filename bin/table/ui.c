/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2018, Free Software Foundation, Inc.

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

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" can be used to view the "
  "information, select columns, or convert tables. The inputs and outputs "
  "can be plain text (with white-space or comma as delimiters), FITS ascii, "
  "or FITS binary tables. The output columns can either be selected by "
  "number (counting from 1), name or using regular expressions. For regular "
  "expressions, enclose the value to the `--column' (`-c') option in "
  "slashes (`\\', as in `-c\\^mag\\'). To print the selected columns on the "
  "command-line, don't specify an output file.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















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


  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        /* Mandatory options. */
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
        case GAL_OPTIONS_KEY_TABLEFORMAT:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;

        /* Options to ignore. */
        case GAL_OPTIONS_KEY_TYPE:
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
  gal_tableintern_check_fits_format(p->cp.output, p->cp.tableformat);

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
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_print_info_exit(struct tableparams *p)
{
  char *tmp;
  int tableformat;
  gal_data_t *allcols;
  gal_list_str_t *lines;
  size_t i, numcols, numrows;

  /* Read the table information for the number of columns and rows. */
  lines=gal_options_check_stdin(p->filename, p->cp.stdintimeout);
  allcols=gal_table_info(p->filename, p->cp.hdu, lines, &numcols,
                         &numrows, &tableformat);
  if(p->filename==NULL) p->filename="Standard-input";
  gal_list_str_free(lines, 1);

  /* If there was no actual data in the file, then inform the user */
  if(allcols==NULL)
    error(EXIT_FAILURE, 0, "%s: no usable data rows", p->filename);


  /* Print the file information. */
  printf("--------\n");
  tmp=gal_fits_name_save_as_string(p->filename, p->cp.hdu);
  printf("%s\n", tmp);
  free(tmp);


  /* Print each column's information. */
  gal_table_print_info(allcols, numcols, numrows);


  /* Free the information from all the columns. */
  for(i=0;i<numcols;++i)
    gal_data_free_contents(&allcols[i]);
  free(allcols);


  /* Free the allocated spaces and exit. Otherwise, add the number of
     columns to the list if the user wanted to print the columns
     (didn't just want their information. */
  ui_free_report(p);
  exit(EXIT_SUCCESS);
}





/* The columns can be given as comma-separated values to one option or
   multiple calls to the column option. Here, we'll check if the input list
   has comma-separated values. If they do then the list will be updated to
   be fully separate. */
static void
ui_columns_prepare(struct tableparams *p)
{
  size_t i;
  char **strarr;
  gal_data_t *strs;
  gal_list_str_t *tmp, *new=NULL;

  /* Go over the whole original list (where each node may have more than
     one value separated by a comma. */
  for(tmp=p->columns;tmp!=NULL;tmp=tmp->next)
    {
      /* Read the different comma-separated strings into an array (within a
         `gal_data_t'). */
      strs=gal_options_parse_csv_strings_raw(tmp->v, NULL, 0);
      strarr=strs->array;

      /* Go over all the elements and add them to the `new' list. */
      for(i=0;i<strs->size;++i)
        {
          gal_list_str_add(&new, strarr[i], 0);
          strarr[i]=NULL;
        }

      /* Clean up. */
      gal_data_free(strs);
    }

  /* Delete the old list. */
  gal_list_str_free(p->columns, 1);

  /* Reverse the new list, then put it into `p->columns'. */
  gal_list_str_reverse(&new);
  p->columns=new;
}





static void
ui_preparations(struct tableparams *p)
{
  gal_list_str_t *lines;
  struct gal_options_common_params *cp=&p->cp;

  /* If there were no columns specified or the user has asked for
     information on the columns, we want the full set of columns. */
  if(p->information)
    ui_print_info_exit(p);

  /* Prepare the column names. */
  ui_columns_prepare(p);

  /* Read in the table columns. */
  lines=gal_options_check_stdin(p->filename, p->cp.stdintimeout);
  p->table=gal_table_read(p->filename, cp->hdu, lines, p->columns,
                          cp->searchin, cp->ignorecase, cp->minmapsize,
                          NULL);
  if(p->filename==NULL) p->filename="Standard-input";
  gal_list_str_free(lines, 1);

  /* If there was no actual data in the file, then inform the user and
     abort. */
  if(p->table==NULL)
    error(EXIT_FAILURE, 0, "%s: no usable data rows (non-commented and "
          "non-blank lines)", p->filename);

  /* Now that the data columns are ready, we can free the string linked
     list. */
  gal_list_str_free(p->columns, 1);
  p->columns=NULL;
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
ui_free_report(struct tableparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
  gal_list_data_free(p->table);
}
