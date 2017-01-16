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
#ifndef ARGS_H
#define ARGS_H





/* Include the common options to all Gnuastro programs. If the program has
   its own HDU parsing procedure, then define the `NOT_COMMON_HDU_PARSER'
   macro before including the `commonopts.h' header. */
#include <commonopts.h>




/* Definition parameters for the argp: */
const char *argp_program_version =
  PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION "\n"
  GAL_STRINGS_COPYRIGHT
  "\n\nWritten/developed by Mohammad Akhlaghi";
const char *argp_program_bug_address=PACKAGE_BUGREPORT;
static char args_doc[] = "ASTRdata";





const char doc[] =
  /* Before the list of options: */
  GAL_STRINGS_TOP_HELP_INFO
  PROGRAM_NAME" can be used to view the information, select columns, or "
  "convert tables. The inputs and outputs can be plain text (with "
  "whitespace or comma as delimiters), FITS ascii, or FITS binary tables. "
  "The output columns can either be selected by number (counting from 1), "
  "name or using regular expressions. For regular expressions, enclose the "
  "value to the `--column' (`-c') option in slashes (`\\', as in "
  "`-c\\^mag\\'). To print the selected columns on the command-line, don't "
  "specify an output file.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b d e f g j k l m n p r u v w x y z
   A B C E F G H J L M O Q R T U W X Y Z  */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_COLUMN_KEY      = 'c',
  ARGS_OPTION_SEARCHIN_KEY    = 's',
  ARGS_OPTION_IGNORECASE_KEY  = 'I',
  ARGS_OPTION_TABLETYPE_KEY   = 't',
  ARGS_OPTION_INFORMATION_KEY = 'i',

  /* Only with long version. */
};

/* Array of acceptable options. */
struct argp_option options[] =
  {
    {
      0, 0, 0, 0,
      "Input:",
      1
    },
    {
      "column",
      ARGS_OPTION_COLUMN_KEY,
      "STR",
      0,
      "Column number (counting from 1) or search string.",
      1,
      NULL, GAL_DATA_TYPE_STRLL
    },
    {
      "searchin",
      ARGS_OPTION_SEARCHIN_KEY,
      "STR",
      0,
      "Search in column `name', `units', or `comments'.",
      1,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "ignorecase",
      ARGS_OPTION_IGNORECASE_KEY,
      0,
      0,
      "Ignore case when matching column information.",
      1,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },




    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "tabletype",
      ARGS_OPTION_TABLETYPE_KEY,
      "STR",
      0,
      "Output table type: `fits-ascii', `fits-binary'.",
      2,
      NULL, GAL_DATA_TYPE_STRING
    },




    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },
    {
      "information",
      ARGS_OPTION_INFORMATION_KEY,
      0,
      0,
      "Only print table and column information.",
      -1,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },


    {0}
  };





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
      if(p->up.filename)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->up.filename=arg;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, options, &p->cp);
    }

  return 0;
}



/* Define the child argp structure. */
static struct argp gal_options_common_child = {gal_commonopts_options,
                                               gal_options_common_argp_parse,
                                               NULL, NULL, NULL, NULL, NULL};

static struct argp_child children[]=
  {
    {&gal_options_common_child, 0, NULL, 0},
    {0, 0, 0, 0}
  };

static struct argp thisargp = {options, parse_opt, args_doc,
                               doc, children, NULL, NULL};
#endif
