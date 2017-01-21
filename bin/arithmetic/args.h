/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef ARGS_H
#define ARGS_H





/* Include necessary headers. */
#define NOT_COMMON_HDU_PARSER 1
#include <commonopts.h>
#include <authors-cite.h>
#include <fixedstringmacros.h>



/* Definition parameters for the argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "ASTRdata or number [ASTRdata] OPERATOR ...";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will do arithmetic "
  "operations on one or multiple images and numbers. Simply put, the name "
  "of the image along with the arithmetic operators and possible numbers "
  "are given as arguments. The extensions of each input image are expected "
  "as options (starting with `hdu') listed below. Please note that currently "
  PROGRAM_NAME" only supports postfix or reverse polish notation. For "
  "example to get the result of `5+6', you should write `5 6 +', or to get "
  "the average of two images, you should write `a.fits b.fits + 2 /' (or "
  "more simply a.fits b.fits average). Please see the manual for more "
  "information. "
  "\n\nThe operators/functions recognized by "PROGRAM_NAME" are: +, -, *, "
  "/, abs, pow, sqrt, log, log10, minvalue, maxvalue, min, max, average, "
  "median, lt, le, gt, ge, eq, ne, and, or, not, isblank, and the full set "
  "of bitwise operators. Please run `info gnuastro \"Arithmetic "
  "operators\"' for detailed information on each operator. Note that "
  "multiplication should be quoted (like \"*\", or '*') to avoid shell "
  "expansion.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b c d e f g i j k l m n p r s t u v w x y z
   A B C E F G H I J L M O Q R T U W X Y Z                */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_HDU_KEY        = 'h',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};

static struct argp_option options[] =
  {
    {
      0, 0, 0, 0,
      "Input:",
      1
    },
    {
      "hdu",
      ARGS_OPTION_HDU_KEY,
      "STR",
      0,
      "Nth call of this option, used for Nth input FITS.",
      1,
      NULL, GAL_DATA_TYPE_STRLL
    },




    {
      0, 0, 0, 0,
      "Output:",
      2
    },




    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },


    {0}
  };





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct imgarithparams *p = state->input;

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
      gal_linkedlist_add_to_stll(&p->tokens, arg, 1);
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, options, &p->cp);
    }

  return 0;
}





/* Define the child argp structure. */
static struct argp
gal_options_common_child = {gal_commonopts_options,
                            gal_options_common_argp_parse,
                            NULL, NULL, NULL, NULL, NULL};

/* Use the child argp structure in list of children (only one for now). */
static struct argp_child
children[]=
{
  {&gal_options_common_child, 0, NULL, 0},
  {0, 0, 0, 0}
};

/* Set all the necessary argp parameters. */
static struct argp
thisargp = {options, parse_opt, args_doc, doc, children, NULL, NULL};
#endif
