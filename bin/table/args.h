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

#include <argp.h>

#include <gnuastro/linkedlist.h>

#include <commonargs.h>
#include <fixedstringmacros.h>










/**************************************************************/
/**************        argp.h definitions       ***************/
/**************************************************************/




/* Definition parameters for the argp: */
const char *argp_program_version=SPACK_STRING"\n"GAL_STRINGS_COPYRIGHT
  "\n\nWritten by Mohammad Akhlaghi";
const char *argp_program_bug_address=PACKAGE_BUGREPORT;
static char args_doc[] = "ASTRdata";





const char doc[] =
  /* Before the list of options: */
  GAL_STRINGS_TOP_HELP_INFO
  SPACK_NAME" prints (in a human readable format) a FITS table or its "
  "information. The output columns can either be selected by number, name "
  "or using regular expressions. The format of their printing can also "
  "be set (based on the type of data in the column).\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b d e f g j k l m n p r u v w x y z
   A B C E F G H J L M O Q R T U W X Y Z

   Number keys used (larger than 1000): 1000

   Options with keys (second structure element) larger than 500 do not
   have a short version.
 */
static struct argp_option options[] =
  {
    {
      0, 0, 0, 0,
      "Input:",
      1
    },
    {
      "column",
      'c',
      "STR",
      0,
      "Column number (counting from 1) or search string.",
      1
    },
    {
      "searchin",
      's',
      "STR",
      0,
      "Search for columns in `name', `units', `comments'.",
      1
    },
    {
      "ignorecase",
      'I',
      0,
      0,
      "Ignore case when matching column names.",
      1
    },




    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "fitstabletype",
      't',
      "STR",
      0,
      "Only `ascii', or `binary' are acceptable.",
      2
    },




    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },
    {
      "information",
      'i',
      0,
      0,
      "Only print table and columns information.",
      -1
    },


    {0}
  };





/* Parse a single option: */
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  /* Save the arguments structure: */
  char *tstring;
  struct tableparams *p = state->input;

  /* Set the pointer to the common parameters for all programs
     here: */
  state->child_inputs[0]=&p->cp;

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

  switch(key)
    {


    /* Input: */
    case 'c':
      gal_checkset_allocate_copy(arg, &tstring);
      gal_linkedlist_add_to_stll(&p->columns, tstring);
      break;

    case 's':
      gal_checkset_allocate_copy(arg, &p->up.searchin);
      p->up.searchinset=1;
      break;

    case 'I':
      p->ignorecase=1;
      p->up.ignorecaseset=1;
      break;


    /* Output: */
    case 't':
      gal_checkset_allocate_copy(arg, &p->up.fitstabletype);
      p->up.fitstabletypeset=1;
      break;


    /* Operating modes: */
    case 'i':
      p->information=1;
      p->up.informationset=1;
      break;


    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      if(p->up.filename)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->up.filename=arg;
      break;


    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
         && p->cp.printparams==0)
        if(state->arg_num==0)
          argp_error(state, "no argument given");
      break;





    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}





/* Specify the children parsers: */
struct argp_child children[]=
  {
    {&commonargp, 0, NULL, 0},
    {0, 0, 0, 0}
  };





/* Basic structure defining the whole argument reading process. */
static struct argp thisargp = {options, parse_opt, args_doc,
                               doc, children, NULL, NULL};

#endif
