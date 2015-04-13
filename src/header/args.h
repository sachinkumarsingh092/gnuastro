/*********************************************************************
Header - View and manipulate a data file header
Header is part of GNU Astronomy Utilities (Gnuastro) package.

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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef ARGS_H
#define ARGS_H

#include <argp.h>

#include "commonargs.h"
#include "linkedlist.h"
#include "fixedstringmacros.h"










/**************************************************************/
/**************        argp.h definitions       ***************/
/**************************************************************/




/* Definition parameters for the argp: */
const char *argp_program_version=SPACK_STRING"\n"COPYRIGHT
  "\n\nWritten by Mohammad Akhlaghi";
const char *argp_program_bug_address=PACKAGE_BUGREPORT;
static char args_doc[] = "ASTRdata";





const char doc[] =
  /* Before the list of options: */
  TOPHELPINFO
  SPACK_NAME" print the header information in any astronomical data file"
  "header. It can also manipulate (add, remove or modify) any of the "
  "existing keywords in a data header. \n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b e f g i j k l m n p s v x y z
   A B C E F G I J L M O R T U W X Y Z

   Number keys used: Nothing!

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
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "delete",
      'd',
      "STR",
      0,
      "Delete a keyword from the header.",
      2
    },
    {
      "rename",
      'r',
      "STR",
      0,
      "Rename keyword, keeping value and comments.",
      2
    },
    {
      "update",
      'u',
      "STR",
      0,
      "Update a keyword value or comments.",
      2
    },
    {
      "write",
      'w',
      "STR",
      0,
      "Write a keyword (with value, comments and units).",
      2
    },
    {
      "history",
      'H',
      "STR",
      0,
      "Add HISTORY keyword, any length is ok.",
      2
    },
    {
      "comment",
      'c',
      "STR",
      0,
      "Add COMMENT keyword, any length is ok.",
      2
    },
    {
      "date",
      't',
      0,
      0,
      "Set the DATE keyword to the current time.",
      2
    },



    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },
    {
      "quitonerror",
      'Q',
      0,
      0,
      "Quit if there is an error on any action.",
      -1
    },


    {0}
  };





/* Parse a single option: */
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  /* Save the arguments structure: */
  struct headerparams *p = state->input;

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
    argp_error(state, "Incorrect use of the equal sign (`=`). For short "
	       "options, `=` should not be used and for long options, "
	       "there should be no space between the option, equal sign "
	       "and value.");

  switch(key)
    {


    /* Input: */



    /* Output: */
    case 'd':
      add_to_stll(&p->delete, arg);
      break;
    case 'r':
      add_to_stll(&p->up.rename, arg);
      break;
    case 'u':
      add_to_stll(&p->up.update, arg);
      break;
    case 'w':
      add_to_stll(&p->up.write, arg);
      break;
    case 'c':
      p->comment=arg;
      break;
    case 'H':
      p->history=arg;
      break;
    case 't':
      p->date=1;
      break;

    /* Operating modes: */
    case 'Q':
      p->quitonerror=1;


    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:

      /* See what type of input value it is and put it in. */
      if( nameisfits(arg) )
        {
          if(p->up.inputname)
            argp_error(state, "Only one input image should be given.");
          else
            p->up.inputname=arg;
	}
      else
        argp_error(state, "%s is not a valid file type.", arg);
      break;





    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
	 && p->cp.printparams==0)
	{
	  if(state->arg_num==0)
	    argp_error(state, "No argument given!");
	  if(p->up.inputname==NULL)
	    argp_error(state, "No input FITS image(s) provided!");
	}
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
