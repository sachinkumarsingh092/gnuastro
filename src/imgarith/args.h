/*********************************************************************
ImageArithmetic - Do arithmetic operations on images.
ImageArithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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
static char args_doc[] = "ASTRdata ASTRdata OPERATOR ...";





const char doc[] =
  /* Before the list of options: */
  TOPHELPINFO
  SPACK_NAME" will do arithmetic operations on one or multiple images.\n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   k l m n p r s t u v w x y z
   A B C E F G H I J L M O Q R T U W X Y Z

   Number keys used: <=500

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
      "mhdu",
      'H',
      "STR",
      0,
      "Mask image header name.",
      1
    },
    {
      "hdu1",
      'a',
      "STR",
      0,
      "2nd image extension name.",
      1
    },
    {
      "hdu2",
      'b',
      "STR",
      0,
      "3rd image extension name.",
      1
    },
    {
      "hdu3",
      'c',
      "STR",
      0,
      "4th image extension name.",
      1
    },
    {
      "hdu4",
      'd',
      "STR",
      0,
      "5th image extension name.",
      1
    },
    {
      "hdu5",
      'e',
      "STR",
      0,
      "6th image extension name.",
      1
    },
    {
      "hdu6",
      'f',
      "STR",
      0,
      "7th image extension name.",
      1
    },
    {
      "hdu7",
      'g',
      "STR",
      0,
      "8th image extension name.",
      1
    },
    {
      "hdu8",
      'i',
      "STR",
      0,
      "9th image extension name.",
      1
    },
    {
      "hdu9",
      'j',
      "STR",
      0,
      "10th image extension name.",
      1
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
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  /* Save the arguments structure: */
  int junkset=0;
  struct imgarithparams *p = state->input;

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
    case 'H':
      allocatecopyset(arg, &p->up.mhdu, &p->up.mhduset);
      break;
    case 'a':
      allocatecopyset(arg, &p->up.hdus[1], &junkset);
      break;
    case 'b':
      allocatecopyset(arg, &p->up.hdus[2], &junkset);
      break;
    case 'c':
      allocatecopyset(arg, &p->up.hdus[3], &junkset);
      break;
    case 'd':
      allocatecopyset(arg, &p->up.hdus[4], &junkset);
      break;
    case 'e':
      allocatecopyset(arg, &p->up.hdus[5], &junkset);
      break;
    case 'f':
      allocatecopyset(arg, &p->up.hdus[6], &junkset);
      break;
    case 'g':
      allocatecopyset(arg, &p->up.hdus[7], &junkset);
      break;
    case 'i':
      allocatecopyset(arg, &p->up.hdus[8], &junkset);
      break;
    case 'j':
      allocatecopyset(arg, &p->up.hdus[9], &junkset);
      break;


    /* Operating modes: */



    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:

      /* Add the argument to the list of tokens: */
      add_to_stll(&p->tokens, arg);
      break;





    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
	 && p->cp.printparams==0)
	{
	  if(state->arg_num==0)
	    argp_error(state, "No argument given!");
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
