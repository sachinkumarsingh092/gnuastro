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
static char args_doc[] = "ASTRdata or number [ASTRdata] OPERATOR ...";





const char doc[] =
  /* Before the list of options: */
  GAL_STRINGS_TOP_HELP_INFO
  SPACK_NAME" will do arithmetic operations on one or multiple images and "
  "numbers. Simply put, the name of the image along with the arithmetic "
  "operators and possible numbers are given as arguments. The extensions of "
  "each input image are expected as options (starting with `hdu') listed "
  "below. Please note that currently "SPACK_NAME" only supports postfix "
  "or reverse polish notation. For example to get the result of `5+6', you "
  "should write `5 6 +', or to get the average of two images, you should "
  "write `a.fits b.fits + 2 /' (or more simply a.fits b.fits average). "
  "Please see the manual for more information. "
  "\n\nThe operators/functions recognized by "SPACK_NAME" are: +, -, *, /, "
  "abs, pow, sqrt, log, log10, minvalue, maxvalue, min, max, average, median, "
  "lt, le, gt, ge, eq, ne, and, or, not, isblank. Please run `info gnuastro "
  "\"Arithmetic operators\"' for detailed information on each operator. Note "
  "that multiplication should be quoted (like \"*\", or '*') to avoid shell "
  "expansion.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b c d e f g i j k l m n p r s t u v w x y z
   A B C E F G H I J L O Q R U W X Y Z

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
      "hdu",
      'h',
      "STR",
      0,
      "Nth call of this option, used for Nth input FITS.",
      1
    },
    {
      "mask",
      'M',
      "STR",
      0,
      "Mask image file name.",
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
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "type",
      'T',
      "STR",
      0,
      "uchar, short, long, longlong, float, double.",
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
  char *tokeephdu;
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
    argp_error(state, "incorrect use of the equal sign (`=`). For short "
               "options, `=` should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  switch(key)
    {

    /* Commandline options don't need to be allocated, since they are
       already in static memory and their pointers will not
       change. They also don't need to be freed for the same
       reason. However, later on, we will also be reading from the
       configuration files and there, we need to allocate space (and
       free it later). So to consistently free all the poped strings,
       we are allocating a copy here too. */
    case 'h':
      gal_checkset_allocate_copy(arg, &tokeephdu);
      gal_linkedlist_add_to_stll(&p->hdus, tokeephdu);
      break;

    /* Input: */
    case 'M':
      gal_checkset_allocate_copy_set(arg, &p->up.maskname,
                                     &p->up.masknameset);
      break;
    case 'H':
      gal_checkset_allocate_copy_set(arg, &p->up.mhdu, &p->up.mhduset);
      break;


    /* Output */
    case 'T':
      gal_checkset_known_types(arg, &p->outtype, NULL, 0);
      p->up.typeset=1;
      break;

    /* Operating modes: */



    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:

      /* Add the argument to the list of tokens: */
      gal_linkedlist_add_to_stll(&p->tokens, arg);
      break;





    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
         && p->cp.printparams==0)
        {
          if(state->arg_num==0)
            argp_error(state, "no argument given");
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
