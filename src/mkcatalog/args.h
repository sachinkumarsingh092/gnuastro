/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

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
  SPACK_NAME" will create a catalog from an input, labeled, and noise "
  "identification images. For each pixel, there should at least be one "
  "labeled image which has the same size as the input but each pixel is "
  "labeled with the ID of the object. A clump labeled image can also be "
  "optionally provided. If the noise mean and standard deviation are "
  "provided, Signal to noise measures will also be output.\n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b d e f g i j k l m n p r u v w x y z
   A B C E F G I J L Q R T U W X Y Z

   Number keys used: <=504

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
      "Mask image header name or number.",
      1
    },
    {
      "objlabs",
      'O',
      "STR",
      0,
      "Image specifying object labeles.",
      1
    },
    {
      "objhdu",
      501,
      "STR",
      0,
      "Object image header name or number.",
      1
    },
    {
      "clumplabs",
      'c',
      "STR",
      0,
      "Image specifying clump labeles.",
      1
    },
    {
      "clumphdu",
      502,
      "STR",
      0,
      "Clumps image header name or number.",
      1
    },
    {
      "sky",
      's',
      "STR",
      0,
      "Sky value image.",
      1
    },
    {
      "skyhdu",
      503,
      "STR",
      0,
      "Sky image header name or number.",
      1
    },
    {
      "std",
      't',
      "STR",
      0,
      "Sky standard deviation image.",
      1
    },
    {
      "stdhdu",
      504,
      "STR",
      0,
      "Sky STD image header name or number.",
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
  struct mkcatalogparams *p = state->input;

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
    case 'M':
      allocatecopyset(arg, &p->up.maskname, &p->up.masknameset);
      break;
    case 'H':
      allocatecopyset(arg, &p->up.mhdu, &p->up.mhduset);
      break;
    case 'O':
      allocatecopyset(arg, &p->up.objlabsname, &p->up.objlabsnameset);
      break;
    case 501:
      allocatecopyset(arg, &p->up.objhdu, &p->up.objhduset);
      break;
    case 'c':
      allocatecopyset(arg, &p->up.clumplabsname,
                      &p->up.clumplabsnameset);
      break;
    case 502:
      allocatecopyset(arg, &p->up.clumphdu, &p->up.clumphduset);
      break;
    case 's':
      allocatecopyset(arg, &p->up.skyname, &p->up.skynameset);
      break;
    case 503:
      allocatecopyset(arg, &p->up.skyhdu, &p->up.skyhduset);
      break;
    case 't':
      allocatecopyset(arg, &p->up.stdname, &p->up.stdnameset);
      break;
    case 504:
      allocatecopyset(arg, &p->up.stdhdu, &p->up.stdhduset);
      break;


    /* Output: */


    /* Operating modes: */


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
