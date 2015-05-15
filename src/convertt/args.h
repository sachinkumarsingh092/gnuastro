/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include "ui.h"
#include "commonargs.h"
#include "fixedstringmacros.h"


/* Definition parameters for the argp: */
const char *argp_program_version=SPACK_STRING"\n"COPYRIGHT
  "\n\nWritten by Mohammad Akhlaghi";
const char *argp_program_bug_address=PACKAGE_BUGREPORT;
static char args_doc[] = "InputFile1 [InputFile2] ... [InputFile4]";





const char doc[] =
  /* Before the list of options: */
  TOPHELPINFO
  SPACK_NAME" will convert any of the known input formats to any other "
  "of the known formats. The output file will have the same number of "
  "pixels.\n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Free letters for options:

   d e f g j k p r s t v y z
   A B E F G I J M O Q R T U W X Y Z

   Free numbers: 503
*/
static struct argp_option options[] =
  {
    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },





    {
      0, 0, 0, 0,
      "Input:",
      1
    },
    {
      "hdu2",
      500,
      "STR",
      0,
      "HDU of second input FITS image.",
      1
    },
    {
      "hdu3",
      501,
      "STR",
      0,
      "HDU of third input FITS image.",
      1
    },
    {
      "hdu4",
      502,
      "STR",
      0,
      "HDU of fourth input FITS image.",
      1
    },




    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "quality",
      'u',
      "INT",
      0,
      "Quality of output JPEG image (1 to 100).",
      2
    },
    {
      "widthincm",
      'w',
      "FLT",
      0,
      "Width in units of centimeters.",
      2
    },
    {
      "borderwidth",
      'b',
      "FLT",
      0,
      "Border width (EPS, PDF) in units of 1/72 inch.",
      2
    },
    {
      "hex",
      'x',
      0,
      0,
      "Hexadecimal encoding in EPS. Default: ASCII85.",
      2
    },




    {
      0, 0, 0, 0,
      "Flux:",
      3
    },
    {
      "fluxlow",
      'L',
      "FLT",
      0,
      "Lower flux truncation value.",
      3
    },
    {
      "fluxhigh",
      'H',
      "FLT",
      0,
      "Higher flux truncation value.",
      3
    },
    {
      "maxbyte",
      'm',
      "INT",
      0,
      "Maximum byte value for all color channels.",
      3
    },
    {
      "flminbyte",
      'i',
      0,
      0,
      "Set value of fluxlow as the minimum byte value.",
      3
    },
    {
      "fhmaxbyte",
      'a',
      0,
      0,
      "Set value of fluxhigh as the maximum byte value.",
      3
    },
    {
      "change",
      'c',
      "STR",
      0,
      "Change pixel values `from_1:to_1,from_2:to_2`.",
      3
    },
    {
      "changeaftertrunc",
      'C',
      0,
      0,
      "First truncate then change pixel values.",
      3
    },
    {
      "log",
      'l',
      0,
      0,
      "Save flux in log scale.",
      3
    },
    {
      "noinvert",
      'n',
      0,
      0,
      "Don't invert the image.",
      3
    },


    {0}
  };




















/* Parse a single option: */
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  int tmp;

  /* Save the arguments structure: */
  struct converttparams *p = state->input;

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

    /* Input:  */
    case 500:
      errno=0;
      p->up.hdu2=malloc(strlen(arg)+1);
      if(p->up.hdu2==NULL) error(EXIT_FAILURE, 0, "Space for hdu2");
      strcpy(p->up.hdu2, arg);
      p->up.hdu2set=1;
      break;
    case 501:
      errno=0;
      p->up.hdu3=malloc(strlen(arg)+1);
      if(p->up.hdu3==NULL) error(EXIT_FAILURE, 0, "Space for hdu3");
      strcpy(p->up.hdu3, arg);
      p->up.hdu3set=1;
      break;
    case 502:
      errno=0;
      p->up.hdu4=malloc(strlen(arg)+1);
      if(p->up.hdu4==NULL) error(EXIT_FAILURE, 0, "Space for hdu4");
      strcpy(p->up.hdu4, arg);
      p->up.hdu4set=1;
      break;


    /* Output: */
    case 'w':
      floatl0(arg, &p->widthincm, "widthincm", key, p->cp.spack, NULL, 0);
      p->up.widthincmset=1;
      break;
    case 'b':
      intelzero(arg, &p->borderwidth, "borderwidth", key,
                p->cp.spack, NULL, 0);
      p->up.borderwidthset=1;
      break;
    case 'u':
      intsmallerequalto(arg, &p->quality, "quality", key,
                        p->cp.spack, NULL, 0, 100);
      if(p->quality<0)
        error(EXIT_FAILURE, 0, "The quality option should be positive.");
      p->up.qualityset=1;
      break;
    case 'x':
      p->hex=1;
      break;



    /* Flux: */
    case 'L':
      anydouble(arg, &p->fluxlow, "fluxlow", key, p->cp.spack, NULL, 0);
      p->up.fluxlowset=1;
      break;
    case 'H':
      anydouble(arg, &p->fluxhigh, "fluxhigh", key, p->cp.spack, NULL, 0);
      p->up.fluxhighset=1;
      break;
    case 'm':
      intsmallerequalto(arg, &tmp, "maxbyte", key,
                        p->cp.spack, NULL, 0, UINT8_MAX);
      if(tmp<0) error(EXIT_FAILURE, 0, "--maxbyte (-m) should be positive.");
      p->maxbyte=tmp;
      p->up.maxbyteset=1;
      break;
    case 'i':
      p->flminbyte=1;
      break;
    case 'a':
      p->fhmaxbyte=1;
      break;
    case 'c':
      p->change=makechangestruct(arg);
      break;
    case 'C':
      p->changeaftertrunc=1;
      break;
    case 'l':
      p->log=1;
      break;
    case 'n':
      p->invert=0;
      break;


    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:
      add_to_stll(&p->inputnames, arg);
      ++p->numinputs;
      break;






    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
	 && p->cp.printparams==0)
	{
	  if(state->arg_num==0)
	    argp_error(state, "No argument given!");
	  if(p->inputnames==NULL)
	    argp_error(state, "No input files provided!");
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
