/*********************************************************************
ImageStatistics - Get general statistics about the image.
ImgeStatistics is part of GNU Astronomy Utilities (Gnuastro) package.

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
  SPACK_NAME" will print the basic statistics of the input image pixel "
  "flux distribution. All blank pixels or pixels specified by a mask "
  "image will be ignored.\n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   c e f j k m s t v w y z
   C E F G I J L O R T W X Y Z

   Number keys used: <=511

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
      "Mask image header name.",
      1
    },
    {
      "ignoremin",
      'r',
      0,
      0,
      "Ignore data with values equal to minimum.",
      1
    },
    {
      "mirrordist",
      'd',
      "FLT",
      0,
      "Distance beyond mirror point. Multiple of std.",
      1
    },





    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "lowerbin",
      'l',
      0,
      0,
      "Interval lower limit for column 1.",
      2
    },
    {
      "onebinvalue",
      'B',
      "FLT",
      0,
      "Shift bins so one bin starts on this value.",
      2
    },
    {
      "noasciihist",
      'A',
      0,
      0,
      "Do not show an ASCII histogram of the data.",
      2
    },
    {
      "checkmode",
      509,
      0,
      0,
      "Mode mirror plot. `_modehist.txt', `_modecfp.txt'",
      2
    },
    {
      "mirrorquant",
      510,
      "FLT",
      0,
      "Mirror quantile. `_mirhist.txt', `_mircfp.txt'.",
      2
    },
    {
      "histrangeformirror",
      511,
      0,
      0,
      "Use input histogram range for mirror plots.",
      2
    },
    {
      "mirrorplotdist",
      503,
      "FLT",
      0,
      "Distance beyond mode to display.",
      2
    },



    {
      0, 0, 0, 0,
      "Histogram (suffix: `_hist.txt'):",
      3
    },
    {
      "nohist",
      500,
      0,
      0,
      "Do not calculate histogram.",
      3
    },
    {
      "normhist",
      501,
      0,
      0,
      "Normalize the  histogram (sum of all bins 1).",
      3
    },
    {
      "maxhistone",
      502,
      0,
      0,
      "Scale such that the maximum bin has value of one.",
      3
    },
    {
      "histnumbins",
      'n',
      "INT",
      0,
      "Number of bins in the histogram.",
      3
    },
    {
      "histmin",
      'i',
      "FLT",
      0,
      "The minimum value for the histogram.",
      3
    },
    {
      "histmax",
      'x',
      "FLT",
      0,
      "The maximum value for the histogram.",
      3
    },
    {
      "histquant",
      'Q',
      "FLT",
      0,
      "Quantile (Q) range. Histogram from Q to 1-Q.",
      3
    },




    {
      0, 0, 0, 0,
      "Cumulative Frequency Plot (suffix: `_cfp.txt'):",
      4
    },
    {
      "nocfp",
      504,
      0,
      0,
      "No Cumulative Frequency Plot.",
      4
    },
    {
      "normcfp",
      505,
      0,
      0,
      "Normalize the CFP (sum of all bins 1).",
      4
    },
    {
      "maxcfpeqmaxhist",
      506,
      0,
      0,
      "Set maximum of CFP to maximum of histogram.",
      4
    },
    {
      "cfpsimhist",
      507,
      0,
      0,
      "Set CFP range and bins similar to histogram.",
      4
    },
    {
      "cfpnum",
      'p',
      "INT",
      0,
      "Number of data points to find CFP.",
      4
    },
    {
      "cfpmin",
      'a',
      "FLT",
      0,
      "Minimum value to use in the CFP.",
      4
    },
    {
      "cfpmax",
      'b',
      "FLT",
      0,
      "Maximum value to use in the CFP.",
      4
    },
    {
      "cfpquant",
      'U',
      "FLT",
      0,
      "Quantile of range: from U to 1-U.",
      4
    },



    {
      0, 0, 0, 0,
      "Sigma clipping:",
      5
    },
    {
      "nosigclip",
      508,
      0,
      0,
      "Do not preform sigma clipping.",
      5
    },
    {
      "sigclipmultip",
      'u',
      "FLT",
      0,
      "Multiple of standard deviation in sigma-clipping.",
      5
    },
    {
      "sigcliptolerance",
      't',
      "FLT",
      0,
      "Difference in STD tolerance to halt iteration.",
      5
    },
    {
      "sigclipnum",
      'g',
      "INT",
      0,
      "Number of times to do sigma clipping.",
      5
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
  struct imgstatparams *p = state->input;

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
    case 'r':
      p->ignoremin=1;
      break;
    case 'd':
      floatl0(arg, &p->mirrordist, "mirrordist", key, SPACK, NULL, 0);
      p->up.mirrordistset=1;
      break;

    /* Output: */
    case 'l':
      p->lowerbin=1;
      break;
    case 'B':
      anyfloat(arg, &p->onebinvalue, "onebinvalue", key, SPACK, NULL, 0);
      p->up.onebinvalueset=1;
      break;
    case 'A':
      p->asciihist=0;
      break;
    case 509:
      p->mhistname="a";
      break;
    case 510:
      floatl0s1(arg, &p->mirror, "mirrorquant", key, SPACK, NULL, 0);
      break;
    case 511:
      p->histrangeformirror=1;
      break;
    case 503:
      floatl0(arg, &p->mirrorplotdist, "mirrorplotdist", key, SPACK, NULL, 0);
      p->up.mirrorplotdistset=1;
      break;

    /* Histogram */
    case 500:
      p->histname=NULL;
      break;
    case 501:
      p->normhist=1;
      break;
    case 502:
      p->maxhistone=1;
      break;
    case 'n':
      sizetlzero(arg, &p->histnumbins, "histnumbins", key, SPACK, NULL, 0);
      p->up.histnumbinsset=1;
      break;
    case 'i':
      anyfloat(arg, &p->histmin, "histmin", key, SPACK, NULL, 0);
      p->up.histminset=1;
      break;
    case 'x':
      anyfloat(arg, &p->histmax, "histmax", key, SPACK, NULL, 0);
      p->up.histmaxset=1;
      break;
    case 'Q':
      floatl0s1(arg, &p->histquant, "histquant", key, SPACK, NULL, 0);
      p->up.histquantset=1;
      break;

    /* Cumulative frequency plot: */
    case 504:
      p->cfpname=NULL;
      break;
    case 505:
      p->normcfp=1;
      break;
    case 506:
      p->maxcfpeqmaxhist=1;
      break;
    case 507:
      p->cfpsimhist=1;
      break;
    case 'p':
      sizetlzero(arg, &p->cfpnum, "cfpnum", key, SPACK, NULL, 0);
      p->up.cfpnumset=1;
      break;
    case 'a':
      anyfloat(arg, &p->cfpmin, "cfpmin", key, SPACK, NULL, 0);
      p->up.cfpminset=1;
      break;
    case 'b':
      anyfloat(arg, &p->cfpmax, "cfpmax", key, SPACK, NULL, 0);
      p->up.cfpmaxset=1;
      break;
    case 'U':
      floatl0s1(arg, &p->cfpquant, "cfpquant", key, SPACK, NULL, 0);
      p->up.cfpquantset=1;
      break;


    /* Sigma clipping: */
    case 508:
      p->sigclip=0;
      break;
    case 'u':
      floatl0(arg, &p->sigclipmultip, "sigclipmultip", key, SPACK, NULL, 0);
      p->up.sigclipmultipset=1;
      break;
    case 't':
      floatl0(arg, &p->sigcliptolerance, "sigcliptolerance", key, SPACK,
              NULL, 0);
      p->up.sigcliptoleranceset=1;
      break;
    case 'g':
      sizetlzero(arg, &p->sigclipnum, "sigclipnum", key, SPACK, NULL, 0);
      p->up.sigclipnumset=1;
      break;

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
