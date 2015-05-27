/*********************************************************************
SubtractSky - Find and subtract the sky value from an image.
SubtractSky is part of GNU Astronomy Utilities (Gnuastro) package.

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
  SPACK_NAME" Finds the sky value over a grid on the input and subtracts "
  "it from the image to give a clear and uniform output. \n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   c e f g i j l m p r v w x y z
   A B C E F G I J O R W X Y Z

   Number keys free: >=508

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
      "kernel",
      'k',
      "STR",
      0,
      "Kernel image file name for convolution.",
      1
    },
    {
      "khdu",
      'U',
      "STR",
      0,
      "Kernel image header name for convolution.",
      1
    },


    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "checksky",
      502,
      0,
      0,
      "Store final sky and its STD in `_sky.fits' file.",
      2
    },
    {
      "checkskystd",
      505,
      0,
      0,
      "Include sky standard deviation in all checks too.",
      2
    },


    {
      0, 0, 0, 0,
      "Mesh grid:",
      3
    },
    {
      "meshsize",
      's',
      "INT",
      0,
      "Size of each mesh (tile) in the grid.",
      3
    },
    {
      "nch1",
      'a',
      "INT",
      0,
      "Number of channels along first FITS axis.",
      3
    },
    {
      "nch2",
      'b',
      "INT",
      0,
      "Number of channels along second FITS axis.",
      3
    },
    {
      "lastmeshfrac",
      'L',
      "INT",
      0,
      "Fraction of last mesh area to add new.",
      3
    },
    {
      "mirrordist",
      'd',
      "FLT",
      0,
      "Distance beyond mirror point. Multiple of std.",
      3
    },
    {
      "minmodeq",
      'Q',
      "FLT",
      0,
      "Minimum acceptable quantile for the mode.",
      3
    },
    {
      "numnearest",
      'n',
      "INT",
      0,
      "Number of nearest neighbors to interpolate.",
      3
    },
    {
      "smoothwidth",
      'T',
      "INT",
      0,
      "Width of smoothing kernel (odd number).",
      3
    },
    {
      "fullconvolution",
      506,
      0,
      0,
      "Ignore channels in imageconvolution.",
      3
    },
    {
      "fullinterpolation",
      503,
      0,
      0,
      "Ignore channels in interpolation.",
      3
    },
    {
      "fullsmooth",
      504,
      0,
      0,
      "Ignore channels in smoothing.",
      3
    },
    {
      "checkmesh",
      500,
      0,
      0,
      "Store mesh IDs in `_mesh.fits' file.",
      3
    },
    {
      "checkinterpolation",
      501,
      0,
      0,
      "Store mesh interpolation in `_interp.fits' file.",
      3
    },
    {
      "checkconvolution",
      507,
      0,
      0,
      "Store convolved image in `_conv.fits' file.",
      3
    },



    {
      0, 0, 0, 0,
      "Statistics:",
      4
    },
    {
      "sigclipmultip",
      'u',
      "FLT",
      0,
      "Multiple of standard deviation in sigma-clipping.",
      4
    },
    {
      "sigcliptolerance",
      't',
      "FLT",
      0,
      "Difference in STD tolerance to halt iteration.",
      4
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
  struct subtractskyparams *p = state->input;

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
      p->up.maskname=arg;
      p->up.masknameset=1;
      break;
    case 'H':
      errno=0;                  /* We want to free it in the end. */
      p->up.mhdu=malloc(strlen(arg)+1);
      if(p->up.mhdu==NULL)
        error(EXIT_FAILURE, errno, "Space for mask HDU.");
      strcpy(p->up.mhdu, arg);
      p->up.mhduset=1;
      break;
    case 'k':
      errno=0;                  /* We want to free it in the end. */
      p->up.kernelname=malloc(strlen(arg)+1);
      if(p->up.kernelname==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for kernel file name.",
              strlen(arg)+1);
      strcpy(p->up.kernelname, arg);
      p->up.kernelnameset=1;
      break;
    case 'U':
      errno=0;                  /* We want to free it in the end. */
      p->up.khdu=malloc(strlen(arg)+1);
      if(p->up.khdu==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for kernel HDU.",
              strlen(arg)+1);
      strcpy(p->up.khdu, arg);
      p->up.khduset=1;
      break;

    /* Output: */
    case 502:
      p->skyname="a";
      break;
    case 505:
      p->checkstd=1;
      break;

    /* Mesh grid: */
    case 's':
      sizetlzero(arg, &p->mp.meshsize, "meshsize", key, SPACK, NULL, 0);
      p->up.meshsizeset=1;
      break;
    case 'a':
      sizetlzero(arg, &p->mp.nch1, "nch1", key, SPACK, NULL, 0);
      p->up.nch1set=1;
      break;
    case 'b':
      sizetlzero(arg, &p->mp.nch2, "nch2", key, SPACK, NULL, 0);
      p->up.nch2set=1;
      break;
    case 'L':
      floatl0s1(arg, &p->mp.lastmeshfrac, "lastmeshfrac", key, SPACK,
                NULL, 0);
      p->up.lastmeshfracset=1;
      break;
    case 'd':
      floatl0(arg, &p->mp.mirrordist, "mirrordist", key, SPACK, NULL, 0);
      p->up.mirrordistset=1;
      break;
    case 'Q':
      floatl0s1(arg, &p->mp.minmodeq, "minmodeq", key, SPACK, NULL, 0);
      p->up.minmodeqset=1;
      break;
    case 'n':
      sizetlzero(arg, &p->mp.numnearest, "numnearest", key, SPACK, NULL, 0);
      p->up.numnearestset=1;
      break;
    case 'T':
      sizetpodd(arg, &p->mp.smoothwidth, "smoothwidth", key, SPACK, NULL, 0);
      p->up.smoothwidthset=1;
      break;
    case 506:
      p->mp.fullconvolution=1;
      break;
    case 503:
      p->mp.fullinterpolation=1;
      break;
    case 504:
      p->mp.fullsmooth=1;
      break;
    case 500:
      p->meshname="a";  /* Just a placeholder! It will be corrected later */
      break;
    case 501:
      p->interpname="a";
      break;
    case 507:
      p->convname="a";
      break;


    /* Statistics */
    case 'u':
      floatl0(arg, &p->sigclipmultip, "sigclipmultip", key, SPACK,
              NULL, 0);
      p->up.sigclipmultipset=1;
      break;
    case 't':
      floatl0s1(arg, &p->sigcliptolerance, "sigcliptolerance", key, SPACK,
              NULL, 0);
      p->up.sigcliptoleranceset=1;
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
