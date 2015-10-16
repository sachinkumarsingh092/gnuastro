/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include "ui.h"
#include "commonargs.h"
#include "fixedstringmacros.h"


/* Definition parameters for the argp: */
const char *argp_program_version=SPACK_STRING"\n"COPYRIGHT
  "\n\nWritten by Mohammad Akhlaghi";
const char *argp_program_bug_address=PACKAGE_BUGREPORT;
static char args_doc[] = "InputFile";





const char doc[] =
  /* Before the list of options: */
  TOPHELPINFO
  SPACK_NAME" will convolve an input image with a given spatial kernel "
  "(image) in the spatial domain (no edge effects) or frequency domain. "
  "The latter suffers from edge effects, but can be much faster.\n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Free letters for options:

   c d e g i j l n r t u w x y z
   A B C E F G H I J O Q R T W X Y Z

   Free numbers: >=504
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
      "Name of kernel for convolution.",
      1
    },
    {
      "khdu",
      'U',
      "STR",
      0,
      "HDU of kernel file.",
      1
    },
    {
      "nokernelflip",
      500,
      0,
      0,
      "Do not flip the kernel image.",
      1
    },
    {
      "nokernelnorm",
      501,
      0,
      0,
      "Do not normalize the kernel image.",
      1
    },




    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "viewfreqsteps",
      'v',
      0,
      0,
      "View the steps in the frequency domain.",
      2
    },



    {
      0, 0, 0, 0,
      "Mesh grid (only for spatial domain):",
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
      "checkmesh",
      503,
      0,
      0,
      "Store mesh IDs in `_mesh.fits' file.",
      3
    },
    {
      "fullconvolution",
      502,
      0,
      0,
      "Ignore channels in imageconvolution.",
      3
    },




    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },
    {
      "spatial",
      'p',
      0,
      0,
      "Spatial domain convolution.",
      -1
    },
    {
      "frequency",
      'f',
      0,
      0,
      "Frequency domain convolution.",
      -1
    },
    {
      "makekernel",
      'm',
      "INT",
      0,
      "Make 2*INT kernel to create input image.",
      -1
    },

    {0}
  };




















/* Parse a single option: */
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{

  /* Save the arguments structure: */
  struct convolveparams *p = state->input;

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

    /* Inputs: */
    case 'M':
      allocatecopyset(arg, &p->up.maskname, &p->up.masknameset);
      break;
    case 'H':
      allocatecopyset(arg, &p->up.mhdu, &p->up.mhduset);
      break;
    case 'k':
      allocatecopyset(arg, &p->up.kernelname, &p->up.kernelnameset);
      break;
    case 'U':
      allocatecopyset(arg, &p->up.khdu, &p->up.khduset);
      break;
    case 500:
      p->kernelflip=0;
      break;
    case 501:
      p->kernelnorm=0;
      break;


    /* Output: */


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
    case 503:
      p->meshname="a";
      break;
    case 502:
      p->mp.fullconvolution=1;
      p->up.fullconvolutionset=1;
      break;


   /* Operating mode: */
    case 'p':
      if(p->up.frequencyset)
	argp_error(state, "Only one of spatial or frequency domain "
                   "convolution modes may be chosen.");
      p->spatial=1;
      p->frequency=0;
      p->up.spatialset=p->up.frequencyset=1;
      break;
    case 'f':
      if(p->up.spatialset)
	argp_error(state, "Only one of spatial or frequency domain "
                   "convolution modes may be chosen.");
      p->spatial=0;
      p->frequency=1;
      p->up.spatialset=p->up.frequencyset=1;
      break;
    case 'v':
      p->viewfreqsteps=1;
      break;
    case 'm':
      intelzero(arg, &p->makekernel, "makekernel", key, SPACK, NULL, 0);
      p->up.makekernelset=1;
      break;




    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:
      if(p->up.inputname)
        argp_error(state, "Only one input file (argument) is required.");
      p->up.inputname=arg;
      break;






    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
	 && p->cp.printparams==0)
	{
	  if(state->arg_num==0)
	    argp_error(state, "No argument given!");
	  if(p->up.inputname==NULL)
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
