/*********************************************************************
MockGals - Create mock galaxies and stars in a noisy image.
MockGals is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef ARGS_H
#define ARGS_H


#include "fixedstringmacros.h"
#include "commonargs.h"







/**************************************************************/
/**************        argp.h definitions       ***************/
/**************************************************************/




/* Definition parameters for the argp: */
const char *argp_program_version=SPACK_STRING"\n"COPYRIGHT
  "\n\nWritten by Mohammad Akhlaghi";
const char *argp_program_bug_address=PACKAGE_BUGREPORT;
static char args_doc[] = "[PSFimage] Catalog";





const char doc[] =
  /* Before the list of options: */
  TOPHELPINFO
  SPACK_NAME" will create a FITS image containing any number of mock "
  "galaxies or stars based on the input catalog. The PSF can either be "
  "given as a FITS file or with Moffat or Gaussian parameters. All the "
  "profiles will be built from the center outwards. First by 10000 "
  "random points, then by integration and finally central pixel "
  "position. The tolerance level specifies when to switch to a less "
  "accurate method.\n"
  MOREHELPINFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/*
   Available letters (-V which is used by GNU is also removed):

   a c d e i j k m n p s u v
   A E F G H I J L M Q R T U W X Y Z
*/
static struct argp_option options[] =
  {
    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },
    {
      "onlypsf",
      'g',
      0,
      0,
      "Only make the PSF and abort.",
      -1
    },


    /* I am forced to put the PSF here, because we want the hdu
       option that comes from commonargs.h. */
    {
      0, 0, 0, 0,
      "PSF:",
      1
    },
    {
      "psffunction",
      'f',
      "STR",
      0,
      "PSF function: `moffat` or `gaussian`.",
      1
    },
    {
      "fwhm",
      'w',
      "FLT",
      0,
      "FWHM of PSF in units of pixels.",
      1
    },
    {
      "moffatbeta",
      'B',
      "FLT",
      0,
      "Moffat function's beta value.",
      1
    },
    {
      "psftrunc",
      'r',
      "FLT",
      0,
      "PSF truncation in units of FWHM/2.",
      1
    },



    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "naxis1",
      'x',
      "INT",
      0,
      "Number of pixels along first FITS axis.",
      2
    },
    {
      "naxis2",
      'y',
      "INT",
      0,
      "Number of pixels along second FITS axis.",
      2
    },
    {
      "noconv",
      'O',
      0,
      0,
      "Save image prior to convolution.",
      2
    },
    {
      "conv",
      'C',
      0,
      0,
      "Save image after convolution, prior to noise.",
      2
    },





    {
      0, 0, 0, 0,
      "Profiles and Noise:",
      3
    },
    {
      "truncation",
      't',
      "FLT",
      0,
      "Profile truncation distance, multiple of radius.",
      3
    },
    {
      "tolerance",
      'l',
      "FLT",
      0,
      "Tolerance to switch to less accurate method.",
      3
    },
    {
      "background",
      'b',
      "FLT",
      0,
      "Image background (amplitude of noise).",
      3
    },
    {
      "zeropoint",
      'z',
      "FLT",
      0,
      "Magnitude zero point.",
      3
    },





    {
      0, 0, 0, 0,
      "Profile catalog (column number, starting from zero):",
      4
    },
    {
      "fcol",
      500,
      "INT",
      0,
      "Function: Sersic (0), Point (3).",
      4
    },
    {
      "xcol",
      501,
      "INT",
      0,
      "Center along first FITS axis (horizontal).",
      4
    },
    {
      "ycol",
      502,
      "INT",
      0,
      "Center along second FITS axis (vertical).",
      4
    },
    {
      "rcol",
      503,
      "INT",
      0,
      "Effective radius in pixels.",
      4
    },
    {
      "ncol",
      504,
      "INT",
      0,
      "Sersic index.",
      4
    },
    {
      "pcol",
      505,
      "INT",
      0,
      "Position angle.",
      4
    },
    {
      "qcol",
      506,
      "INT",
      0,
      "Axis ratio.",
      4
    },
    {
      "mcol",
      507,
      "INT",
      0,
      "Magnitude.",
      4
    },





    {0}
  };



















/* Parse a single option: */
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  /* Save the arguments structure: */
  struct mockgalsparams *p = state->input;

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

    /* Operating modes: */
    case 'g':
      p->up.onlypsf=1;
      break;

    /* Input */
    case 'h':	    /* Because we redefined the explanation. */
      errno=0;
      p->cp.hdu=malloc(strlen(arg)+1);
      if(p->cp.hdu==NULL) error(EXIT_FAILURE, 0, NULL);
      strcpy(p->cp.hdu, arg);
      p->cp.hduset=1;
      break;

    /* Output: */
    case 'x':
      sizetlzero(arg, &p->s1, "naxis1", key, p->cp.spack, NULL, 0);
      p->up.naxis1set=1;
      break;
    case 'y':
      sizetlzero(arg, &p->s0, "naxis2", key, p->cp.spack, NULL, 0);
      p->up.naxis2set=1;
      break;
    case 'O':
      p->noconv=1;
      break;
    case 'C':
      p->conv=1;
      break;

    /* PSF: */
    case 'f':
      if(strcmp(arg, "moffat")==0)
	p->psffunction=1;
      else if(strcmp(arg, "gaussian")==0)
	p->psffunction=2;
      else
	argp_error(state, "The value of the `--psffunction` (`-f`) option "
		   "should be either `moffat` or `gaussian`.");
      p->up.psffunctionset=1;
      break;
    case 'w':
      floatl0(arg, &p->psf_p1, "fwhm", key, p->cp.spack, NULL, 0);
      p->up.fwhmset=1;
      break;
    case 'B':
      floatl0(arg, &p->psf_p2, "moffatbeta", key, p->cp.spack, NULL, 0);
      p->up.moffatbetaset=1;
      break;
    case 'r':
      floatl0(arg, &p->psf_t, "psftrunc", key, p->cp.spack, NULL, 0);
      p->up.psftruncset=1;
      break;

    /* Profiles and noise: */
    case 't':
      floatl0(arg, &p->truncation, "truncation", key, p->cp.spack, NULL, 0);
      p->up.truncationset=1;
      break;
    case 'l':
      floatl0(arg, &p->tolerance, "tolerance", key, p->cp.spack, NULL, 0);
      p->up.toleranceset=1;
      break;
    case 'b':
      floatl0(arg, &p->background, "background", key, p->cp.spack, NULL, 0);
      p->up.backgroundset=1;
      break;
    case 'z':
      floatl0(arg, &p->zeropoint, "zeropoint", key, p->cp.spack, NULL, 0);
      p->up.zeropointset=1;
      break;

   /* Catalog */
    case 500:
      sizetelzero(arg, &p->fcol, "fcol", ' ', p->cp.spack, NULL, 0);
      p->up.fcolset=1;
      break;
    case 501:
      sizetelzero(arg, &p->xcol, "xcol", ' ', p->cp.spack, NULL, 0);
      p->up.xcolset=1;
      break;
    case 502:
      sizetelzero(arg, &p->ycol, "ycol", ' ', p->cp.spack, NULL, 0);
      p->up.ycolset=1;
      break;
    case 503:
      sizetelzero(arg, &p->rcol, "rcol", ' ', p->cp.spack, NULL, 0);
      p->up.rcolset=1;
      break;
    case 504:
      sizetelzero(arg, &p->ncol, "ncol", ' ', p->cp.spack, NULL, 0);
      p->up.ncolset=1;
      break;
    case 505:
      sizetelzero(arg, &p->pcol, "pcol", ' ', p->cp.spack, NULL, 0);
      p->up.pcolset=1;
      break;
    case 506:
      sizetelzero(arg, &p->pcol, "qcol", ' ', p->cp.spack, NULL, 0);
      p->up.qcolset=1;
      break;
    case 507:
      sizetelzero(arg, &p->mcol, "mcol", ' ', p->cp.spack, NULL, 0);
      p->up.mcolset=1;
      break;





    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:

      /* See what type of input value it is and put it in. */
      if( nameisfits(arg) )
	{
	  if(p->up.psfname)
	    argp_error(state, "Only one input FITS image (the PSF) "
		       "should be input. You have given more.");
	  else
	    p->up.psfname=arg;
	}
      else
	{
	  if(p->up.catname)
	    argp_error(state, "Only one catalog file can be given.");
	  else
	    p->up.catname=arg;
	}
      break;






    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
	 && p->cp.printparams==0)
	{
	  if(state->arg_num==0)
	    argp_error(state, "No argument given!");
	  if(p->up.catname==NULL && p->up.onlypsf==0)
	    argp_error(state, "No catalog provided!");
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
