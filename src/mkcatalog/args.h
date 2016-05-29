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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef ARGS_H
#define ARGS_H

#include <argp.h>

#include <gnuastro/commonargs.h>
#include <gnuastro/linkedlist.h>
#include <gnuastro/fixedstringmacros.h>










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
  SPACK_NAME" will create a catalog from an input, labeled, and noise "
  "identification images.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   e f g k l u v w
   F G J L Q R U W X Y Z

   Number keys used: <=533

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
      "Image specifying object labels.",
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
      "Image specifying clump labels.",
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
      "skyfilename",
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
      "stdfilename",
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
      "zeropoint",
      'z',
      "FLT",
      0,
      "Image zeropoint magnitude.",
      1
    },
    {
      "skysubtracted",
      'E',
      0,
      0,
      "Input is already sky subtracted (for S/N).",
      1
    },
    {
      "threshold",
      'T',
      "FLT",
      0,
      "Only values larger than this multiple of STD.",
      1
    },



    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "nsigmag",
      521,
      "FLT",
      0,
      "Multiple of Sky STD to report magnitude of.",
      2
    },
    {
      "intwidth",
      516,
      "INT",
      0,
      "Width of integer columns.",
      2
    },
    {
      "floatwidth",
      517,
      "INT",
      0,
      "Width of floating point columns.",
      2
    },
    {
      "accuwidth",
      518,
      "INT",
      0,
      "Width of more accurate floating point columns.",
      2
    },
    {
      "floatprecision",
      519,
      "INT",
      0,
      "Precision of floating point columns.",
      2
    },
    {
      "accuprecision",
      520,
      "INT",
      0,
      "Precision of more accurate floating pnt. cols.",
      2
    },



    {
      0, 0, 0, 0,
      "Catalog columns:",
      3
    },
    {
      "id",
      'i',
      0,
      0,
      "Overall ID of this object or clump.",
      3
    },
    {
      "hostobjid",
      'j',
      0,
      0,
      "ID of object hosting this clump.",
      3
    },
    {
      "idinhostobj",
      'I',
      0,
      0,
      "ID of clump in host object.",
      3
    },
    {
      "numclumps",
      'C',
      0,
      0,
      "Number of clumps in this object.",
      3
    },
    {
      "area",
      'a',
      0,
      0,
      "Number of pixels.",
      3
    },
    {
      "clumpsarea",
      513,
      0,
      0,
      "Area of clumps in an object.",
      3
    },
    {
      "x",
      'x',
      0,
      0,
      "All obj. flux weighted center (first FITS axis).",
      3
    },
    {
      "y",
      'y',
      0,
      0,
      "All obj. flux weighted center (second FITS axis).",
      3
    },
    {
      "geox",
      522,
      0,
      0,
      "All obj. geometric center (first FITS axis).",
      3
    },
    {
      "geoy",
      523,
      0,
      0,
      "All obj. geometric center (second FITS axis).",
      3
    },
    {
      "clumpsx",
      507,
      0,
      0,
      "Clumps flux weighted center (first FITS axis).",
      3
    },
    {
      "clumpsy",
      508,
      0,
      0,
      "Clumps flux weighted center (second FITS axis).",
      3
    },
    {
      "clumpsgeox",
      524,
      0,
      0,
      "Clumps geometric center (first FITS axis).",
      3
    },
    {
      "clumpsgeoy",
      525,
      0,
      0,
      "Clumps geometric center (second FITS axis).",
      3
    },
    {
      "ra",
      'r',
      0,
      0,
      "All object flux weighted center right ascension.",
      3
    },
    {
      "dec",
      'd',
      0,
      0,
      "All object flux weighted center declination.",
      3
    },
    {
      "geora",
      526,
      0,
      0,
      "All object geometric center right ascension.",
      3
    },
    {
      "geodec",
      527,
      0,
      0,
      "All object geometric center declination.",
      3
    },
    {
      "clumpsra",
      509,
      0,
      0,
      "Clumps flux weighted center right ascension.",
      3
    },
    {
      "clumpsdec",
      510,
      0,
      0,
      "Clumps flux weighted center declination.",
      3
    },
    {
      "clumpsgeora",
      528,
      0,
      0,
      "Clumps geometric center right ascension.",
      3
    },
    {
      "clumpsgeodec",
      529,
      0,
      0,
      "Clumps geometric center declination.",
      3
    },
    {
      "brightness",
      'b',
      0,
      0,
      "Brightness (sum of pixel values).",
      3
    },
    {
      "clumpbrightness",
      511,
      0,
      0,
      "Brightness in clumps of an object.",
      3
    },
    {
      "noriverbrightness",
      533,
      0,
      0,
      "Sky (not river) subtracted clump brightness.",
      3
    },
    {
      "magnitude",
      'm',
      0,
      0,
      "Total magnitude.",
      3
    },
    {
      "clumpsmagnitude",
      512,
      0,
      0,
      "Total magnitude of clumps in this object.",
      3
    },
    {
      "riverave",
      514,
      0,
      0,
      "Average river value surrounding this clump.",
      3
    },
    {
      "rivernum",
      515,
      0,
      0,
      "Number of river pixels surrounding this clump.",
      3
    },
    {
      "sn",
      'n',
      0,
      0,
      "Signal to noise ratio column.",
      3
    },
    {
      "sky",
      505,
      0,
      0,
      "Sky value.",
      3
    },
    {
      "std",
      506,
      0,
      0,
      "Sky standard deviation.",
      3
    },
    {
      "semimajor",
      'A',
      0,
      0,
      "Flux weighted Semi-major axis.",
      3
    },
    {
      "semiminor",
      'B',
      0,
      0,
      "Flux weighted Semi-minor axis.",
      3
    },
    {
      "positionangle",
      'p',
      0,
      0,
      "Flux weighted Position angle.",
      3
    },
    {
      "geosemimajor",
      530,
      0,
      0,
      "Flux weighted Semi-major axis.",
      3
    },
    {
      "geosemiminor",
      531,
      0,
      0,
      "Flux weighted Semi-minor axis.",
      3
    },
    {
      "geopositionangle",
      532,
      0,
      0,
      "Flux weighted Position angle.",
      3
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
      gal_checkset_allocate_copy_set(arg, &p->up.maskname, &p->up.masknameset);
      break;
    case 'H':
      gal_checkset_allocate_copy_set(arg, &p->up.mhdu, &p->up.mhduset);
      break;
    case 'O':
      gal_checkset_allocate_copy_set(arg, &p->up.objlabsname,
                                     &p->up.objlabsnameset);
      break;
    case 501:
      gal_checkset_allocate_copy_set(arg, &p->up.objhdu, &p->up.objhduset);
      break;
    case 'c':
      gal_checkset_allocate_copy_set(arg, &p->up.clumplabsname,
                                     &p->up.clumplabsnameset);
      break;
    case 502:
      gal_checkset_allocate_copy_set(arg, &p->up.clumphdu, &p->up.clumphduset);
      break;
    case 's':
      gal_checkset_allocate_copy_set(arg, &p->up.skyname, &p->up.skynameset);
      break;
    case 503:
      gal_checkset_allocate_copy_set(arg, &p->up.skyhdu, &p->up.skyhduset);
      break;
    case 't':
      gal_checkset_allocate_copy_set(arg, &p->up.stdname, &p->up.stdnameset);
      break;
    case 504:
      gal_checkset_allocate_copy_set(arg, &p->up.stdhdu, &p->up.stdhduset);
      break;
    case 'z':
      gal_checkset_any_float(arg, &p->zeropoint, "zeropoint", key, SPACK,
                             NULL, 0);
      p->up.zeropointset=1;
      break;
    case 'E':
      p->skysubtracted=1;
      p->up.skysubtractedset=1;
      break;
    case 'T':
      gal_checkset_any_double(arg, &p->threshold, "threshold", key, SPACK,
                              NULL, 0);
      p->up.thresholdset=1;
      break;


    /* Output: */
    case 521:
      gal_checkset_any_double(arg, &p->nsigmag, "nsigmag", key, SPACK,
                              NULL, 0);
      p->up.nsigmagset=1;
      break;
    case 516:
      gal_checkset_int_l_zero(arg, &p->intwidth, "intwidth", key, SPACK,
                              NULL, 0);
      p->up.intwidthset=1;
      break;
    case 517:
      gal_checkset_int_l_zero(arg, &p->floatwidth, "floatwidth", key, SPACK,
                              NULL, 0);
      p->up.floatwidthset=1;
      break;
    case 518:
      gal_checkset_int_l_zero(arg, &p->accuwidth, "accuwidth", key, SPACK,
                              NULL, 0);
      p->up.floatwidthset=1;
      break;
    case 519:
      gal_checkset_int_l_zero(arg, &p->floatprecision, "flatprecision",
                              key, SPACK, NULL, 0);
      p->up.floatprecisionset=1;
      break;
    case 520:
      gal_checkset_int_l_zero(arg, &p->accuprecision, "accuprecision",
                              key, SPACK, NULL, 0);
      p->up.accuprecisionset=1;
      break;


    /* Catalog columns: */
    case 'i':
      add_to_sll(&p->allcolsll, CATID);
      p->up.idset=1;
      break;
    case 'j':
      add_to_sll(&p->allcolsll, CATHOSTOBJID);
      p->up.hostobjidset=1;
      break;
    case 'I':
      add_to_sll(&p->allcolsll, CATIDINHOSTOBJ);
      p->up.idinhostobjset=1;
      break;
    case 'C':
      add_to_sll(&p->allcolsll, CATNUMCLUMPS);
      p->up.numclumpsset=1;
      break;
    case 'a':
      add_to_sll(&p->allcolsll, CATAREA);
      p->up.areaset=1;
      break;
    case 513:
      add_to_sll(&p->allcolsll, CATCLUMPSAREA);
      p->up.clumpsareaset=1;
      break;
    case 'x':
      add_to_sll(&p->allcolsll, CATX);
      p->up.xset=1;
      break;
    case 'y':
      add_to_sll(&p->allcolsll, CATY);
      p->up.yset=1;
      break;
    case 522:
      add_to_sll(&p->allcolsll, CATGEOX);
      p->up.geoxset=1;
      break;
    case 523:
      add_to_sll(&p->allcolsll, CATGEOY);
      p->up.geoyset=1;
      break;
    case 507:
      add_to_sll(&p->allcolsll, CATCLUMPSX);
      p->up.clumpsxset=1;
      break;
    case 508:
      add_to_sll(&p->allcolsll, CATCLUMPSY);
      p->up.clumpsyset=1;
      break;
    case 524:
      add_to_sll(&p->allcolsll, CATCLUMPSGEOX);
      p->up.clumpsgeoxset=1;
      break;
    case 525:
      add_to_sll(&p->allcolsll, CATCLUMPSGEOY);
      p->up.clumpsgeoyset=1;
      break;
    case 'r':
      add_to_sll(&p->allcolsll, CATRA);
      p->up.raset=1;
      break;
    case 'd':
      add_to_sll(&p->allcolsll, CATDEC);
      p->up.decset=1;
      break;
    case 526:
      add_to_sll(&p->allcolsll, CATGEORA);
      p->up.georaset=1;
      break;
    case 527:
      add_to_sll(&p->allcolsll, CATGEODEC);
      p->up.geodecset=1;
      break;
    case 509:
      add_to_sll(&p->allcolsll, CATCLUMPSRA);
      p->up.clumpsraset=1;
      break;
    case 510:
      add_to_sll(&p->allcolsll, CATCLUMPSDEC);
      p->up.clumpsdecset=1;
      break;
    case 528:
      add_to_sll(&p->allcolsll, CATCLUMPSGEORA);
      p->up.clumpsgeoraset=1;
      break;
    case 529:
      add_to_sll(&p->allcolsll, CATCLUMPSGEODEC);
      p->up.clumpsgeodecset=1;
      break;
    case 'b':
      add_to_sll(&p->allcolsll, CATBRIGHTNESS);
      p->up.brightnessset=1;
      break;
    case 511:
      add_to_sll(&p->allcolsll, CATCLUMPSBRIGHTNESS);
      p->up.clumpsbrightnessset=1;
      break;
    case 533:
      add_to_sll(&p->allcolsll, CATNORIVERBRIGHTNESS);
      p->up.noriverbrightnessset=1;
      break;
    case 'm':
      add_to_sll(&p->allcolsll, CATMAGNITUDE);
      p->up.magnitudeset=1;
      break;
    case 512:
      add_to_sll(&p->allcolsll, CATCLUMPSMAGNITUDE);
      p->up.clumpsmagnitudeset=1;
      break;
    case 514:
      add_to_sll(&p->allcolsll, CATRIVERAVE);
      p->up.riveraveset=1;
      break;
    case 515:
      add_to_sll(&p->allcolsll, CATRIVERNUM);
      p->up.rivernumset=1;
      break;
    case 'n':
      add_to_sll(&p->allcolsll, CATSN);
      p->up.snset=1;
      break;
    case 505:
      add_to_sll(&p->allcolsll, CATSKY);
      p->up.skyset=1;
      break;
    case 506:
      add_to_sll(&p->allcolsll, CATSTD);
      p->up.stdset=1;
      break;
    case 'A':
      add_to_sll(&p->allcolsll, CATSEMIMAJOR);
      p->up.semimajorset=1;
      break;
    case 'B':
      add_to_sll(&p->allcolsll, CATSEMIMINOR);
      p->up.semiminorset=1;
      break;
    case 'p':
      add_to_sll(&p->allcolsll, CATPOSITIONANGLE);
      p->up.positionangleset=1;
      break;
    case 530:
      add_to_sll(&p->allcolsll, CATGEOSEMIMAJOR);
      p->up.geosemimajorset=1;
      break;
    case 531:
      add_to_sll(&p->allcolsll, CATGEOSEMIMINOR);
      p->up.geosemiminorset=1;
      break;
    case 532:
      add_to_sll(&p->allcolsll, CATGEOPOSITIONANGLE);
      p->up.geopositionangleset=1;
      break;




    /* Operating modes: */


    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:

      /* See what type of input value it is and put it in. */
      if( gal_fitsarray_name_is_fits(arg) )
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
