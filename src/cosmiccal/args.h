/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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
static char args_doc[] = "";





const char doc[] =
  /* Before the list of options: */
  GAL_STRINGS_TOP_HELP_INFO
  SPACK_NAME" will produce cosmological calculations.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   b c d e f g i j k m n p s t u w x y
   A B C E F G I J L M O Q R T U W X Y Z

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
      "redshift",
      'z',
      "FLT",
      0,
      "Redshift of interest.",
      1
    },
    {
      "H0",
      'H',
      "FLT",
      0,
      "Current expansion rate (Hubble constant).",
      1
    },
    {
      "olambda",
      'l',
      "FLT",
      0,
      "Current cosmological cst. dens. per crit. dens.",
      1
    },
    {
      "omatter",
      'm',
      "FLT",
      0,
      "Current matter density per critical density.",
      1
    },
    {
      "oradiation",
      'r',
      "FLT",
      0,
      "Current radiation density per critical density.",
      1
    },



    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "onlyvolume",
      'v',
      0,
      0,
      "Only print comoving volume in Mpc^3",
      2
    },
    {
      "onlyabsmagconv",
      'a',
      0,
      0,
      "Only print conversion to absolute magnitude.",
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
  struct cosmiccalparams *p = state->input;

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


    /* Input: */
    case 'z':
      gal_checkset_double_el_0(arg, &p->redshift, "redshift", key,
                               SPACK, NULL, 0);
      p->up.redshiftset=1;
      break;
    case 'H':
      gal_checkset_double_el_0(arg, &p->H0, "H0", key, SPACK, NULL, 0);
      p->up.H0set=1;
      break;
    case 'l':
      gal_checkset_double_el_0(arg, &p->olambda, "olambda", key,
                               SPACK, NULL, 0);
      p->up.olambdaset=1;
      break;
    case 'm':
      gal_checkset_double_el_0(arg, &p->omatter, "omatter", key,
                               SPACK, NULL, 0);
      p->up.omatterset=1;
      break;
    case 'r':
      gal_checkset_double_el_0(arg, &p->oradiation, "oradiation",
                               key, SPACK, NULL, 0);
      p->up.oradiationset=1;
      break;


    /* Output: */
    case 'v':
      p->onlyvolume=1;
      p->up.onlyvolumeset=1;
      break;
    case 'a':
      p->onlyabsmagconv=1;
      p->up.onlyabsmagconvset=1;
      break;


    /* Operating modes: */


    /* Read the non-option arguments: */
    case ARGP_KEY_ARG:
      argp_error(state, SPACK_NAME" only takes options as "
                 "input, currently no arguments are supported");
      break;





    /* The command line options and arguments are finished. */
    case ARGP_KEY_END:
      /* Currently there are no arguments, or input files, so there is
         no need for such checks.

      if(p->cp.setdirconf==0 && p->cp.setusrconf==0
	 && p->cp.printparams==0)
	{
	  if(state->arg_num==0)
	    argp_error(state, "no argument given");
	  if(p->up.inputname==NULL)
	    argp_error(state, "no input FITS image(s) provided");
	}
      */
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
