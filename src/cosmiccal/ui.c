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
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include <nproc.h>              /* From Gnulib.                     */
#include <gsl/gsl_const_mksa.h>

#include <gnuastro/timing.h>    /* Includes time.h and sys/time.h   */
#include <gnuastro/checkset.h>
#include <gnuastro/txtarrayvv.h>
#include <gnuastro/commonargs.h>
#include <gnuastro/configfiles.h>
#include <gnuastro/fitsarrayvv.h>

#include "main.h"

#include "ui.h"		        /* Needs main.h                   */
#include "args.h"	        /* Needs main.h, includes argp.h. */


/* Set the file names of the places where the default parameters are
   put. */
#define CONFIG_FILE SPACK CONF_POSTFIX
#define SYSCONFIG_FILE SYSCONFIG_DIR "/" CONFIG_FILE
#define USERCONFIG_FILEEND USERCONFIG_DIR CONFIG_FILE
#define CURDIRCONFIG_FILE CURDIRCONFIG_DIR CONFIG_FILE










/**************************************************************/
/**************       Options and parameters    ***************/
/**************************************************************/
void
readconfig(char *filename, struct cosmiccalparams *p)
{
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;
  char key='a';	/* Not used, just a place holder. */

  /* When the file doesn't exist or can't be opened, it is ignored. It
     might be intentional, so there is no error. If a parameter is
     missing, it will be reported after all defaults are read. */
  fp=fopen(filename, "r");
  if (fp==NULL) return;


  /* Allocate some space for `line` with `len` elements so it can
     easily be freed later on. The value of `len` is arbitarary at
     this point, during the run, getline will change it along with the
     pointer to line. */
  errno=0;
  line=malloc(len*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes in readdefaults",
	  len * sizeof *line);

  /* Read the tokens in the file:  */
  while(getline(&line, &len, fp) != -1)
    {
      /* Prepare the "name" and "value" strings, also set lineno. */
      GAL_CONFIGFILES_START_READING_LINE;




      /* Inputs: */
      if(strcmp(name, "redshift")==0)
	{
	  if(up->redshiftset) continue;
          gal_checkset_double_el_0(value, &p->redshift, name, key,
                                      SPACK, filename, lineno);
	  up->redshiftset=1;
	}
      else if(strcmp(name, "H0")==0)
	{
	  if(up->H0set) continue;
          gal_checkset_double_el_0(value, &p->H0, name, key, SPACK,
                                      filename, lineno);
	  up->H0set=1;
	}
      else if(strcmp(name, "olambda")==0)
	{
	  if(up->olambdaset) continue;
          gal_checkset_double_el_0(value, &p->olambda, name, key,
                                      SPACK, filename, lineno);
	  up->olambdaset=1;
	}
      else if(strcmp(name, "omatter")==0)
	{
	  if(up->omatterset) continue;
          gal_checkset_double_el_0(value, &p->omatter, name, key,
                                      SPACK, filename, lineno);
	  up->omatterset=1;
	}
      else if(strcmp(name, "oradiation")==0)
	{
	  if(up->oradiationset) continue;
          gal_checkset_double_el_0(value, &p->oradiation, name, key,
                                   SPACK, filename, lineno);
	  up->oradiationset=1;
	}



      /* Outputs */
      else if(strcmp(name, "onlyvolume")==0)
	{
	  if(up->onlyvolumeset) continue;
          gal_checkset_int_zero_or_one(value, &p->onlyvolume, name, key,
                                       SPACK, filename, lineno);
	  up->onlyvolumeset=1;
	}
      else if(strcmp(name, "onlyabsmagconv")==0)
	{
	  if(up->onlyabsmagconvset) continue;
          gal_checkset_int_zero_or_one(value, &p->onlyabsmagconv, name,
                                       key, SPACK, filename, lineno);
	  up->onlyabsmagconvset=1;
	}



      /* Operating modes: */
      /* Read options common to all programs */
      GAL_CONFIGFILES_READ_COMMONOPTIONS_FROM_CONF


      else
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct cosmiccalparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;


  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input:\n");
  if(up->redshiftset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "redshift", p->redshift);
  if(up->H0set)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "H0", p->H0);


  fprintf(fp, "\n# Current densities per current critical density:\n");
  if(up->olambdaset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "olambda", p->olambda);
  if(up->omatterset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "omatter", p->omatter);
  if(up->oradiationset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "oradiation", p->oradiation);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct cosmiccalparams *p)
{
  struct uiparams *up=&p->up;
  /*struct gal_commonparams *cp=&p->cp;*/

  int intro=0;
  if(up->redshiftset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("redshift");
  if(up->H0set==0)
    GAL_CONFIGFILES_REPORT_NOTSET("H0");
  if(up->olambdaset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("olambda");
  if(up->omatterset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("omatter");
  if(up->oradiationset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("oradiation");


  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/************            Sanity check             *************/
/**************************************************************/
void
sanitycheck(struct cosmiccalparams *p)
{
  int check=p->onlyvolume+p->onlyabsmagconv;

  /* If only one of the single output options are called, then check
     should be 1, if none are called, then it should be zero. However,
     if more than one is called, check will be larger than one. So in
     this case, report an error. */
  if(check>1)
    error(EXIT_FAILURE, 0, "only a single option starting with `--only' "
          "can be called");
}




















/**************************************************************/
/************             Preparations            *************/
/**************************************************************/
void
preparations(struct cosmiccalparams *p)
{
  /* Speed of light: */
  p->c=GSL_CONST_MKSA_SPEED_OF_LIGHT;

  /* The curvature fractional density: */
  p->ocurv=1-(p->olambda+p->omatter+p->oradiation);

  /* Convert H0 from km/sec/Mpc to 1/sec: */
  p->H0s=p->H0/1000/GSL_CONST_MKSA_PARSEC;
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct cosmiccalparams *p)
{
  struct gal_commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = num_processors(NPROC_CURRENT);
  cp->removedirinfo = 1;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Do a sanity check */
  sanitycheck(p);

  /* Make the preparations */
  preparations(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    GAL_CONFIGFILES_REPORT_PARAMETERS_SET;
}
