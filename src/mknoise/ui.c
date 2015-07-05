/*********************************************************************
MakeNoise - Add noise to a dataset.
MakeNoise is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include "timing.h"	/* Includes time.h and sys/time.h   */
#include "checkset.h"
#include "txtarrayvv.h"
#include "commonargs.h"
#include "configfiles.h"
#include "fitsarrayvv.h"

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
readconfig(char *filename, struct mknoiseparams *p)
{
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;
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
      STARTREADINGLINE;




      /* Inputs: */
      if(strcmp(name, "hdu")==0)
        allocatecopyset(value, &cp->hdu, &cp->hduset);

      else if(strcmp(name, "background")==0)
	{
	  if(up->backgroundset) continue;
          anydouble(value, &p->mbackground, value, key, SPACK, filename,
                    lineno);
	  up->backgroundset=1;
	}
      else if(strcmp(name, "zeropoint")==0)
	{
	  if(up->zeropointset) continue;
          anydouble(value, &p->zeropoint, value, key, SPACK, filename,
                    lineno);
	  up->zeropointset=1;
	}
      else if(strcmp(name, "stdadd")==0)
	{
	  if(up->stdaddset) continue;
          doublele0(value, &p->stdadd, value, key, SPACK, filename, lineno);
	  up->stdaddset=1;
	}



      /* Outputs */
      else if(strcmp(name, "output")==0)
        allocatecopyset(value, &cp->output, &cp->outputset);




      /* Operating modes: */
      else if(strcmp(name, "numthreads")==0)
	{
	  if(cp->numthreadsset) continue;
	  sizetlzero(value, &cp->numthreads, name, key, SPACK,
		     filename, lineno);
	  cp->numthreadsset=1;
	}


      else
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct mknoiseparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    PRINTSTINGMAYBEWITHSPACE("hdu", cp->hdu);
  if(up->backgroundset)
    fprintf(fp, CONF_SHOWFMT"%f\n", "background", p->mbackground);
  if(up->zeropointset)
    fprintf(fp, CONF_SHOWFMT"%f\n", "zeropoint", p->zeropoint);
  if(up->stdaddset)
    fprintf(fp, CONF_SHOWFMT"%f\n", "stdadd", p->stdadd);


  fprintf(fp, "\n# Output parameters:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);



  fprintf(fp, "\n# Operating mode:\n");
  /* Number of threads doesn't need to be checked, it is set by
     default */
  fprintf(fp, CONF_SHOWFMT"%lu\n", "numthreads", cp->numthreads);
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct mknoiseparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    REPORT_NOTSET("hdu");

  if(up->backgroundset==0)
    REPORT_NOTSET("background");
  if(up->zeropointset==0)
    REPORT_NOTSET("zeropoint");
  if(up->stdaddset==0)
    REPORT_NOTSET("stdadd");


  END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct mknoiseparams *p)
{
  /* Set the output name: */
  if(p->cp.output)
    checkremovefile(p->cp.output, p->cp.dontdelete);
  else
    automaticoutput(p->up.inputname, "_noised.fits", p->cp.removedirinfo,
                    p->cp.dontdelete, &p->cp.output);

  /* Convert the background value from magnitudes to flux. Note that
     magnitudes are actually calculated from the ratio of brightness,
     not flux. But in the context of MakeNoise where everything is
     done on pixels independently, brightness and flux are the same
     (flux is multiplied by the area of one pixel (=1) to give
     brightness).*/
  p->background=pow(10, (p->zeropoint-p->mbackground)/2.5f);
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
preparearrays(struct mknoiseparams *p)
{
  void *array;

  /* Read in the input image: */
  p->numblank=fitsimgtoarray(p->up.inputname, p->cp.hdu, &p->inputbitpix,
                        &array, &p->is0, &p->is1);
  if(p->inputbitpix==DOUBLE_IMG)
    p->input=array;
  else
    {
      changetype(array, p->inputbitpix, p->is0*p->is1, p->numblank,
                 (void **)&p->input, DOUBLE_IMG);
      free(array);
    }
  readfitswcs(p->up.inputname, p->cp.hdu, &p->nwcs, &p->wcs);

  /* Allocate the random number generator: */
  gsl_rng_env_setup();
  p->rng=gsl_rng_alloc(gsl_rng_default);
  if(p->envseed==0)
    gsl_rng_set(p->rng, timebasedrngseed());
  p->rng_seed=gsl_rng_default_seed;
  strcpy(p->rng_type, gsl_rng_name(p->rng));
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct mknoiseparams *p)
{
  char message[VERBMSGLENGTH_V];
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "Parsing arguments");

  /* Add the user default values and save them if asked. */
  CHECKSETCONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    REPORT_PARAMETERS_SET;

  /* Make the array of input images. */
  preparearrays(p);

  /* Do a sanity check. */
  sanitycheck(p);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      sprintf(message, "Random number generator type: %s",
              gsl_rng_name(p->rng));
      reporttiming(NULL, message, 1);
      if(p->envseed)
        {
          sprintf(message, "Random number generator seed: %lu",
                  gsl_rng_default_seed);
          reporttiming(NULL, message, 1);
        }
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct mknoiseparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->input);
  free(p->cp.hdu);
  free(p->cp.output);

  /* The world coordinate system: */
  if(p->wcs)
    wcsvfree(&p->nwcs, &p->wcs);

  /* Free the random number generator: */
  gsl_rng_free(p->rng);

  /* Print the final message. */
  reporttiming(t1, SPACK_NAME" finished in: ", 0);
}
