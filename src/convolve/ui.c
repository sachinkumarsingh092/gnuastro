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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include "timing.h"	        /* Includes time.h and sys/time.h */
#include "checkset.h"
#include "configfiles.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "ui.h"			/* Needs main.h.                  */
#include "args.h"		/* Needs main.h, includes argp.h. */


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
readconfig(char *filename, struct convolveparams *p)
{
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;
  int zeroorone, spatialset=0, frequencyset=0;
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
	{
	  if(cp->hduset) continue;
	  errno=0;
	  cp->hdu=malloc(strlen(value)+1);
	  if(cp->hdu==NULL)
	    error(EXIT_FAILURE, errno, "Space for hdu");
	  strcpy(cp->hdu, value);
	  cp->hduset=1;
	}
      else if(strcmp(name, "kernel")==0)
	{
	  if(up->kernelnameset) continue;
	  errno=0;
	  up->kernelname=malloc(strlen(value)+1);
	  if(up->kernelname==NULL)
	    error(EXIT_FAILURE, errno, "Space for kernelname");
	  strcpy(up->kernelname, value);
	  up->kernelnameset=1;
	}
      else if(strcmp(name, "khdu")==0)
	{
	  if(up->khduset) continue;
	  errno=0;
	  up->khdu=malloc(strlen(value)+1);
	  if(up->khdu==NULL)
	    error(EXIT_FAILURE, errno, "Space for khdu");
	  strcpy(up->khdu, value);
	  up->khduset=1;
	}



      /* Outputs: */
      else if(strcmp(name, "output")==0)
	{
	  if(cp->outputset) continue;
	  errno=0;
	  cp->output=malloc(strlen(value)+1);
	  if(cp->output==NULL)
	    error(EXIT_FAILURE, errno, "Space for output");
	  strcpy(cp->output, value);
	  cp->outputset=1;
	}





      /* Operating modes: */
      else if(strcmp(name, "spatial")==0)
	{
	  intzeroorone(value, &zeroorone, name, key, SPACK,
		       filename, lineno);
	  if(zeroorone)
	    {
	      spatialset=1;
	      if(frequencyset)
		error_at_line(EXIT_FAILURE, 0, filename, lineno,
			      "Spatial and frequency modes cannot be called "
			      "together. It is ambiguous.");
	      if(up->spatialset==0)
		{
		  p->spatial=1;
		  p->frequency=0;
		  up->spatialset=up->frequencyset=1;
		}
	    }
	}
      else if(strcmp(name, "frequency")==0)
	{
	  intzeroorone(value, &zeroorone, name, key, SPACK,
		       filename, lineno);
	  if(zeroorone)
	    {
	      frequencyset=1;
	      if(spatialset)
		error_at_line(EXIT_FAILURE, 0, filename, lineno,
			      "Spatial and frequency modes cannot be called "
			      "together. It is ambiguous.");
	      if(up->frequencyset==0)
		{
		  p->spatial=0;
		  p->frequency=1;
		  up->spatialset=up->frequencyset=1;
		}
	    }
	}
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
printvalues(FILE *fp, struct convolveparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  fprintf(fp, "\n# Input:\n");
  if(cp->hduset)
    {
      if(stringhasspace(cp->hdu))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "hdu", cp->hdu);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "hdu", cp->hdu);
    }
  if(up->kernelnameset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "kernel", up->kernelname);
  if(up->khduset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "khdu", up->khdu);



  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);



  fprintf(fp, "\n# Operating modes:\n");
  if(up->spatialset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "spatial", p->spatial);
  if(up->frequencyset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "frequency", p->frequency);
  /* Number of threads doesn't need to be checked, it is set by
     default */
  fprintf(fp, CONF_SHOWFMT"%lu\n", "numthreads", p->cp.numthreads);
}





void
checkifset(struct convolveparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;

  if(up->spatialset==0 && up->frequencyset==0)
    REPORT_NOTSET("spatial or frequency");
  if(cp->hduset==0)
    REPORT_NOTSET("hdu");
  if(up->kernelnameset==0)
    REPORT_NOTSET("kernel");
  if(up->khduset==0)
    REPORT_NOTSET("khdu");

  END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/************         Prepare the arrays          *************/
/**************************************************************/
void
preparearrays(struct convolveparams *p)
{
  int bitpix;
  void *array;
  double sum=0.0f;
  size_t numnul, i, size;
  struct uiparams *up=&p->up;
  float *f, *fp, tmp, *kernel;
  struct commonparams *cp=&p->cp;

  /* First read the input image: */
  numnul=fitsimgtoarray(up->inputname, cp->hdu, &bitpix, &array,
                        &p->is0, &p->is1);
  readfitswcs(up->inputname, cp->hdu, &p->nwcs, &p->wcs);
  if(bitpix==FLOAT_IMG)
    p->input=array;
  else
    {
      changetype(array, bitpix, p->is0*p->is1, numnul, (void **)(&p->input),
                 FLOAT_IMG);
      free(array);
    }

  /* Now read the kernel: */
  numnul=fitsimgtoarray(up->kernelname, up->khdu, &bitpix, &array,
                        &p->ks0, &p->ks1);

  /* Convert the kernel to a float image: */
  size=p->ks0*p->ks1;
  if(bitpix==FLOAT_IMG)
    p->kernel=array;
  else
    {
      changetype(array, bitpix, size, numnul, (void **)&p->kernel, FLOAT_IMG);
      free(array);
    }
  kernel=p->kernel;

  /* Convert Nul values to zero: */
  if(numnul)
    { fp=(f=kernel)+size; do *f = isnan(*f) ? 0 : *f; while(++f<fp); }

  /* Normalize the kernel: */
  if(p->kernelflip)
    {
      fp=(f=kernel)+size; do sum+=*f++; while(f<fp);
      fp=(f=kernel)+size; do *f++/=sum; while(f<fp);
    }

  /* Flip the kernel: */
  if(p->spatial && p->kernelflip)
    {
      for(i=0;i<size/2;++i)
        {
          tmp=kernel[i];
          kernel[i]=kernel[size-i-1];
          kernel[size-i-1]=tmp;
        }
    }
}




















/**************************************************************/
/************         Prepare the arrays          *************/
/**************************************************************/
void
sanitycheck(struct convolveparams *p)
{
  /* Check if the kernel has an odd number of pixels: */
  if(p->ks0%2==0 || p->ks1%2==0)
    error(EXIT_FAILURE, 0, "The kernel used for convolution has to have an "
          "odd number of pixels on all sides, %s is %lux%lu.",
          p->up.kernelname, p->ks1, p->ks0);

  /* Check the output file name: */
  if(p->cp.outputset)
    {
      if( dir0file1(p->cp.output, p->cp.dontdelete) == 0 )
        error(EXIT_FAILURE, 0, "Your output name (%s) is a directory.",
              p->cp.output);
    }
  else
    {
      automaticoutput(p->up.inputname, "_convolved.fits",
                      p->cp.removedirinfo, p->cp.dontdelete,
                      &p->cp.output);
      p->cp.outputset=1;
    }
  if(p->frequency && p->viewfreqsteps)
    automaticoutput(p->up.inputname, "_freqsteps.fits",
                    p->cp.removedirinfo, p->cp.dontdelete,
                    &p->up.freqstepsname);
}


















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct convolveparams *p)
{
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;

  /* Convert options: */
  p->kernelflip     = 1;
  p->kernelnorm     = 1;
  p->edgecorrection = 1;

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

  /* Prepare the necessary arrays: */
  preparearrays(p);

  /* Do a sanity check, then remove the possibly existing log file
     created by txttoarray. */
  sanitycheck(p);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf("Convolving %s (hdu: %s)\n"
             " with the kernel %s (hdu: %s).\n"
             " using %lu CPU threads in the %s domain.\n",
             p->up.inputname, cp->hdu, p->up.kernelname, p->up.khdu,
             cp->numthreads, p->spatial?"spatial":"frequency");
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct convolveparams *p, struct timeval *t1)
{
  free(p->input);
  free(p->kernel);
  free(p->cp.hdu);
  free(p->up.khdu);
  free(p->cp.output);
  wcsvfree(&p->nwcs, &p->wcs);

  /* Print the final message. */
  reporttiming(t1, SPACK_NAME" finished in: ", 0);
}
