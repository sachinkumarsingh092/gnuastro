/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include "linkedlist.h"
#include "configfiles.h"
#include "fitsarrayvv.h"

#include "main.h"
#include "args.h"


/* Set the file names of the places where the default parameters are
   put. */
#define CONFIG_FILE SPACK CONF_POSTFIX
#define SYSCONFIG_FILE SYSCONFIG_DIR CONFIG_FILE
#define USERCONFIG_FILEEND USERCONFIG_DIR CONFIG_FILE
#define CURDIRCONFIG_FILE CURDIRCONFIG_DIR CONFIG_FILE





/**************************************************************/
/**************       Options and parameters    ***************/
/**************************************************************/
void
readconfig(char *filename, struct converttparams *p)
{
  int tmp;
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
	{
	  if(cp->hduset) continue;
	  errno=0;
	  cp->hdu=malloc(strlen(value)+1);
	  if(cp->hdu==NULL)
	    error(EXIT_FAILURE, errno, "Space for HDU.");
	  strcpy(cp->hdu, value);
	  cp->hduset=1;
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
      else if(strcmp(name, "quality")==0)
	{
	  if(up->qualityset) continue;
          intsmallerequalto(value, &p->quality, name, key,
                            p->cp.spack, filename, lineno, 100);
          if(p->quality<0)
            error(EXIT_FAILURE, 0, "The quality option should be positive.");
	  up->qualityset=1;
	}
      else if(strcmp(name, "kincmyk")==0)
	{
	  if(up->kincmykset) continue;
          intzeroorone(value, &p->kincmyk, name, key, SPACK,
                       filename, lineno);
	  up->kincmykset=1;
	}
      else if(strcmp(name, "widthincm")==0)
	{
	  if(up->widthincmset) continue;
          floatl0(value, &p->widthincm, name, key, SPACK, filename, lineno);
	  up->widthincmset=1;
	}





      /* Flux: */
      else if(strcmp(name, "fluxlow")==0)
	{
	  if(up->fluxlowset) continue;
          anyfloat(value, &p->fluxlow, name, key, p->cp.spack,
                   filename, lineno);
          up->fluxlowset=1;
	}
      else if(strcmp(name, "fluxhigh")==0)
	{
	  if(up->fluxhighset) continue;
          anyfloat(value, &p->fluxhigh, name, key, p->cp.spack,
                   filename, lineno);
          up->fluxhighset=1;
	}
      else if(strcmp(name, "maxbyte")==0)
	{
	  if(up->maxbyteset) continue;
          intsmallerequalto(value, &tmp, "maxbyte", key,
                            p->cp.spack, NULL, 0, UINT8_MAX);
          if(tmp<0)
            error(EXIT_FAILURE, 0, "--maxbyte (-m) should be positive.");
          p->maxbyte=tmp;
          p->up.maxbyteset=1;
	}



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
setparams(int argc, char *argv[], struct converttparams *p)
{
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;
  p->invert         = 1;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "Parsing arguments");
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct converttparams *p)
{
  free(p->cp.hdu);
  free(p->cp.output);
}
