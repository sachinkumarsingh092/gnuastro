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
#include <math.h>
#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>

#include "timing.h"	/* Includes time.h and sys/time.h   */
#include "txtarrayvv.h"
#include "linkedlist.h"
#include "configfiles.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "eps.h"
#include "jpeg.h"
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
      else if(strcmp(name, "h2")==0)
	{
	  if(up->h2set) continue;
	  errno=0;
	  up->h2=malloc(strlen(value)+1);
	  if(up->h2==NULL)
	    error(EXIT_FAILURE, errno, "Space for h2.");
	  strcpy(up->h2, value);
	  up->h2set=1;
	}
      else if(strcmp(name, "h3")==0)
	{
	  if(up->h3set) continue;
	  errno=0;
	  up->h3=malloc(strlen(value)+1);
	  if(up->h3==NULL)
	    error(EXIT_FAILURE, errno, "Space for h3.");
	  strcpy(up->h3, value);
	  up->h3set=1;
	}
      else if(strcmp(name, "h4")==0)
	{
	  if(up->h4set) continue;
	  errno=0;
	  up->h4=malloc(strlen(value)+1);
	  if(up->h4==NULL)
	    error(EXIT_FAILURE, errno, "Space for h4.");
	  strcpy(up->h4, value);
	  up->h4set=1;
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
printvalues(FILE *fp, struct converttparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    {
      if(stringhasspace(cp->hdu))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "hdu", cp->hdu);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "hdu", cp->hdu);
    }
  if(up->h2set)
    {
      if(stringhasspace(up->h2))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "h2", up->h2);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "h2", up->h2);
    }
  if(up->h3set)
    {
      if(stringhasspace(up->h3))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "h3", up->h3);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "h3", up->h3);
    }
  if(up->h4set)
    {
      if(stringhasspace(up->h4))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "h4", up->h4);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "h4", up->h4);
    }


  fprintf(fp, "\n# Output parameters:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);
  if(up->qualityset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "quality", p->quality);
  if(up->kincmykset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "kincmyk", p->kincmyk);
  if(up->widthincmset)
    fprintf(fp, CONF_SHOWFMT"%.2f\n", "widthincm", p->widthincm);


  fprintf(fp, "\n# Output flux display:\n");
  if(up->fluxlowset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "fluxlow", p->fluxlow);
  if(up->fluxhighset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "fluxhigh", p->fluxhigh);
  if(up->maxbyteset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "maxbyte", (int)(p->maxbyte));
}





/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct converttparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    REPORT_NOTSET("hdu");
  if(up->h2set==0)
    REPORT_NOTSET("h2");
  if(up->h3set==0)
    REPORT_NOTSET("h3");
  if(up->h4set==0)
    REPORT_NOTSET("h4");
  if(cp->outputset==0)
    REPORT_NOTSET("output");
  if(up->qualityset==0)
    REPORT_NOTSET("quality");
  if(up->widthincmset==0)
    REPORT_NOTSET("widthincm");
  if(up->kincmykset==0)
    REPORT_NOTSET("kincmyk");
  if(up->fluxlowset==0)
    REPORT_NOTSET("fluxlow");
  if(up->fluxhighset==0)
    REPORT_NOTSET("fluxhigh");
  if(up->maxbyteset==0)
    REPORT_NOTSET("maxbyte");

  END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct converttparams *p)
{
  /* Make sure there are 1 (for grayscale), 3 (for RGB) or 4 (for
     CMYK) input images. */
  if(p->numinputs!=1 && p->numinputs!=3 && p->numinputs!=4)
    error(EXIT_FAILURE, 0, "The number of input images has to be 1 (for "
          "non image data, grayscale or only K channel in CMYK), 3 (for "
          "RGB) and 4 (for CMYK). You have given %lu input images.",
          p->numinputs);
}




















/**************************************************************/
/***************        Preparations        *******************/
/**************************************************************/
void
preparearrays(struct converttparams *p)
{
  int bitpix;
  void *array;
  double *d, *df;
  size_t i, numnul;
  struct stll *tmp;
  char *hdu, *names[4];

  /* Put the names in the correct order. */
  i=p->numinputs-1;
  for(tmp=p->inputnames; tmp!=NULL; tmp=tmp->next)
    names[i--]=tmp->v;

  p->numch=0;
  for(i=0;i<p->numinputs;++i)
    {
      /* Check if p->numch has not exceeded 4. */
      if(p->numch>=4)
        error(EXIT_FAILURE, 0, "The number of input color channels (not "
              "files) has exceeded 4! Note that one file can contain more "
              "than one color channel.");

      /* FITS: */
      if( nameisfits(names[i]) )
        {
          switch(p->numch) /* Get the HDU value for this channel. */
            {
            case 0: hdu=p->cp.hdu; break;  case 1: hdu=p->up.h2; break;
            case 2: hdu=p->up.h3; break;   case 3: hdu=p->up.h4; break;
            default: error(EXIT_FAILURE, 0, "A bug! In parsing the input "
                           "FITS files, it has gone beyond four! Please "
                           "contact us so we can see what caused this "
                           "problem and fix it.");
            }
          numnul=fitsimgtoarray(names[i], hdu, &bitpix, &array,
                                &p->s0[p->numch],
                                &p->s1[p->numch]);
          changetype(array, bitpix,
                     p->s0[p->numch]*p->s1[p->numch],
                     numnul, (void **)(&p->ch[p->numch]),
                     DOUBLE_IMG);
          free(array);
          ++p->numch;
        }



      /* JPEG: */
      else if ( nameisjpeg(names[i]) )
        preparejpeg(p, names[i]);



      /* EPS:  */
      else if ( nameiseps(names[i]) )
        error(EXIT_FAILURE, 0, "EPS files cannot be used as input. Since "
              "EPS files are not raster graphics, they are only used as "
              "output.");



      /* Text: */
      else
        {
          txttoarray(names[i], &p->ch[p->numch],
                     &p->s0[p->numch], &p->s1[p->numch]);
          df = (d=p->ch[p->numch]) + p->s0[p->numch]*p->s1[p->numch];
          do if(isnan(*d++)) break; while(d<df);
          if(d==df)
            checkremovefile(TXTARRAYVVLOG, 0);
          else
            error(EXIT_FAILURE, 0, "%s contains non-numeric data, see %s.",
                  names[i], TXTARRAYVVLOG);
          ++p->numch;
        }
    }

  /* Check if all the input sources have the same size: */
  if(p->numch>1)
    {
      for(i=1;i<p->numch;++i)
        if(p->s0[i]!=p->s0[0] || p->s1[i]!=p->s1[0])
          break;
      if(i!=p->numch)
        {
          for(i=0;i<p->numch;++i)
            fprintf(stderr, "Channel %lu is %lu x %lu pixels.\n", i,
                    p->s1[i], p->s0[i]);
          error(EXIT_FAILURE, 0, "The input color channels have different "
                "sizes.");
        }
    }

}


















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
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

  /* Add the user default values and save them if asked. */
  CHECKSETCONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    REPORT_PARAMETERS_SET;

  /* Do a sanity check.
  sanitycheck(p);
  */
  /* Prepare the arrays: */
  preparearrays(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct converttparams *p)
{
  size_t i;

  free(p->cp.hdu);
  free(p->cp.output);
  for(i=0;i<4;++i) free(p->ch[i]);
}
