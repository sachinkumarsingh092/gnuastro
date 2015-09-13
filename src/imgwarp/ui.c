/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (Gnuastro) package.

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
readconfig(char *filename, struct imgwarpparams *p)
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





      /* Outputs */
      else if(strcmp(name, "matrix")==0)
        allocatecopyset(value, &up->matrixstring, &up->matrixstringset);

      else if(strcmp(name, "output")==0)
        allocatecopyset(value, &cp->output, &cp->outputset);

      else if(strcmp(name, "maxblankfrac")==0)
        {
	  if(up->maxblankfracset) continue;
	  floatl0s1(value, &p->maxblankfrac, name, key, SPACK,
		       filename, lineno);
	  up->maxblankfracset=1;
        }


      /* Read options common to all programs */
      READ_COMMONOPTIONS_FROM_CONF


      else
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct imgwarpparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    PRINTSTINGMAYBEWITHSPACE("hdu", cp->hdu);


  fprintf(fp, "\n# Output parameters:\n");
  if(up->matrixstringset)
    PRINTSTINGMAYBEWITHSPACE("matrix", up->matrixstring);

  if(cp->outputset)
    PRINTSTINGMAYBEWITHSPACE("output", cp->output);

  if(up->maxblankfracset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "maxblankfrac", p->maxblankfrac);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program) */
  PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct imgwarpparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    REPORT_NOTSET("hdu");
  if(up->matrixstringset==0)
    REPORT_NOTSET("matrix");
  if(up->maxblankfracset==0)
    REPORT_NOTSET("maxblankfrac");

  END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************        Read Matrix         *******************/
/**************************************************************/
void
readmatrixoption(struct imgwarpparams *p)
{
  size_t counter=0;
  char *t, *tailptr;

  errno=0;
  p->matrix=malloc(9*sizeof *p->matrix);
  if(p->matrix==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for temporary array to keep "
          "the matrix option.", 9*sizeof *p->matrix);

  t=p->up.matrixstring;
  while(*t!='\0')
    {
      switch(*t)
        {
        case ' ': case '\t': case ',':
          ++t;
          break;
        default:
          errno=0;
          p->matrix[counter++]=strtod(t, &tailptr);
          if(errno) error(EXIT_FAILURE, errno, "In reading `%s`", t);
          if(tailptr==t)
            error(EXIT_FAILURE, 0, "The provided string `%s' for matrix "
                  "could not be read as a number.", t);
          t=tailptr;
          if(counter>9)       /* Note that it was ++'d! */
            error(EXIT_FAILURE, 0, "There are %lu elements in `%s', there "
                  "should be 4 or 9.", counter, p->up.matrixstring);
          /*printf("%f, %s\n", p->matrix[counter-1], t);*/
        }
    }

  switch(counter)
    {
    case 4:
      p->ms1=p->ms0=2;
      break;
    case 9:
      p->ms1=p->ms0=3;
      break;
    default:
      error(EXIT_FAILURE, 0, "There are %lu numbers in the string `%s'! "
            "It should contain 4 or 9 numbers (for a 2 by 2 or 3 by 3 "
            "matrix).", counter, p->up.matrixstring);
    }
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct imgwarpparams *p)
{
  double *d, *df, *tmp, *m=p->matrix;


  /* Set the output name: */
  if(p->cp.output)
    checkremovefile(p->cp.output, p->cp.dontdelete);
  else
    automaticoutput(p->up.inputname, "_warped.fits", p->cp.removedirinfo,
                    p->cp.dontdelete, &p->cp.output);


  /* Check the size of the input matrix, note that it might only have
     the wrong numbers when it is read from a file. */
  if(p->up.matrixname
     && (p->ms0!=2 || p->ms0!=3 || p->ms1!=2 || p->ms1!=3))
    error(EXIT_FAILURE, 0, "The given matrix in %s has %lu rows and %lu "
          "columns. Its size must be either 2x2 or 3x3.", p->up.matrixname,
          p->ms0, p->ms1);

  /* Check if there are any non-normal numbers in the matrix: */
  df=(d=m)+p->ms0*p->ms1;
  do
    if(!isfinite(*d++))
      error(EXIT_FAILURE, 0, "%f is not a `normal' number!", *(d-1));
  while(d<df);

  /* Check if the determinant is not zero: */
  if( ( p->ms0==2 && (m[0]*m[3] - m[1]*m[2] == 0) )
      ||
      ( p->ms0==3 &&
        (m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7]
         - m[2]*m[4]*m[6] - m[1]*m[3]*m[8] - m[0]*m[5]*m[7] == 0) ) )
      error(EXIT_FAILURE, 0, "The determinant of the given matrix "
            "is zero!");

  /* If the matrix only has two components, then convert it to
     three. */
  if(p->ms0==2)
    {
      errno=0;
      tmp=malloc(9*sizeof *tmp);
      if(tmp==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for 3 by 3 matrix",
              9*sizeof *tmp);
      tmp[0]=m[0];
      tmp[1]=m[1];
      tmp[3]=m[2];
      tmp[4]=m[3];
      tmp[8]=1.0f;
      tmp[7]=tmp[6]=tmp[5]=tmp[2]=0.0f;
      free(p->matrix);
      m=p->matrix=tmp;
    }

  /* Check if the transformation is spatially invariant */
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* It is important that the image names are stored in an array (for
   WCS mode in particular). We do that here. */
void
preparearrays(struct imgwarpparams *p)
{
  void *array;
  size_t numnul;
  double *inv, *m=p->matrix;

  /* Read in the input image: */
  numnul=fitsimgtoarray(p->up.inputname, p->cp.hdu, &p->inputbitpix,
                        &array, &p->is0, &p->is1);
  if(p->inputbitpix==DOUBLE_IMG)
    p->input=array;
  else
    {
      changetype(array, p->inputbitpix, p->is0*p->is1, numnul,
                 (void **)&p->input, DOUBLE_IMG);
      free(array);
    }
  readfitswcs(p->up.inputname, p->cp.hdu, &p->nwcs, &p->wcs);

  /* Make the inverse matrix: */
  errno=0;
  p->inverse=inv=malloc(9*sizeof *inv);
  if(inv==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for the inverse array.",
          9*sizeof *inv);
  inv[0]=m[4]*m[8]-m[5]*m[7];
  inv[1]=m[2]*m[7]-m[1]*m[8];
  inv[2]=m[1]*m[5]-m[2]*m[4];
  inv[3]=m[5]*m[6]-m[3]*m[8];
  inv[4]=m[0]*m[8]-m[2]*m[6];
  inv[5]=m[2]*m[3]-m[0]*m[5];
  inv[6]=m[3]*m[7]-m[4]*m[6];
  inv[7]=m[1]*m[6]-m[0]*m[7];
  inv[8]=m[0]*m[4]-m[1]*m[3];
  /* Just for a test:
  {
    size_t i;
    printf("\nInput matrix:");
    for(i=0;i<9;++i) { if(i%3==0) printf("\n"); printf("%-10.5f", m[i]); }
    printf("\n-----------\n");
    printf("Inverse matrix:");
    for(i=0;i<9;++i) { if(i%3==0) printf("\n"); printf("%-10.5f", inv[i]); }
    printf("\n\n");
  }
  */
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct imgwarpparams *p)
{
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;
  p->correctwcs     = 1;

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

  /* Read catalog if given. */
  if(p->up.matrixname)
    txttoarray(p->up.matrixname, &p->matrix, &p->ms0, &p->ms1);
  else
    readmatrixoption(p);

  /* Do a sanity check. */
  sanitycheck(p);
  checkremovefile(TXTARRAYVVLOG, 0);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf(" Input image: %s\n", p->up.inputname);
      printf(" matrix:"
             "\n\t%.4f   %.4f   %.4f"
             "\n\t%.4f   %.4f   %.4f"
             "\n\t%.4f   %.4f   %.4f\n",
             p->matrix[0], p->matrix[1], p->matrix[2],
             p->matrix[3], p->matrix[4], p->matrix[5],
             p->matrix[6], p->matrix[7], p->matrix[8]);
    }

  /* Make the array of input images. */
  preparearrays(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct imgwarpparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->input);
  free(p->matrix);
  free(p->cp.hdu);
  free(p->inverse);
  free(p->cp.output);

  if(p->wcs)
    wcsvfree(&p->nwcs, &p->wcs);

  /* Print the final message. */
  reporttiming(t1, SPACK_NAME" finished in: ", 0);
}
