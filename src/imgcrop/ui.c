/*********************************************************************
Image Crop - Crop a given size from one or multiple images.
Image Crop is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
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
#include "crop.h"
#include "wcsmode.h"

#include "ui.h"		        /* Needs main.h                   */
#include "args.h"	        /* Needs main.h, includes argp.h. */


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
readconfig(char *filename, struct imgcropparams *p)
{
  FILE *fp;
  int zeroorone;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;
  char key='a';	/* Not used, just a place holder. */
  int imgmodeset=0, wcsmodeset=0; /* For unambiguous default file checking. */

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

      /* Operating modes: */
      if(strcmp(name, "imgmode")==0)
	{
	  intzeroorone(value, &zeroorone, name, key, SPACK,
		       filename, lineno);
	  if(zeroorone)
	    {
	      imgmodeset=1;
	      if(wcsmodeset)
		error_at_line(EXIT_FAILURE, 0, filename, lineno,
			      "Image and WCS modes cannot be called "
			      "together. It is ambiguous.");
	      if(up->imgmodeset==0)
		{
		  p->imgmode=1;
		  p->wcsmode=0;
		  up->imgmodeset=up->wcsmodeset=1;
		}
	    }
	}
      else if(strcmp(name, "wcsmode")==0)
	{
	  intzeroorone(value, &zeroorone, name, key, SPACK,
		       filename, lineno);
	  if(zeroorone)
	    {
	      wcsmodeset=1;
	      if(imgmodeset)
		error_at_line(EXIT_FAILURE, 0, filename, lineno,
			      "Image and WCS modes cannot be called "
			      "together. It is ambiguous.");
	      if(up->wcsmodeset==0)
		{
		  p->imgmode=0;
		  p->wcsmode=1;
		  up->imgmodeset=up->wcsmodeset=1;
		}
	    }
	}
      else if(strcmp(name, "numthreads")==0)
	{
	  if (cp->numthreadsset) continue;
	  sizetlzero(value, &cp->numthreads, name, key, SPACK,
		     filename, lineno);
	  cp->numthreadsset=1;
	}





      /* Inputs: */
      else if(strcmp(name, "hdu")==0)
	{
	  if(cp->hduset) continue;
	  errno=0;
	  cp->hdu=malloc(strlen(value)+1);
	  if(cp->hdu==NULL)
	    error(EXIT_FAILURE, errno, "Space for HDU.");
	  strcpy(cp->hdu, value);
	  cp->hduset=1;
	}
      else if(strcmp(name, "racol")==0)
	{
	  if(up->racolset) continue;
	  sizetelzero(value, &p->racol, name, key, SPACK,
		      filename, lineno);
	  up->racolset=1;
	}
      else if(strcmp(name, "deccol")==0)
	{
	  if(up->deccolset) continue;
	  sizetelzero(value, &p->deccol, name, key, SPACK,
		      filename, lineno);
	  up->deccolset=1;
	}
      else if(strcmp(name, "xcol")==0)
	{
	  if(up->xcolset) continue;
	  sizetelzero(value, &p->xcol, name, key, SPACK,
		      filename, lineno);
	  up->xcolset=1;
	}
      else if(strcmp(name, "ycol")==0)
	{
	  if(up->ycolset) continue;
	  sizetelzero(value, &p->ycol, name, key, SPACK,
		      filename, lineno);
	  up->ycolset=1;
	}
      else if(strcmp(name, "iwidth")==0)
	{
	  if(up->iwidthset) continue;
	  sizetlzero(value, &p->iwidth, name, key, SPACK,
		     filename, lineno);
	  up->iwidthset=1;
	}
      else if(strcmp(name, "wwidth")==0)
	{
	  if(up->wwidthset) continue;
	  doublel0(value, &p->wwidth, name, key, SPACK,
		   filename, lineno);
	  up->wwidthset=1;
	}





      /* Outputs */
      else if(strcmp(name, "checkcenter")==0)
	{
	  if(up->checkcenterset) continue;
	  sizetelzero(value, &p->checkcenter, name, key, SPACK,
		      filename, lineno);
	  up->checkcenterset=1;
	}
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
      else if(strcmp(name, "suffix")==0)
	{
	  if(up->suffixset) continue;
	  errno=0;
	  p->suffix=malloc(strlen(value)+1);
	  if(p->suffix==NULL)
	    error(EXIT_FAILURE, errno, "Space for prefix.");
	  strcpy(p->suffix, value);
	  up->suffixset=1;
	}
      else
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct imgcropparams *p)
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


  fprintf(fp, "\n# Output parameters:\n");
  if(up->checkcenterset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "checkcenter", p->checkcenter);
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);
  if(up->suffixset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "suffix", p->suffix);


  fprintf(fp, "\n# Crop parameters:\n");
  if(up->xcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "xcol", p->xcol);
  if(up->ycolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "ycol", p->ycol);
  if(up->iwidthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "iwidth", p->iwidth);
  if(up->racolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "racol", p->racol);
  if(up->deccolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "deccol", p->deccol);
  if(up->wwidthset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "wwidth", p->wwidth);


  fprintf(fp, "\n# Operating mode:\n");
  if(up->imgmodeset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "imgmode", p->imgmode);
  if(up->wcsmodeset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "wcsmode", p->wcsmode);
  /* Number of threads doesn't need to be checked, it is set by
     default */
  fprintf(fp, CONF_SHOWFMT"%lu\n", "numthreads", cp->numthreads);
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct imgcropparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;
  if(up->imgmodeset==0 && up->wcsmodeset==0)
    REPORT_NOTSET("imgmode or wcsmode");
  if(cp->hduset==0)
    REPORT_NOTSET("hdu");
  if(up->xcolset==0)
    REPORT_NOTSET("xcol");
  if(up->ycolset==0)
    REPORT_NOTSET("ycol");
  if(up->iwidthset==0)
    REPORT_NOTSET("iwidth");
  if(up->racolset==0)
    REPORT_NOTSET("racol");
  if(up->deccolset==0)
    REPORT_NOTSET("deccol");
  if(up->wwidthset==0)
    REPORT_NOTSET("wwidth");
  if(cp->outputset==0)
    REPORT_NOTSET("output");
  if(up->suffixset==0)
    REPORT_NOTSET("suffix");
  if(up->checkcenterset==0)
    REPORT_NOTSET("checkcenter");
  END_OF_NOTSET_REPORT;
}





















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct imgcropparams *p)
{
  int checksum;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;
  char *t1=NULL, *t2=NULL, *t3=NULL, *boxparam;



  /* Check that if multiple threads are beeing requested, CFITSIO is
     configured with the `--enable-reentrant` option. */
  if(cp->numthreads>1 && fits_is_reentrant()==0)
    error(EXIT_FAILURE, 0, "CFITSIO was not configured with the "
	  "`--enable-reentrant` option but you have asked to operate "
	  "on %lu threads. Please configure, make and install CFITSIO "
	  "again with this flag to use multiple threads, run `info %s "
	  "CFITSIO` for more information. Alternatively you can set "
	  "the number of threads to 1 by adding the `--numthreads=1` "
	  "or `-N1` options.", cp->numthreads, PACKAGE);



  /* checkcenter is odd: */
  if(p->checkcenter%2==0)
    p->checkcenter+=1;



  /* Width and checkcenter are odd */
  if(p->iwidth<3)
    error(EXIT_FAILURE, 0, "--iwidth has to be >3 pixels.");
  else if(p->iwidth%2==0)
    p->iwidth+=1;
  if(p->checkcenter<3)
    error(EXIT_FAILURE, 0, "--checkcenter has to be >3 pixels.");
  else if(p->checkcenter%2==0)
    p->checkcenter+=1;



  /* deccol!=racol and xcol!=ycol */
  if(p->racol==p->deccol)
    error(EXIT_FAILURE, 0, "The columns for RA and Dec must be "
	  "different.");
  if(p->xcol==p->ycol)
    error(EXIT_FAILURE, 0, "The columns for x and y must be "
	  "different.");



  /* Make sure that if any of --ra or --dec are given, the other is
     also given. */
  checksum=up->raset+up->decset;
  if(checksum==2) {p->imgmode=0; p->wcsmode=1;}
  else if(checksum)/* is not zero */
    error(EXIT_FAILURE, 0, "The options `--ra` and `--dec` should be "
	  "called together.");



  /* Make sure that if any of --xc or --yc are given, the other is
     also given. */
  checksum=up->xcset+up->ycset;
  if(checksum==2) {p->imgmode=1; p->wcsmode=0;}
  else if(checksum)/* is not zero */
    error(EXIT_FAILURE, 0, "The options `--xc` and `--yc` should be "
	  "called together.");



  /* --section is given, it goes into image mode: */
  if(up->sectionset) {p->imgmode=1; p->wcsmode=0;}



  /* Make sure that the multiple one crop box options have not been
     called together. */
  checksum=up->raset+up->xcset+up->sectionset;
  if(checksum)
    {
      if(up->raset) t1="(`--ra` and `--dec`) ";
      if(up->xcset) t2="(`--xc` and `--yc`) ";
      if(up->sectionset) t3="(`--section) ";

      /* Only one of the three should be called. */
      if(checksum!=1)
	error(EXIT_FAILURE, 0, "There are several ways to specify a crop "
	      "box on the command line, see `--help`. But they should not "
	      "be called together. You have asked for "
	      "%s%s%s simultaneously!", t1?t1:"", t2?t2:"", t3?t3:"");

      /* Check if the value for --output is a file or a directory? */
      p->outnameisfile=nameisawritablefile(cp->output, cp->dontdelete);

      /* When there is only one output, only one thread is needed. */
      cp->numthreads=1;

      /* Not with a catalog. */
      if(up->catname)
	{
	  if(t1) boxparam=t1; else if(t2) boxparam=t2; else boxparam=t3;
	  error(EXIT_FAILURE, 0, "A catalog name (%s) and one box parameters "
		"%s""cannot be given together on the command line. ",
		up->catname, boxparam);
	}
    }
  else
    {
      /* Only one mode. Note that when the box is specified on the
	 command line, in the steps above, we set the image mode or
	 wcs mode.*/
      if(p->imgmode && p->wcsmode)
	error(EXIT_FAILURE, 0, "Only one of imgmode or wcsmode "
	      "must be called. They cannot operate together.");
      else if(p->imgmode==0 && p->wcsmode==0)
	error(EXIT_FAILURE, 0, "At least one of imgmode or "
	      "wcsmode must be called.");

      /* Make sure a catalog is set. */
      if(up->catset)
	{
	  if(p->numimg>1 && p->imgmode)
	    error(EXIT_FAILURE, 0, "In image mode, when a catalog is "
		  "specified, only one image may be provided.");
	}
      else
	error(EXIT_FAILURE, 0, "No catalog. When no crop coordinates "
	      "are specified on the command line, a catalog must be "
	      "provided.");

      /* Make sure the given output is a directory. */
      checkdirwriteaddslash(&cp->output);

      /* Make sure the columns of data are within the catalog range of
	 columns: */
      if(p->imgmode)
	{
	  CHECKCOLINCAT(p->xcol, "xcol");
	  CHECKCOLINCAT(p->ycol, "ycol");
	}
      else
	{
	  CHECKCOLINCAT(p->racol, "racol");
	  CHECKCOLINCAT(p->deccol, "deccol");
	}
    }



  /* If in image mode, there should only be one input image. */
  if(p->imgmode && p->numimg>1)
    error(EXIT_FAILURE, 0, "In image mode, only one input image may be "
	  "specified.");




  /* If we are in WCS mode, noblanks must be off */
  if(p->wcsmode && p->noblank)
    error(EXIT_FAILURE, 0, "`--noblanks` (`-b`) is only for image mode. "
	  "You have called it with WCS mode.");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* It is important that the image names are stored in an array (for
   WCS mode in particular). We do that here. */
void
preparearrays(struct imgcropparams *p)
{
  size_t size, num;
  fitsfile *tmpfits;
  struct timeval t1;
  struct inputimgs *img;
  char msg[VERBMSGLENGTH_V];
  int i, status, firstbitpix=0;

  if(p->cp.verb) gettimeofday(&t1, NULL);

  /* Fill in the WCS information of each image. This is done here
     because WCSLIB is unfortunately not thread-safe when reading the
     WCS information from the FITS files. In cases where the number of
     cropped images are more than the input images, this can also be a
     preformance boost because each image information is only read
     once.

     The images are filled in opposite order because we used a linked
     list to read them in, which is a first in first out structure.*/
  errno=0;
  size=p->numimg*sizeof *p->imgs;
  p->imgs=malloc(size);
  if(p->imgs==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes for p->imgs", size);
  for(i=p->numimg-1;i>=0;--i)
    {
      /* Get the image properties. */
      status=0;
      img=&p->imgs[i];
      pop_from_stll(&p->up.stll, &img->name);
      readfitshdu(img->name, p->cp.hdu, IMAGE_HDU, &tmpfits);
      imgbitpixsize(tmpfits, &p->bitpix, img->naxes);
      readwcs(tmpfits, &img->nwcs, &img->wcs);
      status=wcshdo(WCSHDO_safe, img->wcs, &img->nwcskeys, &img->wcstxt);
      if(status)
	error(EXIT_FAILURE, 0, "wcshdo ERROR %d: %s.", status,
	      wcs_errmsg[status]);
      fits_close_file(tmpfits, &status);
      fitsioerror(status, NULL);

      /* Make sure all the images have the same BITPIX and set the
	 basic BITPIX related parameters. */
      if(firstbitpix==0)
	{
	  firstbitpix=p->bitpix;
	  p->datatype=bitpixtodtype(p->bitpix);
	  p->bitnul=bitpixnull(p->bitpix);
	}
      else if(firstbitpix!=p->bitpix)
	error(EXIT_FAILURE, 0, "%s: BITPIX=%d. Previous images had a "
	      "BITPIX value of %d, For "SPACK_NAME" to work, all images "
	      "must have the same pixel data type.",
	      img->name, p->bitpix, firstbitpix);

      /* In WCS mode, Check resolution and get the first pixel
	 positions. */
      if(p->wcsmode) wcscheckprepare(p, img);
    }

  /* Array of log structures. We will make one more than the needed
     numbers so that we can put a NULL character in the name section
     of it to sign its end (something like a string). This is done so
     we don't have to worry about the length calculation any more! */
  if(p->up.xcset || p->up.sectionset || p->up.raset)
    num=1;
  else
    num=p->cs0;
  errno=0;
  p->log=calloc(num+1, sizeof *p->log);
  if(p->log==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes for p->log",
	  num+1 * sizeof *p->log);

  /* Report timing: */
  if(p->cp.verb)
    {
      sprintf(msg, "Read metadata of %lu images.", p->numimg);
      reporttiming(&t1, msg, 1);
    }
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct imgcropparams *p)
{
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

  /* Read catalog if given. */
  if(p->up.catname)
    txttoarray(p->up.catname, &p->cat, &p->cs0, &p->cs1);

  /* Do a sanity check. */
  sanitycheck(p);
  checkremovefile(TXTARRAYVVLOG, 0);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    printf(SPACK_NAME" started on %s", ctime(&p->rawtime));

  /* Make the array of input images. */
  preparearrays(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct imgcropparams *p, struct timeval *t1)
{
  size_t i;
  int status;

  /* Free the allocated arrays: */
  free(p->cat);
  free(p->cp.hdu);
  free(p->bitnul);
  free(p->suffix);

  /* If these two pointers point to the same place,, that plce will be
     freed below. */
  if(p->log[0].name != p->cp.output)
    free(p->cp.output);

  /* Free the allocated WCS parameters: */
  for(i=0;i<p->numimg;++i)
    {
      free(p->imgs[i].wcstxt);
      status=wcsvfree(&p->imgs[i].nwcs, &p->imgs[i].wcs);
      if(status)
	error(EXIT_FAILURE, 0, "wcsvfree ERROR %d: %s.", status,
	      wcs_errmsg[status]);
    }
  free(p->imgs);

  /* Free the log array: */
  for(i=0;p->log[i].name;++i)
    free(p->log[i].name);
  free(p->log);

  /* Print the final message. */
  reporttiming(t1, SPACK_NAME" finished in: ", 0);
}
