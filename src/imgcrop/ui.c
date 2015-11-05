/*********************************************************************
ImageCrop - Crop a given size from one or multiple images.
ImageCrop is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include "crop.h"
#include "wcsmode.h"

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
readconfig(char *filename, struct imgcropparams *p)
{
  FILE *fp;
  int zeroorone;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  size_t lineno=0, len=200, tmp;
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
	  gal_checkset_int_zero_or_one(value, &zeroorone, name, key, SPACK,
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
	  gal_checkset_int_zero_or_one(value, &zeroorone, name, key, SPACK,
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





      /* Inputs: */
      else if(strcmp(name, "hdu")==0)
        gal_checkset_allocate_copy_set(value, &cp->hdu, &cp->hduset);

      else if(strcmp(name, "racol")==0)
	{
	  if(up->racolset) continue;
	  gal_checkset_sizet_el_zero(value, &p->racol, name, key, SPACK,
		      filename, lineno);
	  up->racolset=1;
	}
      else if(strcmp(name, "deccol")==0)
	{
	  if(up->deccolset) continue;
	  gal_checkset_sizet_el_zero(value, &p->deccol, name, key, SPACK,
		      filename, lineno);
	  up->deccolset=1;
	}
      else if(strcmp(name, "xcol")==0)
	{
	  if(up->xcolset) continue;
	  gal_checkset_sizet_el_zero(value, &p->xcol, name, key, SPACK,
		      filename, lineno);
	  up->xcolset=1;
	}
      else if(strcmp(name, "ycol")==0)
	{
	  if(up->ycolset) continue;
	  gal_checkset_sizet_el_zero(value, &p->ycol, name, key, SPACK,
		      filename, lineno);
	  up->ycolset=1;
	}
      else if(strcmp(name, "iwidth")==0)
	{
	  if(up->iwidthset) continue;
	  gal_checkset_sizet_l_zero(value, &tmp, name, key, SPACK,
		     filename, lineno);
	  p->iwidth[0]=p->iwidth[1]=tmp;
	  up->iwidthset=1;
	}
      else if(strcmp(name, "wwidth")==0)
	{
	  if(up->wwidthset) continue;
	  gal_checkset_double_l_0(value, &p->wwidth, name, key, SPACK,
		   filename, lineno);
	  up->wwidthset=1;
	}
      else if(strcmp(name, "hstartwcs")==0)
	{
	  if(up->hstartwcsset) continue;
	  gal_checkset_sizet_el_zero(value, &p->hstartwcs, name, key, SPACK,
                      filename, lineno);
	  up->hstartwcsset=1;
	}
      else if(strcmp(name, "hendwcs")==0)
	{
	  if(up->hendwcsset) continue;
	  gal_checkset_sizet_el_zero(value, &p->hendwcs, name, key, SPACK,
                      filename, lineno);
	  up->hendwcsset=1;
	}



      /* Outputs */
      else if(strcmp(name, "checkcenter")==0)
	{
	  if(up->checkcenterset) continue;
	  gal_checkset_sizet_el_zero(value, &p->checkcenter, name, key, SPACK,
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
printvalues(FILE *fp, struct imgcropparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    PRINTSTINGMAYBEWITHSPACE("hdu", cp->hdu);


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
    fprintf(fp, CONF_SHOWFMT"%ld\n", "iwidth", p->iwidth[0]);
  if(up->racolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "racol", p->racol);
  if(up->deccolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "deccol", p->deccol);
  if(up->wwidthset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "wwidth", p->wwidth);
  if(up->hstartwcsset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "hstartwcs", p->hstartwcs);
  if(up->hendwcsset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "hendwcs", p->hendwcs);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  PRINT_COMMONOPTIONS;
  if(up->imgmodeset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "imgmode", p->imgmode);
  if(up->wcsmodeset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "wcsmode", p->wcsmode);
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
  char forerr[100];
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;




  /* checkcenter is odd: */
  if(p->checkcenter%2==0)
    p->checkcenter+=1;



  /* Width and checkcenter are odd */
  if(p->iwidth[0]<3)
    error(EXIT_FAILURE, 0, "--iwidth has to be >3 pixels.");
  else if(p->iwidth[0]%2==0)
      p->iwidth[0]+=1;
  p->iwidth[1]=p->iwidth[0];
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
  checksum=up->raset+up->xcset+up->sectionset+up->polygonset;
  if(checksum)
    {
      /* Only one of the three should be called. */
      if(checksum!=1)
        {
          forerr[0]='\0';
          if(up->raset)
            strcat(forerr, "(`--ra' and `--dec'), ");
          if(up->xcset)
            strcat(forerr, "(`--xc' and `--yc'), ");
          if(up->sectionset)
            strcat(forerr, "(`--section'), ");
          if(up->polygonset)
            strcat(forerr, "(`--polygon'), ");
          error(EXIT_FAILURE, 0, "There are several ways to specify a crop "
                "box on the command line, see `--help`. But they should not "
                "be called together. You have asked for %s simultaneously!",
                forerr);
        }

      /* Check if the value for --output is a file or a directory? */
      p->outnameisfile=gal_checkset_dir_0_file_1(cp->output, cp->dontdelete);

      /* When there is only one output, only one thread is needed. */
      cp->numthreads=1;

      /* Not with a catalog. */
      if(up->catname)
	{
          if(up->sectionset) strcpy(forerr, "`--section'");
          if(up->polygonset) strcpy(forerr, "`--polygon'");
          if(up->xcset) strcpy(forerr, "`--xc' and `--yc'");
	  if(up->raset) strcpy(forerr, "`--ra' and `--dec'");
	  error(EXIT_FAILURE, 0, "A catalog name (%s) and command line crop "
		"parameters (%s) cannot be given together.", up->catname,
                forerr);
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
      gal_checkset_check_dir_write_add_slash(&cp->output);

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



  /* Parse the polygon vertices if they are given to make sure that
     there is no error: */
  if(p->up.polygonset)
    {
      polygonparser(p);
      if(p->nvertices<3)
        error(EXIT_FAILURE, 0, "A polygon has to have 3 or more vertices, "
              "you have only given %lu (%s).\n", p->nvertices, p->up.polygon);
      if(p->outpolygon && p->numimg>1)
        error(EXIT_FAILURE, 0, "Currently in WCS mode, outpolygon can only "
              "be set to zero when there is one image, you have given %lu "
              "images. For multiple images the region will be very large. "
              "It is best if you first crop out the larger region you want "
              "into one image, then mask the polygon.", p->numimg);
    }
  else
    p->wpolygon=p->ipolygon=NULL;




  /* Check that if multiple threads are beeing requested, CFITSIO is
     configured with the `--enable-reentrant` option. This is put here
     because the number of threads may change above. */
  if(cp->numthreads>1 && fits_is_reentrant()==0)
    error(EXIT_FAILURE, 0, "CFITSIO was not configured with the "
	  "`--enable-reentrant` option but you have asked to operate "
	  "on %lu threads. Please configure, make and install CFITSIO "
	  "again with this flag to use multiple threads, run `info %s "
	  "CFITSIO` for more information. Alternatively you can set "
	  "the number of threads to 1 by adding the `--numthreads=1` "
	  "or `-N1` options.", cp->numthreads, PACKAGE);
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
      gal_linkedlist_pop_from_stll(&p->up.stll, &img->name);
      gal_fitsarray_read_fits_hdu(img->name, p->cp.hdu, IMAGE_HDU, &tmpfits);
      gal_fitsarray_img_bitpix_size(tmpfits, &p->bitpix, img->naxes);
      gal_fitsarray_read_wcs(tmpfits, &img->nwcs, &img->wcs, p->hstartwcs,
                             p->hendwcs);
      if(img->wcs)
        {
          status=wcshdo(0, img->wcs, &img->nwcskeys, &img->wcstxt);
          if(status)
            error(EXIT_FAILURE, 0, "wcshdo ERROR %d: %s.", status,
                  wcs_errmsg[status]);
        }
      else
        if(p->wcsmode)
          error(EXIT_FAILURE, 0, "The WCS structure of %s (hdu: %s) "
                "image is not recognized. So RA and Dec cannot be used "
                "as input. You can try with pixel coordinates in the "
                "Image Mode (note that the crops will lack WCS "
                "header information).", img->name, p->cp.hdu);
      fits_close_file(tmpfits, &status);
      gal_fitsarray_io_error(status, NULL);

      /* Make sure all the images have the same BITPIX and set the
	 basic BITPIX related parameters. */
      if(firstbitpix==0)
	{
	  firstbitpix=p->bitpix;
	  p->datatype=gal_fitsarray_bitpix_to_dtype(p->bitpix);
	  p->bitnul=gal_fitsarray_bitpix_blank(p->bitpix);
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
  if(p->up.xcset || p->up.sectionset || p->up.raset || p->up.polygonset)
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
      gal_timing_report(&t1, msg, 1);
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
    gal_txtarray_txt_to_array(p->up.catname, &p->cat, &p->cs0, &p->cs1);

  /* If cp->output was not specified on the command line or in any of
     the configuration files, then automatic output should be used, in
     which case, cp->output should be the current directory. */
  if(p->cp.outputset==0)
    {
      p->cp.output=malloc(2+1); /* 2 is length of "./" */
      if(p->cp.output==NULL)
        error(EXIT_FAILURE, errno, "Space for output");
      strcpy(p->cp.output, "./");
      p->cp.outputset=1;
    }

  /* Do a sanity check. */
  sanitycheck(p);
  gal_checkset_check_remove_file(TXTARRAYVVLOG, 0);

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
  free(p->wpolygon);
  free(p->ipolygon);

  /* If these two pointers point to the same place,, that plce will be
     freed below. */
  if(p->log[0].name != p->cp.output)
    free(p->cp.output);

  /* Free the allocated WCS parameters: */
  for(i=0;i<p->numimg;++i)
    if(p->imgs[i].wcs)
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
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in: ", 0);
}
