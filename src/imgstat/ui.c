/*********************************************************************
ImageStatistics - Get general statistics about the image.
ImgeStatistics is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include "timing.h"	/* Includes time.h and sys/time.h   */
#include "checkset.h"
#include "forqsort.h"
#include "arraymanip.h"
#include "statistics.h"
#include "txtarrayvv.h"
#include "commonargs.h"
#include "configfiles.h"
#include "fitsarrayvv.h"
#include "fixedstringmacros.h"

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
readconfig(char *filename, struct imgstatparams *p)
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
      GAL_CONFIGFILES_START_READING_LINE;




      /* Inputs: */
      if(strcmp(name, "hdu")==0)
        gal_checkset_allocate_copy_set(value, &cp->hdu, &cp->hduset);

      else if (strcmp(name, "mask")==0)
        gal_checkset_allocate_copy_set(value, &up->maskname, &up->masknameset);

      else if (strcmp(name, "mhdu")==0)
        gal_checkset_allocate_copy_set(value, &up->mhdu, &up->mhduset);
      else if(strcmp(name, "mirrordist")==0)
	{
	  if(up->mirrordistset) continue;
          gal_checkset_float_l_0(value, &p->mirrordist, name, key, SPACK,
                                 filename, lineno);
          up->mirrordistset=1;
	}



      /* Outputs */
      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output, &cp->outputset);

      else if(strcmp(name, "mirrorplotdist")==0)
	{
	  if(up->mirrorplotdistset) continue;
          gal_checkset_float_l_0(value, &p->mirrorplotdist, name, key, SPACK,
                                 filename, lineno);
          up->mirrorplotdistset=1;
	}
      else if(strcmp(name, "onebinvalue")==0)
	{
	  if(up->onebinvalueset) continue;
          gal_checkset_any_float(value, &p->onebinvalue, name, key, SPACK,
                                 filename, lineno);
          up->onebinvalueset=1;
	}


      /* Histogram: */
      else if(strcmp(name, "histnumbins")==0)
	{
	  if(up->histnumbinsset) continue;
          gal_checkset_sizet_l_zero(value, &p->histnumbins, name, key, SPACK,
                                    filename, lineno);
          up->histnumbinsset=1;
	}
      else if(strcmp(name, "histmin")==0)
	{
	  if(up->histminset) continue;
          gal_checkset_any_float(value, &p->histmin, name, key, SPACK,
                                 filename, lineno);
          up->histminset=1;
	}
      else if(strcmp(name, "histmax")==0)
	{
	  if(up->histmaxset) continue;
          gal_checkset_any_float(value, &p->histmax, name, key, SPACK,
                                 filename, lineno);
          up->histmaxset=1;
	}
      else if(strcmp(name, "histquant")==0)
	{
	  if(up->histquantset) continue;
          gal_checkset_float_l_0_s_1(value, &p->histquant, name, key, SPACK,
                                     filename, lineno);
          up->histquantset=1;
	}


      /* Cumulative Frequency Plot: */
      else if(strcmp(name, "cfpnum")==0)
	{
	  if(up->cfpnumset) continue;
          gal_checkset_sizet_l_zero(value, &p->cfpnum, name, key, SPACK,
                                    filename, lineno);
          up->cfpnumset=1;
	}
      else if(strcmp(name, "cfpmin")==0)
	{
	  if(up->cfpminset) continue;
          gal_checkset_any_float(value, &p->cfpmin, name, key, SPACK,
                                 filename, lineno);
          up->cfpminset=1;
	}
      else if(strcmp(name, "cfpmax")==0)
	{
	  if(up->cfpmaxset) continue;
          gal_checkset_any_float(value, &p->cfpmax, name, key, SPACK,
                                 filename, lineno);
          up->cfpmaxset=1;
	}
      else if(strcmp(name, "cfpquant")==0)
	{
	  if(up->cfpquantset) continue;
          gal_checkset_float_l_0_s_1(value, &p->cfpquant, name, key, SPACK,
                                     filename, lineno);
          up->cfpquantset=1;
	}

      /* Sigma clipping: */
      else if(strcmp(name, "sigclipmultip")==0)
	{
	  if(up->sigclipmultipset) continue;
          gal_checkset_float_l_0(value, &p->sigclipmultip, name, key, SPACK,
                                 filename, lineno);
          up->sigclipmultipset=1;
	}
      else if(strcmp(name, "sigcliptolerance")==0)
	{
	  if(up->sigcliptoleranceset) continue;
          gal_checkset_float_l_0(value, &p->sigcliptolerance, name, key, SPACK,
                                 filename, lineno);
          up->sigcliptoleranceset=1;
	}
      else if(strcmp(name, "sigclipnum")==0)
	{
	  if(up->sigclipnumset) continue;
          gal_checkset_sizet_l_zero(value, &p->sigclipnum, name, key, SPACK,
                                    filename, lineno);
          up->sigclipnumset=1;
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
printvalues(FILE *fp, struct imgstatparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);
  if(up->masknameset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mask", up->maskname);
  if(up->mhdu)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mhdu", up->mhdu);
  if(up->mirrordistset)
    fprintf(fp, CONF_SHOWFMT"%.2f\n", "mirrordist", p->mirrordist);

  /* Output: */
  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);
  if(up->mirrorplotdistset)
    fprintf(fp, CONF_SHOWFMT"%.2f\n", "mirrorplotdist", p->mirrorplotdist);
  if(up->onebinvalueset)
    fprintf(fp, CONF_SHOWFMT"%.5f\n", "onebinvalue", p->onebinvalue);

  /* Histogram: */
  fprintf(fp, "\n# Histogram:\n");
  if(up->histnumbinsset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "histnumbins", p->histnumbins);
  if(up->histminset)
    fprintf(fp, CONF_SHOWFMT"%.5f\n", "histmin", p->histmin);
  if(up->histmaxset)
    fprintf(fp, CONF_SHOWFMT"%.5f\n", "histmax", p->histmax);
  if(up->histquantset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "histquant", p->histquant);

  /* Cumulative frequency plot: */
  fprintf(fp, "\n# Cumulative frequency plot:\n");
  if(up->cfpnumset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "cfpnum", p->cfpnum);
  if(up->cfpminset)
    fprintf(fp, CONF_SHOWFMT"%.5f\n", "cfpmin", p->cfpmin);
  if(up->cfpmaxset)
    fprintf(fp, CONF_SHOWFMT"%.5f\n", "cfpmax", p->cfpmax);
  if(up->cfpquantset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "cfpquant", p->cfpquant);

  /* Sigma clipping: */
  if(up->sigclipmultipset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "sigclipmultip", p->sigclipmultip);
  if(up->sigcliptoleranceset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "sigcliptolerance",
            p->sigcliptolerance);
  if(up->sigclipnumset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "sigclipnum", p->sigclipnum);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. Note that those options which are not mandatory
   must not be listed here. */
void
checkifset(struct imgstatparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");
  if(up->mirrordistset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("mirrordist");
  if(up->mirrorplotdistset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("mirrorplotdist");
  if(up->onebinvalueset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("onebinvalue");
  if(up->histnumbinsset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("histnumbins");
  if(up->cfpnumset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("cfpnum");
  if(up->sigclipmultipset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("sigclipmultip");
  if(up->sigcliptoleranceset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("sigcliptolerance");
  if(up->sigclipnumset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("sigclipnum");


  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct imgstatparams *p)
{
  char *basename;

  /* Set the p->up.maskname accordingly: */
  gal_fitsarray_file_or_ext_name(p->up.inputname, p->cp.hdu, p->up.masknameset,
                                 &p->up.maskname, p->up.mhdu,
                                 p->up.mhduset, "mask");

  /* Set the names of the output files: */
  if(p->cp.outputset) basename=p->cp.output;
  else                basename=p->up.inputname;
  if(p->histname)
    {
      p->histname=NULL;         /* It wasn't allocated. */
      gal_checkset_automatic_output(basename, "_hist.txt", p->cp.removedirinfo,
                                    p->cp.dontdelete, &p->histname);
    }
  if(p->cfpname)
    {
      p->cfpname=NULL;         /* It wasn't allocated. */
      gal_checkset_automatic_output(basename, "_cfp.txt", p->cp.removedirinfo,
                                    p->cp.dontdelete, &p->cfpname);
    }
  if(p->mhistname)              /* The mode mirror distribution will need */
    {                           /* both a histogram and cfp.              */
      p->mcfpname=p->mhistname=NULL;
      gal_checkset_automatic_output(basename, "_modehist.txt",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->mhistname);
      gal_checkset_automatic_output(basename, "_modecfp.txt",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->mcfpname);
    }
  if(isnan(p->mirror)==0)
    {
      p->mirrorhist=p->mirrorcfp=NULL;
      gal_checkset_automatic_output(basename, "_mirrorhist.txt",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->mirrorhist);
      gal_checkset_automatic_output(basename, "_mirrorcfp.txt",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->mirrorcfp);
    }


  /* If the cumulative frequency plot parameters are to depend on the
     histogram, then make sure that the histogram will be created.*/
  if(p->cfpname && p->histname==NULL)
    {
      if(p->cfpsimhist)
	error(EXIT_FAILURE, 0, "Without a histogram, `--cfpsimhist` is "
              "meaningless.");
      if (p->maxcfpeqmaxhist)
	error(EXIT_FAILURE, 0, "Without a histogram, `--maxcfpeqmaxhist` "
		"is meaningless.\n");
    }

  /* Check if `--maxcfpeqmaxhist` and `--normcfp` are not called
     together: */
  if(p->normcfp && p->maxcfpeqmaxhist)
    error(EXIT_FAILURE, 0, "`--normcfp` and `--maxcfpeqmaxhist` "
          "cannot be called together.\n");

  /* Check if `normhist` and `maxhistone` are not called together: */
  if(p->normhist && p->maxhistone)
    error(EXIT_FAILURE, 0, "`--normhist` and `--histnumbins` cannot be "
	      "called together.\n");
}



















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
preparearrays(struct imgstatparams *p)
{
  float min;
  size_t s0, s1;
  int bitpix, anyblank;
  struct uiparams *up=&p->up;

  /* Read the input and mask arrays: */
  gal_fitsarray_file_to_float(up->inputname, up->maskname, p->cp.hdu, up->mhdu,
                              &p->img, &bitpix, &anyblank, &s0, &s1);
  p->size=s0*s1;

  /* If the minimum value is to be used as a mask then do it: */
  if(p->ignoremin)
    {
      gal_statistics_float_min(p->img, p->size, &min);
      gal_arraymanip_freplace_value(p->img, p->size, min, NAN);
    }

  /* Move all the non-nan elements to the start of the array: */
  gal_arraymanip_no_nans(p->img, &p->size);

  /* Make a sorted array for most of the jobs: */
  gal_arraymanip_float_copy(p->img, p->size, &p->sorted);
  qsort(p->sorted, p->size, sizeof *p->sorted, gal_qsort_float_increasing);

  /* Check the given range: */
  if(p->histname || p->asciihist || p->mhistname || p->mirrorhist)
    {
      if(up->histquantset)
        {
          if(p->histquant>=0.5)
            error(EXIT_FAILURE, 0, "The value to `--histquant' (-Q) must "
                  "Be smaller than 0.5, because it sets the lower limit of "
                  "the value range. The higher limit will be 1-Q.");
          p->histmin=
            p->sorted[gal_statistics_index_from_quantile(p->size,
                                                         p->histquant)];
          p->histmax=
            p->sorted[gal_statistics_index_from_quantile(p->size,
                                                         1 - p->histquant)];
        }
      else
        {
          switch(up->histminset+up->histmaxset)
            {
            case 0:
              p->histmin=p->sorted[0];
              p->histmax=p->sorted[p->size-1];
              break;
            case 1:
              error(EXIT_FAILURE, 0, "The options `--histmin' (-i) and "
                    "`--histmax' (-x) should both be specified. You have "
                    "only given the %s."GAL_STRINGS_HOW_TO_CHECK_VALUES,
                    up->histminset==1 ? "former" : "latter");
              break;
            case 2:
              if(p->histmin>=p->histmax)
                error(EXIT_FAILURE, 0, "The value to `--histmin' (-i) (%f) "
                      "is larger or equal to that of `--histmax' (-x) (%f)."
                      GAL_STRINGS_HOW_TO_CHECK_VALUES, p->histmin, p->histmax);
              if(p->histmin>p->sorted[p->size-1] || p->histmax<p->sorted[0])
                error(EXIT_FAILURE, 0, "The range of data is %.5f to %.5f. "
                      "However, you have set `--histmin' (-i) and "
                      "`--histmax' (-x) to %.5f and %.5f respectively. "
                      "They do not overlap!"GAL_STRINGS_HOW_TO_CHECK_VALUES,
                      p->sorted[0], p->sorted[p->size-1], p->histmin,
                      p->histmax);
              break;
            default:
              error(EXIT_FAILURE, 0, "A bug! Please contact us at "
                    PACKAGE_BUGREPORT" So we can solve the problem. the "
                    "value of up->histminset+up->histmaxset is not 0, 1 or "
                    "2.");
            }
        }
    }
  else                          /* For the ascii histogram. */
    {
      p->histmin=p->sorted[0];
      p->histmax=p->sorted[p->size-1];
    }

  if(p->cfpname && p->cfpsimhist==0)
    {
      if(up->cfpquantset)
        {
          if(p->cfpquant>=0.5)
            error(EXIT_FAILURE, 0, "The value to `--cfpquant' (-U) must "
                  "Be smaller than 0.5, because it sets the lower limit of "
                  "the value range. The higher limit will be 1-U.");
          p->cfpmin=p->sorted[gal_statistics_index_from_quantile(p->size,
                                                                 p->cfpquant)];
          p->cfpmax=
            p->sorted[gal_statistics_index_from_quantile(p->size,
                                                         1 - p->cfpquant)];
        }
      else
        {
          switch(up->cfpminset+up->cfpmaxset)
            {
            case 0:
              p->cfpmin=p->sorted[0];
              p->cfpmax=p->sorted[p->size-1];
              break;
            case 1:
              error(EXIT_FAILURE, 0, "The options `--cfpmin' (-a) and "
                    "`--cfpmax' (-b) should both be specified. You have "
                    "only given the %s."GAL_STRINGS_HOW_TO_CHECK_VALUES,
                    up->cfpminset==1 ? "former" : "latter");
              break;
            case 2:
              if(p->cfpmin>p->cfpmax)
                error(EXIT_FAILURE, 0, "The value to `--cfpmin' (-a) (%.f) "
                      "is larger than that of `--cfpmax' (-b) (%f)."
                      GAL_STRINGS_HOW_TO_CHECK_VALUES, p->cfpmin, p->cfpmax);
              if(p->cfpmin>p->sorted[p->size-1] || p->cfpmax<p->sorted[0])
                error(EXIT_FAILURE, 0, "The range of data is %.5f to %.5f. "
                      "However, you have set `--cfpmin' (-a) and "
                      "`--cfpmax' (-b) to %.5f and %.5f respectively. "
                      "They do not overlap!"GAL_STRINGS_HOW_TO_CHECK_VALUES,
                      p->sorted[0], p->sorted[p->size-1], p->cfpmin,
                      p->cfpmax);
              break;
            default:
              error(EXIT_FAILURE, 0, "A bug! Please contact us at "
                    PACKAGE_BUGREPORT" So we can solve the problem. the "
                    "value of up->cfpminset+up->cfpmaxset is not 0, 1 or "
                    "2.");
            }
        }
    }
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct imgstatparams *p)
{
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;

  p->asciihist      = 1;
  p->sigclip        = 1;
  p->mirror         = NAN;
  p->onebinvalue    = NAN;
  p->histname=p->cfpname="a";   /* Will be set later, just a sign that */
                                /* they should be output.              */
  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "Parsing arguments");

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    GAL_CONFIGFILES_REPORT_PARAMETERS_SET;

  /* Do a sanity check. */
  sanitycheck(p);

  /* Make the array of input images. */
  preparearrays(p);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf("  - Input read: %s (hdu: %s)\n", p->up.inputname, p->cp.hdu);
      if(p->up.maskname)
        printf("  - Mask read: %s (hdu: %s)\n", p->up.maskname, p->up.mhdu);
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct imgstatparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->img);
  free(p->sorted);
  free(p->cp.hdu);
  free(p->cfpname);
  free(p->histname);
  free(p->mcfpname);
  free(p->mhistname);
  free(p->cp.output);
  if(p->up.masknameallocated) free(p->up.maskname);

  /* Print the final message. */
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in: ", 0);
}
