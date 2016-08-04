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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include <nproc.h>               /* From Gnulib.                   */

#include <gnuastro/fits.h>
#include <gnuastro/timing.h>     /* Includes time.h and sys/time.h */
#include <gnuastro/checkset.h>
#include <gnuastro/statistics.h>
#include <gnuastro/arraymanip.h>
#include <gnuastro/configfiles.h>

#include "main.h"

#include "ui.h"                  /* Needs main.h.                  */
#include "args.h"                /* Needs main.h, includes argp.h. */


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
  struct gal_commonparams *cp=&p->cp;
  int zeroorone, spatialset=0, frequencyset=0;
  char key='a';        /* Not used, just a place holder. */

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
        gal_checkset_allocate_copy_set(value, &up->maskname,
                                       &up->masknameset);
      else if (strcmp(name, "mhdu")==0)
        gal_checkset_allocate_copy_set(value, &up->mhdu, &up->mhduset);
      else if (strcmp(name, "kernel")==0)
        gal_checkset_allocate_copy_set(value, &up->kernelname,
                                       &up->kernelnameset);
      else if (strcmp(name, "khdu")==0)
        gal_checkset_allocate_copy_set(value, &up->khdu, &up->khduset);



      /* Outputs: */
      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output,
                                       &cp->outputset);




      /* Mesh grid: */
      else if(strcmp(name, "meshsize")==0)
        {
          if(up->meshsizeset) continue;
          gal_checkset_sizet_l_zero(value, &p->mp.meshsize, name, key,
                                    SPACK, filename, lineno);
          up->meshsizeset=1;
        }
      else if(strcmp(name, "nch1")==0)
        {
          if(up->nch1set) continue;
          gal_checkset_sizet_l_zero(value, &p->mp.nch1, name, key, SPACK,
                                    filename, lineno);
          up->nch1set=1;
        }
      else if(strcmp(name, "nch2")==0)
        {
          if(up->nch2set) continue;
          gal_checkset_sizet_l_zero(value, &p->mp.nch2, name, key, SPACK,
                                    filename, lineno);
          up->nch2set=1;
        }
      else if(strcmp(name, "lastmeshfrac")==0)
        {
          if(up->lastmeshfracset) continue;
          gal_checkset_float_l_0_s_1(value, &p->mp.lastmeshfrac, name,
                                     key, SPACK, filename, lineno);
          up->lastmeshfracset=1;
        }
      else if(strcmp(name, "fullconvolution")==0)
        {
          if(up->fullconvolutionset) continue;
          gal_checkset_int_zero_or_one(value, &p->mp.fullconvolution, name,
                                       key, SPACK, filename, lineno);
          up->fullconvolutionset=1;
        }



      /* Operating modes: */
      else if(strcmp(name, "spatial")==0)
        {
          gal_checkset_int_zero_or_one(value, &zeroorone, name, key,
                                       SPACK, filename, lineno);
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
          gal_checkset_int_zero_or_one(value, &zeroorone, name, key,
                                       SPACK, filename, lineno);
          if(zeroorone)
            {
              frequencyset=1;
              if(spatialset)
                error_at_line(EXIT_FAILURE, 0, filename, lineno,
                              "Spatial and frequency modes cannot be "
                              "called together. It is ambiguous.");
              if(up->frequencyset==0)
                {
                  p->spatial=0;
                  p->frequency=1;
                  up->spatialset=up->frequencyset=1;
                }
            }
        }
      else if(strcmp(name, "makekernel")==0)
        {
          if(up->makekernelset) continue;
          gal_checkset_int_el_zero(value, &p->makekernel, name, key,
                                   SPACK, filename, lineno);
          up->makekernelset=1;
        }

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
printvalues(FILE *fp, struct convolveparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_mesh_params *mp=&p->mp;
  struct gal_commonparams *cp=&p->cp;

  fprintf(fp, "\n# Input:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);
  if(up->masknameset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mask", up->maskname);
  if(up->mhduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mhdu", up->mhdu);
  if(up->kernelnameset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("kernel", up->kernelname);
  if(up->khduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("khdu", up->khdu);



  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);


  fprintf(fp, "\n# Mesh grid:\n");
  if(up->meshsizeset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "meshsize", mp->meshsize);
  if(up->nch1set)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "nch1", mp->nch1);
  if(up->nch2set)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "nch2", mp->nch2);
  if(up->lastmeshfracset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "lastmeshfrac", mp->lastmeshfrac);
  if(up->fullconvolutionset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "fullconvolution",
            mp->fullconvolution);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating modes:\n");
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
  if(up->spatialset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "spatial", p->spatial);
  if(up->frequencyset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "frequency", p->frequency);
  if(up->makekernelset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "makekernel", p->makekernel);
}





void
checkifset(struct convolveparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;

  int intro=0;


  /* Input: */
  if(cp->hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");
  if(up->kernelnameset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("kernel");
  if(up->khduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("khdu");


  /* Mesh grid: */
  if(up->meshsizeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("meshsize");
  if(up->nch1set==0)
    GAL_CONFIGFILES_REPORT_NOTSET("nch1");
  if(up->nch2set==0)
    GAL_CONFIGFILES_REPORT_NOTSET("nch2");
  if(up->lastmeshfracset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("lastmeshfrac");
  if(up->fullconvolutionset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("fullconvolution");


  /* Operating mode: */
  if(up->spatialset==0 && up->frequencyset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("spatial or frequency");
  if(up->makekernelset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("makekernel");


  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}



















/**************************************************************/
/************         Prepare the arrays          *************/
/**************************************************************/
void
sanitycheck(struct convolveparams *p)
{
  char *outsuffix = p->makekernel ? "_kernel.fits" : "_convolved.fits";

  /* Make sure the input file exists. */
  gal_checkset_check_file(p->up.inputname);

  /* Set maskname accordingly: */
  gal_fits_file_or_ext_name(p->up.inputname, p->cp.hdu, p->up.masknameset,
                                 &p->up.maskname, p->up.mhdu,
                                 p->up.mhduset, "mask");

  /* Check the output file name: */
  if(p->cp.outputset)
    {
      if( gal_checkset_dir_0_file_1(p->cp.output, p->cp.dontdelete) == 0 )
        error(EXIT_FAILURE, 0, "your output name (%s) is a directory",
              p->cp.output);
    }
  else
    {
      gal_checkset_automatic_output(p->up.inputname, outsuffix,
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->cp.output);
      p->cp.outputset=1;
    }
  if(p->frequency && p->viewfreqsteps)
    gal_checkset_automatic_output(p->up.inputname, "_freqsteps.fits",
                    p->cp.removedirinfo, p->cp.dontdelete,
                    &p->up.freqstepsname);

  /* Check output names: */
  if(p->meshname)
    {
      p->meshname=NULL;           /* Was not allocated before!  */
      gal_checkset_automatic_output(p->up.inputname, "_mesh.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->meshname);
    }

  /* makekernel can only operate in frequency mode: */
  if(p->makekernel && p->spatial)
    error(EXIT_FAILURE, 0, "`--makekernel' (`-m') can only be defined "
          "in frequency domain operations, not spatial domain. You can "
          "either explicitly call for frequency mode on the command line "
          "or change your nearest configuration file so Convolve defaults "
          "to the frequency domain. Please see the Gnuastro for more "
          "information. On the command line, the manual can be seen with "
          "the `info gnuastro' command");
}




















/**************************************************************/
/************         Prepare the arrays          *************/
/**************************************************************/
void
preparearrays(struct convolveparams *p)
{
  size_t i, size;
  int bitpix, anyblank;
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;
  float *f, *fp, tmp, *kernel, sum;

  /* First read the input image: */
  gal_fits_file_to_float(up->inputname, up->maskname, cp->hdu, up->mhdu,
                              &p->input, &bitpix, &p->anyblank,
                              &p->is0, &p->is1);
  gal_fits_read_wcs(up->inputname, cp->hdu, 0, 0, &p->nwcs, &p->wcs);
  if(p->frequency && p->anyblank)
    fprintf(stderr, "\n----------------------------------------\n"
            "######## %s WARNING ########\n"
            "There are blank (masked) pixels in %s (hdu: %s) and you "
            "have asked for frequency domain convolution.%s All the "
            "convolved pixels will become blank. Only spatial domain "
            "convolution can account for blank (masked) pixels in the "
            "input data.\n"
            "----------------------------------------\n\n",
            SPACK_NAME, up->inputname, cp->hdu, up->maskname ?
            "" : " Even though you have not provided any mask image, "
            "these are the blank pixels in the input image, see the `Blank "
            "pixels' section of the Gnuastro manual for more information.");


  /* Read the file specified by --kernel. If makekernel is specified,
     then this is actually the low resolution image. */
  if(p->makekernel)
    {
      /* Read in the kernel array: */
      gal_fits_file_to_float(up->kernelname, NULL, up->khdu, NULL,
                                  &p->kernel, &bitpix, &anyblank,
                                  &p->ks0, &p->ks1);
      if(p->ks0!=p->is0 || p->ks1!=p->is1)
        error(EXIT_FAILURE, 0, "with the `--makekernel' (`-m') option, "
              "the input image and the image specified with the kernel "
              "option should have the same size. The lower resolution "
              "input image (%s) has %lux%lu pixels while the higher "
              "resolution image (%s) specified with the kernel option has "
              "%lux%lu pixels", up->inputname, p->is1, p->is0,
              up->kernelname, p->ks1, p->ks0);

      /* Divide both images by their sum so their lowest frequency
         becomes 1 (and their division would be meaningful!).*/
      size=p->is0*p->is1;
      sum=gal_statistics_float_sum(p->input, size);
      gal_arraymanip_fmultip_const(p->input, size, 1/sum);
      sum=gal_statistics_float_sum(p->kernel, size);
      gal_arraymanip_fmultip_const(p->kernel, size, 1/sum);
    }

  /* Read the kernel. If there is anything particular to Convolve, then
     don't use the standard kernel reading function in fits.c. Otherwise
     just use the same one that all programs use. The standard one is
     faster because it mixes the NaN conversion and also the normalization
     into one loop.*/
  else
    {
      if(p->kernelnorm==0 || p->kernelflip==0)
        {
          /* Read in the kernel array: */
          gal_fits_file_to_float(up->kernelname, NULL, up->khdu, NULL,
                                      &p->kernel, &bitpix, &anyblank,
                                      &p->ks0, &p->ks1);
          size=p->ks0*p->ks1;
          kernel=p->kernel;

          if(p->ks0%2==0 || p->ks1%2==0)
            error(EXIT_FAILURE, 0, "the kernel image has to have an odd "
                  "number of pixels on both sides (there has to be on pixel "
                  "in the center). %s (hdu: %s) is %lu by %lu",
                  p->up.kernelname, p->up.khdu, p->ks1, p->ks0);

          /* Convert all the NaN pixels to zero if the kernel contains
             blank pixels. */
          if(anyblank)
            { fp=(f=kernel)+size; do if(isnan(*f)) *f=0.0f; while(++f<fp); }

          /* Normalize the kernel: */
          if(p->kernelnorm)
            {
              sum=gal_statistics_float_sum(kernel, size);
              gal_arraymanip_fmultip_const(kernel, size, 1/sum);
            }

          /* Flip the kernel: */
          if(p->spatial && p->kernelflip)
            for(i=0;i<size/2;++i)
              {
                tmp=kernel[i];
                kernel[i]=kernel[size-i-1];
                kernel[size-i-1]=tmp;
              }
        }
      else
        gal_fits_prep_float_kernel(up->kernelname, up->khdu, &p->kernel,
                                        &p->ks0, &p->ks1);
    }
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct convolveparams *p)
{
  struct gal_commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = num_processors(NPROC_CURRENT);
  cp->removedirinfo = 1;

  /* Set non-zero options: */
  p->kernelflip     = 1;
  p->kernelnorm     = 1;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    GAL_CONFIGFILES_REPORT_PARAMETERS_SET;

  /* Do a sanity check, then remove the possibly existing log file
     created by gal_txtarray_txt_to_array. */
  sanitycheck(p);

  /* Prepare the necessary arrays: */
  preparearrays(p);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf("  - Using %lu CPU threads.\n", p->cp.numthreads);
      printf("  - Input: %s (hdu: %s)\n", p->up.inputname, p->cp.hdu);
      if(p->up.maskname)
        printf("  - Mask: %s (hdu: %s)\n", p->up.maskname, p->up.mhdu);
      printf("  - Kernel: %s (hdu: %s)\n", p->up.kernelname, p->up.khdu);
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
  free(p->meshname);
  free(p->cp.output);
  wcsvfree(&p->nwcs, &p->wcs);

  /* Free the mask image name. Note that p->up.inputname was not
     allocated, but given to the program by the operating system. */
  if(p->up.maskname && p->up.maskname!=p->up.inputname)
    free(p->up.maskname);

  /* Print the final message. */
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in: ", 0);
}
