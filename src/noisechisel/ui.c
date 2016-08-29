/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <nproc.h>               /* From Gnulib.                   */

#include <gnuastro/fits.h>
#include <gnuastro/timing.h>     /* includes time.h and sys/time.h */
#include <gnuastro/checkset.h>
#include <gnuastro/txtarray.h>
#include <gnuastro/commonargs.h>
#include <gnuastro/arraymanip.h>
#include <gnuastro/statistics.h>
#include <gnuastro/configfiles.h>

#include "main.h"

#include "ui.h"                  /* Needs main.h                   */
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
readconfig(char *filename, struct noisechiselparams *p)
{
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;
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

      else if(strcmp(name, "mask")==0)
        gal_checkset_allocate_copy_set(value, &up->maskname,
                                       &up->masknameset);

      else if(strcmp(name, "mhdu")==0)
        gal_checkset_allocate_copy_set(value, &up->mhdu, &up->mhduset);

      else if(strcmp(name, "kernel")==0)
        gal_checkset_allocate_copy_set(value, &up->kernelname,
                                       &up->kernelnameset);

      else if(strcmp(name, "khdu")==0)
        gal_checkset_allocate_copy_set(value, &up->khdu, &up->khduset);

      else if(strcmp(name, "skysubtracted")==0)
        {
          if(up->skysubtractedset) continue;
          gal_checkset_int_zero_or_one(value, &p->skysubtracted, name,
                                       key, SPACK, filename, lineno);
          up->skysubtractedset=1;
        }
      else if(strcmp(name, "minbfrac")==0)
        {
          if(up->minbfracset) continue;
          gal_checkset_float_l_0_s_1(value, &p->minbfrac, name,
                                     key, SPACK, filename, lineno);
          up->minbfracset=1;
        }
      else if(strcmp(name, "minnumfalse")==0)
        {
          if(up->minnumfalseset) continue;
          gal_checkset_sizet_l_zero(value, &p->minnumfalse, name,
                                    key, SPACK, filename, lineno);
          up->minnumfalseset=1;
        }



      /* Outputs */
      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output, &cp->outputset);

      else if(strcmp(name, "grownclumps")==0)
        {
          if(up->grownclumpsset) continue;
          gal_checkset_int_zero_or_one(value, &p->grownclumps, name,
                                       key, SPACK, filename, lineno);
          up->grownclumpsset=1;
        }


      /* Mesh grid: */
      else if(strcmp(name, "smeshsize")==0)
        {
          if(up->smeshsizeset) continue;
          gal_checkset_sizet_l_zero(value, &p->smp.meshsize, name,
                                    key, SPACK, filename, lineno);
          up->smeshsizeset=1;
        }
      else if(strcmp(name, "lmeshsize")==0)
        {
          if(up->lmeshsizeset) continue;
          gal_checkset_sizet_l_zero(value, &p->lmp.meshsize, name,
                                    key, SPACK, filename, lineno);
          up->lmeshsizeset=1;
        }
      else if(strcmp(name, "nch1")==0)
        {
          if(up->nch1set) continue;
          gal_checkset_sizet_l_zero(value, &p->smp.nch1, name, key,
                                    SPACK, filename, lineno);
          up->nch1set=1;
        }
      else if(strcmp(name, "nch2")==0)
        {
          if(up->nch2set) continue;
          gal_checkset_sizet_l_zero(value, &p->smp.nch2, name, key,
                                    SPACK, filename, lineno);
          up->nch2set=1;
        }
      else if(strcmp(name, "lastmeshfrac")==0)
        {
          if(up->lastmeshfracset) continue;
          gal_checkset_float_l_0_s_1(value, &p->smp.lastmeshfrac, name,
                                     key, SPACK, filename, lineno);
          up->lastmeshfracset=1;
        }
      else if(strcmp(name, "mirrordist")==0)
        {
          if(up->mirrordistset) continue;
          gal_checkset_float_l_0(value, &p->smp.mirrordist, name,
                                 key, SPACK, filename, lineno);
          up->mirrordistset=1;
        }
      else if(strcmp(name, "minmodeq")==0)
        {
          if(up->minmodeqset) continue;
          gal_checkset_float_l_0_s_1(value, &p->smp.minmodeq, name,
                                     key, SPACK, filename, lineno);
          up->minmodeqset=1;
        }
      else if(strcmp(name, "numnearest")==0)
        {
          if(up->numnearestset) continue;
          gal_checkset_sizet_l_zero(value, &p->smp.numnearest, name,
                                    key, SPACK, filename, lineno);
          up->numnearestset=1;
        }
      else if(strcmp(name, "smoothwidth")==0)
        {
          if(up->smoothwidthset) continue;
          gal_checkset_sizet_p_odd(value, &p->smp.smoothwidth, name,
                                   key, SPACK, filename, lineno);
          up->smoothwidthset=1;
        }
      else if(strcmp(name, "fullconvolution")==0)
        {
          if(up->fullconvolutionset) continue;
          gal_checkset_int_zero_or_one(value, &p->smp.fullconvolution,
                                       name, key, SPACK, filename, lineno);
          up->fullconvolutionset=1;
        }
      else if(strcmp(name, "fullinterpolation")==0)
        {
          if(up->fullinterpolationset) continue;
          gal_checkset_int_zero_or_one(value, &p->smp.fullinterpolation,
                                       name, key, SPACK, filename, lineno);
          up->fullinterpolationset=1;
        }
      else if(strcmp(name, "fullsmooth")==0)
        {
          if(up->fullsmoothset) continue;
          gal_checkset_int_zero_or_one(value, &p->smp.fullsmooth, name, key,
                                       SPACK, filename, lineno);
          up->fullsmoothset=1;
        }


      /* Detection: */
      else if(strcmp(name, "qthresh")==0)
        {
          if(up->qthreshset) continue;
          gal_checkset_float_l_0_s_1(value, &p->qthresh, name, key,
                                     SPACK, filename, lineno);
          up->qthreshset=1;
        }
      else if(strcmp(name, "erode")==0)
        {
          if(up->erodeset) continue;
          gal_checkset_sizet_el_zero(value, &p->erode, name, key,
                                     SPACK, filename, lineno);
          up->erodeset=1;
        }
      else if(strcmp(name, "erodengb")==0)
        {
          if(up->erodengbset) continue;
          gal_checkset_int_4_or_8(value, &p->erodengb, name, key,
                                  SPACK, filename, lineno);
          up->erodengbset=1;
        }
      else if(strcmp(name, "noerodequant")==0)
        {
          if(up->noerodequantset) continue;
          gal_checkset_float_l_0_s_1(value, &p->noerodequant, name, key,
                                     SPACK, filename, lineno);
          up->noerodequantset=1;
        }
      else if(strcmp(name, "opening")==0)
        {
          if(up->openingset) continue;
          gal_checkset_sizet_el_zero(value, &p->opening, name, key,
                                     SPACK, filename, lineno);
          up->openingset=1;
        }
      else if(strcmp(name, "openingngb")==0)
        {
          if(up->openingngbset) continue;
          gal_checkset_int_4_or_8(value, &p->openingngb, name,
                                  key, SPACK, filename, lineno);
          up->openingngbset=1;
        }
      else if(strcmp(name, "sigclipmultip")==0)
        {
          if(up->sigclipmultipset) continue;
          gal_checkset_float_l_0(value, &p->sigclipmultip, name,
                                 key, SPACK, filename, lineno);
          up->sigclipmultipset=1;
        }
      else if(strcmp(name, "sigcliptolerance")==0)
        {
          if(up->sigcliptoleranceset) continue;
          gal_checkset_float_l_0_s_1(value, &p->sigcliptolerance, name,
                                     key, SPACK, filename, lineno);
          up->sigcliptoleranceset=1;
        }
      else if(strcmp(name, "dthresh")==0)
        {
          if(up->dthreshset) continue;
          gal_checkset_any_float(value, &p->dthresh, name, key,
                                 SPACK, filename, lineno);
          up->dthreshset=1;
        }
      else if(strcmp(name, "detsnminarea")==0)
        {
          if(up->detsnminareaset) continue;
          gal_checkset_sizet_l_zero(value, &p->detsnminarea, name,
                                    key, SPACK, filename, lineno);
          up->detsnminareaset=1;
        }
      else if(strcmp(name, "detsnhistnbins")==0)
        {
          if(up->detsnhistnbinsset) continue;
          gal_checkset_sizet_el_zero(value, &p->detsnhistnbins, name,
                                     key, SPACK, filename, lineno);
          up->detsnhistnbinsset=1;
        }
      else if(strcmp(name, "detquant")==0)
        {
          if(up->detquantset) continue;
          gal_checkset_float_l_0_s_1(value, &p->detquant, name, key,
                                     SPACK, filename, lineno);
          up->detquantset=1;
        }
      else if(strcmp(name, "dilate")==0)
        {
          if(up->dilateset) continue;
          gal_checkset_sizet_el_zero(value, &p->dilate, name, key,
                                     SPACK, filename, lineno);
          up->dilateset=1;
        }


      /* Segmentation: */
      else if(strcmp(name, "segsnminarea")==0)
        {
          if(up->segsnminareaset) continue;
          gal_checkset_sizet_l_zero(value, &p->segsnminarea, name,
                                    key, SPACK, filename, lineno);
          up->segsnminareaset=1;
        }
      else if(strcmp(name, "keepmaxnearriver")==0)
        {
          if(up->keepmaxnearriverset) continue;
          gal_checkset_int_zero_or_one(value, &p->keepmaxnearriver, name,
                                       key, SPACK, filename, lineno);
          up->keepmaxnearriverset=1;
        }
      else if(strcmp(name, "segquant")==0)
        {
          if(up->segquantset) continue;
          gal_checkset_float_l_0_s_1(value, &p->segquant, name, key,
                                     SPACK, filename, lineno);
          up->segquantset=1;
        }
      else if(strcmp(name, "clumpsnhistnbins")==0)
        {
          if(up->clumpsnhistnbinsset) continue;
          gal_checkset_sizet_el_zero(value, &p->clumpsnhistnbins, name,
                                     key, SPACK,filename, lineno);
          up->clumpsnhistnbinsset=1;
        }
      else if(strcmp(name, "gthresh")==0)
        {
          if(up->gthreshset) continue;
          gal_checkset_any_float(value, &p->gthresh, name, key, SPACK,
                                 filename, lineno);
          up->gthreshset=1;
        }
      else if(strcmp(name, "minriverlength")==0)
        {
          if(up->minriverlengthset) continue;
          gal_checkset_sizet_l_zero(value, &p->minriverlength, name,
                                    key, SPACK, filename, lineno);
          up->minriverlengthset=1;
        }
      else if(strcmp(name, "objbordersn")==0)
        {
          if(up->objbordersnset) continue;
          gal_checkset_float_l_0(value, &p->objbordersn, name, key,
                                 SPACK, filename, lineno);
          up->objbordersnset=1;
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
printvalues(FILE *fp, struct noisechiselparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;
  struct gal_mesh_params *smp=&p->smp, *lmp=&p->lmp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);
  if(up->masknameset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mask", up->maskname);
  if(up->mhdu)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mhdu", up->mhdu);
  if(up->kernelnameset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("kernel", up->kernelname);
  if(up->khdu)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("khdu", up->khdu);
  if(up->skysubtractedset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "skysubtracted", p->skysubtracted);
  if(up->minbfracset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "minbfrac", p->minbfrac);
  if(up->minnumfalseset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "minnumfalse", p->minnumfalse);


  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);
  if(up->grownclumpsset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "grownclumps", p->grownclumps);


  fprintf(fp, "\n# Mesh grid:\n");
  if(up->smeshsizeset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "smeshsize", smp->meshsize);
  if(up->lmeshsizeset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "lmeshsize", lmp->meshsize);
  if(up->nch1set)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "nch1", smp->nch1);
  if(up->nch2set)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "nch2", smp->nch2);
  if(up->lastmeshfracset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "lastmeshfrac", smp->lastmeshfrac);
  if(up->mirrordistset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "mirrordist", smp->mirrordist);
  if(up->minmodeqset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "minmodeq", smp->minmodeq);
  if(up->numnearestset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "numnearest", smp->numnearest);
  if(up->smoothwidthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "smoothwidth", smp->smoothwidth);
  if(up->fullconvolutionset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "fullconvolution",
            smp->fullconvolution);
  if(up->fullinterpolationset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "fullinterpolation",
            smp->fullinterpolation);
  if(up->fullsmoothset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "fullsmooth", smp->fullsmooth);


  fprintf(fp, "\n# Detection:\n");
  if(up->qthreshset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "qthresh", p->qthresh);
  if(up->erodeset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "erode", p->erode);
  if(up->erodengbset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "erodengb", p->erodengb);
  if(up->noerodequantset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "noerodequant", p->noerodequant);
  if(up->openingset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "opening", p->opening);
  if(up->openingngbset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "openingngb", p->openingngb);
  if(up->sigclipmultipset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "sigclipmultip", p->sigclipmultip);
  if(up->sigcliptoleranceset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "sigcliptolerance",
            p->sigcliptolerance);
  if(up->dthreshset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "dthresh", p->dthresh);
  if(up->detsnminareaset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "detsnminarea", p->detsnminarea);
  if(up->detsnhistnbinsset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "detsnhistnbins", p->detsnhistnbins);
  if(up->detquantset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "detquant", p->detquant);
  if(up->dilateset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "dilate", p->dilate);


  fprintf(fp, "\n# Segmentation:\n");
  if(up->segsnminareaset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "segsnminarea", p->segsnminarea);
  if(up->keepmaxnearriverset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "keepmaxnearriver", p->keepmaxnearriver);
  if(up->segquantset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "segquant", p->segquant);
  if(up->clumpsnhistnbinsset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "clumpsnhistnbins", p->clumpsnhistnbins);
  if(up->gthreshset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "gthresh", p->gthresh);
  if(up->minriverlengthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "minriverlength", p->minriverlength);
  if(up->objbordersnset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "objbordersn", p->objbordersn);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct noisechiselparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");
  if(up->khduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("khdu");
  if(up->skysubtractedset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("skysubtracted");
  if(up->minbfracset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("minbfrac");
  if(up->minnumfalseset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("minnumfalse");

  /* Output */
  if(up->grownclumpsset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("grownclumps");

  /* Mesh grid: */
  if(up->smeshsizeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("smeshsize");
  if(up->lmeshsizeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("lmeshsize");
  if(up->nch1set==0)
    GAL_CONFIGFILES_REPORT_NOTSET("nch1");
  if(up->nch2set==0)
    GAL_CONFIGFILES_REPORT_NOTSET("nch2");
  if(up->lastmeshfracset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("lastmeshfrac");
  if(up->mirrordistset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("mirrordist");
  if(up->minmodeqset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("minmodeq");
  if(up->numnearestset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("numnearest");
  if(up->smoothwidthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("smoothwidth");
  if(up->fullconvolutionset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("fullconvolution");
  if(up->fullinterpolationset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("fullinterpolation");
  if(up->fullsmoothset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("fullsmooth");

  /* Detection: */
  if(up->qthreshset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("qthresh");
  if(up->erodeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("erode");
  if(up->erodengbset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("erodengb");
  if(up->noerodequantset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("noerodequant");
  if(up->openingset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("opening");
  if(up->openingngbset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("openingngb");
  if(up->sigclipmultipset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("sigclipmultip");
  if(up->sigcliptoleranceset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("sigcliptolerance");
  if(up->dthreshset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("dthresh");
  if(up->detsnminareaset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("detsnminarea");
  if(up->detsnhistnbinsset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("detsnhistnbins");
  if(up->detquantset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("detquant");
  if(up->dilateset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("dilate");

  /* Segmentation: */
  if(up->segsnminareaset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("segsnminarea");
  if(up->keepmaxnearriverset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("keepmaxnearriver");
  if(up->segquantset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("segquant");
  if(up->clumpsnhistnbinsset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("clumpsnhistnbins");
  if(up->gthreshset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("gthresh");
  if(up->minriverlengthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("minriverlength");
  if(up->objbordersnset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("objbordersn");

  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct noisechiselparams *p)
{
  struct gal_mesh_params *smp=&p->smp;

  /* Make sure the input file exists. */
  gal_checkset_check_file(p->up.inputname);

  /* Make sure that the noerode quantile is larger than qthresh. */
  if( p->noerodequant <= p->qthresh)
    error(EXIT_FAILURE, 0, "The quantile for no erosion (`--noerodequant') "
          "must be larger than the base quantile threshold (`--qthresh', "
          "or `-t'). You have provided %.4f and %.4f for the former and "
          "latter, respectively.", p->noerodequant, p->qthresh);

  /* Set the maskname and mask hdu accordingly: */
  gal_fits_file_or_ext_name(p->up.inputname, p->cp.hdu, p->up.masknameset,
                                 &p->up.maskname, p->up.mhdu, p->up.mhduset,
                                 "mask");

  /* Set the output name: */
  if(p->cp.output)
    {
      gal_checkset_check_remove_file(p->cp.output, p->cp.dontdelete);

      /* When the output name is given (possibly with directory
         information), the user certainly wants the directory
         information, if they have bothered to include it. */
      p->cp.removedirinfo=0;

    }
  else
    gal_checkset_automatic_output(p->up.inputname, "_labeled.fits",
                                  p->cp.removedirinfo, p->cp.dontdelete,
                                  &p->cp.output);

  /* Set the check image names: */
  if(p->meshname)
    {
      p->meshname=NULL;         /* Was not allocated before!  */
      gal_checkset_automatic_output(p->cp.output, "_meshs.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->meshname);
    }
  if(p->threshname)
    {
      p->threshname=NULL;
      gal_checkset_automatic_output(p->cp.output, "_thresh.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->threshname);
    }
  if(p->detectionname)
    {
      p->detectionname=NULL;
      gal_checkset_automatic_output(p->cp.output, "_det.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->detectionname);
    }
  if(p->detectionskyname)
    {
      p->detectionskyname=NULL;
      gal_checkset_automatic_output(p->cp.output, "_detsky.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->detectionskyname);
    }
  if(p->detsnhistnbins)
    {
      p->detectionsnhist=NULL;
      gal_checkset_automatic_output(p->cp.output, "_detsn.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->detectionsnhist);
    }
  if(p->skyname)
    {
      p->skyname=NULL;
      gal_checkset_automatic_output(p->cp.output, "_sky.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->skyname);
    }
  if(p->segmentationname)
    {
      p->segmentationname=NULL;
      gal_checkset_automatic_output(p->cp.output, "_seg.fits",
                                    p->cp.removedirinfo,
                                    p->cp.dontdelete, &p->segmentationname);
    }
  if(p->clumpsnhistnbins)
    {
      p->clumpsnhist=NULL;
      gal_checkset_automatic_output(p->cp.output, "_clumpsn.txt",
                                    p->cp.removedirinfo,
                                    p->cp.dontdelete, &p->clumpsnhist);
    }
  if(p->maskdetname)
    {
      p->maskdetname=NULL;
      gal_checkset_automatic_output(p->cp.output, "_maskdet.fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &p->maskdetname);
    }

  /* Other checks: */
  if(smp->numnearest<GAL_MESH_MIN_ACCEPTABLE_NEAREST)
    error(EXIT_FAILURE, 0, "the smallest possible number for `--numnearest' "
          "(`-n') is %d. You have asked for: %lu",
          GAL_MESH_MIN_ACCEPTABLE_NEAREST, smp->numnearest);
}



















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* The default PSF. It was created by saving the following commands in
   a script and running it. The crop is because the first and last
   rows of all PSFs made by MakeProfiles is blank (zero). You can keep
   the spaces when copying and pasting ;-). Just make it executable
   and run it.

   set -o errexit           # Stop if a program returns false.
   echo "0    0.0    0.0   2   2   0   0   1   1   5" > tmp.txt
   export GSL_RNG_TYPE=ranlxs2
   export GSL_RNG_SEED=1
   astmkprof tmp.txt --oversample=1 --envseed --numrandom=10000 \
             --tolerance=0.01
   astimgcrop 0.fits --section=2:*,2:* --zeroisnotblank --output=fwhm2.fits
   astconvertt fwhm2.fits --output=fwhm2.txt
   rm 0.fits tmp.fits *.log tmp.txt
*/
size_t defaultkernel_s0=11;
size_t defaultkernel_s1=11;
float defaultkernel[121]=
  {
    0, 0, 0, 0, 0, 2.58073e-08, 0, 0, 0, 0, 0,

    0, 0, 2.90237e-08, 6.79851e-07, 4.4435e-06, 8.31499e-06,
    4.50166e-06, 6.97185e-07, 3.00904e-08, 0, 0,

    0, 2.87873e-08, 2.48435e-06, 5.81339e-05, 0.000379508, 0.000709334,
    0.000383714, 5.94125e-05, 2.56498e-06, 3.00032e-08, 0,

    0, 6.70501e-07, 5.77826e-05, 0.00134992, 0.00879665, 0.0164126,
    0.00886609, 0.00137174, 5.92134e-05, 6.92853e-07, 0,

    0, 4.3798e-06, 0.000376616, 0.00877689, 0.0570404, 0.106142, 0.0572108,
    0.00883846, 0.000381257, 4.46059e-06, 0,

    2.54661e-08, 8.24845e-06, 0.00070725, 0.0164287, 0.10639, 0.19727,
    0.106003, 0.0163402, 0.000703951, 8.23152e-06, 2.55057e-08,

    0, 4.5229e-06, 0.000386632, 0.00894947, 0.0577282, 0.106614, 0.0570877,
    0.00877699, 0.000377496, 4.41036e-06, 0,

    0, 7.1169e-07, 6.0678e-05, 0.00140013, 0.00899917, 0.0165582, 0.00883658,
    0.00135509, 5.81823e-05, 6.79067e-07, 0,

    0, 3.12002e-08, 2.65502e-06, 6.11192e-05, 0.000391739, 0.000718637,
    0.000382453, 5.85194e-05, 2.50864e-06, 2.9249e-08, 0,

    0, 0, 3.14197e-08, 7.22146e-07, 4.61954e-06, 8.45613e-06, 4.49082e-06,
    6.85919e-07, 2.9364e-08, 0, 0,

    0, 0, 0, 0, 0, 2.63305e-08, 0, 0, 0, 0, 0
  };





void
preparearrays(struct noisechiselparams *p)
{
  struct gal_mesh_params *smp=&p->smp, *lmp=&p->lmp;

  long *meshindexs;
  float *f, *ff, *fp;
  size_t *relngb=p->relngb, s0, s1;

  /* Read the input image in. Note that the pointer to the image is
     also kept in p->img. Since some of the mesh operations should be
     done on the convolved image and some on the actual image, we will
     need to change the mesh's img value some times and the p->img
     will be used to keep its actual value. */
  gal_fits_file_to_float(p->up.inputname, p->up.maskname, p->cp.hdu,
                              p->up.mhdu, (float **)&smp->img, &p->bitpix,
                              &p->anyblank, &smp->s0, &smp->s1);
  gal_fits_read_wcs(p->up.inputname, p->cp.hdu, 0, 0, &p->nwcs, &p->wcs);
  s0=smp->s0; s1=smp->s1;

  /* make sure the channel sizes fit the channel sizes. */
  if( s0%smp->nch2 || s1%smp->nch1 )
    error(EXIT_FAILURE, 0, "the input image size (%lu x %lu) is not an "
          "exact multiple of the number of the given channels (%lu, %lu) "
          "in the respective axis", s1, s0, smp->nch1, smp->nch2);

  /* p->imgss (image-sky-subtracted) is the sky subtracted input
     image. For both the removal of false detections and also the
     segmentation, it is important to use the sky subtracted image,
     not the actual input image for some operations. In both cases the
     input image is used for Signal to noise ratio measurements (where
     subtracting the sky will add noise).

       False detection removal: The sky subtracted image will be used
           for thresholding over the detected and undetected regions.

       Segmentation: The sky subtracted image will also be used for
           generating the catalog.

     Since both operations involve the sky subtracted input image
     (with different sky values) one array playing the role can really
     help in following the code. It will also help in the memory usage
     of the program. This array is allocated in the beginning and
     freed in the end and used throughout. A huge chunk of memory
     doesn't have to be allocated and de-allocated on every step.  */
  errno=0; p->imgss=malloc(s0*s1*sizeof *p->imgss);
  if(p->imgss==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->imgss in preparearrays "
          "(ui.c)", s0*s1*sizeof *p->imgss);

  /* Read the kernel: */
  if(p->up.kernelnameset)
    gal_fits_prep_float_kernel(p->up.kernelname, p->up.khdu, &smp->kernel,
                                    &smp->ks0, &smp->ks1);
  else
    {
      errno=0;
      smp->ks0=defaultkernel_s0;
      smp->ks1=defaultkernel_s1;
      smp->kernel=malloc(smp->ks0*smp->ks1*sizeof *smp->kernel);
      if(smp->kernel==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for default kernel",
              smp->ks0*smp->ks1);
      ff=defaultkernel;
      fp=(f=smp->kernel)+smp->ks0*smp->ks1;
      do *f=*ff++; while(++f<fp);
    }

  /* Allocate the other necessary arrays: */
  errno=0; p->byt=malloc(s0*s1*sizeof *p->byt);
  if(p->byt==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->byt (ui.c)",
          s0*s1*sizeof *p->byt);
  errno=0; p->olab=malloc(s0*s1*sizeof *p->olab);
  if(p->olab==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->olab (ui.c)",
          s0*s1*sizeof *p->olab);
  errno=0; p->clab=malloc(s0*s1*sizeof *p->clab);
  if(p->clab==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->clab (ui.c)",
          s0*s1*sizeof *p->clab);

  /* This ngb array is used to keep the relative indexs of the
     neighbors of a pixel. It is used in the over-segmentation step
     (clumps.c). The labelings are such that the first four elements
     are the four-connected ones and the second four are
     8-connected. There is no problem with the negative values that
     are stored as size_t (which is an unsigned type): they will be
     added with positive values during the processing to give correct
     values. */
  relngb[4]=    s1-1;    relngb[0]=    s1;    relngb[5]=    s1+1;
  relngb[1]=      -1;                         relngb[2]=       1;
  relngb[6]= -1*s1-1;    relngb[3]= -1*s1;    relngb[7]= -1*s1+1;

  /* Set the parameters for both mesh grids. */
  lmp->s0=smp->s0;
  lmp->s1=smp->s1;
  lmp->ks0=smp->ks0;
  lmp->ks1=smp->ks1;
  lmp->nch1=smp->nch1;
  lmp->nch2=smp->nch2;
  lmp->kernel=smp->kernel;
  lmp->img=p->img=smp->img;
  lmp->params=smp->params=p;
  lmp->minmodeq=smp->minmodeq;
  lmp->mirrordist=smp->mirrordist;
  lmp->fullsmooth=smp->fullsmooth;
  lmp->numnearest=smp->numnearest;
  lmp->smoothwidth=smp->smoothwidth;
  lmp->lastmeshfrac=smp->lastmeshfrac;
  lmp->meshbasedcheck=smp->meshbasedcheck;
  lmp->interponlyblank=smp->interponlyblank;
  lmp->fullinterpolation=smp->fullinterpolation;
  lmp->numthreads=smp->numthreads=p->cp.numthreads;


  /* Prepare the mesh structures. */
  gal_mesh_make_mesh(smp);
  gal_mesh_make_mesh(lmp);
  if(p->meshname)
    {
      gal_fits_array_to_file(p->meshname, "Input", FLOAT_IMG,
                             smp->img, s0, s1, p->anyblank, p->wcs,
                             NULL, SPACK_STRING);
      gal_check_mesh_id(smp, &meshindexs);
      gal_fits_array_to_file(p->meshname, "SmallMeshIndexs",
                             LONG_IMG, meshindexs, s0, s1, 0, p->wcs,
                             NULL, SPACK_STRING);
      free(meshindexs);
      gal_check_mesh_id(lmp, &meshindexs);
      gal_fits_array_to_file(p->meshname, "LargeMeshIndexs", LONG_IMG,
                             meshindexs, s0, s1, 0, p->wcs,
                             NULL, SPACK_STRING);
      free(meshindexs);
    }
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct noisechiselparams *p)
{
  struct gal_commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = num_processors(NPROC_CURRENT);
  cp->removedirinfo = 1;

  /* NoiseChisel parameter initializations. */
  p->detsnhistnbins=p->clumpsnhistnbins=0;

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

  /* Do a sanity check. */
  sanitycheck(p);

  /* Make the array of input images. */
  preparearrays(p);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf("  - Using %lu CPU thread%s\n", p->cp.numthreads,
             p->cp.numthreads==1 ? "." : "s.");
      printf("  - Input: %s (hdu: %s)\n", p->up.inputname, p->cp.hdu);
      if(p->up.maskname)
        printf("  - Mask: %s (hdu: %s)\n", p->up.maskname, p->up.mhdu);
      if(p->up.kernelnameset)
        printf("  - Kernel: %s (hdu: %s)\n", p->up.kernelname,
               p->up.khdu);
      else
        printf("  - Kernel: FWHM=2 pixel Gaussian.\n");
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct noisechiselparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->img);
  free(p->byt);
  free(p->olab);
  free(p->clab);
  free(p->imgss);
  free(p->cp.hdu);
  free(p->up.mhdu);
  free(p->up.khdu);
  free(p->cp.output);
  free(p->smp.kernel);
  free(p->up.kernelname);

  /* Free the mask image name. Note that p->up.inputname was not
     allocated, but given to the program by the operating system. */
  if(p->up.maskname && p->up.maskname!=p->up.inputname)
    free(p->up.maskname);

  /* Free all the allocated names. Note that detsnhist */
  free(p->skyname);
  free(p->meshname);
  free(p->threshname);
  free(p->maskdetname);
  free(p->detectionname);
  free(p->segmentationname);
  free(p->detectionskyname);

  /* Free the WCS structure: */
  if(p->wcs)
    wcsvfree(&p->nwcs, &p->wcs);

  /* Print the final message. */
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in", 0);
}
