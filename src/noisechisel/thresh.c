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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include "mode.h"
#include "timing.h"
#include "forqsort.h"
#include "checkset.h"
#include "statistics.h"
#include "spatialconvolve.h"

#include "main.h"
#include "noisechisel.h"










/*********************************************************************/
/***************        Quantile threshold       *********************/
/*********************************************************************/
void *
qthreshonmesh(void *inparam)
{
  struct gal_mesh_thread_params *mtp=(struct gal_mesh_thread_params *)inparam;
  struct gal_mesh_params *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  float *mponeforall=mp->oneforall;
  float *oneforall=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  float qthresh=p->qthresh;
  size_t modeindex=(size_t)(-1);
  size_t s0, s1, ind, row, num, start, is1=mp->s1;
  size_t i, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  float *f, *img, *imgend, *inimg=p->conv, modesym=0.0f;
  float mirrordist=mp->mirrordist, minmodeq=mp->minmodeq;

  /* Start this thread's work: */
  for(i=0;indexs[i]!=GAL_THREADS_NON_THRD_INDEX;++i)
    {
      /* Prepare the values: */
      num=row=0;
      f=oneforall;
      ind=indexs[i];
      start=mp->start[ind];
      s0=mp->ts0[mp->types[ind]];
      s1=mp->ts1[mp->types[ind]];

      /* Copy all the non-NaN pixels images pixels of this mesh into
         the mesh array. Note that currently, the spatial positioning
         of the pixels is irrelevant, so we only keep those that are
         non-NaN.*/
      do
        {
          imgend=(img = inimg + start + row++ * is1 ) + s1;
          do
            if(!isnan(*img))
            {
              ++num;
              *f++ = *img;
            }
          while(++img<imgend);
        }
      while(row<s0);

      /* Do the desired operation on the mesh: */
      if(num)
        {
          qsort(oneforall, num, sizeof *oneforall, gal_qsort_float_increasing);
          gal_mode_index_in_sorted(oneforall, num, mirrordist, &modeindex,
                                   &modesym);
          if( modesym>GAL_MODE_SYM_GOOD && (float)modeindex/(float)num>minmodeq)
            mp->garray1[ind]=
              oneforall[gal_statistics_index_from_quantile(num, qthresh)];
        }
    }

  /* Free alltype and if multiple threads were used, wait until all
     other threads finish. */
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





/* The threshold values are stored in the garray1 array of the
   gal_mesh_params structure. This function is basically identical to the
   gal_mesh_check_garray function in mesh.c, only in the middle, it does
   something else.

   See the comments there to better understand the process.
*/
void
applythreshold(struct noisechiselparams *p)
{
  struct gal_mesh_params *mp=&p->smp; /* `mp' instead of `smp' so you can try */
                                      /* with p->lmp if you like.             */
  unsigned char *b, *byt=p->byt;
  size_t gs0=mp->gs0, gs1=mp->gs1, nch1=mp->nch1;
  size_t s0, s1, fs1=mp->gs1*mp->nch1, chid, inchid;
  size_t *types=mp->types, *ts0=mp->ts0, *ts1=mp->ts1;
  size_t i, f0, f1, nmeshi=mp->nmeshi, nmeshc=mp->nmeshc;
  size_t row, startind, *start=mp->start, meshid, is1=mp->s1;
  float *f, *fp, thresh, *garray1=mp->garray1, *array=p->conv;

    /* Fill the array: */
  for(i=0;i<nmeshi;++i)
    {
      /* Set the proper meshid. */
      if(garray1==mp->fgarray1)
        {
          f0=i/fs1;
          f1=i%fs1;
          inchid = (f0%gs0) * gs1  + f1%gs1;
          chid   = (f0/gs0) * nch1 + f1/gs1;
          meshid = chid * nmeshc + inchid;
        }
      else meshid=i;

      /* Fill the output array with the value in this mesh: */
      row=0;
      thresh=garray1[i];
      s0=ts0[types[meshid]];
      s1=ts1[types[meshid]];
      startind=start[meshid];
      do
        {
          b = byt + startind + row * is1;
          fp= ( f = array + startind + row * is1 ) + s1;
          do *b++ = isnan(*f) || *f>thresh ? 1 : 0; while(++f<fp);
          ++row;
        }
      while(row<s0);
    }
}





void
findapplyqthreshold(struct noisechiselparams *p)
{
  struct gal_mesh_params *mp=&p->smp;

  /* Find the threshold on each mesh: */
  gal_mesh_operate_on_mesh(mp, qthreshonmesh, sizeof(float), 0, 1);
  if(p->threshname)
    gal_mesh_value_file(mp, p->threshname, "Quantile values", NULL,
                        p->wcs, SPACK_STRING);

  gal_mesh_interpolate(mp, "Interpolating quantile threshold");
  if(p->threshname)
    gal_mesh_value_file(mp, p->threshname, "Interpolated", NULL,
                        p->wcs, SPACK_STRING);

  gal_mesh_smooth(mp);
  if(p->threshname)
    gal_mesh_value_file(mp, p->threshname, "smoothed", NULL,
                        p->wcs, SPACK_STRING);

  /* Apply the threshold on all the pixels: */
  applythreshold(p);
}




















/*********************************************************************/
/**************      Average and STD threshold     *******************/
/*********************************************************************/
/* This is very similar to the gal_mesh_check_garray function. The sky and its
   Standard deviation are stored in the garray1 and garray2 arrays of
   smp gal_mesh_params structure. */
void
applydetectionthresholdskysub(struct noisechiselparams *p)
{
  struct gal_mesh_params *smp=&p->smp;

  size_t is0=smp->s0;
  float *f, *in, sky, std;
  float dthresh=p->dthresh;
  unsigned char *b, *bf, *dbyt;
  size_t gid, row, start, chbasedid, *types=smp->types;
  size_t s0, s1, is1=smp->s1, *ts0=smp->ts0, *ts1=smp->ts1;

  /* Allocate the array to keep the threshold value: */
  errno=0; dbyt=p->dbyt=malloc(is0*is1*sizeof *p->dbyt);
  if(dbyt==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for dbyt in "
          "applydetectionthreshold (detection.c)", is0*is1*sizeof *dbyt);


  /* Apply the threshold */
  for(gid=0;gid<smp->nmeshi;++gid)
    {
      /* Get the meshid from i: */
      chbasedid=gal_mesh_ch_based_id_from_gid(smp, gid);

      /* Fill the output array with the value in this mesh: */
      row=0;
      sky = smp->garray1[gid];
      std = smp->garray2[gid];
      s0=ts0[types[chbasedid]];
      s1=ts1[types[chbasedid]];
      start=smp->start[chbasedid];
      do
        {
          in = p->img + start + row*is1;
          f = p->imgss + start + row*is1;
          bf= ( b = dbyt + start + row++ * is1) + s1;
          do
            {
              /* The threshold is always very low. So for the majority
                 of non-NaN pixels in the image, the condition above
                 will be true. If we come over a NaN pixel, then by
                 definition of NaN, all conditionals will fail.

                 The benefit of this method of testing is that if an
                 image doesn't have any NaN pixels, only the pixels
                 below the threshold have to be checked for a NaN
                 which are by definition a very small fraction of the
                 total pixels. And if there are NaN pixels in the
                 image, they will be checked. */
              *b = ( (*f++=*in-sky) > dthresh*std
                     ? 1
                     : isnan(*in) ? GAL_FITSARRAY_BYTE_BLANK : 0 );
              ++in;
            }
          while(++b<bf);
        }
      while(row<s0);
    }
}




















/*********************************************************************/
/***************      S/N Quantile threshold     *********************/
/*********************************************************************/
void
snthresh(struct noisechiselparams *p, float *sntable, size_t size,
         int det0seg1)
{
  double sn;
  char report[200], cline[1000];
  char *job = det0seg1 ? "Clump" : "Detection";
  float quant = det0seg1 ? p->segquant : p->detquant;
  char *name = det0seg1 ? "clumps" : "pseudo-detections";
  char *histname= det0seg1 ? p->clumpsnhist : p->detectionsnhist;
  size_t snhistnbins= det0seg1 ? p->clumpsnhistnbins : p->detsnhistnbins;


  /* Check if the number is acceptable to the user. */
  if(size<p->minnumfalse)
    error(EXIT_FAILURE, 0, "There are only %lu %s in the sky region of "
          "the image. This is smaller than the minimum number you "
          "specified: %lu. You can decrease this minimum with the "
          "`--minnumfalse' (`-F') option or you can decrease the other "
          "parameters that determine the %s. See the GNU Astronomy "
          "Utilities manual (section on NoiseChisel) or Akhlaghi and "
          "Ichikawa (2015) for more information.", size, name,
          p->minnumfalse, name);


  /* Sort the signal to noise ratios and remove their outliers */
  qsort(sntable, size, sizeof *sntable, gal_qsort_float_increasing);

  /* The removal of outliers was useful when the S/N was calculated
     separately on each mesh (thus there were very few
     points). However, now that the S/N is calculated over the full
     image, the number of pseudo-detections and clumps is so high that
     these outliers will not play a significant role in the S/N
     threshold (unless you set unreasonably high quantiles!). So This
     is commented now. */
  /*gal_statistics_remove_outliers_flat_cdf(sntable, &size);*/


  /* Store the SN value. */
  sn=sntable[gal_statistics_index_from_quantile(size, quant)];
  if(p->cp.verb)
    {
      sprintf(report, "%s S/N: %.3f (%.3f quantile of %lu %s).",
              job, sn, quant, size, name);
      gal_timing_report(NULL, report, 2);
    }

  /* Put the S/N value in its proper place. */
  if(det0seg1) p->clumpsn = sn;
  else         p->detsn   = sn;


  /* If the user has asked for it, make the histogram of the S/N
     distribution. */
  if(snhistnbins)
    {
      /* For a check:
         if(ind!=0) continue;
         ff=(f=sntable)+numlabs; do printf("%f\n", *f++); while(f<ff);
      */

      /* histname has to be set to NULL so gal_checkset_automatic_output can
         safey free it. */
      sprintf(cline, "# %s\n# %s started on %s"
              "# Input: %s (hdu: %s)\n"
              "# S/N distribution histogram of %lu sky %s.\n"
              "# The %.3f quantile has an S/N of %.4f.",
              SPACK_STRING, SPACK_NAME, ctime(&p->rawtime),
              p->up.inputname, p->cp.hdu, size, name, quant, sn);
      gal_statistics_save_hist(sntable, size, snhistnbins, histname, cline);
      free(histname);
    }
}
