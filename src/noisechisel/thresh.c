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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include "mode.h"
#include "timing.h"
#include "forqsort.h"
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
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparam;
  struct meshparams *mp=mtp->mp;
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
  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
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
      qsort(oneforall, num, sizeof *oneforall, floatincreasing);
      modeindexinsorted(oneforall, num, mirrordist, &modeindex, &modesym);
      if( modesym>MODESYMGOOD && (float)modeindex/(float)num>minmodeq)
        mp->garray1[ind]=oneforall[indexfromquantile(num, qthresh)];
    }

  /* Free alltype and if multiple threads were used, wait until all
     other threads finish. */
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





/* The threshold values are stored in the garray1 array of the
   meshparams structure. This function is basically identical to the
   checkgarray function in mesh.c, only in the middle, it does
   something else.

   See the comments there to better understand the process.
*/
void
applythreshold(struct noisechiselparams *p)
{
  struct meshparams *mp=&p->smp; /* `mp' instead of `smp' so you can try */
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
  float *thresh=NULL;
  struct meshparams *mp=&p->smp; /* `mp' instead of `smp' for you to */
  size_t  s0=mp->s0,  s1=mp->s1; /* try this function with p->lmp if */
                                 /* you like to see the effect. */


  /* Find the threshold on each mesh: */
  operateonmesh(mp, qthreshonmesh, sizeof(float), 0, 1);
  if(p->threshname)
    {
      checkgarray(mp, &thresh, NULL);
      arraytofitsimg(p->threshname, "QuantileThreshold", FLOAT_IMG,
                     thresh, s0, s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      free(thresh);
    }

  meshinterpolate(mp, "Interpolating quantile threshold");
  if(p->threshname)
    {
      checkgarray(mp, &thresh, NULL);
      arraytofitsimg(p->threshname, "Interpolated", FLOAT_IMG,
                     thresh, s0, s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      free(thresh);
    }

  meshsmooth(mp);
  if(p->threshname)
    {
      checkgarray(mp, &thresh, NULL);
      arraytofitsimg(p->threshname, "smoothed", FLOAT_IMG,
                     thresh, s0, s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      free(thresh);
    }


  /* Apply the threshold on all the pixels: */
  applythreshold(p);
}




















/*********************************************************************/
/**************      Average and STD threshold     *******************/
/*********************************************************************/
/* This is very similar to the checkgarray function. The sky and its
   Standard deviation are stored in the garray1 and garray2 arrays of
   smp meshparams structure. */
void
applydetectionthresholdskysub(struct noisechiselparams *p)
{
  struct meshparams *smp=&p->smp;

  size_t is0=smp->s0;
  unsigned char *b, *dbyt;
  float dthresh=p->dthresh;
  float *f, *fp, *in, sky, std;
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
      chbasedid=chbasedidfromgid(smp, gid);

      /* Fill the output array with the value in this mesh: */
      row=0;
      sky = smp->garray1[gid];
      std = smp->garray2[gid];
      s0=ts0[types[chbasedid]];
      s1=ts1[types[chbasedid]];
      start=smp->start[chbasedid];
      do
        {
          b = dbyt + start + row*is1;
          in = p->img + start + row*is1;
          fp= ( f = p->imgss + start + row++ * is1 ) + s1;
          do
            {
              /*
                The basic idea behind choosing this comparison is that
                we want it to work for NaN values too.

                (From the GNU C library manual) "In comparison
                operations, ... NaN is `unordered': it is not equal
                to, greater than, or less than anything, _including
                itself_."

                In this context, we want NaN pixels to be treated like
                those that are above the threshold (a region that is
                only NaN will be removed later because it has no
                flux). So if the checking condition below fails, we
                want *b=1. To make this work when *f!=NAN also, we
                check if the flux is below the threshold (so when it
                fails, *b=1). In this manner we don't have to add an
                `isnan' check and make this a tiny bit faster.
              */
              *f = *in++ - sky;
              *b++ = *f<dthresh*std ? 0 : 1;
            }
          while(++f<fp);
        }
      while(row<s0);
    }
}
