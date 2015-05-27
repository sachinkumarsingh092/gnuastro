/*********************************************************************
SubtractSky - Find and subtract the sky value from an image.
SubtractSky is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include "mesh.h"
#include "mode.h"
#include "timing.h"
#include "forqsort.h"
#include "statistics.h"
#include "arraymanip.h"
#include "fitsarrayvv.h"

#include "main.h"


void *
avestdonthread(void *inparam)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparam;
  struct meshparams *mp=mtp->mp;
  struct subtractskyparams *p=(struct subtractskyparams *)mp->params;

  float *mponeforall=mp->oneforall;
  float *oneforall=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  int setnan;
  float *cofa;                   /* convolved-oneforall */
  size_t modeindex=(size_t)(-1);
  size_t s0, s1, ind, row, num, start, is1=mp->s1;
  size_t i, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  float mirrordist=mp->mirrordist, minmodeq=mp->minmodeq;
  float ave, med, std, *sorted, sigclipmultip=p->sigclipmultip;
  float *imgend, *inimg=mp->img, *inconv=p->conv, modesym=0.0f;
  float *f, *c=NULL, *conv, *img, sigcliptolerance=p->sigcliptolerance;

  /* Allocate the oneforall array for the convolved image, since we
     only want to use those meshs whose convolved image mode is larger
     than minmodeq. */
  if(p->conv!=mp->img)
    {
      errno=0; cofa=malloc(mp->maxs0*mp->maxs1*sizeof *oneforall);
      if(cofa==NULL)
        error(EXIT_FAILURE, errno, "Unable to allocate %lu bytes for"
              "cofa in avestdonthread of subtractsky.c.",
              mp->maxs0*mp->maxs1*sizeof *mp->img);
    }
  else cofa=NULL;

  /* Start this thread's work: */
  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
    {
      /* Prepare the values: */
      setnan=0;
      num=row=0;
      f=oneforall;
      ind=indexs[i];
      start=mp->start[ind];
      s0=mp->ts0[mp->types[ind]];
      s1=mp->ts1[mp->types[ind]];
      if(cofa) c=sorted=cofa; else sorted=oneforall;

      /* Copy all the non-NaN pixels images pixels of this mesh into
         the mesh array. Note that currently, the spatial positioning
         of the pixels is irrelevant, so we only keep those that are
         non-NaN. Recall that both the convolved an unconvolved image
         have the same NaN pixels.*/
      do
        {
          if(cofa) conv = inconv + start + row * is1;
          imgend=(img = inimg + start + row++ * is1 ) + s1;
          do
            {
              if(!isnan(*img))
                {
                  ++num;
                  *f++ = *img;
                  if(cofa) *c++=*conv;
                }
              if(cofa) ++conv;
            }
          while(++img<imgend);
        }
      while(row<s0);

      /* Do the desired operation on the mesh: */
      qsort(sorted, num, sizeof *oneforall, floatincreasing);
      modeindexinsorted(sorted, num, mirrordist, &modeindex, &modesym);
      if( modesym>MODESYMGOOD && (float)modeindex/(float)num>minmodeq )
        {
          /* If cofa was defined, then oneforall was not sorted. */
          if(cofa)
            qsort(oneforall, num, sizeof *oneforall, floatincreasing);

          /* Do sigma-clipping and save the result if it is
             accurate. */
          if(sigmaclip_converge(oneforall, 1, num, sigclipmultip,
                                sigcliptolerance, &ave, &med, &std, 0))
            {
              mp->cgarray1[ind]=ave;
              if(mp->cgarray2) mp->cgarray2[ind]=std;
            }
          else setnan=1;
        }
      else setnan=1;

      /* Set this mesh should not be used: */
      if(setnan)
        {
          mp->cgarray1[ind]=NAN;
          if(mp->cgarray2) mp->cgarray2[ind]=NAN;
        }
    }

  /* Free any allocated space and if multiple threads were used, wait
     until all other threads finish. */
  if(p->conv!=mp->img) free(cofa);
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





void
subtractsky(struct subtractskyparams *p)
{
  struct meshparams *mp=&p->mp;

  long *meshindexs;
  struct timeval t1;
  int checkstd=p->checkstd;
  size_t s0=mp->s0, s1=mp->s1;
  float *sky, *std, *skysubtracted;


  /* Prepare the mesh array. */
  gettimeofday(&t1, NULL);
  makemesh(mp);
  if(p->meshname)
    {
      checkmeshid(mp, &meshindexs);
      arraytofitsimg(p->meshname, "Input", FLOAT_IMG, p->mp.img, s0, s1,
                     p->numblank, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(p->meshname, "MeshIndexs", LONG_IMG, meshindexs,
                     s0, s1, 0, p->wcs, NULL, SPACK_STRING);
      free(meshindexs);
    }
  if(p->cp.verb) reporttiming(&t1, "Mesh grid ready.", 1);



  /* Convolve the image if the user has asked for it: */
  if(p->up.kernelnameset)
    {
      spatialconvolveonmesh(mp, &p->conv);
      if(p->convname)
        {
          arraytofitsimg(p->convname, "Input", FLOAT_IMG, p->mp.img, s0, s1,
                         p->numblank, p->wcs, NULL, SPACK_STRING);
          arraytofitsimg(p->convname, "Input", FLOAT_IMG, p->conv, s0, s1,
                         p->numblank, p->wcs, NULL, SPACK_STRING);
        }
    }
  else p->conv=p->mp.img;
  if(p->cp.verb)
    reporttiming(&t1, "Input image convolved with kernel.", 1);


  /* Find the sky value and its standard deviation on each mesh. */
  operateonmesh(mp, avestdonthread, sizeof(float), checkstd);
  if(p->interpname)
    {
      checkgarray(mp, &sky, &std);
      arraytofitsimg(p->interpname, "Sky", FLOAT_IMG, sky, s0, s1,
                     p->numblank, p->wcs, NULL, SPACK_STRING);
      free(sky);
      if(checkstd)
        {
          arraytofitsimg(p->interpname, "SkySTD", FLOAT_IMG, std, s0, s1,
                         p->numblank, p->wcs, NULL, SPACK_STRING);
          free(std);
        }
    }
  if(p->cp.verb)
    reporttiming(&t1, "Sky and its STD found on some meshes.", 1);


  /* Interpolate over the meshs to fill all the blank ones in both the
     sky and the standard deviation arrays: */
  meshinterpolate(mp);
  if(p->interpname)
    {
      checkgarray(mp, &sky, &std);
      arraytofitsimg(p->interpname, "Sky", FLOAT_IMG, sky, s0, s1, 0,
                     p->wcs, NULL, SPACK_STRING);
      free(sky);
      if(checkstd)
        {
          arraytofitsimg(p->interpname, "SkySTD", FLOAT_IMG, std, s0, s1, 0,
                         p->wcs, NULL, SPACK_STRING);
          free(std);
        }
    }
  if(p->cp.verb)
    reporttiming(&t1, "All blank meshs filled (interplated).", 1);



  /* Smooth the interpolated array:  */
  if(mp->smoothwidth>1)
    {
      meshsmooth(mp);
      if(p->cp.verb)
        reporttiming(&t1, "Mesh grid smoothed.", 1);
    }


  /* Make the sky array and save it if the user has asked for it: */
  checkgarray(mp, &sky, &std);
  if(p->skyname)
    {
      arraytofitsimg(p->skyname ,"Sky", FLOAT_IMG, sky, s0, s1, 0,
                     p->wcs, NULL, SPACK_STRING);
      if(checkstd)
        arraytofitsimg(p->skyname, "SkySTD", FLOAT_IMG, std, s0, s1, 0,
                       p->wcs, NULL, SPACK_STRING);
    }


  /* Subtract the sky value */
  fmultipconst(sky, s0*s1, -1.0f);
  skysubtracted=fsumarrays(mp->img, sky, s0*s1);
  arraytofitsimg(p->cp.output ,"SkySubtracted", FLOAT_IMG, skysubtracted,
                 s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);

  /* Clean up: */
  free(sky);
  freemesh(mp);
  free(skysubtracted);
  if(checkstd) free(std);
  if(p->up.kernelnameset) free(p->conv);
}
