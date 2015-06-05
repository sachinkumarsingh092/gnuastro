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

#include "mesh.h"
#include "forqsort.h"
#include "statistics.h"

#include "main.h"

#include "sky.h"

void *
avestdonthread(void *inparam)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparam;
  struct meshparams *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  float *mponeforall=mp->oneforall;
  float *oneforall=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  unsigned char *byt, *inbyt=p->byt;
  size_t s0, s1, ind, start, is1=mp->s1;
  float *f, *img, *imgend, *inimg=p->img;
  float ave, med, std, minbfrac=p->minbfrac;
  size_t i, num, row, *indexs=&mp->indexs[mtp->id*mp->thrdcols];

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
         non-NaN. Recall that both the convolved an unconvolved image
         have the same NaN pixels.*/
      do
        {
          byt = inbyt + start + row*is1;
          imgend=(img = inimg + start + row++ * is1 ) + s1;
          do
            /* Only input pixels that have byt==0 and are not NaN. */
            if(*byt++==0 && isnan(*img)==0)
              {
                ++num;
                *f++ = *img;
              }
          while(++img<imgend);
        }
      while(row<s0);

      /* Do the desired operation on the mesh, all the meshs were
         initialized to NaN, so if they don't fit the criteria, they
         can just be ignored. */
      if( (float)num/(float)(s0*s1)>minbfrac )
        {
          /* Sort the array of values: */
          qsort(oneforall, num, sizeof *oneforall, floatincreasing);

          /* Do sigma-clipping and save the result if it is
             accurate. */
          if(sigmaclip_converge(oneforall, 1, num, p->sigclipmultip,
                                p->sigcliptolerance, &ave, &med, &std, 0))
            {
              mp->garray1[ind]=ave;
              mp->garray2[ind]=std;
            }
        }
    }

  /* Free any allocated space and if multiple threads were used, wait
     until all other threads finish. */
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





/* Using the smaller mesh and the p->byt array, find the average and
   standard deviation of the undetected pixels and put them in the
   smp->garray1 and smp->garray2 arrays. This function will be used
   multiple times, the outputs for each should be different. So it
   takes the second argument as the name.*/
void
findavestdongrid(struct noisechiselparams *p, char *outname)
{
  struct meshparams *smp=&p->smp;

  float *sky, *std;
  size_t s0=smp->s0, s1=smp->s1;


  /* Find the average and standard deviation */
  operateonmesh(smp, avestdonthread, sizeof(float), 1, 1);
  if(outname)
    {
      checkgarray(smp, &sky, &std);
      arraytofitsimg(outname, "Detected", BYTE_IMG, p->byt, s0, s1,
                     0, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(outname, "Sky", FLOAT_IMG, sky, s0, s1,
                     p->numblank, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(outname, "SkySTD", FLOAT_IMG, std, s0, s1,
                     p->numblank, p->wcs, NULL, SPACK_STRING);
      free(sky);
      free(std);
    }


  /* In case the image is in electrons or counts per second the
     standard deviation of the noise will become smaller than
     unity. You have to find the minimum STD value (which is always
     positive) for later corrections. */
  floatmin(smp->garray2, smp->nmeshi, &p->cpscorr);
  if(p->cpscorr>1) p->cpscorr=1.0f;

  /* Interpolate over the meshs to fill all the blank ones in both the
     sky and the standard deviation arrays: */
  meshinterpolate(smp, "Interpolating sky value and its standard deviation");
  if(outname)
    {
      checkgarray(smp, &sky, &std);
      arraytofitsimg(outname, "SkyInterpolated", FLOAT_IMG, sky,
                     s0, s1, 0, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(outname, "SkySTDInterpolated", FLOAT_IMG, std,
                     s0, s1, 0, p->wcs, NULL, SPACK_STRING);
      free(sky);
      free(std);
    }


  /* Smooth the interpolated array:  */
  if(smp->smoothwidth>1)
    {
      meshsmooth(smp);
      if(outname)
        {
          checkgarray(smp, &sky, &std);
          arraytofitsimg(outname,"SkySmoothed", FLOAT_IMG, sky,
                         s0, s1, 0, p->wcs, NULL, SPACK_STRING);
          arraytofitsimg(outname, "SkySTDSmoothed", FLOAT_IMG, std,
                         s0, s1, 0, p->wcs, NULL, SPACK_STRING);
          free(sky);
          free(std);
        }
    }
}





/* Using the p->byt array find the sky value on the input and
   convolved images. Then subtract the sky value from both and save
   the standard deviation for every pixel in p->std. */
void
findsubtractskyimgconv(struct noisechiselparams *p)
{
  struct meshparams *smp=&p->smp;

  float *tmpg1, *tmpg2;
  float *f, *s, *fp, sky, std, *tmpimg, *img=p->img;
  size_t i, s0, s1, row, start, meshid, is1=smp->s1;
  float *c, csky, *tmpcg1, *tmpcg2, *tmpfg1, *tmpfg2, *convsky;


  /* Allocate an array for the standard deviation, we will need it for
     every pixel. */
  errno=0; p->std=malloc(smp->s0*smp->s1*sizeof *p->std);
  if(p->std==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->std in findsubtractsky "
          "(sky.c)", smp->s0*smp->s1*sizeof *p->std);



  /* Replace the necessary arrays to find the sky value on the
     convolved image. */
  tmpimg=p->img;          /* Keep backup pointers to the main arrays.  */
  tmpg1=smp->garray1;           tmpg2=smp->garray2;
  tmpcg1=smp->cgarray1;         tmpcg2=smp->cgarray2;
  tmpfg1=smp->fgarray1;         tmpfg2=smp->fgarray2;

  p->img=smp->img=p->conv;  /* Prepare for working on convolved image. */
  smp->cgarray1=smp->cgarray2=NULL;
  smp->fgarray1=smp->fgarray2=NULL;
  findavestdongrid(p, NULL);

  convsky=smp->garray1;                /* Keep garray1, free the rest. */
  if(smp->garray1==smp->cgarray1)  free(smp->fgarray1);
  else                             free(smp->cgarray1);
  free(smp->cgarray2);             free(smp->fgarray2);

  p->img=smp->img=tmpimg;          /* Set back to their previous state */
  smp->garray1=tmpg1;           smp->garray2=tmpg2;
  smp->cgarray1=tmpcg1;         smp->cgarray2=tmpcg2;
  smp->fgarray1=tmpfg1;         smp->fgarray2=tmpfg2;


  /* Find the sky value and its STD on the input image and put it in
     the p->smp.garray1 and p->smp.garray2 arrays respectively.

     VERY IMPORTANT:
     ###############

     The sky value for the input image should be found after the
     convolved image. This is because findavestdongrid, will also set
     the p->cpscorr value and we want that from the input image, not
     the convolved image.
  */
  findavestdongrid(p, p->skyname);


  /* Apply the threshold */
  for(i=0;i<smp->nmeshi;++i)
    {
      /* Get the meshid from i: */
      meshid=setmeshid(smp, i);

      /* Subtract the sky for each pixel. */
      row=0;
      csky = convsky[i];
      sky = smp->garray1[i];
      std = smp->garray2[i];
      start=smp->start[meshid];
      s0=smp->ts0[smp->types[meshid]];
      s1=smp->ts1[smp->types[meshid]];
      do
        {
          s = p->std + start + row*is1;
          c = p->conv + start + row*is1;
          fp= ( f = img + start + row++ * is1 ) + s1;
          do
            {
              *s++  = std;      /* Fill STD array.              */
              *f++ -= sky;      /* Subtract sky from image.     */
              *c++ -= csky;     /* Subtract sky from convolved. */
            }
          while(f<fp);
        }
      while(row<s0);
    }

  /* Write the output as the last frame of the sky subtracted
     image if the user wants to check it: */
  if(p->skyname)
    arraytofitsimg(p->skyname, "SkySubtracted", FLOAT_IMG, p->img,
                   p->smp.s0, p->smp.s1, 0, p->wcs, NULL, SPACK_STRING);

  /* Clean up: */
  free(convsky);
}
