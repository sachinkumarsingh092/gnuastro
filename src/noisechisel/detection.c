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
#include "timing.h"
#include "forqsort.h"
#include "checkset.h"
#include "arraymanip.h"
#include "statistics.h"
#include "astrthreads.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "label.h"
#include "binary.h"
#include "detection.h"










/*********************************************************************/
/**************    Average and STD of undetected   *******************/
/*********************************************************************/
void *
avestdonthread(void *inparam)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparam;
  struct meshparams *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  float *mponeforall=mp->oneforall;
  float *oneforall=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  int setnan;
  unsigned char *byt, *inbyt=p->byt;
  size_t s0, s1, ind, start, is1=mp->s1;
  float *f, *img, *imgend, *inimg=mp->img;
  float ave, med, std, minbfrac=p->minbfrac;
  size_t i, num, row, *indexs=&mp->indexs[mtp->id*mp->thrdcols];

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

      /* Do the desired operation on the mesh: */
      if( (float)num/(float)(s0*s1)>minbfrac )
        {
          /* Sort the array of values: */
          qsort(oneforall, num, sizeof *oneforall, floatincreasing);

          /* Do sigma-clipping and save the result if it is
             accurate. */
          if(sigmaclip_converge(oneforall, 1, num, p->sigclipmultip,
                                p->sigcliptolerance, &ave, &med, &std, 0))
            {
              mp->cgarray1[ind]=ave;
              mp->cgarray2[ind]=std;
            }
          else setnan=1;
        }
      else setnan=1;

      /* Set this mesh should not be used: */
      if(setnan)
        {
          mp->cgarray1[ind]=NAN;
          mp->cgarray2[ind]=NAN;
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
findavestdongrid(struct noisechiselparams *p, char *outname, int i0f1)
{
  struct meshparams *smp=&p->smp;

  float *sky, *std;
  int verb=p->cp.verb;
  size_t s0=smp->s0, s1=smp->s1;


  /* Find the average and standard deviation */
  operateonmesh(smp, avestdonthread, sizeof(float), 1);
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
  if(verb)
    {
      if(i0f1)
        reporttiming(NULL, "Sky and its STD found on some meshes.", 2);
      else
        reporttiming(NULL, "Initial sky and its STD found on some meshes.",
                     2);
    }


  /* In case the image is in electrons or counts per second the
     standard deviation of the noise will become smaller than
     unity. You have to find the minimum STD value (which is always
     positive) for later corrections. */
  floatmin(smp->garray2, smp->nmeshi, &p->cpscorr);
  if(p->cpscorr>1) p->cpscorr=1.0f;


  /* Interpolate over the meshs to fill all the blank ones in both the
     sky and the standard deviation arrays: */
  meshinterpolate(smp);
  if(outname)
    {
      checkgarray(smp, &sky, &std);
      arraytofitsimg(outname, "Sky", FLOAT_IMG, sky, s0, s1, 0,
                     p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(outname, "SkySTD", FLOAT_IMG, std, s0, s1, 0,
                     p->wcs, NULL, SPACK_STRING);
      free(sky);
      free(std);
    }
  if(verb)
    reporttiming(NULL, "All blank meshs filled (interplated).", 2);


  /* Smooth the interpolated array:  */
  if(smp->smoothwidth>1)
    {
      meshsmooth(smp);
      if(outname)
        {
          checkgarray(smp, &sky, &std);
          arraytofitsimg(outname,"Sky", FLOAT_IMG, sky, s0, s1, 0,
                         p->wcs, NULL, SPACK_STRING);
          arraytofitsimg(outname, "SkySTD", FLOAT_IMG, std, s0, s1, 0,
                         p->wcs, NULL, SPACK_STRING);
        }
      if(verb)
        reporttiming(NULL, "Grid smoothed.", 2);
    }
}






/* This is very similar to the checkgarray function. The sky and its
   Standard deviation are stored in the garray1 and garray2 arrays of
   smp meshparams structure. */
void
applydetectionthresholdskysub(struct noisechiselparams *p)
{
  unsigned char *b, *dbyt;
  float dthresh=p->dthresh;
  struct meshparams *smp=&p->smp;
  float *f, *fp, sky, std, *img=p->img;
  size_t is0=smp->s0, gs0=smp->gs0, gs1=smp->gs1;
  size_t i, row, start, meshid, *types=smp->types;
  size_t f0, f1, fs1=smp->gs1*smp->nch1, chid, inchid;
  size_t s0, s1, is1=smp->s1, *ts0=smp->ts0, *ts1=smp->ts1;

  /* Allocate the array to keep the threshold value: */
  errno=0; dbyt=p->dbyt=malloc(is0*is1*sizeof *p->dbyt);
  if(dbyt==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for dbyt in "
          "applydetectionthreshold (detection.c)", is0*is1*sizeof *dbyt);

  /* Apply the threshold */
  for(i=0;i<smp->nmeshi;++i)
    {
      /* Similar to checkgarray (mesh.c) */
      if(smp->garray1==smp->fgarray1)
        {
          f0=i/fs1;
          f1=i%fs1;
          inchid = (f0%gs0)*gs1       + f1%gs1;
          chid   = (f0/gs0)*smp->nch1 + f1/gs1;
          meshid = chid * smp->nmeshc + inchid;
        }
      else meshid=i;

      /* Fill the output array with the value in this mesh: */
      row=0;
      s0=ts0[types[meshid]];
      s1=ts1[types[meshid]];
      sky = smp->garray1[i];
      std = smp->garray2[i];
      start=smp->start[meshid];
      do
        {
          b = dbyt + start + row*is1;
          fp= ( f = img + start + row++ * is1 ) + s1;
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
              *f-=sky; /* Sky subtracted for later, p->img is a copy */
              *b++ = *f<dthresh*std ? 0 : 1;
            }
          while(++f<fp);
        }
      while(row<s0);
    }
}



















/****************************************************************
 *****************     Multi-threaded detection     *************
 ****************************************************************/
/* The p->byt points to the array keeping the positions of initial
   detections. p->dbyt points to the array keeping the thresholded
   values for removal of false detections. We want to copy the above
   threshold (with value 1) parts of p->dbyt to `out' which are either
   foreground or background in p->byt (depending on if the process is
   done over the detections or over the noise). */
void
bytpartfromlarge(struct noisechiselparams *p, unsigned char *out,
               size_t start, size_t s0, size_t s1)
{
  size_t r=0, is1=p->smp.s1;
  unsigned char *b, *d, *bf, b0f1=p->b0f1;

  do
    {
      d = p->dbyt + start + r*is1;
      bf = ( b = p->byt + start + r++ * is1 ) + s1;
      do
        {
          *out++ = *b==b0f1 ? *d : 0;
          ++d;
        }
      while(++b<bf);
    }
  while(r<s0);
}





/* Copy a region of a smaller array into a larger array. *in is the
   smaller array and *out is the larger array. */
void
bytparttolarge(struct noisechiselparams *p, unsigned char *in,
               size_t start, size_t s0, size_t s1)
{
  unsigned char *d, *df;
  size_t r=0, is1=p->smp.s1;

  do
    {
      df = ( d = p->dbyt + start + r++ * is1 ) + s1;
      do *d=*in++; while(++d<df);
    }
  while(r<s0);
}





void
floatpartfromlarge(struct noisechiselparams *p, float *img, float *std,
                   size_t start, size_t s0, size_t s1)
{
  float *f, *ff, *s;
  size_t r=0, is1=p->smp.s1;

  do
    {
      s = p->std + start + r*is1;
      ff = ( f = p->img + start + r++ * is1 ) + s1;
      do { *img++ = *f; *std++=*s; ++s;} while(++f<ff);
    }
  while(r<s0);
}





void*
snthreshonmesh(void *inparams)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparams;
  struct meshparams *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  unsigned char *mponeforall=mp->oneforall;
  unsigned char *thisbyt=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  long *thislab;
  size_t i, s0, s1;
  unsigned char *dbyt=p->dbyt;
  char *detectionname=p->detectionname;
  float *thisimg, *thisstd, minbfrac=p->minbfrac;
  size_t ind, size, numcon, startind, is1=mp->s1;
  size_t nf, nb, *indexs=&mp->indexs[mtp->id*mp->thrdcols];


  /* Allocate all the necessary arrays: */
  errno=0; thislab=malloc(mp->maxs0*mp->maxs1*sizeof *thislab);
  if(thislab==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for thislab (detection.c)",
          mp->maxs0*mp->maxs1*sizeof *thislab);
  errno=0; thisimg=malloc(mp->maxs0*mp->maxs1*sizeof *thisimg);
  if(thisimg==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for thisimg (detection.c)",
          mp->maxs0*mp->maxs1*sizeof *thisimg);
  errno=0; thisstd=malloc(mp->maxs0*mp->maxs1*sizeof *thisstd);
  if(thisstd==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for thisstd (detection.c)",
          mp->maxs0*mp->maxs1*sizeof *thisstd);


  /* Do the job for each mesh: */
  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
    {
      /* Set the necesary parameters: */
      ind=indexs[i];
      s0=mp->ts0[mp->types[ind]];
      s1=mp->ts1[mp->types[ind]];
      startind=mp->start[ind];
      size=s0*s1;

      /* If we are on "noise" mode, check to see if we have enough
	 blank area for getting the background noise statistics. */
      if(p->b0f1==0)
	{
	  count_f_b_onregion(p->byt, startind, s0, s1, is1, &nf, &nb);
	  if( (float)nb < (float)(size)*minbfrac )
	    {
              mp->cgarray1[ind]=NAN;
	      if(detectionname)
		ucharinitonregion(dbyt, 0, startind, s0, s1, is1);
	      continue;
	    }
	}

      /* Copy this mesh into a separate array. In the meantime, only
         bring in the desired pixels. If the user wanted to see the
         steps, then save the result to the p->dbyt and come leave the
         work for this mesh if this is the first step. */
      bytpartfromlarge(p, thisbyt, startind, s0, s1);
      if(detectionname && p->stepnum==1)
	{ bytparttolarge(p, thisbyt, startind, s0, s1); continue; }


      /* Find the connected components */
      fillboundedholes(thisbyt, s0, s1);
      if(detectionname && p->stepnum==2)
	{ bytparttolarge(p, thisbyt, startind, s0, s1); continue; }


      /* Open the image once: */
      opening(thisbyt, s0, s1, 1, 4);
      if(detectionname && p->stepnum==3)
	{ bytparttolarge(p, thisbyt, startind, s0, s1); continue; }


      /* Label the opened regions, remove the small area regions, and
         relabel all the ones that are larger.  */
      numcon=BF_concmp(thisbyt, thislab, s0, s1, 4);
      removesmallarea_relabel(thislab, thisbyt, size, &numcon,
                              p->detsnminarea);
      if(p->b0f1 && numcon<100)
        {
          if(detectionname)
            ucharinitonregion(dbyt, 0, startind, s0, s1, is1);
          continue;
        }
      if(detectionname && p->stepnum==4)
	{ bytparttolarge(p, thisbyt, startind, s0, s1); continue; }


      /* Make a copy of the input image and the standard deviation on
         each pixel to calculate the signal to noise ratio. */
      floatpartfromlarge(p, thisimg, thisstd, startind, s0, s1);

    }

  /* Free any allocated space and if multiple threads were used, wait
     until all other threads finish. */
  free(thisstd);
  free(thisimg);
  free(thislab);
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





/* Doing object detection is done separately on each mesh both for the
   first part of finding the SN limit and the second part of applying
   it and removing false detections. Therefore, if the user asks to
   view the steps, the mesh process has to be stopped at differnt
   times for all the meshs and re-done for each step. This function
   does that job. */
void
snthreshongrid(struct noisechiselparams *p)
{
  struct meshparams *lmp=&p->lmp;

  char *extname=NULL;
  unsigned char *tmp;
  size_t s0=lmp->s0, s1=lmp->s1;


  /* If lmp->garray1 (the signal-to-noise limit on each mesh) is
     allocated, then it shows that we are working on the data pixels,
     otherwise we are working on the noise pixels. */
  p->b0f1 = lmp->garray1 ? 1 : 0;


  /* Find the SN threshold and save the steps in a FITS file: */
  if(p->detectionname)
    {
      p->stepnum=1;
      ucharcopy(p->dbyt, s0*s1, &tmp); /* Backup of p->dbyt in tmp */
      while( (p->b0f1==0 && p->stepnum<5)           /* Undetected. */
             || (p->b0f1 && p->stepnum<6))          /* Detected.   */
        {
          free(p->dbyt);
          ucharcopy(tmp, s0*s1, &p->dbyt);
          operateonmesh(lmp, snthreshonmesh, sizeof(unsigned char), 0);
          switch(p->stepnum)
            {
            case 1:
              extname = p->b0f1 ? "ThresholdDetections" : "ThresholdNoise";
              break;
            case 2:
              extname="HolesFilled";
              break;
            case 3:
              extname="Opened";
              break;
            case 4:
              extname="SmallRemoved";
              break;
            default:
              error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so "
                    "we can fix the problem. For some reason p->stepnum "
                    "(detection.c) has a value of %d which is not "
                    "recognized.", PACKAGE_BUGREPORT, p->stepnum);
            }
          arraytofitsimg(p->detectionname, extname, BYTE_IMG, p->dbyt,
                         s0, s1, 0, p->wcs, NULL, SPACK_STRING);
          ++(p->stepnum);
        }
      free(tmp);
    }
  else
    operateonmesh(lmp, snthreshonmesh, sizeof(unsigned char), 0);

}



















/*********************************************************************/
/**************       Main detection function      *******************/
/*********************************************************************/
void
detectonmesh(struct noisechiselparams *p, size_t *numlabs)
{
  struct meshparams *lmp=&p->lmp;

  int verb=p->cp.verb;
  size_t s0=lmp->s0, s1=lmp->s1;
  char *detectionname=p->detectionname;
  float *imgcopy, *ave, *inputimage=p->img;


  /* In case you want to view the steps, put the correct arrays in the
     FITS file. */
  if(detectionname)
    {
      arraytofitsimg(detectionname, "Input", FLOAT_IMG, p->img,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(detectionname, "InitialDetections", LONG_IMG, p->olab,
                     s0, s1, 0, p->wcs, NULL, SPACK_STRING);
    }


  /* Find the average and STD of the undetected pixels for each
     mesh for the initial detection threshold: */
  findavestdongrid(p, p->detectionskyname, 0);


  /* Put a copy of the input array in imgcopy to operate on here and
     point p->img to that copy, not the original image. For the
     removal of false detections, it is important that the initial sky
     value be subtracted. However we don't want to touch the input
     image, since each subtraction adds noise. A copy of the pointer
     to the input array is kept in `inputimage'. After the job is
     done, this temporary copy will be freed and p->img will point to
     the actual input image again. */
  floatcopy(p->img, s0*s1, &imgcopy);
  p->img=imgcopy;


  /* Apply the false detection removal threshold to the image. */
  applydetectionthresholdskysub(p);
  if(p->detectionname)
    arraytofitsimg(detectionname, "InitalSkySubtracted", FLOAT_IMG, p->img,
                   s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
  if(verb)
    reporttiming(NULL, "Initial sky value threshold applied.", 2);

  /* Save the standard deviation image. Since we have subtracted the
     average, we only need the standard deviation. */
  checkgarray(&p->smp, &ave, &p->std);
  free(ave);

  /* Find the Signal to noise ratio threshold on the grid: */
  snthreshongrid(p);

  /* Point the input image to its correct place: */
  free(p->img);
  p->img=inputimage;

  /* Clean up: */
  free(p->std);
  free(p->dbyt);
}
