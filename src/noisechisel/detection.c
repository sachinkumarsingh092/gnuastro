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
#include "thresh.h"
#include "binary.h"
#include "detection.h"
















/****************************************************************
 ************           Initial detection            ************
 ****************************************************************/
void
initialdetection(struct noisechiselparams *p)
{
  int verb=p->cp.verb;
  char report[VERBMSGLENGTHS2_V];
  size_t i, s0=p->smp.s0, s1=p->smp.s1;
  char *initdetectionname=p->initdetectionname;

  /* Find the threshold and apply it: */
  findapplyqthreshold(p);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Thresholded", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "%.2f quantile threshold found and applied.",
              p->qthresh);
      reporttiming(NULL, report, 2);
    }


  /* Erode the thresholded image: */
  if(p->erodengb==4)
    for(i=0;i<p->numerosion;++i)
      dilate0_erode1_4con(p->byt, s0, s1, 1);
  else
    for(i=0;i<p->numerosion;++i)
      dilate0_erode1_8con(p->byt, s0, s1, 1);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Eroded", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Eroded %lu times (%s connectivity).", p->numerosion,
              p->erodengb==4 ? "4" : "8");
      reporttiming(NULL, report, 2);
    }



  /* Do the opening: */
  opening(p->byt, s0, s1, p->opening, p->openingngb);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Opened", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Opened (depth: %lu, %s connectivity).",
              p->opening, p->openingngb==4 ? "4" : "8");
      reporttiming(NULL, report, 2);
    }



  /* Label the connected regions. Note that p->olab was allocated in
     ui.c and will be freed there. */
  p->numobjects=BF_concmp(p->byt, p->olab, s0, s1, 4);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Labeled", LONG_IMG, p->olab,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
}




















/****************************************************************
 ************        Remove false detections         ************
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




/*
void
removefalsedetections(unsigned char *byt, long *lab, size_t size,
                      size_t *numlabs)
{

  long *newlabs;
  size_t i, curlab=1;

  errno=0; newlabs=calloc(*numlabs, sizeof *newlabs);
  if(newlabs==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for numlabs in "
          "removefalsedetections (detections.c)", *numlabs*sizeof *newlabs);

  for(i=1;i<*numlabs;++i)
    if(vartable[i]>minvar)
      newlabs[i]=curlab++;

  if(byt)
    for(i=0;i<size;++i)
      byt[i] = (in[i]=newlabs[in[i]]) > 0;
  else
    for(i=0;i<size;++i)
      in[i]=newlabs[in[i]];

  free(newlabs);
}
*/




void*
snthreshonmesh(void *inparams)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparams;
  struct meshparams *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  unsigned char *mponeforall=mp->oneforall;
  unsigned char *thisbyt=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  long *thislab;
  char cline[1000];
  size_t detsnhistnbins=p->detsnhistnbins;
  unsigned char b0f1=p->b0f1, *dbyt=p->dbyt;
  size_t i, s0, s1, minnumfalse=p->minnumfalse;
  size_t ind, size, numlabs, startind, is1=mp->s1;
  float *f, *ff, *to, *sntable, minbfrac=p->minbfrac;
  size_t nf, nb, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  char suffix[50], *histname, *detectionname=p->detectionname;


  /* Allocate all the necessary arrays: */
  errno=0; thislab=malloc(mp->maxs0*mp->maxs1*sizeof *thislab);
  if(thislab==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for thislab (detection.c)",
          mp->maxs0*mp->maxs1*sizeof *thislab);


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
      if(b0f1==0)
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


      /* Find the connected components. Detections that are smaller
         than detsnminarea and the meshs that don't have enough noise
         detections will be removed in next step when the user doesn't
         want to view the steps. However, when the user does want to
         view the steps, they will be small detections and the meshs
         will be removed for them to visually understand the condition
         better. */
      numlabs=BF_concmp(thisbyt, thislab, s0, s1, 4);
      if(detectionname && p->stepnum==4)
	{
          removesmallarea_relabel(thislab, thisbyt, size, &numlabs,
                                  p->detsnminarea);
          if(b0f1==0 && numlabs<minnumfalse)
            ucharinitonregion(dbyt, 0, startind, s0, s1, is1);
          bytparttolarge(p, thisbyt, startind, s0, s1); continue;
        }


      /* Find the signal to noise of all the detections in this mesh
         and allocate an array to keep the S/N values in *sntable. If
         we are on the noise, then check to see if we have enough
         points or not. */
      detlabelsn(p, thislab, numlabs, startind, s0, s1, &sntable);
      if(b0f1)                  /* Initial detections. */
        {
          /*
          removefalsedetections(thisbyt, thislab, s0*s1, &numlabs);
          */
        }
      else                      /* Noise (background). */
        {
          /* We are in the noise region. From this point on, the order
             of elements in sntable doesn't matter. So put all the
             non-zero elements contiguously beside each other at the
             start of the array. Also check if the number is not
             enough, then go onto the next mesh. */
          to=sntable; ff=(f=sntable)+numlabs;
          do if(*f>0) *to++=*f; while(++f<ff);
          numlabs=to-sntable;
          if(numlabs<minnumfalse) {free(sntable); continue;}

          /* Sort the signal to noise array for the proper values,
             then remove any possible outliers using the cumulative
             frequency plot of the distribution: */
          qsort(sntable, numlabs, sizeof *sntable, floatincreasing);
          removeoutliers_flatcdf(sntable, &numlabs);
          if(numlabs<minnumfalse) {free(sntable); continue;}

          /* Put the signal to noise ratio quantile into the
             mp->garray1 array. Note that since garray1 was
             initialized to NaN, when a mesh doesn't reach this
             point, it will be NaN. */
          mp->cgarray1[ind]=sntable[indexfromquantile(numlabs, p->detquant)];

          /* If the user has asked for it, make the histogram of the
             S/N distribution. */
          if(detsnhistnbins)
            {
              /* For a check:
              if(ind!=0) continue;
              ff=(f=sntable)+numlabs; do printf("%f\n", *f++); while(f<ff);
              */

              /* histname has to be set to NULL so automaticoutput can
                 safey free it. */
              histname=NULL;
              sprintf(suffix, "_%lu_detsn.txt", ind);
              sprintf(cline, "# %s\n# %s started on %s"
                      "# Input: %s (hdu: %s)\n"
                      "# Histogram for S/N distribution of false "
                      "detections.\n"
                      "# On large mesh id %lu.\n"
                      "# The %.3f quantile has a value of %.4f on "
                      "this bin.", SPACK_STRING, SPACK_NAME,
                      ctime(&p->rawtime), p->up.inputname, p->cp.hdu,
                      ind, p->detquant, mp->cgarray1[ind]);
              automaticoutput(p->up.inputname, suffix, p->cp.removedirinfo,
                              p->cp.dontdelete, &histname);
              savehist(sntable, numlabs, detsnhistnbins,
                       histname, cline);
              free(histname);
            }
        }

      /* Clean up: */
      free(sntable);
    }

  /* Free any allocated space and if multiple threads were used, wait
     until all other threads finish. */
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
      while( (p->b0f1==0 && p->stepnum<6)           /* Undetected. */
             || (p->b0f1 && p->stepnum<7))          /* Detected.   */
        {
          free(p->dbyt);    /* Free the old, p->dbyt, put the original */
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
              extname=NULL;
              break;
            }
          if(extname)
            arraytofitsimg(p->detectionname, extname, BYTE_IMG, p->dbyt,
                           s0, s1, 0, p->wcs, NULL, SPACK_STRING);
          else
            break;
          ++(p->stepnum);
        }
      free(p->dbyt);
      p->dbyt=tmp;
    }
  else
    operateonmesh(lmp, snthreshonmesh, sizeof(unsigned char), 0);

}





void
findsnthreshongrid(struct noisechiselparams *p)
{
  struct meshparams *lmp=&p->lmp;

  float *sn;
  size_t s0=lmp->s0, s1=lmp->s1;


  /* Find the Signal to noise ratio threshold for the good meshs. */
  snthreshongrid(p);
  if(p->detectionsnname)
    {
      checkgarray(lmp, &sn, NULL);
      arraytofitsimg(p->detectionsnname, "S/N", FLOAT_IMG, sn, s0, s1,
                     0, p->wcs, NULL, SPACK_STRING);
      free(sn);
    }


  /* Interpolate over the meshs to fill all the blank ones in both the
     sky and the standard deviation arrays: */
  meshinterpolate(lmp);
  if(p->detectionsnname)
    {
      checkgarray(lmp, &sn, NULL);
      arraytofitsimg(p->detectionsnname, "Interpolated", FLOAT_IMG, sn,
                     s0, s1, 0, p->wcs, NULL, SPACK_STRING);
      free(sn);
    }


  /* Smooth the interpolated array:  */
  if(lmp->smoothwidth>1)
    {
      meshsmooth(lmp);
      if(p->detectionsnname)
        {
          checkgarray(lmp, &sn, NULL);
          arraytofitsimg(p->detectionsnname,"Smoothed", FLOAT_IMG, sn,
                         s0, s1, 0, p->wcs, NULL, SPACK_STRING);
          free(sn);
        }
    }
}


















/*********************************************************************/
/**************       Main detection function      *******************/
/*********************************************************************/
void
onlytruedetections(struct noisechiselparams *p)
{
  struct meshparams *lmp=&p->lmp;

  float snave;
  int verb=p->cp.verb;
  size_t s0=lmp->s0, s1=lmp->s1;
  char report[VERBMSGLENGTHS2_V];
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
  if(verb)
    reporttiming(NULL, "Initial sky and its STD found.", 2);


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

  /* Find the Signal to noise ratio threshold on the grid and keep it
     in one array. */
  findsnthreshongrid(p);
  if(verb)
    {
      snave=floataverage(lmp->garray1, lmp->nmeshi);
      sprintf(report, "S/N limit found on larger grid (Average: %.3f).",
              snave);
      reporttiming(NULL, report, 2);
    }

  /* Apply the SN threshold to all the detections. */
  snthreshongrid(p);

  /* Point the input image to its correct place: */
  free(p->img);
  p->img=inputimage;

  /* Clean up: */
  free(p->std); p->std=NULL;
  free(p->dbyt);
}
