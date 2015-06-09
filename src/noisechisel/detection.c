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

#include "sky.h"
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
  char *detectionname=p->detectionname;

  /* Find the threshold and apply it: */
  findapplyqthreshold(p);
  if(detectionname)
    arraytofitsimg(detectionname, "Thresholded", BYTE_IMG, p->byt,
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
  if(detectionname)
    arraytofitsimg(detectionname, "Eroded", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Eroded %lu times (%s connectivity).", p->numerosion,
              p->erodengb==4 ? "4" : "8");
      reporttiming(NULL, report, 2);
    }



  /* Do the opening: */
  opening(p->byt, s0, s1, p->opening, p->openingngb);
  if(detectionname)
    arraytofitsimg(detectionname, "Opened", BYTE_IMG, p->byt,
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
  if(detectionname)
    arraytofitsimg(detectionname, "Labeled", LONG_IMG, p->olab,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
}




















/****************************************************************
 *********    Signal to noise ratio calculation      ************
 ****************************************************************/
/* Calculate the Signal to noise ratio for all the labels in
   `labinmesh'. The sky subtracted input image and the standard
   deviations image for each pixel are over the full image: p->img and
   p->std. However, the array keeping the labels of each object is
   only over the mesh under consideration.
*/
void
detlabelsn(struct noisechiselparams *p, long *labinmesh, size_t *numlabs,
           size_t start, size_t s0, size_t s1, float **outsntable)
{
  float *imgss=p->imgss;
  double *fluxs, err, *xys;
  struct meshparams *smp=&p->smp;
  size_t i, ind, r=0, is1=p->smp.s1, counter=0;
  float *f, *ff, *s, ave, *sntable, cpscorr=p->cpscorr;
  long *areas, *l=labinmesh, detsnminarea=p->detsnminarea;


  /* Allocate the necessary arrays to keep the label information: */
  errno=0; sntable=*outsntable=calloc(*numlabs, sizeof *sntable);
  if(sntable==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for sntable in labelsn (label.c)",
          *numlabs*sizeof *sntable);
  errno=0; fluxs=calloc(*numlabs, sizeof *fluxs);
  if(fluxs==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for fluxs in labelsn (label.c)",
          *numlabs*sizeof *fluxs);
  errno=0; xys=calloc(*numlabs*2, sizeof *xys);
  if(xys==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for xys in labelsn (label.c)",
          *numlabs*sizeof *xys);
  errno=0; areas=calloc(*numlabs, sizeof *areas);
  if(areas==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for areas in labelsn (label.c)",
          *numlabs*sizeof *areas);


  /* For a check:
  arraytofitsimg("test.fits", "Labs", LONG_IMG, labinmesh,
                 s0, s1, 0, NULL, NULL, SPACK_STRING);
  */


  /* Go over the arrays and fill in the information. Note that since
     labinmesh is over the mesh and not the full image, setting it
     once in the definition was enough. */
  do
    {
      ff = ( f = imgss + start + r++ * is1 ) + s1;
      do
        {
          if(*l && !isnan(*f))     /* There is a label on this mesh. */
            {
              ++areas[*l];
              fluxs[*l]   += *f;
              xys[*l*2]   += (double)((f-imgss)/is1) * *f;
              xys[*l*2+1] += (double)((f-imgss)%is1) * *f;
            }
          ++l; ++s;
        }
      while(++f<ff);
    }
  while(r<s0);


  /* calculate the signal to noise for successful detections: */
  for(i=1;i<*numlabs;++i)
    if(areas[i]>detsnminarea && (ave=fluxs[i]/areas[i]) >0 )
      {
        /* Find the flux weighted center of this object in each
           dimension to find the mesh it belongs to and get the
           standard deviation at this point. Note that the standard
           deviation on the grid was stored in smp->garray2. The error
           should then be taken to the power of two and if the sky is
           subtracted, a 2 should be multiplied to it.*/
        err = smp->garray2[imgxytomeshid(smp, xys[i*2]/fluxs[i],
                                         xys[i*2+1]/fluxs[i])];
        err *= p->skysubtracted ? err : 2.0f*err;

        /* Set the index in the sntable to store the Signal to noise
           ratio. When we are dealing with the noise, we only want the
           non-zero signal to noise values, so we will just use a
           counter. But for initial detections, it is very important
           that their Signal to noise ratio be placed in the same
           index as their label. */
        ind = p->b0f1 ? i : counter++;
        sntable[ind]=sqrt( (float)(areas[i])/cpscorr ) * ave / sqrt(ave+err);
      }
  if(p->b0f1==0) *numlabs=counter;


  /* For a check:
  for(i=0;i<*numlabs;++i)
    printf("%lu: %-10ld %-10.3f %-10.3f %-10.3f\n",
           i, areas[i], fluxs[i], stds[i]/areas[i], sntable[i]);
  */


  /* Clean up */
  free(xys);
  free(fluxs);
  free(areas);
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





void
removefalsedetections(unsigned char *byt, long *lab, size_t size,
                      size_t numlabs, float *sntable, float minsn)
{

  long *newlabs;
  size_t i, curlab=1;

  errno=0; newlabs=calloc(numlabs, sizeof *newlabs);
  if(newlabs==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for numlabs in "
          "removefalsedetections (detections.c)", numlabs*sizeof *newlabs);

  for(i=1;i<numlabs;++i)
    if(sntable[i]>minsn)
      newlabs[i]=curlab++;
  /*
  for(i=1;i<numlabs;++i)
    printf("   %lu: %.3f %lu\n", i, sntable[i], newlabs[i]);
  */
  for(i=0;i<size;++i)
    byt[i] = newlabs[lab[i]] > 0;

  free(newlabs);
}





void*
detsnthreshonmesh(void *inparams)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparams;
  struct meshparams *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  unsigned char *mponeforall=mp->oneforall;
  unsigned char *thisbyt=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  long *thislab;
  char cline[1000];
  float *sntable, minbfrac=p->minbfrac;
  unsigned char b0f1=p->b0f1, *dbyt=p->dbyt;
  size_t i, s0, s1, minnumfalse=p->minnumfalse;
  size_t gid, detsnhistnbins=p->detsnhistnbins;
  size_t ind, size, numlabs, startind, is1=mp->s1;
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
         and allocate an array to keep the S/N values in *sntable. */
      detlabelsn(p, thislab, &numlabs, startind, s0, s1, &sntable);
      if(b0f1)                  /* Initial detections. */
        {
          /* The index is based on the initial channel-based ID (since
             it was used in mp->start). So if the interpolation and/or
             smoothing were done on the full array, then there is
             going to be a problem and we need to use gidfromchbasedid
             to get the appropriate garray-based ID.*/
          gid=gidfromchbasedid(mp, ind);

          /* A simple sanity check: */
          if(isnan(mp->garray1[gid]))
            error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                  "can fix the problem. For some reason, the minimum Signal "
                  "to noise ratio for mesh number %lu is a NaN!",
                  PACKAGE_BUGREPORT, gid);

          /* Remove the false detections and copy thisbyt into
             p->dbyt. */
          removefalsedetections(thisbyt, thislab, s0*s1, numlabs,
                                sntable, mp->garray1[gid]);
          bytparttolarge(p, thisbyt, startind, s0, s1);
        }
      else                      /* Noise (background). */
        {
          /* We are in the noise region. Note that in detlabelsn, the
             numlabs was changed to only keep those that have the good
             conditions. The order of object indexs is no longer
             important. So in detlabelsn the labels were moved and the
             non-zero elements were put contiguously beside each other
             at the start of the array. These were only done for the
             noisy regions, not detections. */
          if(numlabs<minnumfalse) {free(sntable); continue;}

          /* Sort the signal to noise array for the proper values,
             then remove any possible outliers using the cumulative
             frequency plot of the distribution: */
          qsort(sntable, numlabs, sizeof *sntable, floatincreasing);
          removeoutliers_flatcdf(sntable, &numlabs);
          if(numlabs<minnumfalse) {free(sntable); continue;}

          /* Put the signal to noise ratio quantile into the
             mp->garray1 array. Note that since garray1 was
             initialized to NaN, when a mesh doesn't reach this point,
             it will be NaN. */
          mp->garray1[ind]=sntable[indexfromquantile(numlabs, p->detquant)];

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
                      ind, p->detquant, mp->garray1[ind]);
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
detsnthreshongrid(struct noisechiselparams *p)
{
  struct meshparams *lmp=&p->lmp;

  int initialize;
  char *extname=NULL;
  unsigned char *tmp;
  size_t s0=lmp->s0, s1=lmp->s1;


  /* If lmp->garray1 (the signal-to-noise limit on each mesh) is
     allocated, then it shows that we are working on the data pixels,
     otherwise we are working on the noise pixels. */
  p->b0f1 = lmp->garray1 ? 1 : 0;
  initialize = p->b0f1 ? 0 : 1;


  /* Find the SN threshold and save the steps in a FITS file: */
  if(p->detectionname)
    {
      p->stepnum=1;
      ucharcopy(p->dbyt, s0*s1, &tmp); /* Backup of p->dbyt in tmp */
      while(p->stepnum<6)
        {
          free(p->dbyt);    /* Free the old, p->dbyt, put the original */
          ucharcopy(tmp, s0*s1, &p->dbyt);
          operateonmesh(lmp, detsnthreshonmesh, sizeof(unsigned char), 0,
                        initialize);
          switch(p->stepnum)
            {
            case 1:
              extname = p->b0f1 ? "ThresholdDetections" : "ThresholdNoise";
              break;
            case 2:
              extname="HolesFilled"; break;
            case 3:
              extname="Opened"; break;
            case 4:
              extname="SmallRemoved"; break;
            default:
              extname = p->b0f1 ? "True" : NULL; break;
            }
          if(extname)
            arraytofitsimg(p->detectionname, extname, BYTE_IMG, p->dbyt,
                           s0, s1, 0, p->wcs, NULL, SPACK_STRING);
          else
            break;
          ++p->stepnum;
        }
      if(p->b0f1)
        free(tmp);
      else
        {
          free(p->dbyt);
          p->dbyt=tmp;
        }
    }
  else
    operateonmesh(lmp, detsnthreshonmesh, sizeof(unsigned char), 0,
                  initialize);

}





/* The p->lmp.garray1 has been filled with the successful mesh
   elements. The ones with an unsuccessful value are blank. This
   function will first interpolate the grid to have no blank
   elements. Then it will smooth it. It is used both by the detection
   process and the segmentation process to find an S/N value for
   each grid element.*/
void
findsnthreshongrid(struct meshparams *lmp, char *filename,
                   char *comment, struct wcsprm *wcs)
{
  float *sn;
  size_t s0=lmp->s0, s1=lmp->s1;


  /* Find the Signal to noise ratio threshold for the good meshs. */
  if(filename)
    {
      checkgarray(lmp, &sn, NULL);
      arraytofitsimg(filename, "S/N", FLOAT_IMG, sn, s0, s1,
                     0, wcs, NULL, SPACK_STRING);
      free(sn);
    }


  /* Interpolate over the meshs to fill all the blank ones in both the
     sky and the standard deviation arrays: */
  meshinterpolate(lmp, comment);
  if(filename)
    {
      checkgarray(lmp, &sn, NULL);
      arraytofitsimg(filename, "Interpolated", FLOAT_IMG, sn,
                     s0, s1, 0, wcs, NULL, SPACK_STRING);
      free(sn);
    }


  /* Smooth the interpolated array:  */
  if(lmp->smoothwidth>1)
    {
      meshsmooth(lmp);
      if(filename)
        {
          checkgarray(lmp, &sn, NULL);
          arraytofitsimg(filename, "Smoothed", FLOAT_IMG, sn,
                         s0, s1, 0, wcs, NULL, SPACK_STRING);
          free(sn);
        }
    }
}





void
dbytolaboverlap(struct noisechiselparams *p)
{
  long *tokeep, *olab=p->olab;
  size_t numobjects=p->numobjects;
  unsigned char *dbyt=p->dbyt, *byt=p->byt;
  size_t i, size=p->smp.s0*p->smp.s1, curlab=1;

  /* Allocate array to keep the overlapping labels. */
  errno=0; tokeep=calloc(numobjects, sizeof*tokeep);
  if(tokeep==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for tokeep in dbytolaboverlap "
          "(detection.c)", numobjects*sizeof*tokeep);

  /* Note that the zeroth element of tokeep can also be non zero, this
     is because the holes of the labeled regions might be filled
     during filling the holes, but have not been filled in the
     original labeled array. They are not important so you can just
     ignore them.

     The basic idea is that after the first pixel of a label overlaps
     with dbyt[i], we don't need to check the rest of that object's
     pixels. tokeep is only binary: 0 or 1.
  */
  for(i=0;i<size;++i)
    tokeep[ olab[i] ] =
      tokeep[ olab[i] ]		/* Check if this label is to be kept.    */
      ? 1			/* It has, its all we need!              */
      : dbyt[i]; 		/* It hasn't, check if it should be kept.*/
  tokeep[0]=0;

  /* Set the new labels: */
  for(i=1;i<numobjects;++i)
    if(tokeep[i])
      tokeep[i]=curlab++;

  /* Replace the byt and olab values with their proper values. If the
     user doesn't want to dilate, then change the labels in `lab'
     too. Otherwise, you don't have to worry about the label
     array. After dilation a new labeling will be done and the whole
     labeled array will be re-filled.*/
  if(p->dilate)
    for(i=0;i<size;++i)
      byt[i] = tokeep[ olab[i] ] > 0;
  else
    for(i=0;i<size;++i)
      byt[i] = (olab[i] = tokeep[ olab[i] ]) > 0;

  /* Note that until here, numdetected is the number of labels
     assigned. But since later we will want the labels to be
     indexes in an array, numdetected has to be added with one. */
  p->numobjects=curlab;

  free(tokeep);
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
  char report[VERBMSGLENGTHS2_V];
  char *detectionname=p->detectionname;
  size_t s0=lmp->s0, s1=lmp->s1, numobjects=p->numobjects;


  /* Find the average and STD of the undetected pixels for each
     mesh for the initial detection threshold: */
  findavestdongrid(p, p->detectionskyname);


  /* Apply the false detection removal threshold to the image. */
  applydetectionthresholdskysub(p);
  if(p->detectionname)
    arraytofitsimg(detectionname, "InitalSkySubtracted", FLOAT_IMG,
                   p->imgss, s0, s1, p->numblank, p->wcs, NULL,
                   SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Initial sky threshold (%.3f sigma) applied.",
              p->dthresh);
      reporttiming(NULL, report, 2);
    }


  /* Find the Signal to noise ratio threshold on the grid and keep it
     in one array. */
  detsnthreshongrid(p);
  findsnthreshongrid(&p->lmp, p->detectionsnname, "Interpolating the "
                     "DETECTION Signal to noise ratio threshold", p->wcs);
  if(verb)
    {
      snave=floataverage(lmp->garray1, lmp->nmeshi);
      sprintf(report, "Detection S/N limit found (Average: %.3f).",
              snave);
      reporttiming(NULL, report, 2);
    }


  /* Apply the SN threshold to all the detections. */
  detsnthreshongrid(p);
  dbytolaboverlap(p);
  if(detectionname)
    arraytofitsimg(detectionname, "TrueDetections", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {            /* p->numobjects changed in dbytlaboverlap. */
      sprintf(report, "%lu false detections removed.",
              numobjects-p->numobjects);
      reporttiming(NULL, report, 2);
    }


  /* Clean up: */
  free(p->dbyt);
}
