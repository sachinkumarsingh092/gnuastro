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

#include <gnuastro/fits.h>
#include <gnuastro/mesh.h>
#include <gnuastro/qsort.h>
#include <gnuastro/timing.h>
#include <gnuastro/threads.h>
#include <gnuastro/checkset.h>
#include <gnuastro/arraymanip.h>
#include <gnuastro/statistics.h>

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
  char report[GAL_TIMING_VERB_MSG_LENGTHS_2_V];
  size_t i, s0=p->smp.s0, s1=p->smp.s1;
  char *detectionname=p->detectionname;

  /* Find the threshold and apply it, if there are blank pixels in the
     image, set those pixels to the binary blank value too: */
  findapplyqthreshold(p);
  if(p->anyblank)  setbytblank(p->img, p->byt, s0*s1);
  if(detectionname)
    gal_fits_array_to_file(detectionname, "Thresholded", BYTE_IMG,
                           p->byt, s0, s1, p->anyblank, p->wcs, NULL,
                           SPACK_STRING);
  if(verb)
    {
      sprintf(report, "%.2f quantile threshold found and applied.",
              p->qthresh);
      gal_timing_report(NULL, report, 2);
    }



  /* Erode the thresholded image: */
  if(p->erodengb==4)
    for(i=0;i<p->erode;++i)
      dilate0_erode1_4con(p->byt, s0, s1, 1);
  else
    for(i=0;i<p->erode;++i)
      dilate0_erode1_8con(p->byt, s0, s1, 1);
  if(detectionname)
    gal_fits_array_to_file(detectionname, "Eroded", BYTE_IMG, p->byt,
                           s0, s1, p->anyblank, p->wcs, NULL,
                           SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Eroded %lu times (%s connectivity).",
              p->erode, p->erodengb==4 ? "4" : "8");
      gal_timing_report(NULL, report, 2);
    }



  /* Do the opening: */
  opening(p->byt, s0, s1, p->opening, p->openingngb);
  if(detectionname)
    gal_fits_array_to_file(detectionname, "Opened", BYTE_IMG, p->byt,
                           s0, s1, p->anyblank, p->wcs, NULL,
                           SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Opened (depth: %lu, %s connectivity).",
              p->opening, p->openingngb==4 ? "4" : "8");
      gal_timing_report(NULL, report, 2);
    }



  /* Label the connected regions. Note that p->olab was allocated in
     ui.c and will be freed there. */
  p->numobjects=BF_concmp(p->byt, p->olab, s0, s1, p->anyblank, 4);
  if(detectionname)
    gal_fits_array_to_file(detectionname, "Labeled", LONG_IMG, p->olab,
                           s0, s1, p->anyblank, p->wcs, NULL,
                           SPACK_STRING);
}




















/****************************************************************
 *********    Signal to noise ratio calculation      ************
 ****************************************************************/
/* Calculate the Signal to noise ratio for all the labels in p->clab
   which is keeping all the pseudo detections' labels.  However, the
   array keeping the labels of each object is only over the mesh under
   consideration. */
void
detlabelsn(struct noisechiselparams *p, size_t *numlabs, float **outsntable)
{
  float *imgss=p->imgss;
  struct gal_mesh_params *smp=&p->smp;
  double *brightnesses, err, *xys;
  float *f, *ff, ave, *sntable, cpscorr=p->cpscorr;
  size_t i, ind, xyscol=3, is1=p->smp.s1, counter=0;
  unsigned char *coversdet=NULL, b0f1=p->b0f1, *b=p->byt;
  long *areas, *lab=p->clab, *lf, detsnminarea=p->detsnminarea;


  /* Allocate the necessary arrays to keep the label information: */
  errno=0; sntable=*outsntable=calloc(*numlabs, sizeof *sntable);
  if(sntable==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for sntable in detlabelsn "
          "(detection.c)", *numlabs*sizeof *sntable);
  errno=0; brightnesses=calloc(*numlabs, sizeof *brightnesses);
  if(brightnesses==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for brightnesses in detlabelsn "
          "(detection.c)", *numlabs*sizeof *brightnesses);
  errno=0; xys=calloc(*numlabs*xyscol, sizeof *xys);
  if(xys==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for xys in detlabelsn "
          "(detection.c)", *numlabs*sizeof *xys);
  errno=0; areas=calloc(*numlabs, sizeof *areas);
  if(areas==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for areas in detlabelsn "
          "(detection.c)", *numlabs*sizeof *areas);
  if(b0f1==0)
    {
      errno=0; coversdet=calloc(*numlabs, sizeof *coversdet);
      if(coversdet==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for coversdet in detlabelsn "
              "(detection.c)", *numlabs*sizeof *coversdet);
    }


  /* Go over the area and fill in the necessary information for
     calculating the Signal to noise ratio. */
  ff = (f=imgss) + p->smp.s0*p->smp.s1;
  do
    {
      /* See label.h for this macro. It checks speciie*/
      if(ISINDEXABLELABEL)
        {
          /* This test is done first so if the sky region
             pseudo-detection lies over a detected object it doesn't
             waste our time for the next calculations. */
          if(b0f1==0)
            {
              if(coversdet[*lab]) {++lab; ++b; continue;}
              else if(*b==1)
                {coversdet[*lab]=1; areas[*lab]=0; ++lab; ++b; continue;}
            }

          if(isnan(*f))
            error(EXIT_FAILURE, 0, "*f was nan");

          ++areas[*lab];
          brightnesses[*lab]     += *f;
          if(f-imgss>0)
            {
              xys[*lab*xyscol+2] += f-imgss;
              xys[*lab*xyscol  ] += (double)((f-imgss)/is1) * *f;
              xys[*lab*xyscol+1] += (double)((f-imgss)%is1) * *f;
            }
        }
      ++lab;
      if(b0f1==0) ++b;
    }
  while(++f<ff);


  /* If the user wants to see the steps (on the background), remove
     all the pseudo-detections that will not be used in the histogram
     calcluation. */
  if(p->detectionname && b0f1==0)
    {
      b=p->dbyt;
      lf=(lab=p->clab)+p->smp.s0*p->smp.s1;
      do
        {
          if( ISINDEXABLELABEL
              && (areas[*lab]<detsnminarea || brightnesses[*lab]<0) )
              *b=0;
          ++b;
        }
      while(++lab<lf);
      gal_fits_array_to_file(p->detectionname, "For S/N", BYTE_IMG,
                             p->dbyt, p->smp.s0, p->smp.s1,
                             p->anyblank, p->wcs, NULL, SPACK_STRING);
    }



  /* calculate the signal to noise for successful detections: */
  for(i=1;i<*numlabs;++i)
    {
      ave=brightnesses[i]/areas[i];
      if(areas[i]>detsnminarea && ave>0.0f && xys[i*xyscol+2]>0.0f)
        {
          /* Find the flux weighted center of this object in each
             dimension to find the mesh it belongs to and get the
             standard deviation at this point. Note that the standard
             deviation on the grid was stored in smp->garray2. The error
             should then be taken to the power of two and if the sky is
             subtracted, a 2 should be multiplied to it.*/
          err = smp->garray2[
                 gal_mesh_img_xy_to_mesh_id(smp,
                                    xys[i*xyscol]/xys[i*xyscol+2],
                                    xys[i*xyscol+1]/xys[i*xyscol+2]) ];
          err *= p->skysubtracted ? err : 2.0f*err;

          /* Set the index in the sntable to store the Signal to noise
             ratio. When we are dealing with the noise, we only want the
             non-zero signal to noise values, so we will just use a
             counter. But for initial detections, it is very important
             that their Signal to noise ratio be placed in the same
             index as their label. */
          ind = p->b0f1 ? i : counter++;
          sntable[ind] = ( sqrt( (float)(areas[i])/cpscorr )
                           * ave / sqrt(ave+err) );
        }
    }


  /* In background mode, set the number of labels to the number of
     acceptable S/N values measured. */
  if(p->b0f1==0) *numlabs=counter;


  /* For a check. IMPORTANT NOTE: in background mode, the sntable is
     not ordered like the areas, brightnesses and coversdet
     arrays. Since indexes are irrelevant for this mode, they are just
     put contigusly in memory for the next step.
  if(b0f1)
    {
      for(i=0;i<*numlabs;++i)
        printf("%lu: %-10ld %-10.3f %-10.3f\n",
               i, areas[i], brightnesses[i], sntable[i]);
    }
  */

  /* Clean up */
  free(xys);
  free(areas);
  free(brightnesses);
  if(b0f1==0) free(coversdet);
}





/* Apply the S/N threshold on the pseudo-detections. */
void
applydetsn(struct noisechiselparams *p, float *sntable, size_t numpseudo)
{
  size_t i, curlab=1;
  unsigned char *b=p->dbyt;
  long *newlabs, *lab, *lf, *clab=p->clab;


  errno=0; newlabs=calloc(numpseudo, sizeof *newlabs);
  if(newlabs==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for numlabs in "
          "removefalsedetections (detections.c)", numpseudo*sizeof *newlabs);


  for(i=1;i<numpseudo;++i)
    if(sntable[i] > p->detsn)
      newlabs[i]=curlab++;


  lf= (lab=clab) + p->smp.s0*p->smp.s1;
  do
    {
      if(*lab!=GAL_FITS_LONG_BLANK)
        *b = newlabs[*lab] > 0;
      ++b;
    }
  while(++lab<lf);


  if(p->detectionname)
    gal_fits_array_to_file(p->detectionname, "True pseudo-detections",
                           BYTE_IMG, p->dbyt, p->smp.s0, p->smp.s1,
                           p->anyblank, p->wcs, NULL, SPACK_STRING);

  free(newlabs);
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
          *out++ = *b==b0f1 ? *d : (*b==GAL_FITS_BYTE_BLANK
                                    ? GAL_FITS_BYTE_BLANK : 0);
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






void*
detectpseudos(void *inparams)
{
  struct gal_mesh_thread_params *mtp=(struct gal_mesh_thread_params *)inparams;
  struct gal_mesh_params *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  unsigned char *mponeforall=mp->oneforall;
  unsigned char *thisbyt=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  int anyblank;
  size_t i, s0, s1;
  unsigned char *b, *bf;
  size_t ind, startind;
  char *detectionname=p->detectionname;
  size_t *indexs=&mp->indexs[mtp->id*mp->thrdcols];


  /* Do the job for each mesh: */
  for(i=0;indexs[i]!=GAL_THREADS_NON_THRD_INDEX;++i)
    {
      /* Set the necesary parameters: */
      ind=indexs[i];
      s0=mp->ts0[mp->types[ind]];
      s1=mp->ts1[mp->types[ind]];
      startind=mp->start[ind];
      anyblank=0;

      /* Copy this mesh into a separate array. In the meantime, only
         bring in the desired pixels. If the user wanted to see the
         steps, then save the result to the p->dbyt and come leave the
         work for this mesh if this is the first step. */
      bytpartfromlarge(p, thisbyt, startind, s0, s1);
      if(detectionname && p->stepnum==1)
	{ bytparttolarge(p, thisbyt, startind, s0, s1); continue; }


      /* If there were NaN pixels in the image, make sure if this
         mesh has any or not. */
      if(p->anyblank)
        {
          bf=(b=thisbyt)+s0*s1;
          do if(*b++==GAL_FITS_BYTE_BLANK) { anyblank=1; break; }
          while(b<bf);
        }


      /* Fill the bounded holes. */
      fillboundedholes(thisbyt, s0, s1, anyblank);
      if(detectionname && p->stepnum==2)
	{ bytparttolarge(p, thisbyt, startind, s0, s1); continue; }

      /* Open the image once: */
      opening(thisbyt, s0, s1, 1, 4);


      /* Put this region back into p->dbyte. */
      bytparttolarge(p, thisbyt, startind, s0, s1);
    }


  /* If multiple threads were used, wait until all other threads
     finish. */
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);


  return NULL;
}





/* Doing object detection is done separately on each mesh both for the
   first part of finding the SN limit and the second part of applying
   it and removing false detections. Therefore, if the user asks to
   view the steps, the mesh process has to be stopped at differnt
   times for all the meshs and re-done for each step. This function
   does that job.

   Here we temporarily (clab will be cleaned when it is needed)
   operate on the p->clab array which is finally used to keep the
   clump labels.
*/
void
detsnthresh(struct noisechiselparams *p)
{
  struct gal_mesh_params *lmp=&p->lmp;

  float *sntable;
  static int b0f1;
  char *extname=NULL;
  unsigned char *tmp, *originaldbyt;
  size_t numpseudo, s0=lmp->s0, s1=lmp->s1;

  /* General preparations: */
  p->b0f1 = b0f1++;

  /* All the operations below happen on p->dbyt, but we need it two
     times (once for the sky region and once for the detected
     region). So the first time we operate on it,  */
  if(p->b0f1==0)
    gal_arraymanip_uchar_copy(p->dbyt, s0*s1, &originaldbyt);

  /* Find the psudo-detections: */
  if(p->detectionname)
    {
      p->stepnum=1;
      /* Backup of p->dbyt in tmp */
      gal_arraymanip_uchar_copy(p->dbyt, s0*s1, &tmp);
      while(p->stepnum<4)
        {
          free(p->dbyt);    /* Free the old, p->dbyt, put the original */
          gal_arraymanip_uchar_copy(tmp, s0*s1, &p->dbyt);
          gal_mesh_operate_on_mesh(lmp, detectpseudos, sizeof(unsigned char),
                                   0, 0);
          switch(p->stepnum)
            {
            case 1:
              extname = p->b0f1 ? "ThresholdDetections" : "ThresholdNoise";
              break;
            case 2:
              extname="HolesFilled"; break;
            case 3:
              extname="Opened"; break;
            default:
              error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we "
                    "can find the solution. For some reason, the value to "
                    "p->stepnum (%d) is not recognized in detsnthresh "
                    "(detection.c)", PACKAGE_BUGREPORT, p->stepnum);
            }
          gal_fits_array_to_file(p->detectionname, extname, BYTE_IMG,
                                 p->dbyt, s0, s1, p->anyblank, p->wcs,
                                 NULL, SPACK_STRING);
          ++p->stepnum;
        }
      free(tmp);
    }
  else
    gal_mesh_operate_on_mesh(lmp, detectpseudos, sizeof(unsigned char), 0, 0);


  /* Find the connected components and the signal to noise ratio of
     all the detections. Note that when dealing with the background
     pixels in p->byt, the order of sntable is no longer valid, since
     we don't care any more, we just want their sorted values for the
     quantile. */
  numpseudo=BF_concmp(p->dbyt, p->clab, s0, s1, p->anyblank, 4);
  if(p->detectionname)
    gal_fits_array_to_file(p->detectionname, "Labeled", LONG_IMG,
                           p->clab, s0, s1, p->anyblank, p->wcs,
                           NULL, SPACK_STRING);
  detlabelsn(p, &numpseudo, &sntable);



  /* Find or apply the threshold depending on if we are dealing with
     the background or foreground. */
  if(p->b0f1)
    applydetsn(p, sntable, numpseudo);
  else
    {
      snthresh(p, sntable, numpseudo, 0);
      free(p->dbyt);
      p->dbyt=originaldbyt;
    }

  /* Clean up: */
  free(sntable);
}





void
dbytolaboverlap(struct noisechiselparams *p)
{
  unsigned char *byt;
  long *lf, *tokeep, *lab;
  size_t numobjects=p->numobjects;
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
     pixels. At this point, tokeep is only binary: 0 or 1.
  */
  byt=p->dbyt; lf=(lab=p->olab)+size;
  do
    {
      if(ISINDEXABLELABEL)
        {
          tokeep[ *lab ] =
            tokeep[ *lab ]	/* Check if this label is to be kept.    */
            ? 1			/* It has, its all we need!              */
            : *byt; 		/* It hasn't, check if it should be kept.*/
        }
      ++byt;
    }
  while(++lab<lf);
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
  byt=p->byt; lf=(lab=p->olab)+size;
  if(p->dilate)
    do
      {
        if(*lab!=GAL_FITS_LONG_BLANK)
          *byt=tokeep[*lab]>0;
        ++byt;
      }
    while(++lab<lf);
  else
    do
      {
        if(*lab!=GAL_FITS_LONG_BLANK)
          *byt = ( *lab = tokeep[*lab] ) > 0;
        ++byt;
      }
    while(++lab<lf);


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
  struct gal_mesh_params *lmp=&p->lmp;

  int verb=p->cp.verb;
  char report[GAL_TIMING_VERB_MSG_LENGTHS_2_V];
  char *detectionname=p->detectionname;
  size_t s0=lmp->s0, s1=lmp->s1, numobjects=p->numobjects;


  /* Find the average and STD of the undetected pixels for each
     mesh for the initial detection threshold: */
  findavestdongrid(p, p->detectionskyname);


  /* Apply the false detection removal threshold to the image. */
  applydetectionthresholdskysub(p);
  if(detectionname)
    gal_fits_array_to_file(detectionname, "InitalSkySubtracted",
                           FLOAT_IMG, p->imgss, s0, s1, p->anyblank,
                           p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Initial sky threshold (%.3f sigma) applied.",
              p->dthresh);
      gal_timing_report(NULL, report, 2);
    }


  /* Find the Signal to noise ratio threshold on the grid and keep it
     in one array. */
  detsnthresh(p);

  /* Apply the SN threshold to all the detections. */
  detsnthresh(p);

  /* Select the true detections: */
  dbytolaboverlap(p);
  if(detectionname)
    gal_fits_array_to_file(detectionname, "TrueDetections",
                           BYTE_IMG, p->byt, s0, s1, p->anyblank,
                           p->wcs, NULL, SPACK_STRING);
  if(verb)
    {            /* p->numobjects changed in dbytlaboverlap. */
      sprintf(report, "%lu false detections removed.",
              numobjects-p->numobjects);
      gal_timing_report(NULL, report, 2);
    }


  /* Clean up: */
  free(p->dbyt);
}
