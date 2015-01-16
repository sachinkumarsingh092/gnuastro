/********************************************************************
mkprof (MakeProfiles) - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include "box.h"
#include "timing.h"
#include "checkset.h"
#include "statistics.h"
#include "arraymanip.h"
#include "astrthreads.h"
#include "fitsarrayvv.h"

#include "main.h"
#include "mkprof.h"
#include "oneprofile.h"


















/**************************************************************/
/************        builtqueue linked list       *************/
/**************************************************************/
/* Add an empty element to the existing builtqueue. */
void
builtqueue_addempty(struct builtqueue **bq)
{
  struct builtqueue *tbq;

  /* Allocate the element. */
  errno=0;
  tbq=malloc(sizeof *tbq);
  if(tbq==NULL)
    error(EXIT_FAILURE, 0, "%lu byte element in builtqueue_addempty.",
	  sizeof *tbq);
  tbq->img=NULL;

  /* Set its next element to the input bq and re-set the input bq. */
  tbq->next=*bq;
  *bq=tbq;
}




















/**************************************************************/
/************           Save individual           *************/
/**************************************************************/
#define NUMBERNAMESTRLEN 100
void
saveindividual(struct mkonthread *mkp)
{
  size_t len;
  char *jobname;
  char *outname, *outdir=mkp->p->cp.output;

  /* Save the array to an image. */
  if(mkp->p->dir0file1==0)
    len=strlen(outdir)+NUMBERNAMESTRLEN;
  else
    {
      outdir="";
      len=NUMBERNAMESTRLEN;
    }

  /* Allocate the space for the name of this file. */
  errno=0;
  outname=malloc(len*sizeof *outname);
  if(outname==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for name of object in "
	  "row %lu of %s", len*sizeof *outname, mkp->ibq->id,
	  mkp->p->up.catname);

  /* Write the name and save the FITS image: */
  sprintf(outname, "%s%lu.fits", outdir, mkp->ibq->id+1);
  checkremovefile(outname, mkp->p->cp.dontdelete);
  arraytofitsimg(outname, "MockImg", FLOAT_IMG, mkp->ibq->img,
		 mkp->p->oversample*mkp->width[1],
		 mkp->p->oversample*mkp->width[0],
		 NULL, SPACK_STRING);

  /* Report if in verbose mode. */
  if(mkp->p->cp.verb)
    {
      errno=0;
      jobname=malloc(len+100*sizeof *jobname);
      if(jobname==NULL)	error(EXIT_FAILURE, errno, "jobname in mkprof.c");
      sprintf(jobname, "%s created.", outname);
      reporttiming(NULL, jobname, 2);
      free(jobname);
    }
  free(outname);

  /* Free the image array. */
  free(mkp->ibq->img);
  mkp->ibq->img=NULL;
}



















/**************************************************************/
/************            The builders             *************/
/**************************************************************/
/* Build the profiles that are indexed in the indexs array of the
   mkonthread structure that was assigned to it.

   See the explanation above overlap (/lib/box.c) for a complete
   explanation of fpixel_i, lpixel_i, fpixel_o and lpixel_o.

   =========================================================
   About the Central x and y of each profile:

   The user has asked for the profile to be built on the coordinates
   (real numbers) of `x` and `y` in an output image in the FITS
   format. We are building the full image for each galaxy separately
   in an array with an odd number of sides which maybe oversampled.

   In the FITS format, the pixel centers have an integer value. So for
   example in 1D, a pixel whose center value is 10.00 covers the area
   of: [9.5,10.5). We want the fractional part of `x` (don't forget,
   this example is 1D) to be in the central pixel of this separate
   array (with odd sides) that we will be building.

   The result of this convention is that in 1D, a continuous space
   pixel with a fractional value of 0.1 is going to be after the
   central pixel's center, while one with a fractional value of 0.9
   will be before it. In this manner, later on, when we want to find
   the overlap between this array and the output array, if we have a
   fractional value >=0.5, we will just shift the integer part of the
   central pixel by one and completely ignore the fractional part.
*/
void *
build(void *inparam)
{
  struct mkonthread *mkp=(struct mkonthread *)inparam;
  struct mkprofparams *p=mkp->p;

  size_t i;
  int lockresult;
  double *cat, pixfrac, junk;
  struct builtqueue *fbq=NULL;
  pthread_mutex_t *qlock=&p->qlock;
  pthread_cond_t *qready=&p->qready;
  long os=p->oversample, lpixel_o[2];

  /* Make each profile that was specified for this thread. */
  for(i=0;mkp->indexs[i]!=NONTHRDINDEX;++i)
    {
      /* Create a new builtqueue element with all the information. fbq
	 will be used when we want to add ibq to p->bq. It is defined
	 so we don't have to waste time traversing the ibq. Its
	 characteristic compared to the other elements of ibq is that
	 fbq->next==NULL. So to add ibq to p->bq, we just have to set
	 fbq->next=p->bq and then set p->bq to ibq.*/
      builtqueue_addempty(&mkp->ibq);
      mkp->ibq->id=mkp->indexs[i];
      if(fbq==NULL) fbq=mkp->ibq;
      cat=&p->cat[mkp->ibq->id*p->cs1];


      /* Prepare the parameter for building the profile.*/
      setprofparams(mkp);


      /* Find the bounding box size. */
      if((int)cat[p->fcol]==POINTCODE)
	mkp->width[0]=mkp->width[1]=1;
      else
	ellipseinbox(mkp->truncr, mkp->q*mkp->truncr,
		     cat[p->pcol]*DEGREESTORADIANS, mkp->width);


      /* Based on the width, set the central positions of this profile
	 in the array that will be built, see the comments above for a
	 full explanation.*/
      pixfrac = modf(fabs(cat[p->xcol]), &junk);
      mkp->yc = ( os * (mkp->width[0]/2 + pixfrac)
		  + (pixfrac<0.50f ? os/2 : -1*os/2-1) );
      mkp->yc = round(mkp->yc*100)/100;

      pixfrac = modf(fabs(cat[p->ycol]), &junk);;
      mkp->xc = ( os*(mkp->width[1]/2 + pixfrac)
		  + (pixfrac<0.5f ? os/2 : -1*os/2-1) );
      mkp->xc = round(mkp->xc*100)/100;

      /* Get the overlapping pixels using the starting points. */
      borderfromcenter(cat[p->xcol], cat[p->ycol], mkp->width,
		       mkp->ibq->fpixel_i, mkp->ibq->lpixel_i);
      if(p->individual
	 || overlap(p->naxes, mkp->ibq->fpixel_i, mkp->ibq->lpixel_i,
		    mkp->ibq->fpixel_o, lpixel_o))
	makeoneprofile(mkp);


      /* For individually built images. */
      if( p->individual || (mkp->ispsf && p->psfinimg==0) )
	saveindividual(mkp);

      /* Add ibq to bq. */
      if(p->cp.numthreads>1)
	{
	  /* Try locking the mutex so no thread can change the value
	     of p->bq. If you can lock it, then put the internal
	     builtqueue on top of the general builtqueue. If you
	     can't, continue adding to the internal builtqueue (make
	     the next profiles) until you find a chance to lock the
	     mutex. */
	  lockresult=pthread_mutex_trylock(qlock);
	  if(lockresult==0)		/* Mutex was successfully locked. */
	    {
	      /* Add this internal queue to system queue. */
	      fbq->next=p->bq;
	      p->bq=mkp->ibq;

	      /* If the list was empty when you locked the mutex, then
		 either `write` is waiting behind a condition variable
		 for you to fill it up or not (either it hasn't got to
		 setting the condition variable yet (this function
		 locked the mutex before `write`) or it just got the
		 list to be made and is busy writing the arrays in the
		 output). In either case, pthread_cond_signal will
		 work. */
	      if(fbq->next==NULL)
		pthread_cond_signal(qready);
	      pthread_mutex_unlock(qlock);

	      /* Finally set both the internal queue and the first
		 internal queue element to NULL.*/
	      fbq=NULL;
	      mkp->ibq=NULL;
	    }
	}
    }

  /* Wait until all other threads finish. */
  if(p->cp.numthreads==1)
    p->bq=mkp->ibq;
  else
    pthread_barrier_wait(mkp->b);

  return NULL;
}




















/**************************************************************/
/************              The writer             *************/
/**************************************************************/
void
write(struct mkprofparams *p)
{
  size_t complete=0;
  struct builtqueue *ibq=NULL, *tbq;

  /* Allocate the output array. */

  /* Write each image into the output array. */
  while(complete<p->cs0)
    {
      /* Set ibq. */
      if(ibq==NULL)
	{
	  if(p->cp.numthreads==1)
	    ibq=p->bq;
	  else
	    {
	      pthread_mutex_lock(&p->qlock);
	      while(p->bq==NULL)
		pthread_cond_wait(&p->qready, &p->qlock);
	      ibq=p->bq;
	      p->bq=NULL;
	      pthread_mutex_unlock(&p->qlock);
	    }
	}

      /* Write the array pointed to by ibq into the output image and
	 also fill the log. */

      /* Free the array and the queue element and change it to the
	 next one and increment complete. Note that there is no
	 problem to free a NULL pointer (when the built array didn't
	 overlap). */
      free(ibq->img);
      tbq=ibq->next;
      free(ibq);
      ibq=tbq;
      ++complete;
    }

  /* Write the final array to the final FITS image. */
}




















/**************************************************************/
/************           Outside function          *************/
/**************************************************************/
void
mkprof(struct mkprofparams *p)
{
  int err;
  size_t nt;
  pthread_t t;		 /* Thread id not used, all are saved here. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct mkonthread *mkp;
  size_t i, *indexs, thrdcols;

  /* Set the number of threads. In the multi-threaded case, we want to
     use numthreads-1 to built and one (this main one) to write. */
  nt = p->cp.numthreads==1 ? 1 : p->cp.numthreads-1;

  /* Allocate the arrays to keep the thread and parameters for each
     thread. Note that we only want nt-1 threads to do the
     building. */
  errno=0;
  mkp=malloc(nt*sizeof *mkp);
  if(mkp==NULL)
    error(EXIT_FAILURE, errno,
	  "%lu bytes in mkprof (mkprof.c) for mkp", (nt-1)*sizeof *mkp);

  /* Distribute the different profiles for different threads. Note
     that one thread is left out for writing, while nt-1 are left
     for building. */
  distinthreads(p->cs0, nt, &indexs, &thrdcols);

  /* Build the profiles: */
  if(nt==1)
    {
      mkp[0].p=p;
      mkp[0].indexs=indexs;
      build(&mkp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that nt was set to the number
	 of builder threads. This main thread will write the results,
	 so we need nt+1 threads to finish. */
      attrbarrierinit(&attr, &b, nt+1);

      /* Initialize the condition variable and mutex. */
      err=pthread_mutex_init(&p->qlock, NULL);
      if(err) error(EXIT_FAILURE, 0, "Mutex not initialized.");
      err=pthread_cond_init(&p->qready, NULL);
      if(err) error(EXIT_FAILURE, 0, "Condition variable not initialized.");

      /* Spin off the threads: */
      for(i=0;i<nt;++i)
	if(indexs[i*thrdcols]!=NONTHRDINDEX)
	  {
	    mkp[i].p=p;
	    mkp[i].b=&b;
	    mkp[i].ibq=NULL;
	    mkp[i].indexs=&indexs[i*thrdcols];
	    err=pthread_create(&t, &attr, build, &mkp[i]);
	    if(err)
	      error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }
    }


  /* Write the created arrays into the image. */
  write(p);


  /* If numthreads>1, then wait for all the jobs to finish and destroy
     the attribute and barrier. */
  if(nt>1)
    {
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
      pthread_cond_destroy(&p->qready);
      pthread_mutex_destroy(&p->qlock);
    }
  free(mkp);
}
