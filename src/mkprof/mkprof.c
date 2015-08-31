/********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include "txtarrayvv.h"
#include "astrthreads.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "mkprof.h"             /* Needs main.h astrthreads.h */
#include "oneprofile.h"

















/**************************************************************/
/************        Prepare WCS parameters       *************/
/**************************************************************/
void
preparewcs(struct mkprofparams *p)
{
  int status;
  struct wcsprm wcs;
  long os=p->oversample;

  /* Initialize the structure (allocate all the arrays). */
  wcs.flag=-1;
  if( (status=wcsini(1, 2, &wcs)) )
    error(EXIT_FAILURE, 0, "wcsinit error %d: %s.",
	  status, wcs_errmsg[status]);

  /* Correct the CRPIX values. */
  p->crpix[0]=p->crpix[0]*os+p->shift[0]-os/2;
  p->crpix[1]=p->crpix[1]*os+p->shift[1]-os/2;

  /* Fill in all the important input array values. */
  wcs.equinox=2000.0f;
  wcs.crpix[0]=p->crpix[0];
  wcs.crpix[1]=p->crpix[1];
  wcs.crval[0]=p->crval[0];
  wcs.crval[1]=p->crval[1];
  wcs.pc[0]=-1.0f*p->resolution/3600/p->oversample;
  wcs.pc[3]=p->resolution/3600/p->oversample;
  wcs.pc[1]=wcs.pc[2]=0.0f;
  wcs.cdelt[0]=wcs.cdelt[1]=1.0f;
  strcpy(wcs.cunit[0], "deg");
  strcpy(wcs.cunit[1], "deg");
  strcpy(wcs.ctype[0], "RA---TAN");
  strcpy(wcs.ctype[1], "DEC--TAN");

  /* Set up the wcs structure: */
  if( (status=wcsset(&wcs)) )
    error(EXIT_FAILURE, 0, "wcsset error %d: %s.", status,
	  wcs_errmsg[status]);

  /* Write the WCS structure to a header string. */
  if( (status=wcshdo(WCSHDO_safe, &wcs, &p->wcsnkeyrec, &p->wcsheader)) )
    error(EXIT_FAILURE, 0, "wcshdo error %d: %s.", status,
	  wcs_errmsg[status]);

  /* Free the allocated spaces by wcsini/wcsset: */
  if( (status=wcsfree(&wcs)) )
    error(EXIT_FAILURE, 0, "wcsfree error %d: %s.", status,
	  wcs_errmsg[status]);
}




















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

  /* Initialize some of the values. */
  tbq->img=NULL;
  tbq->numaccu=0;
  tbq->accufrac=0.0f;
  tbq->indivcreated=0;

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
  struct mkprofparams *p=mkp->p;

  size_t len;
  char *jobname;
  double crpix[2];
  long os=p->oversample;
  struct builtqueue *ibq=mkp->ibq;
  char *outname, *outdir=p->cp.output;

  /* Save the array to an image. */
  if(p->dir0file1==0)
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
	  "row %lu of %s", len*sizeof *outname, ibq->id,
	  p->up.catname);

  /* Save the correct CRPIX values: */
  crpix[0] = p->crpix[0] - os*(mkp->fpixel_i[0]-1);
  crpix[1] = p->crpix[1] - os*(mkp->fpixel_i[1]-1);

  /* Write the name and save the FITS image: */
  sprintf(outname, "%s%lu.fits", outdir, ibq->id);
  checkremovefile(outname, p->cp.dontdelete);

  /* Change NaN values to 0.0f: */
  freplacevalue(ibq->img, mkp->width[1]*mkp->width[0], NAN, 0.0f);
  if(p->setconsttonan)
    freplacevalue(ibq->img, mkp->width[1]*mkp->width[0], CONSTFORNAN, NAN);

  /* Write the array to file (A separately built PSF doesn't need WCS
     coordinates): */
  if(ibq->ispsf && p->psfinimg==0)
    arraytofitsimg(outname, "MockImg", FLOAT_IMG, ibq->img,
		   mkp->width[1], mkp->width[0], 0, NULL, NULL,
                   SPACK_STRING);
  else
    atofcorrectwcs(outname, "MockImg", FLOAT_IMG, ibq->img,
		   mkp->width[1], mkp->width[0], p->wcsheader,
		   p->wcsnkeyrec, crpix, SPACK_STRING);
  ibq->indivcreated=1;

  /* Change 0.0f values to NAN: */
  if(p->setconsttonan)
    freplacevalue(ibq->img, mkp->width[1]*mkp->width[0], NAN, CONSTFORNAN);
  freplacevalue(ibq->img, mkp->width[1]*mkp->width[0], 0.0f, NAN);

  /* Report if in verbose mode. */
  if(p->cp.verb)
    {
      errno=0;
      jobname=malloc((len+100)*sizeof *jobname);
      if(jobname==NULL)	error(EXIT_FAILURE, errno, "jobname in mkprof.c");
      sprintf(jobname, "%s created.", outname);
      reporttiming(NULL, jobname, 2);
      free(jobname);
    }
  free(outname);
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
  double *cat;
  int lockresult;
  long lpixel_o[2];
  pthread_mutex_t *qlock=&p->qlock;
  struct builtqueue *ibq, *fbq=NULL;
  pthread_cond_t *qready=&p->qready;

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
      ibq=mkp->ibq;
      ibq->id=mkp->indexs[i];
      if(fbq==NULL) fbq=ibq;
      cat=&p->cat[ibq->id*p->cs1];


      /* Write the necessary parameters for this profile into mkp.*/
      setprofparams(mkp);


      /* Find the bounding box size (NOT oversampled). */
      if((int)cat[p->fcol]==POINTCODE)
	mkp->width[0]=mkp->width[1]=1;
      else
	ellipseinbox(mkp->truncr, mkp->q*mkp->truncr,
		     cat[p->pcol]*DEGREESTORADIANS, mkp->width);


      /* Get the overlapping pixels using the starting points (NOT
	 oversampled). */
      borderfromcenter(cat[p->xcol], cat[p->ycol], mkp->width,
		       ibq->fpixel_i, ibq->lpixel_i);
      mkp->fpixel_i[0]=ibq->fpixel_i[0];
      mkp->fpixel_i[1]=ibq->fpixel_i[1];
      ibq->overlaps = overlap(mkp->onaxes, ibq->fpixel_i, ibq->lpixel_i,
			      ibq->fpixel_o, lpixel_o);


      /* Build the profile if necessary, After this, the width is
	 oversampled. */
      if(ibq->overlaps || p->individual || (ibq->ispsf && p->psfinimg==0))
	{
          /* Put a copy of the main random number generator for this
             thread to use for this profile. */
          gsl_rng_memcpy(mkp->rng, p->rng);

          /* Set the seed of the random number generator if the
             environment is not to be used. */
          if(mkp->p->envseed==0)
            gsl_rng_set(mkp->rng, timebasedrngseed());

          /* Make the profile */
	  makeoneprofile(mkp);
	  if( p->individual || (ibq->ispsf && p->psfinimg==0))
	    {
	      saveindividual(mkp);
	      if(ibq->ispsf && p->psfinimg==0)
		ibq->overlaps=0;
	    }
	}

      /* Add ibq to bq if you can lock the mutex. */
      if(p->cp.numthreads>1)
	{
	  /* Try locking the mutex so no thread can change the value
	     of p->bq. If you can lock it, then put the internal
	     builtqueue on top of the general builtqueue. If you
	     can't, continue adding to the internal builtqueue (make
	     the next profiles) until you find a chance to lock the
	     mutex. */
	  lockresult=pthread_mutex_trylock(qlock);
	  if(lockresult==0)     /* Mutex was successfully locked. */
	    {
	      /* Add this internal queue to system queue. */
	      fbq->next=p->bq;
	      p->bq=ibq;

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
          /* The mutex couldn't be locked and there are no more
             objects for this thread to build (giving a chance for
             this thread to add up its built profiles). So we have to
             lock the mutex to pass on this built structure to the
             builtqueue. */
          else if (mkp->indexs[i+1]==NONTHRDINDEX)
            {
	      pthread_mutex_lock(qlock);
	      fbq->next=p->bq;
	      p->bq=ibq;
              pthread_cond_signal(qready);
	      pthread_mutex_unlock(qlock);
            }
	}
    }

  /* Free the allocated space for this thread and wait until all other
     threads finish. */
  gsl_rng_free(mkp->rng);
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
writelog(struct mkprofparams *p)
{
  char comments[1000];
  int space[]={6, 10, 15}, prec[]={3, 6};
  int int_cols[]={0, 2, 4, -1}, accu_cols[]={-1};

  sprintf(comments, "# Log file for "SPACK_STRING".\n"
	  "# Run on %s"
	  "# Column 0: Row number in catalog (starting from zero).\n"
	  "# Column 1: Overlap magnitude with final image "
	  "(zeropoint: %.3f).\n"
	  "# Column 2: Number of Monte Carlo integration pixels.\n"
	  "# Column 3: Fraction of brightness in Monte Carlo "
          "integrated pixels.\n"
	  "# Column 4: An individual image was created.\n",
	  ctime(&p->rawtime), p->zeropoint);


  arraytotxt(p->log, p->cs0, LOGNUMCOLS, comments, int_cols,
	     accu_cols, space, prec, 'f', LOGFILENAME);
}





void
write(struct mkprofparams *p)
{
  void *array;
  char *jobname;
  double sum, *log;
  struct timeval t1;
  long os=p->oversample;
  int replace=p->replace;
  size_t complete=0, cs0=p->cs0;
  struct builtqueue *ibq=NULL, *tbq;
  int verb=p->cp.verb, bitpix=p->bitpix;
  float *out, *to, *from, *colend, *rowend;
  size_t i, j, iw, jw, ii, jj, w=p->naxes[0], ow;


  /* Allocate the output array. */
  if(p->up.backname)
    out=p->out;
  else
    {
      errno=0;
      out=calloc(p->naxes[0]*p->naxes[1], sizeof *out);
      if(out==NULL)
        error(EXIT_FAILURE, 0, "%lu bytes for output image",
              p->naxes[0]*p->naxes[1]*sizeof *out);
    }


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
      sum=0.0f;

      /* Write the array pointed to by ibq into the output image. Note
	 that the FITS and C arrays have opposite axis orders and FITS
	 counting starts from 1, not zero. Also fpixel is the first
	 (inclusive) pixel and so is lpixel (it is inclusive). */
      if(ibq->overlaps && p->nomerged==0)
	{
	  /* Set the starting and ending points in the complete image. */
	  i  = os * (ibq->fpixel_i[1]-1);
	  j  = os * (ibq->fpixel_i[0]-1);

	  /* Set the starting and ending points in the overlapping
	     image. Note that oversampling has already been taken
	     into account in ibq->width. */
	  ow = ibq->imgwidth;
	  ii = os * (ibq->fpixel_o[1]-1);
	  jj = os * (ibq->fpixel_o[0]-1);

	  /* Find the width of the overlapping region: */
	  iw = os*(ibq->lpixel_i[1]-ibq->fpixel_i[1]+1);
	  jw = os*(ibq->lpixel_i[0]-ibq->fpixel_i[0]+1);

	  /* Write the overlap to the actual image. Instead of writing
	     two for loops and summing all the row and column indexs
	     for every pixel and each image, we use pointer arithmetic
	     which is much more efficient. Just think of one pointer
	     that is advancing over the final image (*to) and one that
	     is advancing over the overlap image (*from). Since we
	     know the images overlap, iw and jw are both smaller than
	     the two image number of columns and number of rows, so
	     w-jw and ow-jw will always be positive. */
	  to=out+i*w+j;
	  from=ibq->img+ii*ow+jj;
	  rowend=to+iw*w;
          do
            {
              colend=to+jw;
              do
                {
                  if(!isnan(*from))
                    {
                      *from = p->setconsttonan ? NAN : *from;
                      sum+=*from;
                      *to = replace ? *from : *to+*from;
                    }
                  ++from;
                }
              while(++to<colend);
              to+=w-jw; from+=ow-jw;	     /* Go to next row. */
            }
          while(to<rowend);
	}

      /* Fill the log array. */
      log=&p->log[ibq->id*LOGNUMCOLS];
      log[0] = ibq->id;
      log[1] = sum>0.0f ? -2.5f*log10(sum)+p->zeropoint : NAN;
      log[2] = ibq->numaccu;
      log[3] = ibq->accufrac;
      log[4] = ibq->indivcreated;

      /* Report if in verbose mode. */
      ++complete;
      if(verb && p->nomerged==0)
	{
	  errno=0;
	  jobname=malloc(100*sizeof *jobname);
	  if(jobname==NULL)
	    error(EXIT_FAILURE, errno, "jobname in mkprof.c");
	  sprintf(jobname, "Row %lu complete, %lu left to go.",
		  ibq->id, cs0-complete);
	  reporttiming(NULL, jobname, 2);
	  free(jobname);
	}

      /* Free the array and the queue element and change it to the
	 next one and increment complete. Note that there is no
	 problem to free a NULL pointer (when the built array didn't
	 overlap). */
      free(ibq->img);
      tbq=ibq->next;
      free(ibq);
      ibq=tbq;
    }

  /* Write the final array to the final FITS image. */
  if(p->nomerged==0)
    {
      if(verb) gettimeofday(&t1, NULL);
      if(p->up.backname)
        {
          if(bitpix==FLOAT_IMG) array=out;
          else changetype(p->out, FLOAT_IMG, p->naxes[1]*p->naxes[0],
                          p->anyblank, &array, bitpix);
          arraytofitsimg(p->mergedimgname, "MockImg on back", bitpix,
                         array, p->naxes[1], p->naxes[0], p->anyblank,
                         p->wcs, NULL, SPACK_STRING);
          if(bitpix!=FLOAT_IMG) free(array);
        }
      else
        atofcorrectwcs(p->mergedimgname, "MockImg", FLOAT_IMG, out,
                       p->naxes[1], p->naxes[0], p->wcsheader,
                       p->wcsnkeyrec, NULL, SPACK_STRING);
      if(verb)
	{
	  errno=0;
	  jobname=malloc((strlen(p->mergedimgname)+100)*sizeof *jobname);
	  if(jobname==NULL)
	    error(EXIT_FAILURE, errno, "final report in mkprof.c");
	  sprintf(jobname, "%s created.", p->mergedimgname);
	  reporttiming(&t1, jobname, 1);
	  free(jobname);
	}
    }

  free(out);
}




















/**************************************************************/
/************           Outside function          *************/
/**************************************************************/
void
mkprof(struct mkprofparams *p)
{
  int err;
  pthread_t t;		 /* Thread id not used, all are saved here. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct mkonthread *mkp;
  size_t i, *indexs, thrdcols;
  size_t nt=p->cp.numthreads, nb;
  long onaxes[2], os=p->oversample;

  /* Get the WCS header strings ready to put into the FITS
     image(s). */
  if(p->up.backname==NULL)
    preparewcs(p);

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

  /* onaxes are sides of the image without over-sampling. */
  onaxes[0] = (p->naxes[0]-2*p->shift[0])/os + 2*p->shift[0]/os;
  onaxes[1] = (p->naxes[1]-2*p->shift[1])/os + 2*p->shift[1]/os;

  /* Build the profiles: */
  if(nt==1)
    {
      mkp[0].p=p;
      mkp[0].onaxes=onaxes;
      mkp[0].indexs=indexs;
      mkp[0].rng=gsl_rng_clone(p->rng);
      build(&mkp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this main thread will
	 also have to be kept behind the barrier, so we need nt+1
	 barrier stops. */
      if(p->cs0<nt) nb=p->cs0+1;
      else nb=nt+1;
      attrbarrierinit(&attr, &b, nb);

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
	    mkp[i].onaxes=onaxes;
            mkp[i].rng=gsl_rng_clone(p->rng);
	    mkp[i].indexs=&indexs[i*thrdcols];
	    err=pthread_create(&t, &attr, build, &mkp[i]);
	    if(err)
	      error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }
    }

  /* Write the created arrays into the image. */
  write(p);
  writelog(p);

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

  /* Free the allocated spaces. */
  free(mkp);
  free(indexs);
  free(p->wcsheader);
}
