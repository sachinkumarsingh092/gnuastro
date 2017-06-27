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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/box.h>
#include <gnuastro/git.h>
#include <gnuastro/fits.h>
#include <gnuastro/threads.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "mkprof.h"             /* Needs main.h astrthreads.h */
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
    error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for `tbq'",
          __func__, sizeof *tbq);

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

  double crpix[2];
  gal_data_t *data;
  long os=p->oversample;
  struct builtqueue *ibq=mkp->ibq;
  char *filename, *jobname, *outdir=p->outdir;

  /* Note that `width' is in FITS format, not C. */
  size_t dsize[2]={mkp->width[1], mkp->width[0]};

  /* Write the name and remove a similarly named file. */
  asprintf(&filename, "%s%zu_%s", outdir, ibq->id, p->basename);
  gal_checkset_check_remove_file(filename, 0, p->cp.dontdelete);

  /* Put the array into a data structure */
  data=gal_data_alloc(ibq->img, GAL_TYPE_FLOAT32, 2, dsize, NULL, 0,
                      p->cp.minmapsize, "MockImage", "Brightness", NULL);

  /* Write the array to file (a separately built PSF doesn't need WCS
     coordinates). */
  if(ibq->ispsf && p->psfinimg==0)
    gal_fits_img_write(data, filename, NULL, PROGRAM_STRING);
  else
    {
      /* Save the correct CRPIX values: */
      crpix[0] = p->crpix[0] - os*(mkp->fpixel_i[0]-1);
      crpix[1] = p->crpix[1] - os*(mkp->fpixel_i[1]-1);

      /* Write the image. */
      gal_fits_img_write_corr_wcs_str(data, filename, p->wcsheader,
                                      p->wcsnkeyrec, crpix, NULL,
                                      PROGRAM_STRING);
    }
  ibq->indivcreated=1;

  /* Report if in verbose mode. */
  if(!p->cp.quiet)
    {
      asprintf(&jobname, "%s created.", filename);
      gal_timing_report(NULL, jobname, 2);
      free(jobname);
    }
  free(filename);
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
mkprof_build(void *inparam)
{
  struct mkonthread *mkp=(struct mkonthread *)inparam;
  struct mkprofparams *p=mkp->p;

  size_t i, id;
  int lockresult;
  double center[2];
  long lpixel_o[2];
  pthread_mutex_t *qlock=&p->qlock;
  struct builtqueue *ibq, *fbq=NULL;
  pthread_cond_t *qready=&p->qready;

  /* Make each profile that was specified for this thread. */
  for(i=0; mkp->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* Create a new builtqueue element with all the information. fbq
         will be used when we want to add ibq to p->bq. It is defined
         so we don't have to waste time traversing the ibq. Its
         characteristic compared to the other elements of ibq is that
         fbq->next==NULL. So to add ibq to p->bq, we just have to set
         fbq->next=p->bq and then set p->bq to ibq.*/
      builtqueue_addempty(&mkp->ibq);
      ibq=mkp->ibq;
      id=ibq->id=mkp->indexs[i];
      if(fbq==NULL) fbq=ibq;


      /* Write the necessary parameters for this profile into mkp.*/
      oneprof_set_prof_params(mkp);


      /* Find the bounding box size (NOT oversampled). */
      if( p->f[id] == PROFILE_POINT )
        mkp->width[0]=mkp->width[1]=1;
      else
        gal_box_ellipse_in_box(mkp->truncr, mkp->q*mkp->truncr,
                               p->p[id]*DEGREESTORADIANS, mkp->width);


      /* Get the overlapping pixels using the starting points (NOT
         oversampled). */
      center[0]=p->x[id];
      center[1]=p->y[id];
      gal_box_border_from_center(center, 2, mkp->width,
                                 ibq->fpixel_i, ibq->lpixel_i);
      mkp->fpixel_i[0]=ibq->fpixel_i[0];
      mkp->fpixel_i[1]=ibq->fpixel_i[1];
      ibq->overlaps = gal_box_overlap(mkp->onaxes, ibq->fpixel_i,
                                      ibq->lpixel_i, ibq->fpixel_o,
                                      lpixel_o, 2);


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
            gsl_rng_set(mkp->rng, gal_timing_time_based_rng_seed());

          /* Make the profile */
          oneprofile_make(mkp);
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
                 either `mkprof_write` is waiting behind a condition
                 variable for you to fill it up or not (either it hasn't
                 got to setting the condition variable yet (this function
                 locked the mutex before `mkprof_write`) or it just got the
                 list to be made and is busy writing the arrays in the
                 output). In either case, pthread_cond_signal will work. */
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
          else if (mkp->indexs[i+1]==GAL_BLANK_SIZE_T)
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
mkprof_write(struct mkprofparams *p)
{
  double sum;
  char *jobname;
  struct timeval t1;
  long os=p->oversample;
  int replace=p->replace;
  gal_data_t *out=p->out, *log;
  struct builtqueue *ibq=NULL, *tbq;
  float *to, *from, *colend, *rowend;
  size_t complete=0, num=p->num, clog;
  size_t i, j, iw, jw, ii, jj, w=p->naxes[0], ow;


  /* Write each image into the output array. */
  while(complete<p->num)
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
      if(ibq->overlaps && out->array)
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
          to = (float *)(out->array) + i*w+j;
          from=ibq->img+ii*ow+jj;
          rowend=to+iw*w;
          do
            {
              /* Go over all the pixels in this row and write this profile
                 into the final output array. Just note that when
                 replacing, we don't want to replace those pixels that have
                 a zero value, since no profile existed there. */
              colend=to+jw;
              do
                {
                  *to  = ( replace
                           ? ( *from==0.0f ? *to : *from )
                           :  *to + *from );
                  sum += *from;
                  ++from;
                }
              while(++to<colend);

              /* Go to the next row. */
              to   += w-jw;
              from += ow-jw;
            }
          while(to<rowend);
        }

      /* Fill the log array. */
      if(p->cp.log)
        {
          clog=0;
          for(log=p->log; log!=NULL; log=log->next)
            switch(++clog)
              {
              case 5:
                ((unsigned char *)(log->array))[ibq->id] = ibq->indivcreated;
                break;
              case 4:
                ((float *)(log->array))[ibq->id] = ibq->accufrac;
                break;
              case 3:
                ((unsigned long *)(log->array))[ibq->id]=ibq->numaccu;
                break;
              case 2:
                ((float *)(log->array))[ibq->id] =
                  sum>0.0f ? -2.5f*log10(sum)+p->zeropoint : NAN;
                break;
              case 1:
                ((unsigned long *)(log->array))[ibq->id]=ibq->id+1;
                break;
              }
        }

      /* Report if in verbose mode. */
      ++complete;
      if(!p->cp.quiet)
        {
          asprintf(&jobname, "row %zu complete, %zu left to go",
                   ibq->id, num-complete);
          gal_timing_report(NULL, jobname, 2);
          free(jobname);
        }

      /* Free the array and the queue element and change it to the next one
         and increment complete. Note that there is no problem to free a
         NULL pointer (when the built array didn't overlap). */
      free(ibq->img);
      tbq=ibq->next;
      free(ibq);
      ibq=tbq;
    }

  /* Write the final array to the output FITS image if a merged image is to
     be created. */
  if(out->array)
    {
      /* Get the current time for verbose output. */
      if(!p->cp.quiet) gettimeofday(&t1, NULL);

      /* Write the final image into a FITS file with the requested type. */
      gal_fits_img_write_to_type(out, p->mergedimgname, NULL,
                                 PROGRAM_STRING, p->cp.type);

      /* Clean up */
      gal_data_free(out);

      /* In verbose mode, print the information. */
      if(!p->cp.quiet)
        {
          asprintf(&jobname, "%s created.", p->mergedimgname);
          gal_timing_report(&t1, jobname, 1);
          free(jobname);
        }
    }

  /* Even with no merged image, there might still be pointers in `out' that
     need to be freed. */
  else
    gal_data_free(p->out);
}




















/**************************************************************/
/************           Outside function          *************/
/**************************************************************/
void
mkprof(struct mkprofparams *p)
{
  int err;
  char *tmp;
  pthread_t t;            /* Thread id not used, all are saved here. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct mkonthread *mkp;
  size_t i, *indexs, thrdcols;
  gal_list_str_t *comments=NULL;
  size_t nt=p->cp.numthreads, nb;
  long onaxes[2], os=p->oversample;


  /* Allocate the arrays to keep the thread and parameters for each
     thread. Note that we only want nt-1 threads to do the
     building. */
  errno=0;
  mkp=malloc(nt*sizeof *mkp);
  if(mkp==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `mkp'",
          __func__, (nt-1)*sizeof *mkp);

  /* Distribute the different profiles for different threads. Note
     that one thread is left out for writing, while nt-1 are left
     for building. */
  gal_threads_dist_in_threads(p->num, nt, &indexs, &thrdcols);

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
      mkprof_build(&mkp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this main thread will
         also have to be kept behind the barrier, so we need nt+1
         barrier stops. */
      if(p->num<nt) nb=p->num+1;
      else nb=nt+1;
      gal_threads_attr_barrier_init(&attr, &b, nb);

      /* Initialize the condition variable and mutex. */
      err=pthread_mutex_init(&p->qlock, NULL);
      if(err) error(EXIT_FAILURE, 0, "%s: mutex not initialized", __func__);
      err=pthread_cond_init(&p->qready, NULL);
      if(err) error(EXIT_FAILURE, 0, "%s: condition variable not initialized",
                    __func__);

      /* Spin off the threads: */
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=GAL_BLANK_SIZE_T)
          {
            mkp[i].p=p;
            mkp[i].b=&b;
            mkp[i].ibq=NULL;
            mkp[i].onaxes=onaxes;
            mkp[i].rng=gsl_rng_clone(p->rng);
            mkp[i].indexs=&indexs[i*thrdcols];
            err=pthread_create(&t, &attr, mkprof_build, &mkp[i]);
            if(err)
              error(EXIT_FAILURE, 0, "%s: can't create thread %zu",
                    __func__, i);
          }
    }

  /* Write the created arrays into the image. */
  mkprof_write(p);

  /* Write the log file. */
  if(p->cp.log)
    {
      asprintf(&tmp, "Zeropoint: %g", p->zeropoint);
      gal_list_str_add(&comments, tmp, 0);
      gal_table_write_log(p->log, PROGRAM_STRING, &p->rawtime, comments,
                          LOGFILENAME, p->cp.dontdelete, p->cp.quiet);
      gal_list_str_free(comments, 1);
    }

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
}
