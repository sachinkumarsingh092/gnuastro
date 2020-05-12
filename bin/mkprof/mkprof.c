/********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "ui.h"
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
    error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for 'tbq'",
          __func__, sizeof *tbq);

  /* Initialize the values (same order as in structure definition). */
  tbq->id           = GAL_BLANK_SIZE_T;
  tbq->ispsf        = 0;
  tbq->overlaps     = 0;
  tbq->image        = NULL;
  tbq->overlap_i    = NULL;
  tbq->overlap_m    = NULL;
  tbq->func         = PROFILE_MAXIMUM_CODE;
  tbq->indivcreated = 0;
  tbq->numaccu      = 0;
  tbq->accufrac     = 0.0f;

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

  double *crpix;
  long os=p->oversample;
  gal_fits_list_key_t *keys=NULL;
  struct builtqueue *ibq=mkp->ibq;
  size_t i, ndim=p->ndim, id=mkp->ibq->id;
  char *filename, *jobname, *outdir=p->outdir;


  /* Write the name and remove a similarly named file when the '--kernel'
     option wasn't called. If '--kernel' is called, then we should just use
     the final merged filename. */
  if(p->kernel)
    filename=p->mergedimgname;
  else
    {
      if( asprintf(&filename, "%s%zu_%s", outdir, ibq->id, p->basename)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_checkset_writable_remove(filename, 0, p->cp.dontdelete);
    }


  /* Write the array to the file (a separately built PSF doesn't need WCS
     coordinates). */
  if(ibq->ispsf && p->psfinimg==0)
    gal_fits_img_write(ibq->image, filename, NULL, PROGRAM_NAME);
  else
    {
      /* Allocate space for the corrected crpix and fill it in. Both
         'crpix' and 'fpixel_i' are in FITS order. */
      crpix=gal_pointer_allocate(GAL_TYPE_FLOAT64, ndim, 0, __func__,
                                 "crpix");
      for(i=0;i<ndim;++i)
        crpix[i] = ((double *)(p->crpix->array))[i] - os*(mkp->fpixel_i[i]-1);

      /* Write the image. */
      gal_fits_img_write_corr_wcs_str(ibq->image, filename, p->wcsheader,
                                      p->wcsnkeyrec, crpix, NULL,
                                      PROGRAM_NAME);
    }
  ibq->indivcreated=1;


  /* Write profile settings into the FITS file. */
  gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "PROFILE", 0,
                        ui_profile_name_write(mkp->func), 0,
                        "Radial function", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT64, "XCENTER", 0,
                        &p->x[id], 0, "Center of profile in catalog "
                        "(FITS axis 1)", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT64, "YCENTER", 0,
                        &p->y[id], 0, "Center of profile in catalog "
                        "(FITS axis 2)", 0, NULL);
  if(ndim==3)
    gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT64, "ZCENTER", 0,
                          &p->z[id], 0, "Center of profile in catalog "
                          "(FITS axis 3)", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "RADIUS", 0,
                        &p->r[id], 0, "Radial parameter in catalog",
                        0, NULL);
  if( mkp->func==PROFILE_SERSIC || mkp->func==PROFILE_MOFFAT )
    gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "PINDEX", 0,
                          &p->r[id], 0, "Index (Sersic or Moffat) of profile"
                          "in catalog", 0, NULL);
  if(ndim==2)
    {
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "PA_DEG", 0,
                            &p->p1[id], 0, "Position angle of profile in "
                            "catalog", 0, "deg");
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "AXISRATIO", 0,
                            &p->q1[id], 0, "Axis ratio of profile in catalog",
                            0, NULL);
    }
  else
    {
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "PA1_DEG", 0,
                            &p->p1[id], 0, "First X-Z-X Euler angle in 3D", 0,
                            "deg");
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "PA2_DEG", 0,
                            &p->p2[id], 0, "Second X-Z-X Euler angle in 3D", 0,
                            "deg");
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "PA3_DEG", 0,
                            &p->p3[id], 0, "Third X-Z-X Euler angle in 3D", 0,
                            "deg");
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "AXISRATIO1", 0,
                            &p->q1[id], 0, "Axis ratio along second dim",
                            0, NULL);
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "AXISRATIO2", 0,
                            &p->q2[id], 0, "Axis ratio along third dim",
                            0, NULL);
    }
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "MAGNITUDE", 0,
                        &p->m[id], 0, "Magnitude of profile in catalog",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "TRUNCATION", 0,
                        &p->t[id], 0, "Truncation of profile in catalog",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "RNGNAME", 0,
                        (void *)(p->rng_name), 0,
                        "Name of random number generator", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_ULONG, "RNGSEED", 0,
                        &mkp->rng_seed, 0, "Seed of random number generator",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_SIZE_T, "NUMRANDOM", 0,
                        &p->numrandom, 0,
                        "Number of random points in central pixels", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "TOLERANCE", 0,
                        &p->tolerance, 0,
                        "Tolerance level to stop random integration",
                        0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_STRING, "MODE", 0,
                        p->mode==MKPROF_MODE_IMG?"img":"wcs", 0,
                        "Coordinates in image or WCS units", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "OVERSAMPLE", 0,
                        &p->oversample, 0, "Oversampling factor", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "TUNITINP", 0,
                        &p->tunitinp, 0, "Truncation is in units of pixels, "
                        "not radius", 0, NULL);
  if( !isnan(p->zeropoint) )
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "ZEROPOINT", 0,
                            &p->zeropoint, 0, "Zeropoint magnitude", 0, NULL);
  if( mkp->func==PROFILE_CIRCUMFERENCE )
      gal_fits_key_list_add(&keys, GAL_TYPE_FLOAT32, "CIRCUMWIDTH", 0,
                            &p->circumwidth, 0, "Width of circumference "
                            "(inward) profiles", 0, NULL);
  if( mkp->func==PROFILE_FLAT || mkp->func==PROFILE_CIRCUMFERENCE )
      gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "MFORFLATPIX", 0,
                            &p->mforflatpix, 0, "Magnitude is flat pixel "
                            "value", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "MCOLISBRIGHTNESS", 0,
                        &p->mcolisbrightness, 0, "Catalog's magnitude is "
                        "actually brightness", 0, NULL);
  gal_fits_key_list_add(&keys, GAL_TYPE_UINT8, "MAGATPEAK", 0,
                        &p->magatpeak, 0, "Magnitude is for peak pixel, "
                        "not full profile", 0, NULL);

  gal_fits_key_list_reverse(&keys);
  gal_fits_key_write_config(&keys, "Profile configuration", "PROFILE-CONFIG",
                            filename, "0");


  /* Report if in verbose mode. */
  if(!p->cp.quiet)
    {
      if( asprintf(&jobname, "%s created.", filename)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(NULL, jobname, 2);
      free(jobname);
    }


  /* Clean up. */
  if(p->kernel==NULL) free(filename);
}



















/**************************************************************/
/************            The builders             *************/
/**************************************************************/
/* High-level function to built a single profile and prepare it for the
   next steps. */
static void
mkprof_build_single(struct mkonthread *mkp, long *fpixel_i, long *lpixel_i,
                    long *fpixel_o)
{
  struct mkprofparams *p = mkp->p;
  struct builtqueue *ibq = mkp->ibq;

  void *ptr;
  int needs_crop=0;
  size_t i, ind, fits_i, ndim=p->ndim;
  size_t start_indiv[3], start_mrg[3], dsize[3], os=p->oversample;

  /* Make a copy of the main random number generator to use for this
     profile (in this thread). */
  gsl_rng_memcpy(mkp->rng, p->rng);

  /* Set the seed of the random number generator if the
     environment is not to be used. */
  if(mkp->p->envseed)
    mkp->rng_seed=mkp->p->rng_seed;
  else
    {
      mkp->rng_seed=gal_timing_time_based_rng_seed();
      gsl_rng_set(mkp->rng, mkp->rng_seed);
    }

  /* Make the profile */
  oneprofile_make(mkp);

  /* Build an individual image if necessary. */
  if( p->individual || (ibq->ispsf && p->psfinimg==0))
    {
      saveindividual(mkp);
      if(ibq->ispsf && p->psfinimg==0)
        ibq->overlaps=0;
    }

  /* If we want a merged image, then a tile needs to be defined over the
     individual profile array and the output merged array to define the
     overlapping region. */
  if(p->out)
    {
      /* Note that 'fpixel_i' and 'lpixel_o' were in the un-oversampled
         image, they are also in the FITS coordinates. */
      for(i=0;i<ndim;++i)
        {
          /* Set the start and width of the overlap. */
          fits_i = ndim-i-1;
          start_indiv[i] = os * (fpixel_o[fits_i] - 1);
          start_mrg[i]   = os * (fpixel_i[fits_i] - 1);
          dsize[i]       = os * (lpixel_i[fits_i] - fpixel_i[fits_i] + 1);

          /* Check if we need to crop the individual image or not. */
          if(dsize[i] != ibq->image->dsize[i]) needs_crop=1;
        }


      /* Define the individual overlap tile. */
      if(needs_crop)
        {
          /* If a crop is needed, set the starting pointer. */
          ind=gal_dimension_coord_to_index(ndim, ibq->image->dsize,
                                           start_indiv);
          ptr=gal_pointer_increment(ibq->image->array, ind, ibq->image->type);
        }
      else ptr=ibq->image->array;
      ibq->overlap_i=gal_data_alloc(ptr, ibq->image->type, ndim, dsize, NULL,
                                    0, -1, 1, NULL, NULL, NULL);
      ibq->overlap_i->block=ibq->image;


      /* Define the merged overlap tile. */
      ind=gal_dimension_coord_to_index(ndim, p->out->dsize, start_mrg);
      ptr=gal_pointer_increment(p->out->array, ind, p->out->type);
      ibq->overlap_m=gal_data_alloc(ptr, p->out->type, ndim, dsize, NULL,
                                    0, -1, 1, NULL, NULL, NULL);
      ibq->overlap_m->block=p->out;
    }
}





/* The profile has been built, now add it to the queue of profiles that
   must be written into the final merged image. */
static void
mkprof_add_built_to_write_queue(struct mkonthread *mkp,
                                struct builtqueue *ibq,
                                struct builtqueue **fbq, size_t counter)
{
  struct mkprofparams *p = mkp->p;

  int lockresult;
  pthread_mutex_t *qlock=&p->qlock;
  pthread_cond_t *qready=&p->qready;

  /* Try locking the mutex so no thread can change the value of p->bq. If
     you can lock it, then put the internal builtqueue on top of the
     general builtqueue. If you can't, continue adding to the internal
     builtqueue (make the next profiles) until you find a chance to lock
     the mutex. */
  lockresult=pthread_mutex_trylock(qlock);
  if(lockresult==0)     /* Mutex was successfully locked. */
    {
      /* Add this internal queue to system queue. */
      (*fbq)->next=p->bq;
      p->bq=ibq;

      /* If the list was empty when you locked the mutex, then either
         'mkprof_write' is waiting behind a condition variable for you to
         fill it up or not (either it hasn't got to setting the condition
         variable yet (this function locked the mutex before
         'mkprof_write') or it just got the list to be made and is busy
         writing the arrays in the output). In either case,
         pthread_cond_signal will work. */
      if((*fbq)->next==NULL)
        pthread_cond_signal(qready);
      pthread_mutex_unlock(qlock);

      /* Finally set both the internal queue and the first internal queue
         element to NULL.*/
      (*fbq)=NULL;
      mkp->ibq=NULL;
    }

  /* The mutex couldn't be locked and there are no more objects for this
     thread to build (giving a chance for this thread to add up its built
     profiles). So we have to lock the mutex to pass on this built
     structure to the builtqueue. */
  else if (mkp->indexs[counter+1]==GAL_BLANK_SIZE_T)
    {
      pthread_mutex_lock(qlock);
      (*fbq)->next=p->bq;
      p->bq=ibq;
      pthread_cond_signal(qready);
      pthread_mutex_unlock(qlock);
    }
}





/* Build the profiles that are indexed in the indexs array of the
   mkonthread structure that was assigned to it.

   See the explanation above overlap (/lib/box.c) for a complete
   explanation of fpixel_i, lpixel_i, fpixel_o and lpixel_o.

   =========================================================
   About the Central x and y of each profile:

   The user has asked for the profile to be built on the coordinates
   (real numbers) of 'x' and 'y' in an output image in the FITS
   format. We are building the full image for each galaxy separately
   in an array with an odd number of sides which maybe oversampled.

   In the FITS format, the pixel centers have an integer value. So for
   example in 1D, a pixel whose center value is 10.00 covers the area
   of: [9.5,10.5). We want the fractional part of 'x' (don't forget,
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
static void *
mkprof_build(void *inparam)
{
  struct mkonthread *mkp=(struct mkonthread *)inparam;
  struct mkprofparams *p=mkp->p;
  size_t i, id, ndim=p->ndim;
  struct builtqueue *ibq, *fbq=NULL;
  double center[3], semiaxes[3], euler_deg[3];
  long fpixel_i[3], lpixel_i[3], fpixel_o[3], lpixel_o[3];


  /* Make each profile that was specified for this thread. */
  for(i=0; mkp->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* Create a new builtqueue element with all the information. 'fbq'
         will be used when we want to add 'ibq' to 'p->bq'. It is defined
         so we don't have to waste time traversing the 'ibq'. Its
         characteristic compared to the other elements of 'ibq' is that
         'fbq->next==NULL'. So to add ibq to p->bq, we just have to set
         'fbq->next=p->bq' and then set 'p->bq' to 'ibq'.*/
      builtqueue_addempty(&mkp->ibq);
      ibq=mkp->ibq;
      id=ibq->id=mkp->indexs[i];
      if(fbq==NULL) fbq=ibq;


      /* Write the necessary parameters for this profile into 'mkp'.*/
      oneprofile_set_prof_params(mkp);


      /* Find the bounding box size (NOT oversampled). */
      if( p->f[id] == PROFILE_POINT )
        mkp->width[0]=mkp->width[1]=1;
      else
        switch(ndim)
          {
          case 2:
            gal_box_bound_ellipse(mkp->truncr, mkp->q[0]*mkp->truncr,
                                  p->p1[id], mkp->width);
            break;

          case 3:
            euler_deg[0] = p->p1[id];
            euler_deg[1] = p->p2[id];
            euler_deg[2] = p->p3[id];
            semiaxes[0]  = mkp->truncr;
            semiaxes[1]  = mkp->truncr * mkp->q[0];
            semiaxes[2]  = mkp->truncr * mkp->q[1];
            gal_box_bound_ellipsoid(semiaxes, euler_deg, mkp->width);
            break;

          default:
            error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to "
                  "address the issue. %zu is not recognized for 'ndim'",
                  __func__, PACKAGE_BUGREPORT, ndim);
          }


      /* Get the overlapping pixels using the starting points (NOT
         oversampled). */
      if(p->out)
        {
          center[0]=p->x[id];
          center[1]=p->y[id];
          if(ndim==3) center[2]=p->z[id];
          gal_box_border_from_center(center, ndim, mkp->width, fpixel_i,
                                     lpixel_i);
          memcpy(mkp->fpixel_i, fpixel_i, ndim*sizeof *fpixel_i);
          ibq->overlaps = gal_box_overlap(mkp->onaxes, fpixel_i, lpixel_i,
                                          fpixel_o, lpixel_o, ndim);
        }


      /* Build the profile if necessary. */
      if(ibq->overlaps || p->individual || (ibq->ispsf && p->psfinimg==0))
        mkprof_build_single(mkp, fpixel_i, lpixel_i, fpixel_o);


      /* Add this profile to the list of profiles that must be written onto
         the final merged image with another thread. */
      if(p->cp.numthreads>1)
        mkprof_add_built_to_write_queue(mkp, ibq, &fbq, i);
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
static void
mkprof_write(struct mkprofparams *p)
{
  double sum;
  char *jobname;
  struct timeval t1;
  gal_data_t *out=p->out, *log;
  struct builtqueue *ibq=NULL, *tbq;
  size_t complete=0, num=p->num, clog;

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


      /* During the build process, we also defined the overlap tiles of
         both the individual array and the final merged array, here we will
         use those to put the required profile pixels into the final
         array. */
      if(ibq->overlaps && out)
        GAL_TILE_PO_OISET(float,float,ibq->overlap_i,ibq->overlap_m,1,0, {
            *o  = p->replace ? ( *i==0.0f ? *o : *i ) :  (*i + *o);
            sum += *i;
          });


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
      if(!p->cp.quiet && p->num>1)
        {
          if( asprintf(&jobname, "row %zu complete, %zu left to go",
                       ibq->id+1, num-complete)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(NULL, jobname, 2);
          free(jobname);
        }


      /* Free the array and the queue element and change it to the next one
         and increment complete. Note that there is no problem to free a
         NULL pointer (when the built array didn't overlap). */
      gal_data_free(ibq->overlap_i);
      gal_data_free(ibq->overlap_m);
      gal_data_free(ibq->image);
      tbq=ibq->next;
      free(ibq);
      ibq=tbq;
    }


  /* Write the final array to the output FITS image if a merged image is to
     be created. */
  if(out)
    {
      /* Get the current time for verbose output. */
      if(!p->cp.quiet) gettimeofday(&t1, NULL);

      /* Write the final image into a FITS file with the requested
         type. Until now, we were using 'p->wcs' for the WCS, but from now
         on, will put it in 'out' to also free it while freeing 'out'. */
      out->wcs=p->wcs;
      gal_fits_img_write_to_type(out, p->mergedimgname, NULL,
                                 PROGRAM_NAME, p->cp.type);
      p->wcs=NULL;

      /* Clean up */
      gal_data_free(out);

      /* Write the configuration keywords. */
      gal_fits_key_write_filename("input", p->catname, &p->cp.okeys, 1);
      gal_fits_key_write_config(&p->cp.okeys, "MakeProfiles configuration",
                                "MKPROF-CONFIG", p->mergedimgname, "0");

      /* In verbose mode, print the information. */
      if(!p->cp.quiet)
        {
          if( asprintf(&jobname, "%s created.", p->mergedimgname)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(&t1, jobname, 1);
          free(jobname);
        }
    }
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
  gal_list_str_t *comments=NULL;
  size_t i, fi, *indexs, thrdcols;
  long *onaxes=NULL, os=p->oversample;
  size_t nb, ndim=p->ndim, nt=p->cp.numthreads;

  /* Allocate the arrays to keep the thread and parameters for each
     thread. Note that we only want nt-1 threads to do the
     building. */
  errno=0;
  mkp=malloc(nt*sizeof *mkp);
  if(mkp==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'mkp'",
          __func__, (nt-1)*sizeof *mkp);


  /* Distribute the different profiles for different threads. Note
     that one thread is left out for writing, while nt-1 are left
     for building. */
  gal_threads_dist_in_threads(p->num, nt, &indexs, &thrdcols);


  /* 'onaxes' are size of the merged output image without over-sampling or
     shifting in FITS order. When no output merged image is needed, we can
     ignore it. */
  if(p->out)
    {
      onaxes=gal_pointer_allocate(GAL_TYPE_LONG, ndim, 0, __func__, "onaxes");
      for(fi=0; fi < ndim; ++fi)
        {
          i=ndim-fi-1;
          onaxes[fi] = ( ( p->dsize[i] - 2 * p->shift[i] ) / os
                         + 2 * p->shift[i]/os );
        }
    }


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
      if( asprintf(&tmp, "Zeropoint: %g", p->zeropoint)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_str_add(&comments, tmp, 0);
      gal_checkset_writable_remove(LOGFILENAME, 0, p->cp.dontdelete);
      gal_table_write_log(p->log, PROGRAM_STRING, &p->rawtime, comments,
                          LOGFILENAME, p->cp.quiet);
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


  /* Clean up. */
  free(mkp);
  free(indexs);
  if(onaxes) free(onaxes);
}
