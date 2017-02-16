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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <sys/time.h>            /* generate random seed */
#include <gsl/gsl_rng.h>         /* used in setrandoms   */
#include <gsl/gsl_randist.h>     /* To make noise.       */
#include <gsl/gsl_integration.h> /* gsl_integration_qng  */

#include <gnuastro/fits.h>
#include <gnuastro/array.h>
#include <gnuastro/statistics.h>
#include <gnuastro/linkedlist.h>

#include <timing.h>
#include <neighbors.h>

#include "main.h"

#include "mkprof.h"              /* Needs main.h and astrthreads.h */
#include "profiles.h"
#include "oneprofile.h"






/****************************************************************
 **************        Elliptical radius       ******************
 ****************************************************************/
/* Convert cartesian coordinates to the rotated elliptical radius. */
void
r_el(struct mkonthread *mkp)
{
  double c=mkp->c, s=mkp->s, q=mkp->q, x=mkp->x, y=mkp->y;
  mkp->r = sqrt( (x*c+y*s)*(x*c+y*s) + ((y*c-x*s)*(y*c-x*s)/q/q) );
}





/* Calculate the cercular distance of a pixel to the profile center. */
float
r_circle(size_t p, struct mkonthread *mkp)
{
  double x, y;

  x = p/mkp->width[0];   /* Note that width[0] is the First FITS */
  y = p%mkp->width[0];   /* axis, not first C axis.              */

  return sqrt( (x-mkp->xc)*(x-mkp->xc) + (y-mkp->yc)*(y-mkp->yc) );
}



















/****************************************************************
 **************          Random points         ******************
 ****************************************************************/
/* Fill pixel with random values */
float
randompoints(struct mkonthread *mkp)
{
  double xrange, yrange, sum=0.0f;
  size_t i, numrandom=mkp->p->numrandom;

  /* Set the range of the x and y: */
  xrange=mkp->xh-mkp->xl;
  yrange=mkp->yh-mkp->yl;

  /* Find the sum of the profile on the random positions */
  for(i=0;i<numrandom;++i)
    {
      mkp->x = mkp->xl + gsl_rng_uniform(mkp->rng)*xrange;
      mkp->y = mkp->yl + gsl_rng_uniform(mkp->rng)*yrange;
      r_el(mkp);
      sum+=mkp->profile(mkp);
    }

  return sum/numrandom;
}




















/****************************************************************
 *****************      2D integration       ********************
 ****************************************************************/
/* This function is used in the integration of a profile. It
   assumes a fixed y and integrates over a range of x values.  */
double
twod_over_x(double x, void *params)
{
  struct mkonthread *mkp=(struct mkonthread *) params;

  mkp->x=x;
  r_el(mkp);
  return mkp->profile(mkp);
}





/* Find the 2d integration over the region. */
double
twod_over_xy(double y, void *params)
{
  gsl_function F;
  static double abserr;
  static size_t neval=0;
  double epsabs=0, epsrel=EPSREL_FOR_INTEG, result;
  struct mkonthread *mkp=(struct mkonthread *) params;

  F.function = &twod_over_x;
  F.params = params;

  mkp->y=y;
  gsl_integration_qng(&F, mkp->xl, mkp->xh, epsabs, epsrel,
                      &result, &abserr, &neval);
  return result;
}




/* 2D integration of a profile.*/
double
integ2d(struct mkonthread *mkp)
{
  gsl_function F;
  static double abserr;
  static size_t neval=0;
  double epsabs=0, epsrel=EPSREL_FOR_INTEG, result;

  F.function = &twod_over_xy;
  F.params = mkp;
  gsl_integration_qng(&F, mkp->yl, mkp->yh, epsabs,
                      epsrel, &result, &abserr, &neval);
  return result;
}




















/**************************************************************/
/************       Pixel by pixel building       *************/
/*********        Positions are in C not FITS         *********/
/**************************************************************/
static void
makepixbypix(struct mkonthread *mkp)
{
  float circ_r;
  struct gal_linkedlist_sll *Q=NULL;
  unsigned char *byt;
  float *img=mkp->ibq->img;
  int use_rand_points=1, ispeak=1;
  struct builtqueue *ibq=mkp->ibq;
  size_t is1=mkp->width[0], is0=mkp->width[1];
  size_t p, ngb[4], *ind=&p, numngb, *n, *nf, x, y;
  double tolerance=mkp->p->tolerance, pixfrac, junk;
  double (*profile)(struct mkonthread *)=mkp->profile;
  double xc=mkp->xc, yc=mkp->yc, os=mkp->p->oversample;
  double truncr=mkp->truncr, approx, hp=0.5f/mkp->p->oversample;

  /* lQ: Largest. sQ: Smallest in queue */
  struct gal_linkedlist_tosll *lQ=NULL, *sQ;

  /* Find the nearest pixel to the profile center and add it to the
     queue. */
  pixfrac = modf(mkp->xc, &junk);
  x=(long)mkp->xc + ( pixfrac<0.5f ? 0 : 1 );
  pixfrac = modf(mkp->yc, &junk);
  y=(long)mkp->yc + ( pixfrac<0.5f ? 0 : 1 );
  p=x*mkp->width[0]+y;

  /* If this is a point source, just fill that one pixel and go. */
  if(mkp->func==PROFILE_POINT)
    { img[p]=1; return; }

  /* Allocate the byt array to not repeat completed pixels. */
  errno=0;
  byt=calloc(is0*is1, sizeof *byt);
  if(byt==NULL)
    error(EXIT_FAILURE, 0, "%zu bytes for map of object in row %zu of "
          "data in %s", is0*is1*sizeof *byt, ibq->id,
          mkp->p->catname);

  /* Start the queue: */
  byt[p]=1;
  gal_linkedlist_add_to_tosll_end( &lQ, &sQ, p, r_circle(p, mkp) );

  /* If random points are necessary, then do it: */
  if(mkp->func==PROFILE_SERSIC || mkp->func==PROFILE_MOFFAT
     || mkp->func==PROFILE_GAUSSIAN)
    {
      while(sQ)
        {
          /* In case you want to see the status of the twosided ordered
             queue, increasing and decreasing side by side, uncomment this
             line. Note that there will be a lot of lines printed! */
          /*print_tossll(lQ, sQ);*/

          /* Pop the pixel from the queue and check if it is within the
             truncation radius. Note that `xc` and `p` both belong to the
             over sampled image. But all the profile parameters are in the
             non-oversampled image. So we divide the distance by os
             (p->oversample in double type) */
          gal_linkedlist_pop_from_tosll_start(&lQ, &sQ,
                                              ind, &circ_r); /* ind=&p */
          mkp->x=(p/is1-xc)/os;
          mkp->y=(p%is1-yc)/os;
          r_el(mkp);
          if(mkp->r>truncr) continue;

          /* Find the value for this pixel: */
          mkp->xl=mkp->x-hp;
          mkp->xh=mkp->x+hp;
          mkp->yl=mkp->y-hp;
          mkp->yh=mkp->y+hp;
          /*
            printf("Center (%zu, %zu). r: %.4f. x: [%.4f--%.4f], "
                   "y: [%.4f, %.4f]\n", p%is1+1, p/is1+1, mkp->r, mkp->xl,
                   mkp->xh, mkp->yl, mkp->yh);
          */
          /* Find the random points and profile center. */
          img[p]=randompoints(mkp);
          approx=profile(mkp);
          if (fabs(img[p]-approx)/img[p] < tolerance)
            use_rand_points=0;

          /* Save the peak flux if this is the first pixel: */
          if(ispeak) { mkp->peakflux=img[p]; ispeak=0; }

          /* For the log file: */
          ++ibq->numaccu;
          ibq->accufrac+=img[p];

          /*
            printf("\tac: %f, ap: %f, frac: %f\n", img[p], approx,
            fabs(img[p]-approx)/img[p]);
            gal_fits_array_to_file("tmp.fits", "", FLOAT_IMG, img, is0,
                                   is1, NULL, SPACK_STRING);
          */

          /* Go over the neighbours and add them to queue of elements
             to check. */
          GAL_NEIGHBORS_FILL_4_ALLIMG;
          nf=(n=ngb)+numngb;
          do
            if(byt[*n]==0)
              {
                byt[*n]=1;
                gal_linkedlist_add_to_tosll_end( &lQ, &sQ, *n,
                                                 r_circle(*n, mkp) );
              }
          while(++n<nf);

          if(use_rand_points==0) break;
        }
    }


  /* All the pixels that required integration or random points are now
     done, so we don't need an ordered array any more. */
  gal_linkedlist_tosll_into_sll(lQ, &Q);


  /* Order doesn't matter any more, add all the pixels you find. */
  while(Q)
    {
      gal_linkedlist_pop_from_sll(&Q, ind);        /* ind=&p */
      mkp->x=(p/is1-xc)/os;
      mkp->y=(p%is1-yc)/os;
      r_el(mkp);
      if(mkp->r>truncr)
        {
          /* For the circumference, if the profile is too elongated
             and circumwidth is too small, then some parts of the
             circumference will not be shown without this condition. */
          if(mkp->func==PROFILE_CIRCUMFERENCE) img[p]=profile(mkp);
          continue;
        }

      /* Find the value for this pixel: */
      img[p]=profile(mkp);

      /* Save the peak flux if this is the first pixel: */
      if(ispeak) { mkp->peakflux=img[p]; ispeak=0; }

      /*
      gal_fits_array_to_file("tmp.fits", "", FLOAT_IMG, img, is0, is1,
                              NULL, SPACK_STRING);
      */
      /* Go over the neighbours and add them to queue of elements
         to check. */
      GAL_NEIGHBORS_FILL_4_ALLIMG;
      nf=(n=ngb)+numngb;
      do
        if(byt[*n]==0)
          {
            byt[*n]=1;
            gal_linkedlist_add_to_sll(&Q, *n);
          }
      while(++n<nf);
    }
  free(byt);
}



















/**************************************************************/
/************        Set profile parameters       *************/
/**************************************************************/
int
oneprofile_ispsf(int fcode)
{
  return fcode==PROFILE_MOFFAT || fcode==PROFILE_GAUSSIAN;
}





/* About the shifts on the X column and y column:*/
void
oneprof_set_prof_params(struct mkonthread *mkp)
{
  struct mkprofparams *p=mkp->p;

  double sigma;
  int tp=p->tunitinp;
  size_t id=mkp->ibq->id;

  /* Fill in the profile independant parameters. */
  p->x[id]       += p->shift[0]/p->oversample; /* Shifts were multiplied by */
  p->y[id]       += p->shift[1]/p->oversample; /* `p->oversample' before.   */
  mkp->c          = cos( (90-p->p[id]) * DEGREESTORADIANS );
  mkp->s          = sin( (90-p->p[id]) * DEGREESTORADIANS );
  mkp->q          = p->q[id];
  mkp->brightness = pow( 10, (p->zeropoint - p->m[id]) / 2.5f );
  mkp->ibq->ispsf = oneprofile_ispsf(p->f[id]);
  mkp->func       = mkp->ibq->func = p->f[id];


  /* Fill the profile dependent parameters. */
  switch (mkp->func)
    {
    case PROFILE_SERSIC:
      mkp->correction       = 1;
      mkp->profile          = &Sersic;
      mkp->sersic_re        = p->r[id];
      mkp->sersic_inv_n     = 1.0f/p->n[id];
      mkp->sersic_nb        = -1.0f*sersic_b(p->n[id]);
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id];
      break;



    case PROFILE_MOFFAT:
      mkp->correction       = 1;
      mkp->profile          = &Moffat;
      mkp->moffat_nb        = -1.0f*p->n[id];
      mkp->moffat_alphasq   = moffat_alpha(p->r[id], p->n[id]);
      mkp->moffat_alphasq  *= mkp->moffat_alphasq;
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id]/2;
      if(p->psfinimg==0 && p->individual==0)
        {
          mkp->brightness   = 1.0f; /* When the PSF is a separate image, */
          p->x[id]          = 0.0f; /* it should be centered and have a  */
          p->y[id]          = 0.0f; /* total brightness of 1.0f. */
        }
      break;



    case PROFILE_GAUSSIAN:
      mkp->correction       = 1;
      mkp->profile          = &Gaussian;
      sigma                 = p->r[id]/2.35482f;
      mkp->gaussian_c       = -1.0f/(2.0f*sigma*sigma);
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id]/2;
      if(p->psfinimg==0 && p->individual==0)
        {
          mkp->brightness   = 1.0f; /* Same as the explanations for    */
          p->x[id]          = 0.0f; /* The Moffat profile. */
          p->y[id]          = 0.0f;
        }
      break;



    case PROFILE_POINT:
      mkp->correction       = 1;
      mkp->fixedvalue       = 1.0f;
      mkp->profile          = &Flat;
      break;



    case PROFILE_FLAT:
      mkp->profile          = &Flat;
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id];
      if(p->mforflatpix)
        {
          mkp->correction   = 0;
          mkp->fixedvalue   = p->m[id];
        }
      else
        {
          mkp->correction   = 1;
          mkp->fixedvalue   = 1.0f;
        }
      break;



    case PROFILE_CIRCUMFERENCE:
      mkp->profile          = &Circumference;
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id];
      mkp->intruncr         = mkp->truncr - p->circumwidth;
      if(p->mforflatpix)
        {
          mkp->correction   = 0;
          mkp->fixedvalue   = p->m[id];
        }
      else
        {
          mkp->correction   = 1;
          mkp->fixedvalue   = 1.0f;
        }
      if(mkp->intruncr<0.0f)
        mkp->intruncr       = 0.0f;
      break;



    default:
      error(EXIT_FAILURE, 0, "a bug in setprofparams (oneprofile.c)! "
            "The profile code is not recognized. This should have been "
            "seen and reported prior to this step. Please contact us so "
            "we can correct this");
    }
}




















/**************************************************************/
/************          Outside functions          *************/
/**************************************************************/
void
oneprofile_make(struct mkonthread *mkp)
{
  struct mkprofparams *p=mkp->p;

  float *f, *ff;
  long os=p->oversample;
  double sum, pixfrac, intpart;
  size_t size, id=mkp->ibq->id;


  /* Find the profile center (see comments above
     `mkprof_build'). mkp->width is in the non-oversampled scale.*/
  pixfrac = modf(fabs(p->x[id]), &intpart);
  mkp->yc = ( os * (mkp->width[0]/2 + pixfrac)
              + (pixfrac<0.50f ? os/2 : -1*os/2-1) );
  mkp->yc = round(mkp->yc*100)/100;

  pixfrac = modf(fabs(p->y[id]), &intpart);
  mkp->xc = ( os*(mkp->width[1]/2 + pixfrac)
              + (pixfrac<0.5f ? os/2 : -1*os/2-1) );
  mkp->xc = round(mkp->xc*100)/100;


  /* From this point on, the widths are the actual pixel
     widths (with onversampling). */
  mkp->width[0] *= os;
  mkp->width[1] *= os;
  mkp->ibq->imgwidth=mkp->width[0];


  /* Allocate and clear the array for this one profile. */
  errno=0;
  size=mkp->width[0]*mkp->width[1];
  mkp->ibq->img=calloc(size, sizeof *mkp->ibq->img);
  if(mkp->ibq->img==NULL)
    error(EXIT_FAILURE, 0, "%zu bytes for object in row %zu of data in %s",
          size*sizeof *mkp->ibq->img, mkp->ibq->id, p->catname);


  /* Build the profile in the image. */
  makepixbypix(mkp);


  /* Correct the sum of pixels in the profile so it has the fixed total
     magnitude or pixel value, mkp->correction was set in
     setprofparams. Note that the profiles were not normalized during the
     building.*/
  if(mkp->correction)
    {
      /* First get the sum of all the pixels in the profile. */
      sum=0.0f; ff=(f=mkp->ibq->img)+size; do sum+=*f++; while(f<ff);

      /* Correct the fraction of brightness that was calculated
         accurately (not using the pixel center). */
      mkp->ibq->accufrac /= sum;

      /* Correct all the profile pixels. */
      if(p->magatpeak)
        {
          ff=(f=mkp->ibq->img)+size;
          do *f++ *= mkp->brightness/mkp->peakflux; while(f<ff);
        }
      else
        {
          ff=(f=mkp->ibq->img)+size;
          do *f++ *= mkp->brightness/sum; while(f<ff);
        }
    }
}
