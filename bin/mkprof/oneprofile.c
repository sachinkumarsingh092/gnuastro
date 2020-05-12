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
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "mkprof.h"              /* Needs main.h and astrthreads.h */
#include "profiles.h"
#include "oneprofile.h"






/****************************************************************
 **************          Radial distance       ******************
 ****************************************************************/
/* Set the center position of the profile in the oversampled image. Note
   that 'mkp->width' is in the non-oversampled scale. IMPORTANT: the
   ordering is in FITS coordinate order. */
static void
oneprofile_center_oversampled(struct mkonthread *mkp)
{
  struct mkprofparams *p=mkp->p;

  long os=p->oversample;
  double *dim, r=1000000;
  size_t i, id=mkp->ibq->id;
  double val, pixfrac, intpart;

  for(i=0;i<p->ndim;++i)
    {
      dim = i==0 ? p->x : (i==1 ? p->y : p->z);
      pixfrac = modf(fabs(dim[id]), &intpart);
      val     = ( os*(mkp->width[i]/2 + pixfrac)
                  + (pixfrac<0.5f ? os/2 : -1*os/2-1) );
      mkp->center[i] = round(val*r)/r;
    }
}





static void
oneprofile_set_coord(struct mkonthread *mkp, size_t index)
{
  size_t i, coord_c[3];
  uint8_t os=mkp->p->oversample;
  size_t ndim=mkp->ibq->image->ndim, *dsize=mkp->ibq->image->dsize;

  /* Get the coordinates in C order. */
  gal_dimension_index_to_coord(index, ndim, dsize, coord_c);

  /* Convert these coordinates to one where the profile center is at the
     center and the image is not over-sampled. Note that only 'coord_c' is
     in C order.*/
  for(i=0;i<ndim;++i)
    mkp->coord[i] = ( coord_c[ndim-i-1] - mkp->center[i] )/os;
}





/* Convert cartesian coordinates to the rotated elliptical radius. See the
   "Defining an ellipse and ellipsoid" section of the book for the full
   derivation. */
static void
oneprofile_r_el(struct mkonthread *mkp)
{
  double Xr, Yr, Zr;                   /* Rotated x, y, z. */
  double q1=mkp->q[0],   q2=mkp->q[1];
  double c1=mkp->c[0],   s1=mkp->s[0];
  double c2=mkp->c[1],   s2=mkp->s[1];
  double c3=mkp->c[2],   s3=mkp->s[2];
  double x=mkp->coord[0], y=mkp->coord[1], z=mkp->coord[2];

  switch(mkp->p->ndim)
    {
    case 2:
      /* The parenthesis aren't necessary, but help in readability and
         avoiding human induced bugs. */
      Xr = x * ( c1       )     +   y * ( s1 );
      Yr = x * ( -1.0f*s1 )     +   y * ( c1 );
      mkp->r = sqrt( Xr*Xr + Yr*Yr/q1/q1 );
      break;

    case 3:
      Xr = x*(  c3*c1   - s3*c2*s1 ) + y*( c3*s1   + s3*c2*c1) + z*( s3*s2 );
      Yr = x*( -1*s3*c1 - c3*c2*s1 ) + y*(-1*s3*s1 + c3*c2*c1) + z*( c3*s2 );
      Zr = x*(  s1*s2              ) + y*(-1*s2*c1           ) + z*( c2    );
      mkp->r = sqrt( Xr*Xr + Yr*Yr/q1/q1 + Zr*Zr/q2/q2 );
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
            "the problem. The value %zu is not recognized for "
            "'mkp->p->ndim'", __func__, PACKAGE_BUGREPORT, mkp->p->ndim);
    }
}





/* Calculate the circular/spherical distance of a pixel to the profile
   center. This is just used to add pixels in the stack. Later, when the
   pixels are popped from the stack, the elliptical radius will be used to
   give them a value.*/
static float
oneprofile_r_circle(size_t index, struct mkonthread *mkp)
{

  size_t i, c[3];
  double d, sum=0.0f;
  size_t ndim=mkp->ibq->image->ndim, *dsize=mkp->ibq->image->dsize;

  /* Convert the index into a coordinate. */
  gal_dimension_index_to_coord(index, ndim, dsize, c);

  /* Find the distance to the center along each dimension (in FITS
     order). */
  for(i=0;i<ndim;++i)
    {
      d = c[ndim-i-1] - mkp->center[i];
      sum += d*d;
    }

  /* Return the distance. */
  return sqrt(sum);
}



















/****************************************************************
 **************          Random points         ******************
 ****************************************************************/
/* Fill pixel with random values */
float
oneprofile_randompoints(struct mkonthread *mkp)
{
  double r_before=mkp->r;
  double range[3], sum=0.0f;
  size_t i, j, numrandom=mkp->p->numrandom, ndim=mkp->p->ndim;
  double coord_before[3]={mkp->coord[0], mkp->coord[1], mkp->coord[2]};

  /* Set the range in each dimension. */
  for(i=0;i<ndim;++i)
    range[i] = mkp->higher[i] - mkp->lower[i];

  /* Find the sum of the profile on the random positions. */
  for(i=0;i<numrandom;++i)
    {
      for(j=0;j<ndim;++j)
        mkp->coord[j] = mkp->lower[j] + gsl_rng_uniform(mkp->rng) * range[j];
      oneprofile_r_el(mkp);
      sum+=mkp->profile(mkp);
    }

  /* Reset the original distance and coordinate of the pixel and return the
     average random value. The resetting is mostly redundant (only useful
     in checks), but since it has a very negligible cost (compared to the
     random checks above) cost, its good to reset it to help in debugging
     when necessary (avoid confusion when un-commenting the checks in
     'oneprofile_pix_by_pix'). */
  mkp->r=r_before;
  mkp->coord[0]=coord_before[0];
  mkp->coord[1]=coord_before[1];
  mkp->coord[2]=coord_before[2];
  return sum/numrandom;
}




















/****************************************************************
 *****************      2D integration       ********************
 ****************************************************************/
/* This is an old implementation which we are not using now. But it is kept
   here in case it can be useful */
#if 0
double
twod_over_x(double x, void *params)
{
  struct mkonthread *mkp=(struct mkonthread *) params;

  mkp->x=x;
  oneprofile_r_el(mkp);
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
#endif


















/**************************************************************/
/************       Pixel by pixel building       *************/
/*********        Positions are in C not FITS         *********/
/**************************************************************/
/* 'oneprofile_center_oversampled' stored the center of the profile in
   floating point coordinates. This function will convert that into a
   pixel index. */
static size_t
oneprofile_center_pix_index(struct mkonthread *mkp)
{
  double pixfrac, intpart;
  size_t *dsize=mkp->ibq->image->dsize;
  size_t i, coord[3], ndim=mkp->p->ndim;

  /* Find the coordinates of the center point. Note 'mkp->center' is in
     FITS coordinates, while coord must be in C coordinates (to be used in
     'gal_dimension_coord_to_index'). */
  for(i=0;i<ndim;++i)
    {
      pixfrac = modf(mkp->center[i], &intpart);
      coord[ndim-i-1] = (long)(mkp->center[i]) + ( pixfrac<0.5f ? 0 : 1 );
    }

  /* Retun the pixel index of this coordinate. */
  return gal_dimension_coord_to_index(ndim, dsize, coord);
}





static void
oneprofile_pix_by_pix(struct mkonthread *mkp)
{
  struct builtqueue *ibq=mkp->ibq;
  size_t ndim=ibq->image->ndim, *dsize=ibq->image->dsize;

  uint8_t *byt;
  gal_list_sizet_t *Q=NULL;
  int use_rand_points=1, ispeak=1;
  double tolerance=mkp->p->tolerance;
  float circ_r, *array=mkp->ibq->image->array;
  double (*profile)(struct mkonthread *)=mkp->profile;
  double truncr=mkp->truncr, approx, hp=0.5f/mkp->p->oversample;
  size_t i, p, *dinc=gal_dimension_increment(ndim, dsize);

  /* lQ: Largest. sQ: Smallest in queue */
  gal_list_dosizet_t *lQ=NULL, *sQ;

  /* Find the nearest pixel to the profile center and add it to the
     queue. */
  p=oneprofile_center_pix_index(mkp);

  /* If this is a point source, just fill that one pixel and leave this
     function. */
  if(mkp->func==PROFILE_POINT)
    { array[p]=1; return; }

  /* Allocate the 'byt' array. It is used as a flag to make sure that we
     don't re-calculate the profile value on a pixel more than once. */
  byt = gal_pointer_allocate(GAL_TYPE_UINT8,
                             gal_dimension_total_size(ndim, dsize), 1,
                             __func__, "byt");

  /* Start the queue: */
  byt[p]=1;
  gal_list_dosizet_add( &lQ, &sQ, p, oneprofile_r_circle(p, mkp) );

  /* If random points are necessary, then do it: */
  switch(mkp->func)
    {
    case PROFILE_SERSIC:
    case PROFILE_MOFFAT:
    case PROFILE_GAUSSIAN:
      while(sQ)
        {
          /* In case you want to see the status of the twosided ordered
             queue, increasing and decreasing side by side, uncomment this
             line. Note that there will be a lot of lines printed! */
          /*print_tossll(lQ, sQ);*/

          /* Pop a pixel from the queue, convert its index into coordinates
             and use them to estimate the elliptical radius of the
             pixel. If the pixel is outside the truncation radius, ignore
             it. */
          p=gal_list_dosizet_pop_smallest(&lQ, &sQ, &circ_r);
          oneprofile_set_coord(mkp, p);
          oneprofile_r_el(mkp);
          if(mkp->r > truncr) continue;

          /* Set the range for this pixel. */
          for(i=0;i<ndim;++i)
            {
              mkp->lower[i]  = mkp->coord[i] - hp;
              mkp->higher[i] = mkp->coord[i] + hp;
            }

          /* Find the random points and profile center. */
          array[p]=oneprofile_randompoints(mkp);
          approx=profile(mkp);
          if (fabs(array[p]-approx)/array[p] < tolerance)
            use_rand_points=0;

          /* For a check:
          printf("coord: %g, %g\n", mkp->coord[0], mkp->coord[1]);
          printf("r_rand: %g (rand: %g, center: %g)\n\n", mkp->r, array[p],
                 approx);
          */

          /* Save the peak flux if this is the first pixel: */
          if(ispeak) { mkp->peakflux=array[p]; ispeak=0; }

          /* For the log file: */
          ++ibq->numaccu;
          ibq->accufrac+=array[p];

          /* Go over the neighbors and add them to queue of elements to
             check if they haven't been done already. */
          GAL_DIMENSION_NEIGHBOR_OP(p, ndim, dsize, 1, dinc,
            {
              if(byt[nind]==0)
                {
                  byt[nind]=1;
                  gal_list_dosizet_add( &lQ, &sQ, nind,
                                        oneprofile_r_circle(nind, mkp) );
                }
            } );

          if(use_rand_points==0) break;
        }
    }


  /* All the pixels that required integration or random points are now
     done, so we don't need an ordered array any more. */
  gal_list_dosizet_to_sizet(lQ, &Q);


  /* Order doesn't matter any more, add all the pixels you find. */
  while(Q)
    {
      p=gal_list_sizet_pop(&Q);
      oneprofile_set_coord(mkp, p);
      oneprofile_r_el(mkp);

      /* See if this pixel's radial distance is larger than the truncation
         radius. If so, then don't add its neighbors to the queue and
         continue to the next pixel in the queue. */
      if(mkp->r>truncr)
        {
          /* For the circumference, if the profile is too elongated
             and circumwidth is too small, then some parts of the
             circumference will not be shown without this condition. */
          if(mkp->func==PROFILE_CIRCUMFERENCE) array[p]=profile(mkp);
          continue;
        }

      /* Find the value for this pixel: */
      array[p]=profile(mkp);

      /* For a check:
      printf("r_center: %g\n", mkp->r);
      */

      /* Save the peak flux if this is the first pixel: */
      if(ispeak) { mkp->peakflux=array[p]; ispeak=0; }

      /* Go over the neighbours and add them to queue of elements
         to check. */
      GAL_DIMENSION_NEIGHBOR_OP(p, ndim, dsize, 1, dinc,
        {
          if(byt[nind]==0)
            {
              byt[nind]=1;
              gal_list_sizet_add(&Q, nind);
            }
        } );
    }

  /* Clean up. */
  free(byt);
  free(dinc);
}



















/**************************************************************/
/************        Set profile parameters       *************/
/**************************************************************/
int
oneprofile_ispsf(uint8_t fcode)
{
  return fcode==PROFILE_MOFFAT || fcode==PROFILE_GAUSSIAN;
}





/* Prepare all the parameters for any type of profile. */
void
oneprofile_set_prof_params(struct mkonthread *mkp)
{
  struct mkprofparams *p=mkp->p;

  double sigma;
  int tp=p->tunitinp;
  size_t id=mkp->ibq->id, ndim=p->ndim;

  /* Fill the most basic profile agnostic parameters. */
  mkp->brightness = ( p->mcolisbrightness
                      ? p->m[id]
                      : pow( 10, (p->zeropoint - p->m[id]) / 2.5f ) );
  mkp->ibq->ispsf = p->kernel ? 1 : oneprofile_ispsf(p->f[id]);
  mkp->func       = mkp->ibq->func = p->f[id];


  /* Fill in the dimension-dependent parameters. */
  switch(ndim)
    {
    case 2:
      /* Shifts were already multiplied with oversample. Just note that
         p->x and p->y are in the FITS ordering, while p->shift is in C
         ordering. */
      mkp->q[0]       = p->q1[id];
      p->x[id]       += p->shift[1]/p->oversample;
      p->y[id]       += p->shift[0]/p->oversample;
      mkp->c[0]       = cos( p->p1[id] * DEGREESTORADIANS );
      mkp->s[0]       = sin( p->p1[id] * DEGREESTORADIANS );
      break;

    case 3:
      /* See comments for 2D. */
      mkp->q[0]       = p->q1[id];
      mkp->q[1]       = p->q2[id];
      p->x[id]       += p->shift[2]/p->oversample;
      p->y[id]       += p->shift[1]/p->oversample;
      p->z[id]       += p->shift[0]/p->oversample;
      mkp->c[0]       = cos( p->p1[id] * DEGREESTORADIANS );
      mkp->s[0]       = sin( p->p1[id] * DEGREESTORADIANS );
      mkp->c[1]       = cos( p->p2[id] * DEGREESTORADIANS );
      mkp->s[1]       = sin( p->p2[id] * DEGREESTORADIANS );
      mkp->c[2]       = cos( p->p3[id] * DEGREESTORADIANS );
      mkp->s[2]       = sin( p->p3[id] * DEGREESTORADIANS );
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
            "address the problem. The value '%zu' is not recognized for "
            "'ndim'", __func__, PACKAGE_BUGREPORT, ndim);
    }


  /* Fill the profile-dependent parameters. */
  switch (mkp->func)
    {
    case PROFILE_SERSIC:
      mkp->correction       = 1;
      mkp->profile          = &profiles_sersic;
      mkp->sersic_re        = p->r[id];
      mkp->sersic_inv_n     = 1.0f/p->n[id];
      mkp->sersic_nb        = -1.0f*profiles_sersic_b(p->n[id]);
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id];
      break;



    case PROFILE_MOFFAT:
      mkp->correction       = 1;
      mkp->profile          = &profiles_moffat;
      mkp->moffat_nb        = -1.0f*p->n[id];
      mkp->moffat_alphasq   = profiles_moffat_alpha(p->r[id], p->n[id]);
      mkp->moffat_alphasq  *= mkp->moffat_alphasq;
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id]/2;
      if(p->psfinimg==0 && p->individual==0)
        {
          mkp->brightness   = 1.0f; /* When the PSF is a separate image, */
          p->x[id]          = 0.0f; /* it should be centered and have a  */
          p->y[id]          = 0.0f; /* total brightness of 1.0f. */
          if(ndim==3)
            p->z[id]        = 0.0f;
        }
      break;



    case PROFILE_GAUSSIAN:
      mkp->correction       = 1;
      mkp->profile          = &profiles_gaussian;
      sigma                 = p->r[id]/2.35482f;
      mkp->gaussian_c       = -1.0f/(2.0f*sigma*sigma);
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id]/2;
      if(p->psfinimg==0 && p->individual==0)
        {
          mkp->brightness   = 1.0f; /* Same as the explanations for    */
          p->x[id]          = 0.0f; /* The Moffat profile. */
          p->y[id]          = 0.0f;
          if(ndim==3)
            p->z[id]        = 0.0f;
        }
      break;



    case PROFILE_POINT:
      mkp->correction       = 1;
      mkp->fixedvalue       = 1.0f;
      mkp->profile          = &profiles_flat;
      break;



    case PROFILE_FLAT:
      mkp->profile          = &profiles_flat;
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
      mkp->profile          = &profiles_circumference;
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



    case PROFILE_DISTANCE:
      mkp->profile          = profiles_radial_distance;
      mkp->truncr           = tp ? p->t[id] : p->t[id]*p->r[id];
      mkp->correction       = 0;
      break;



    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us so we can "
            "correct this problem. The profile code %u is not recognized.",
            __func__, mkp->func);
    }
}




















/**************************************************************/
/************          Outside functions          *************/
/**************************************************************/
void
oneprofile_make(struct mkonthread *mkp)
{
  struct mkprofparams *p=mkp->p;

  double sum;
  float *f, *ff;
  size_t i, dsize[3], ndim=p->ndim;


  /* Find the profile center in the over-sampled image in C
     coordinates. IMPORTANT: width must not be oversampled.*/
  oneprofile_center_oversampled(mkp);


  /* From this point on, the widths will be in the actual pixel widths
     (with oversampling). */
  for(i=0;i<ndim;++i)
    {
      mkp->width[i]  *= p->oversample;
      dsize[ndim-i-1] = mkp->width[i];
    }


  /* Allocate and clear the array for this one profile. */
  mkp->ibq->image=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, ndim, dsize,
                                 NULL, 1, p->cp.minmapsize, p->cp.quietmmap,
                                 "MOCK", "Brightness", NULL);


  /* Build the profile in the image. */
  oneprofile_pix_by_pix(mkp);


  /* Correct the sum of pixels in the profile so it has the fixed total
     magnitude or pixel value, mkp->correction was set in
     setprofparams. Note that the profiles were not normalized during the
     building.*/
  if(mkp->correction)
    {
      /* First get the sum of all the pixels in the profile. */
      ff=(f=mkp->ibq->image->array) + mkp->ibq->image->size;
      sum=0.0f; do sum+=*f++; while(f<ff);

      /* Correct the fraction of brightness that was calculated
         accurately (not using the pixel center). */
      mkp->ibq->accufrac /= sum;

      /* Correct all the profile pixels. Note that ideally, if a user wants
         a NaN valued profile, they should use the 'flat' profile with
         '--mforflatpix', which won't need this correction. However, it
         might happen that they forget the later, or the catalog might be
         generated by a script that gives a NaN value for the magnitude
         with any kind of profile. In such cases if we don't check the NaN
         value, then the whole profile's box is going to be NaN values,
         which is inconvenient and with the simple check here we can avoid
         it (only have the profile's pixels set to NaN. */
      ff = (f=mkp->ibq->image->array) + mkp->ibq->image->size;
      if(isnan(mkp->brightness))
        do *f = *f ? NAN : *f ; while(++f<ff);
      else
        {
          if(p->magatpeak)
            do *f++ *= mkp->brightness/mkp->peakflux; while(f<ff);
          else
            do *f++ *= mkp->brightness/sum; while(f<ff);
        }
    }
}
