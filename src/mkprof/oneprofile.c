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
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <sys/time.h>		 /* generate random seed */
#include <gsl/gsl_rng.h>	 /* used in setrandoms   */
#include <gsl/gsl_randist.h>	 /* To make noise.       */
#include <gsl/gsl_integration.h> /* gsl_integration_qng  */

#include "neighbors.h"
#include "arraymanip.h"
#include "statistics.h"
#include "linkedlist.h"
#include "fitsarrayvv.h"

#include "main.h"
#include "mkprof.h"
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





float
r_circle(size_t p, struct mkonthread *mkp)
{
  double x, y;

  x = p/mkp->width[0];		/* Note that width[0] is the First FITS */
  y = p%mkp->width[0];		/* axis, not first C axis.              */

  return sqrt( (x-mkp->xc)*(x-mkp->xc) + (y-mkp->yc)*(y-mkp->yc) );
}



















/****************************************************************
 **************          Random points         ******************
 ****************************************************************/
unsigned long int
random_seed()
{
  struct timeval tv;
  gettimeofday(&tv,0);
  return(tv.tv_sec + tv.tv_usec);
}





/* Fill pixel with random values */
float
randompoints(struct mkonthread *mkp)
{
  gsl_rng * rng;
  const gsl_rng_type *T;
  double xrange, yrange, sum=0.0f;
  size_t i, numrandom=mkp->p->numrandom;

  xrange=mkp->xh-mkp->xl;
  yrange=mkp->yh-mkp->yl;

  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  gsl_rng_set(rng,random_seed());

  for(i=0;i<numrandom;++i)
    {
      mkp->x = mkp->xl + gsl_rng_uniform(rng)*xrange;
      mkp->y = mkp->yl + gsl_rng_uniform(rng)*yrange;
      r_el(mkp);
      sum+=mkp->profile(mkp);
    }

  gsl_rng_free(rng);

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
void
makepixbypix(struct mkonthread *mkp)
{
  float circ_r;
  struct sll *Q=NULL;
  unsigned char *byt;
  int use_rand_points=1;
  float *img=mkp->ibq->img;
  size_t is1=mkp->width[0], is0=mkp->width[1];
  size_t p, ngb[4], *ind=&p, numngb, *n, *nf, x, y;
  double tolerance=mkp->p->tolerance, pixfrac, junk;
  double (*profile)(struct mkonthread *)=mkp->profile;
  double xc=mkp->xc, yc=mkp->yc, os=mkp->p->oversample;
  double truncr=mkp->truncr, approx, hp=mkp->p->halfpixel;
  struct tosll *lQ=NULL, *sQ;	/* lQ: Largest. sQ: Smallest in queue */

  /* Find the nearest pixel to the profile center and add it to the
     queue. */
  pixfrac = modf(mkp->xc, &junk);
  x=(long)mkp->xc + ( pixfrac<0.5f ? 0 : 1 );
  pixfrac = modf(mkp->yc, &junk);
  y=(long)mkp->yc + ( pixfrac<0.5f ? 0 : 1 );
  p=x*mkp->width[0]+y;

  /* If this is a point source, just fill that one pixel and go. */
  if(mkp->type==POINTCODE)
    { img[p]=1; return; }

  /* Allocate the byt array to not repeat completed pixels. */
  errno=0;
  byt=calloc(is0*is1, sizeof *byt);
  if(byt==NULL)
    error(EXIT_FAILURE, 0, "%lu bytes for map of object in row %lu of "
	  "data in %s.", is0*is1*sizeof *byt, mkp->ibq->id,
	  mkp->p->up.catname);

  /* Start filling the queue. */
  byt[p]=1;
  add_to_tosll_end( &lQ, &sQ, p, r_circle(p, mkp) );
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
      pop_from_tosll_start(&lQ, &sQ, ind, &circ_r); /* ind=&p */
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
      printf("Center (%lu, %lu). r: %.4f. x: [%.4f--%.4f], y: [%.4f, %.4f]\n",
	     p%is1+1, p/is1+1, mkp->r, mkp->xl, mkp->xh, mkp->yl, mkp->yh);
      */
      /* Find the random points and profile center. */
      img[p]=randompoints(mkp);
      approx=profile(mkp);
      if (fabs(img[p]-approx)/img[p] < tolerance)
	use_rand_points=0;
      /*
      printf("\tac: %f, ap: %f, frac: %f\n", img[p], approx,
	     fabs(img[p]-approx)/img[p]);
      arraytofitsimg("tmp.fits", "", FLOAT_IMG, img, is0, is1,
		     NULL, SPACK_STRING);
      */

      /* Go over the neighbours and add them to queue of elements
	 to check. */
      FILL_NGB_4_ALLIMG;
      nf=(n=ngb)+numngb;
      do
	if(byt[*n]==0)
	  {
	    byt[*n]=1;
	    add_to_tosll_end( &lQ, &sQ, *n, r_circle(*n, mkp) );
	  }
      while(++n<nf);

      if(use_rand_points==0) break;
    }



  /* All the pixels that required integration or random points are now
     done, so we don't need an ordered array any more. */
  tosll_into_sll(lQ, &Q);


  /* Order doesn't matter any more, add all the pixels you find. */
  while(Q)
    {
      pop_from_sll(&Q, ind);	/* ind=&p */
      mkp->x=(p/is1-xc)/os;
      mkp->y=(p%is1-yc)/os;
      r_el(mkp);
      if(mkp->r>truncr) continue;

      /* Find the value for this pixel: */
      img[p]=profile(mkp);

      /*
      arraytofitsimg("tmp.fits", "", FLOAT_IMG, img, is0, is1,
		     NULL, SPACK_STRING);
      */
      /* Go over the neighbours and add them to queue of elements
	 to check. */
      FILL_NGB_4_ALLIMG;
      nf=(n=ngb)+numngb;
      do
	if(byt[*n]==0)
	  {
	    byt[*n]=1;
	    add_to_sll(&Q, *n);
	  }
      while(++n<nf);
    }
  free(byt);
}



















/**************************************************************/
/************        Set profile parameters       *************/
/**************************************************************/
void
setprofparams(struct mkonthread *mkp)
{
  struct mkprofparams *p=mkp->p;

  double *cat, sigma;
  int tp=p->tunitinp;
  size_t rcol=p->rcol, tcol=p->tcol;

  /* Fill in the profile independant parameters. */
  cat=&p->cat[mkp->ibq->id*p->cs1];
  cat[p->xcol] += p->shift[0];
  cat[p->ycol] += p->shift[1];
  mkp->c        = cos((90-cat[p->pcol])*DEGREESTORADIANS);
  mkp->s        = sin((90-cat[p->pcol])*DEGREESTORADIANS);
  mkp->q        = cat[p->qcol];
  mkp->totflux  = pow( 10, (p->zeropoint - cat[p->mcol]) / 2.5f );
  mkp->type     = cat[p->fcol];

  /* Fill the profile dependent parameters. */
  switch (mkp->type)
    {
    case SERSICCODE:
      mkp->ispsf            = 0;
      mkp->profile          = &Sersic;
      mkp->sersic_re        = cat[rcol];
      mkp->sersic_inv_n     = 1.0f/cat[p->ncol];
      mkp->sersic_nb        = -1.0f*sersic_b(cat[p->ncol]);
      mkp->truncr           = tp ? cat[tcol] : cat[tcol]*cat[rcol];
      break;

    case MOFFATCODE:
      mkp->ispsf            = 1;
      mkp->profile          = &Moffat;
      mkp->moffat_nb        = -1.0f*cat[p->ncol];
      mkp->moffat_alphasq   = moffat_alpha(cat[rcol], cat[p->ncol]);
      mkp->moffat_alphasq  *= mkp->moffat_alphasq;
      mkp->truncr           = tp ? cat[tcol] : cat[tcol]*cat[rcol]/2;
      if(p->psfinimg==0 && p->individual==0)
	{ mkp->totflux=1.0f; cat[p->xcol]=0.0f; cat[p->ycol]=0.0f; }
      break;

    case GAUSSIANCODE:
      mkp->ispsf            = 1;
      mkp->profile          = &Gaussian;
      sigma                 = cat[rcol]/2.35482f;
      mkp->gaussian_c       = -1.0f/(2.0f*sigma*sigma);
      mkp->truncr           = tp ? cat[tcol] : cat[tcol]*cat[rcol]/2;
      if(p->psfinimg==0 && p->individual==0)
	{ mkp->totflux=1.0f; cat[p->xcol]=0.0f; cat[p->ycol]=0.0f; }
      break;

    case POINTCODE:
      mkp->ispsf            = 0;
      mkp->profile          = &Point;
      mkp->point_v          = 1;
      break;

    default:
      error(EXIT_FAILURE, 0, "A bug in setprofparams (oneprofile.c)! "
	    "The profile code is not recognized. This should have been "
	    "seen and reported prior to this step. Please contact us so "
	    "we can correct this.");
    }

}




















/**************************************************************/
/************          Outside functions          *************/
/**************************************************************/
void
makeoneprofile(struct mkonthread *mkp)
{
  float sum;
  size_t size;


  /* From this point on, the widths are the actual pixel widths (with
     onversampling). */
  mkp->width[0] *= mkp->p->oversample;
  mkp->width[1] *= mkp->p->oversample;


  /* Allocate the array and build it. */
  errno=0;
  size=mkp->width[0]*mkp->width[1];
  mkp->ibq->img=calloc(size, sizeof *mkp->ibq->img);
  if(mkp->ibq->img==NULL)
    error(EXIT_FAILURE, 0, "%lu bytes for object in row %lu of data in %s.",
	  size*sizeof *mkp->ibq->img, mkp->ibq->id, mkp->p->up.catname);


  /* Build the profile in the image. */
  makepixbypix(mkp);


  /* Correct the sum of pixels in it. */
  sum=floatsum(mkp->ibq->img, size);
  fmultipconst(mkp->ibq->img, size, mkp->totflux/sum);


  /* Return the widths back to the un-oversampled size: */
  mkp->width[0] /= mkp->p->oversample;
  mkp->width[1] /= mkp->p->oversample;
}
