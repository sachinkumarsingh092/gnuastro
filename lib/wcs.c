/*********************************************************************
Functions to that only use WCSLIB functionality.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include <gsl/gsl_linalg.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>



/**************************************************************/
/**********              Utilities                 ************/
/**************************************************************/
double *
gal_wcs_array_from_wcsprm(struct wcsprm *wcs)
{
  double *out;
  size_t i, j, size=wcs->naxis*wcs->naxis;

  /* Allocate the necessary array. */
  errno=0;
  out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `out' in "
          "`gal_wcs_array_from_wcsprm'", size*sizeof *out);

  /* Fill in the array. */
  if(wcs->altlin |= 1)          /* Has a PCi_j array. */
    {
      for(i=0;i<wcs->naxis;++i)
        for(j=0;j<wcs->naxis;++j)
          out[i*wcs->naxis+j] = wcs->cdelt[i] * wcs->pc[i*wcs->naxis+j];
    }
  else if(wcs->altlin |= 2)     /* Has CDi_j array */
    {
      for(i=0;i<size;++i)
        out[i]=wcs->cd[i];
    }
  else
    error(EXIT_FAILURE, 0, "currently, `gal_wcs_pixel_scale_deg' only "
          "recognizes PCi_ja and CDi_ja keywords");

  /* Return the result */
  return out;
}




















/**************************************************************/
/**********              XY to RADEC               ************/
/**************************************************************/
/* Use the X and Y columns in a larger array to fill the RA and Dec columns
   in that same array. `xy' points to the first element in the X column and
   `radec' points to the first element in the RA column. The columns for Y
   and Dec have to be immediately after X and RA.

   It appears that WCSLIB can only deal with static allocation. At least in
   its tests it only uses static allocation. I tried dynamic allocation,
   but it didn't work. So I can't use the vector functionalities of WCSLIB
   and have to translate each point separately.
*/
void
gal_wcs_xy_array_to_radec(struct wcsprm *wcs, double *xy, double *radec,
                          size_t number, size_t stride)
{
  size_t i;
  double imgcrd[2], phi, theta;
  int status=0, ncoord=1, nelem=2;

  for(i=0;i<number;++i)
    {
      if(isnan(xy[i*stride]) || isnan(xy[i*stride+1]))
        radec[i*stride]=radec[i*stride+1]=NAN;
      else
        {
          wcsp2s(wcs, ncoord, nelem, xy+i*stride, imgcrd, &phi,
                 &theta, radec+i*stride, &status);
          if(status)
            error(EXIT_FAILURE, 0, "wcsp2s ERROR %d: %s", status,
                  wcs_errmsg[status]);

          /* For a check:
             printf("(%f, %f) --> (%f, %f)\n", xy[i*stride], xy[i*stride+1],
                    radec[i*stride], radec[i*stride+1]);
          */
        }
    }
}





/* Similar to the gal_wcs_xyarray_to_radec, but the reverse: to convert an
   array of RA-Dec to X-Y. */
void
gal_wcs_radec_array_to_xy(struct wcsprm *wcs, double *radec, double *xy,
                          size_t number, size_t stride)
{
  size_t i;
  double imgcrd[2], phi, theta;
  int status=0, ncoord=1, nelem=2;

  for(i=0;i<number;++i)
    {
      if(isnan(radec[i*stride]) || isnan(radec[i*stride+1]))
        radec[i*stride]=radec[i*stride+1]=NAN;
      else
        {
          wcss2p(wcs, ncoord, nelem, radec+i*stride, &phi, &theta,
                 imgcrd, xy+i*stride, &status);
          if(status)
            error(EXIT_FAILURE, 0, "wcss2p ERROR %d: %s", status,
                  wcs_errmsg[status]);

          /* For a check:
          printf("(%f, %f) --> (%f, %f)\n", radec[i*stride], radec[i*stride+1],
                 xy[i*stride], xy[i*stride+1]);
          */
        }
    }
}




/* The distance (along a great circle) on a sphere between two points
   is calculated here. Since the pixel sides are usually very small,
   we won't be using the direct formula:

   cos(distance)=sin(d1)*sin(d2)+cos(d1)*cos(d2)*cos(r1-r2)

   We will be using the haversine formula which better considering
   floating point errors (from Wikipedia:)

   sin^2(distance)/2=sin^2( (d1-d2)/2 )+cos(d1)*cos(d2)*sin^2( (r1-r2)/2 )

   Inputs and outputs are all in degrees.
*/
double
gal_wcs_angular_distance_deg(double r1, double d1, double r2, double d2)
{
  /* Convert degrees to radians. */
  double r1r=r1*M_PI/180, d1r=d1*M_PI/180;
  double r2r=r2*M_PI/180, d2r=d2*M_PI/180;

  /* To make things easier to read: */
  double a=sin( (d1r-d2r)/2 );
  double b=sin( (r1r-r2r)/2 );

  /* Return the result: */
  return 2*asin( sqrt( a*a + cos(d1r)*cos(d2r)*b*b) ) * 180/M_PI;
}




/* Return the pixel scale of the image in both dimentions in degrees. */
double *
gal_wcs_pixel_scale_deg(struct wcsprm *wcs)
{
  gsl_vector S;
  gsl_matrix A, V;
  size_t n=wcs->naxis;
  double *a, *v, *pixscale;

  /* Allocate space for the `v' and `pixscale' arrays. */
  errno=0; v=malloc(n*n*sizeof *v);
  if(v==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `v' in "
          "`gal_wcs_pixel_scale_deg'", n*n*sizeof *v);
  errno=0; pixscale=malloc(n*sizeof *pixscale);
  if(pixscale==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `pixscale' in "
          "`gal_wcs_pixel_scale_deg'", n*sizeof *pixscale);

  /* Write the full matrix into an array, irrespective of what type it was
     stored in the wcsprm structure (`PCi_j' style or `CDi_j' style). */
  a=gal_wcs_array_from_wcsprm(wcs);

  /* Fill in the necessary GSL vector and Matrix structures. */
  S.size=n;     S.stride=1;                S.data=pixscale;
  V.size1=n;    V.size2=n;    V.tda=n;     V.data=v;
  A.size1=n;    A.size2=n;    A.tda=n;     A.data=a;

  /* Run GSL's Singular Value Decomposition, using one-sided Jacobi
     orthogonalization which computes the singular (scale) values to a
     higher relative accuracy.*/
  gsl_linalg_SV_decomp_jacobi(&A, &V, &S);

  /* Clean up and return. */
  free(a);
  free(v);
  return pixscale;
}





/* Report the arcsec^2 area of the pixels in the image based on the
   WCS information in that image. We first use the angular distance of
   two edges of one pixel in radians. Then the radians are multiplied
   to give stradians and finally, the stradians are converted to
   arcsec^2. */
double
gal_wcs_pixel_area_arcsec2(struct wcsprm *wcs)
{
  double out;
  double *pixscale;

  /* A small sanity check. Later, when higher dimensions are necessary, we
     can find which ones correlate to RA and Dec and use them to find the
     pixel area in arcsec^2. */
  if(wcs->naxis!=2)
    error(EXIT_FAILURE, 0, "`gal_wcs_pixel_area_arcsec2' can currently "
          "calculate the area only when the image has 2 dimensions.");

  /* Get the pixel scales along each axis in degrees, then multiply. */
  pixscale=gal_wcs_pixel_scale_deg(wcs);

  /* Clean up and return the result. */
  out = pixscale[0] * pixscale[1] * 3600.0f * 3600.0f;
  free(pixscale);
  return out;
}
