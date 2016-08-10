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

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>





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
                          size_t number, size_t width)
{
  size_t i;
  double imgcrd[2], phi, theta;
  int stat, status=0, ncoord=1, nelem=2;

  for(i=0;i<number;++i)
    {
      if(isnan(xy[i*width]) || isnan(xy[i*width+1]))
        radec[i*width]=radec[i*width+1]=NAN;
      else
        {
          status=wcsp2s(wcs, ncoord, nelem, xy+i*width, imgcrd, &phi,
                        &theta, radec+i*width, &stat);
          if(status)
            error(EXIT_FAILURE, 0, "wcsp2s ERROR %d: %s", status,
                  wcs_errmsg[status]);

          /* For a check:
             printf("(%f, %f) --> (%f, %f)\n", xy[i*width], xy[i*width+1],
                    radec[i*width], radec[i*width+1]);
          */
        }
    }
}





/* Similar to the gal_wcs_xyarray_to_radec, but the reverse: to convert an
   array of RA-Dec to X-Y. */
void
gal_wcs_radec_array_to_xy(struct wcsprm *wcs, double *radec, double *xy,
                          size_t number, size_t width)
{
  size_t i;
  double imgcrd[2], phi, theta;
  int stat, status=0, ncoord=1, nelem=2;

  for(i=0;i<number;++i)
    {
      if(isnan(radec[i*width]) || isnan(radec[i*width+1]))
        radec[i*width]=radec[i*width+1]=NAN;
      else
        {
          status=wcss2p(wcs, ncoord, nelem, radec+i*width, &phi, &theta,
                        imgcrd, xy+i*width, &stat);
          if(status)
            error(EXIT_FAILURE, 0, "wcss2p ERROR %d: %s", status,
                  wcs_errmsg[status]);

          /* For a check:
          printf("(%f, %f) --> (%f, %f)\n", radec[i*width], radec[i*width+1],
                 xy[i*width], xy[i*width+1]);
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
*/
double
angulardistance(double r1, double d1, double r2, double d2)
{
  double a=sin( (d1-d2)/2 );
  double b=sin( (r1-r2)/2 );

  return 2*asin( sqrt( a*a + cos(d1)*cos(d2)*b*b) );
}





/* Report the arcsec^2 area of the pixels in the image based on the
   WCS information in that image. We first use the angular distance of
   two edges of one pixel in radians. Then the radians are multiplied
   to give stradians and finally, the stradians are converted to
   arcsec^2. */
double
gal_wcs_pixel_area_arcsec2(struct wcsprm *wcs)
{
  double xy[]={0,0,1,0,0,1};
  double st, *d, *df, radec[6];

  /* Get the RA and Dec of the bottom left, bottom right and top left
     sides of the first pixel in the image. */
  gal_wcs_xy_array_to_radec(wcs, xy, radec, 3, 2);

  /* Covert the RA and dec values to radians for easy calculation: */
  df=(d=radec)+6; do *d++ *= M_PI/180.0f; while(d<df);

  /* For a check:
  printf("\n\nAlong first axis: %g\nAlong second axis: %g\n\n",
         ( angulardistance(radec[0], radec[1], radec[2], radec[3])
           *180/M_PI*3600 ),
         ( angulardistance(radec[0], radec[1], radec[4], radec[5])
           *180/M_PI*3600 ) );
  */

  /* Get the area in stradians. */
  st= ( angulardistance(radec[0], radec[1], radec[2], radec[3]) *
        angulardistance(radec[0], radec[1], radec[4], radec[5]) );

  /* Convert the stradians to arcsec^2:

     1deg^2 = (180/PI)^2 * 1stradian.
     1arcsec^2 = (3600*3600) * 1degree^2
   */
  return st*180.0f*180.0f*3600.0f*3600.0f/(M_PI*M_PI);
}
