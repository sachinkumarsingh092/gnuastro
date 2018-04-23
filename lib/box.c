/*********************************************************************
Box -- Define bounding and overlapping boxes.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2018, Free Software Foundation, Inc.

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

#include <gnuastro/box.h>



/*                        IMPORTANT NOTE:
         All the axises are based on the FITS standard, not C.
*/




/* Any ellipse can be enclosed into a rectangular box. The purpose of
   this function is to give the height and width of that box. The
   logic behind it is this: All the points on the circumference of an
   ellipse that is aligned on the x axis can be written as:

   (acos(t),bsin(t)) where 0<t<2\pi.               (1)

   But when we rotate the ellipse by \theta, the points can be
   characterized by:

   ( acos(t)cos(\theta)+bsin(t)sin(\theta),        (2)
    -acos(t)sin(\theta)+bsin(t)cos(\theta) )

   To find the maximum and minimum points of this function you just
   have to take the derivative of each with respect to "t" and set it
   to zero. This will give you the "t" that maximizes both x and the
   "t" that maximizes y.  Once you do that, you will get:

   For x: tan(t)=(b/a)tan(\theta)                  (3)
   For y: tan(t)=(-b/a)cot(\theta)

   Once you find the "t", put it in (2) for the respective coordinate
   and you will find the distance (about the center of the ellipse
   that encloses the whole ellipse. */
void
gal_box_bound_ellipse_extent(double a, double b, double theta_deg,
                             double *extent)
{
  double t_r=theta_deg*M_PI/180;
  double ct=cos(t_r), st=sin(t_r);
  double t_x=atan(b/a*tan(t_r)), t_y=atan(-1.0f*b/a/tan(t_r));

  /* Calculate the maxima along each direction. */
  extent[0] = fabs( a*cos(t_x)*ct    + b*sin(t_x)*st );
  extent[1] = fabs( -1*a*cos(t_y)*st + b*sin(t_y)*ct );
}





void
gal_box_bound_ellipse(double a, double b, double theta_deg, long *width)
{
  double extent[2];

  /* Find the extent of the ellipse. */
  gal_box_bound_ellipse_extent(a, b, theta_deg, extent);

  /* max_x and max_y are calculated from the center of the ellipse. We
     want the final height and width of the box enclosing the
     ellipse. So we have to multiply them by two, then take one from
     them (for the center). */
  width[0] = 2 * ( (long)extent[0] ) + 1;
  width[1] = 2 * ( (long)extent[1] ) + 1;
}





/* We have the central pixel and box size of the crop box, find the
   starting (first) and ending (last) pixels: */
void
gal_box_border_from_center(double *center, size_t ndim, long *width,
                           long *fpixel, long *lpixel)
{
  size_t i;
  long tmp;
  double intpart;

  /* Go over the dimensions. */
  for(i=0;i<ndim;++i)
    {
      /* When the decimal fraction of the center (in floating point) is
         larger than 0.5, then it is must be incremented. */
      tmp = center[i] + ( fabs(modf(center[i], &intpart))>=0.5 ? 1 : 0 );

      /* Set the first and last pixels in this dimension. */
      fpixel[i]=tmp-width[i]/2;
      lpixel[i]=tmp+width[i]/2;
    }

  /* For a check:
  if(ndim==2)
    {
      printf("center: %g, %g\n", center[0], center[1]);
      printf("fpixel: %ld, %ld\n", fpixel[0], fpixel[1]);
      printf("lpixel: %ld, %ld\n", lpixel[0], lpixel[1]);
    }
  else if(ndim==3)
    {
      printf("center: %g, %g, %g\n", center[0], center[1], center[2]);
      printf("fpixel: %ld, %ld, %ld\n", fpixel[0], fpixel[1], fpixel[2]);
      printf("lpixel: %ld, %ld, %ld\n", lpixel[0], lpixel[1], lpixel[2]);
    }
  */
}






/* Problem to solve: We have set the first and last pixels of a box in an
   input image (fpixel_i[2] and lpixel_i[2]). But those first and last
   pixels don't necessarily lie within the image's boundaries. They can be
   outside of it or patially overlap with it (see examples below). The job
   of this function is to corret for such situations and find the starting
   and ending points of any overlap.

   It is assumed that your output (overlap) image's first pixel lies right
   ontop of the fpixel_i[0] in the input image. But since fpixel_i might be
   outside of the image, in this function we find the fpixel_o[2] and
   lpixel_o[2] in the overlap image coordinates that overlap with the input
   image. So the values of all four points might change after this
   function.

   Before:
   =======

   fpixel_o and lpixel_o are not shown. fpixel_o and lpixel_o point
   to the same place as fpixel_i and lpixel_i, But in the coordinates
   of the overlap image.

                                 -----------------lpixel_i
                                 |  overlap      |
                                 |   image       |
                                 |               |
           ----------------------|------         |
           |                     |     |         |
           |            fpixel_i -----------------
           |                           |
           |                           |
           |                           |
           |      Input image          |
           -----------------------------

   After
   =====

                                 -----------------
                                 |  overlap      |
                                 |   image       |
                                 |               |
           ----------------------|------lpixel_i |
           |                     |     |         |
           |            fpixel_i -----------------
           |                           |
           |                           |
           |                           |
           |      Input image          |
           -----------------------------


   So, in short the arrays we are dealing with here are:

   fpixel_i:    Coordinates of the first pixel in input image.
   lpixel_i:    Coordinates of the last pixel in input image.
   fpixel_o:    Coordinates of the first pixel in overlap image.
   lpixel_o:    Coordinates of the last pixel in overlap image.

   NOTES:
   ======

     - lpixel is the last pixel in the image (not outside of it).

     - The coordinates are in the FITS format.

   Return value:
   =============
   1: There is an overlap
   0: There is no overlap
*/
int
gal_box_overlap(long *naxes, long *fpixel_i, long *lpixel_i,
                long *fpixel_o, long *lpixel_o, size_t ndim)
{
  size_t i;
  long width;

  /* In case you want to see how things are going:
  printf("\n\nImage size: [");
  for(i=0;i<ndim;++i) {printf("%ld, ", naxes[i]);} printf("\b\b]\n");
  printf("fpixel_i -- lpixel_i: (");
  for(i=0;i<ndim;++i) {printf("%ld, ", fpixel_i[i]);} printf("\b\b) -- (");
  for(i=0;i<ndim;++i) {printf("%ld, ", lpixel_i[i]);} printf("\b\b)\n");
  */

  /* Do the calculations for each dimension: */
  for(i=0;i<ndim;++i)
    {
      /* Set the width and initialize the overlap values. */
      fpixel_o[i] = 1;
      lpixel_o[i] = width = lpixel_i[i] - fpixel_i[i] + 1;


      /* Check the four corners to see if they should be adjusted: To
         understand the first, look at this, suppose | separates pixels and
         the * shows where the image actually begins.

             |-2|-1| 0* 1| 2| 3| 4|        ( input image )
             *1 | 2| 3| 4| 5| 6| 7|        ( crop image  )

         So when fpixel_i is negative, e.g., fpixel_i=-2, then the index of
         the pixel in the cropped image we want to begin with, that
         corresponds to 1 in the survey image is:

             fpixel_o = 4 = 2 + -1*fpixel_i
      */
      if(fpixel_i[i]<1)
        {
          fpixel_o[i] = -1*fpixel_i[i]+2;
          fpixel_i[i] = 1;
        }


      /*The same principle applies to the end of an image. Take "s" is the
        maximum size along a specific axis in the survey image and "c" is
        the size along the same axis on the cropped image.  Assume the the
        cropped region's last pixel in that axis will be 2 pixels larger
        than s:

             |s-1|   s* s+1| s+2|           (survey image)
             |c-3| c-2| c-1|   c*           (crop image)

        So you see that if the outer pixel is n pixels away then in the
        cropped image we should only fill upto c-n.*/
      if(lpixel_i[i]>naxes[i])
        {
          lpixel_o[i] = width - (lpixel_i[i]-naxes[i]);
          lpixel_i[i] = naxes[i];
        }
    }

  /* In case you want to see the results.
  printf("\nAfter correction:\n");
  printf("Input image: (");
  for(i=0;i<ndim;++i) {printf("%ld, ", fpixel_i[i]);} printf("\b\b) -- (");
  for(i=0;i<ndim;++i) {printf("%ld, ", lpixel_i[i]);} printf("\b\b)\n");
  printf("output image: (");
  for(i=0;i<ndim;++i) {printf("%ld, ", fpixel_o[i]);} printf("\b\b) -- (");
  for(i=0;i<ndim;++i) {printf("%ld, ", lpixel_o[i]);} printf("\b\b)\n");
  */


  /* If the first image pixel in every dimension is larger than the input
     image's width or the last image pixel is before the input image, then
     there is no overlap and we must return 0. */
  for(i=0;i<ndim;++i)
    if( fpixel_i[i]>naxes[i] || lpixel_i[i]<1 )
      return 0;


  /* The first and last image pixels are within the image, so there is
     overlap. */
  return 1;
}
