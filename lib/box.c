/*********************************************************************
overlap -- Find the overlap of a region and an image.
This is part of GNU Astronomy Utilities (gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "box.h"



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
   that encloses the whole ellipse.
*/
void
ellipseinbox(double a, double b, double theta_rad, long *width)
{
  double t_x, t_y, max_x, max_y;
  t_x=atan(b/a*tan(theta_rad));
  t_y=atan(-1*b/a/tan(theta_rad));

  max_x=a*cos(t_x)*cos(theta_rad)+b*sin(t_x)*sin(theta_rad);
  max_y=-1*a*cos(t_y)*sin(theta_rad)+b*sin(t_y)*cos(theta_rad);

  /* max_x and max_y are calculated from the center of the ellipse. We
     want the final height and width of the box enclosing the
     ellipse. So we have to multiply them by two, then take one from
     them (for the center).
  */
  width[0]=2*( (size_t)fabs(max_x)+1 ) + 1;
  width[1]=2*( (size_t)fabs(max_y)+1 ) + 1;
}





/* We have the central pixel and box size of the crop box, find the
   starting and ending pixels: */
void
borderfromcenter(double xc, double yc, long *width,
		 long *fpixel, long *lpixel)
{
  long lxc, lyc;
  double intpart;

  /* Round the double values in a the long values: */
  lxc=xc;
  if (fabs(modf(xc, &intpart))>=0.5)
    ++lxc;
  lyc=yc;
  if (fabs(modf(yc, &intpart))>=0.5)
    ++lyc;

  /* Set the initial values for the actual image: */
  fpixel[0]=lxc-width[0]/2;      fpixel[1]=lyc-width[1]/2;
  lpixel[0]=lxc+width[0]/2;      lpixel[1]=lyc+width[1]/2;
  /*
  printf("\n\nCenter is on: %ld, %ld\n", lxc, lyc);
  printf("Starting and ending pixels: (%ld, %ld) -- (%ld, %ld)\n\n\n",
	 fpixel[0], fpixel[1], lpixel[0], lpixel[1]);
  */
}





/* Problem to solve: We have set the first and last pixels in an input
   image (fpixel_i[2] and lpixel_i[2]). But those first and last
   pixels don't necessarily have to lie within the image, they can be
   outside of it or patially overlap with it. So the job of this
   function is to correct for such situations and find the starting
   and ending points of any overlap.

   It is assumed that your output (overlap) image's first pixel lies
   right ontop of the fpixel_i[0] in the input image. But since
   fpixel_i might be outside of the image, in this function we find
   the fpixel_o[2] and lpixel_o[2] in the overlap image coordinates
   that overlap with the input image. So the values of all four points
   might change after this function.

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


   So in short the arrays we are dealing with here are:

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
overlap(long *naxes, long *fpixel_i, long *lpixel_i,
	long *fpixel_o, long *lpixel_o)
{
  long width[2];

  /* In case you want to see how things are going:
  printf("\n\nImage size: [%ld,  %ld]\n", naxes[0], naxes[1]);
  printf("fpixel_i -- lpixel_i: (%ld, %ld) -- (%ld, %ld)\n", fpixel_i[0],
         fpixel_i[1], lpixel_i[0], lpixel_i[1]);
  */

  width[0]=lpixel_i[0]-fpixel_i[0]+1;
  width[1]=lpixel_i[1]-fpixel_i[1]+1;

  /* Set the initial values for the cropped array: */
  fpixel_o[0]=1;          fpixel_o[1]=1;
  lpixel_o[0]=width[0];   lpixel_o[1]=width[1];


  /* Check the four corners to see if they should be adjusted: To
     understand the first, look at this, suppose | separates pixels
     and the * shows where the image actually begins.

     |-2|-1| 0* 1| 2| 3| 4|        (survey image)
     *1 | 2| 3| 4| 5| 6| 7|        (crop image)

     So when fpixel_i is negative, e.g., fpixel_i=-2, then the index
     of the pixel in the cropped image we want to begin with, that
     corresponds to 1 in the survey image is: fpixel_o= 4 = 2 +
     -1*fpixel_i.*/
  if(fpixel_i[0]<1)
    {
      fpixel_o[0]=-1*fpixel_i[0]+2;
      fpixel_i[0]=1;
    }
  if(fpixel_i[1]<1)
    {
      fpixel_o[1]=-1*fpixel_i[1]+2;
      fpixel_i[1]=1;
    }


  /*The same principle applies to the end of an image. Take "s"
    is the maximum size along a specific axis in the survey image
    and "c" is the size along the same axis on the cropped image.
    Assume the the cropped region's last pixel in that axis will
    be 2 pixels larger than s:

    |s-1|   s* s+1| s+2|           (survey image)
    |c-3| c-2| c-1|   c*           (crop image)

    So you see that if the outer pixel is n pixels away then in
    the cropped image we should only fill upto c-n.*/
  if (lpixel_i[0]>naxes[0])
    {
      lpixel_o[0]=width[0]-(lpixel_i[0]-naxes[0]);
      lpixel_i[0]=naxes[0];
    }
  if (lpixel_i[1]>naxes[1])
    {
      lpixel_o[1]=width[1]-(lpixel_i[1]-naxes[1]);
      lpixel_i[1]=naxes[1];
    }

  /* In case you wish to see the results.
  printf("\nAfter correction:\n");
  printf("Input image: (%ld, %ld) -- (%ld, %ld)\n", fpixel_i[0],
	 fpixel_i[1], lpixel_i[0], lpixel_i[1]);
  printf("output image:(%ld, %ld) -- (%ld, %ld)\n\n\n", fpixel_o[0],
	 fpixel_o[1], lpixel_o[0], lpixel_o[1]);
  */

  if(fpixel_i[0]>naxes[0] || fpixel_i[1]>naxes[1]
     || lpixel_i[0]<1 || lpixel_i[1]<1)
    return 0;
  else
    return 1;
}
