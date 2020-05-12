/*********************************************************************
Box -- Define bounding and overlapping boxes.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro/box.h>



/*                        IMPORTANT NOTE:
         All the axises are based on the FITS standard, not C.
*/




/* Any ellipse can be enclosed into a rectangular box. The purpose of this
   function is to give the height and width of that box. The logic behind
   it is this: All the points on the circumference of an ellipse that is
   aligned on the x axis can be written as:

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





/* Find the bounding box of an ellipsoid. An ellipsoid is defined by its
   three axises: the first (a) must be the major axis, the other two must
   be smaller than 'a' but no particular relation between them is
   assumed. We will define the orientation of the ellipsoid from its major
   axis and use the "Proper Euler angles" (ZXZ order) to define that
   orientation.

   The following derivation is taken from:

      https://tavianator.com/exact-bounding-boxes-for-spheres-ellipsoids/

   For the definition of the Euler angles (and the ZXZ order) see Wikipedia
   in the link below. In short: first rotate around the Z axis, then the
   (rotated) first axis, then the (rotated) Z axis.

      https://en.wikipedia.org/wiki/Euler_angles

   Defining the general point 'p' as (the transpose of 'p' with 'p^T' and
   its inverse with 'p^-1'):

           | x |
       p = | y |
           | z |
           | 1 |

   A sphere satisfies the following homogenous coordinate matrix and
   general quadratic surface:

           | 1  0  0  0  |
       S = | 0  1  0  0  |       -->    p^T * S * p = 0
           | 0  0  1  0  |
           | 0  0  0  -1 |

   The conversion from a sphere to an arbitrarily oriented ellipsoid will
   be done by rotating the three semi-axes lengths and merging the three
   vertical vectors into one matrix. The rotation can be written as the
   following combined affine transformation (only in the three main axises
   (since we aren't dealing with translation here). Here, we'll call
   'sin(alpha)' as 's1', 'cos(beta)' as 'c2' and 'sin(gamma)' as 's3'.

                   Rotate by Euler angles
                ----------------------------
        | c1  -s1  0 |   | 1   0    0  |   | c3  -s3   0 |
    R = | s1   c1  0 | * | 0   c2  -s2 | * | s3   c3   0 |
        | 0    0   1 |   | 0   s2   c2 |   | 0    0    1 |

   Then 'M' (rotation and scaling to obtain ellipsoid from sphere) will be:

            | a |          | 0 |          | 0 |              | A1 B1 C1 |
    A = R * | 0 |, B = R * | b |, C = R * | 0 |    -->   M = | A2 B2 C2 |
            | 0 |          | 0 |          | c |              | A3 B3 C3 |

   The final result is:

        |  a*c1*c3-a*s1*c2*s3   -b*c1*s3-b*s1*c2*c3    c*s1*s2  |
    M = |  a*s1*c3+a*c1*c2*s3   -b*s1*s3+b*c1*c2*c3   -c*c1*s2  |
        |  a*s2*s3               b*s2*c3               c*c2     |

   The quadratic surface for this ellipsoid can now be written as:

     (M^-1 * p)^T * S * (M^-1 * p) = 0  --> p^T * (M^-T * S * M^-1) * p = 0

   Writing Q = M^-T * S * M^-1, we get: 'p^T * Q * p = 0'. Now, we define a
   plane with a horizontal vector 'u = [a b c d ]', such that 'u.p=0'. For
   a point on the ellipsoid (at 'p') we have: 'u^T=p^T * Q'. This is
   because: 'u.p = u^T * p = p^T * Q * p = 0' (as we showed above).

   Tangent planes will have the following useful property:

        u^T * Q^-1 * u = p^T * Q * Q^-1 * Q * p = p^T * Q * p = 0

   Now, the particular plane that is perpendicular to the X axis has the
   general form: 'u = [ 1 0 0 -x ]'. So, defining 'R = Q^1', and using the
   property above for tangential planes, we can find the X axis position.

   However, getting to 'R' from 'M' as described above is not easy. So,
   taking the following considerations into account, we can derive the
   final values: [from that webpage] "Several details of the problem can
   make computing the planes more efficient. The first is that 'S' is
   involutory, meaning 'S^-1 = S'. This means that the product 'M * S^-1'
   can be computed implicitly: it is simply 'M' with its last column
   negated. The last column of 'R = M * S^-1 * MT' is the same, because the
   last column of 'M^T' is '[ 0 0 0 1 ]'. In particular, 'R[4,4]=-1'.

   Not all values of RR are used; in fact, only values from the last column
   and the diagonal appear in the formulae. We know the last column
   implicitly, and the diagonal has the formula: "

      R[i,i] = (sum_{j=1 to 3} M[i,j]^2 - M_[i,j]

   So the bounding box lengths along each dimension are the
   following. Recall that in homogenous coordinates, the last column is for
   translation. So in the case of this function all the 'M[i,4]' values are
   zero.

      x = M[1,4] \pm sqrt( M[1,1]^2 + M[1,2]^2 + M[1,3]^2 )
      y = M[2,4] \pm sqrt( M[2,1]^2 + M[2,2]^2 + M[2,3]^2 )
      z = M[3,4] \pm sqrt( M[3,1]^2 + M[3,2]^2 + M[3,3]^2 )        */
void
gal_box_bound_ellipsoid_extent(double *semiaxes, double *euler_deg,
                               double *extent)
{
  size_t i, j;
  double a=semiaxes[0], b=semiaxes[1], c=semiaxes[2];
  double c1=cos(euler_deg[0]*M_PI/180), s1=sin(euler_deg[0]*M_PI/180);
  double c2=cos(euler_deg[1]*M_PI/180), s2=sin(euler_deg[1]*M_PI/180);
  double c3=cos(euler_deg[2]*M_PI/180), s3=sin(euler_deg[2]*M_PI/180);
  double R[9]={ a*c1*c3-a*s1*c2*s3, -1.0f*b*c1*s3-b*s1*c2*c3,   c*s1*s2,
                a*s1*c3+a*c1*c2*s3, -1.0f*b*s1*s3+b*c1*c2*c3,  -1.0f*c*c1*s2,
                a*s2*s3,             b*s2*c3,                   c*c2     };

  /* Sanity check. */
  if(b>a || c>a)
    error(EXIT_FAILURE, 0, "%s: the second and third semi-axes lengths "
          "(%g, %g respectively) must both be smaller or equal to the "
          "first (%g)", __func__, b, c, a);

  /* Find the width along each dimension. */
  for(i=0;i<3;++i)
    {
      extent[i]=0.0;
      for(j=0;j<3;++j)
        extent[i] += R[i*3+j] * R[i*3+j];
      extent[i] = sqrt(extent[i]);
    }

  /* For a check:
  printf("Axies: %g, %g, %g\n", a, b, c);
  printf("Euler angles: %g, %g, %g\n", euler_deg[0], euler_deg[1],
         euler_deg[2]);
  printf("c1: %f, s1: %f\nc2: %f, s2: %f\nc3: %f, s3: %f\n",
         c1, s1, c2, s2, c3, s3);
  PRINT3BY3("R", R);
  printf("extent: %ld, %ld, %ld\n", extent[0], extent[1], extent[2]);
  exit(0);
  */
}





/* Using 'gal_box_bound_ellipsoid_extent', find the integer width of a box
   that contains the ellipsoid. */
#define PRINT3BY3(C, A){                                                \
    printf("%s: | %-15g%-15g%-15g |\n"                                  \
           "   | %-15g%-15g%-15g |\n"                                   \
           "   | %-15g%-15g%-15g |\n\n", (C), (A)[0], (A)[1], (A)[2],   \
           (A)[3], (A)[4], (A)[5], (A)[6], (A)[7], (A)[8]);             \
  }
void
gal_box_bound_ellipsoid(double *semiaxes, double *euler_deg, long *width)
{
  size_t i;
  double extent[3];

  /* Find the extent of the ellipsoid in each axis. */
  gal_box_bound_ellipsoid_extent(semiaxes, euler_deg, extent);

  /* Find the width along each dimension. */
  for(i=0;i<3;++i)
    width[i] = 2 * ( (long)extent[i] ) + 1;
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
   of this function is to correct for such situations and find the starting
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
          /* Along any dimension, if 'lpixel_i' is also smaller than 1,
             then there is no overlap. */
          if(lpixel_i[i]<1) return 0;

          /* Correct the coordinates. */
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
          /* Along any dimension, if 'fpixel_i' is larger than the image
             size, there is no overlap. */
          if(fpixel_i[i]>naxes[i]) return 0;

          /* Correct the coordinates. */
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
