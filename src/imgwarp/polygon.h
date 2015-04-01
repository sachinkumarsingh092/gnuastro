/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (gnuastro) package.

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
#ifndef POLYGON_H
#define POLYGON_H





#define MAXPOLYGONCORNERS 50
#define ROUNDERR          1e-5





/***************************************************************/
/**************     Function declarations     ******************/
/***************************************************************/
void
orderedpolygoncorners(double *in, size_t n, size_t *ordinds);

double
polygonarea(double *v, size_t n);

int
pinpolygon(double *v, double *p, size_t n);

int
ppropinpolygon(double *v, double *p, size_t n);

void
polygonclip(double *s, size_t n, double *c, size_t m,
            double *o, size_t *numcrn);


















/***************************************************************/
/**************            MACROS             ******************/
/***************************************************************/
/* The cross product of two points from the center. */
#define crossproduct(A, B) ( (A)[0]*(B)[1] - (B)[0]*(A)[1] )




/* Find the cross product (2*area) between three points. Each point is
   assumed to be a pointer that has atleast two values within it. */
#define tricrossproduct(A, B, C)                  \
  ( ( (B)[0]-(A)[0] ) * ( (C)[1]-(A)[1] ) -       \
    ( (C)[0]-(A)[0] ) * ( (B)[1]-(A)[1] ) )       \





/* We have the line A-B. We want to see if C is to the left of this
   line or to its right. This function will return 1 if it is to the
   left. It uses the basic property of vector multiplication: If the
   three points are anti-clockwise (the point is to the left), then
   the vector multiplication is positive, if it is negative, then it
   is clockwise (c is to the right).

   Ofcourse it is very important that A be below or equal to B in both
   the X and Y directions. The rounding error might give
   -0.0000000000001 (I didn't count the number of zeros!!) instead of
   zero for the area. Zero would indicate that they are on the same
   line in this case this should give a true result.
*/
#define pleftofline(A, B, C)                            \
  ( tricrossproduct((A), (B), (C)) > -ROUNDERR ) /* >= 0 */




/* See if the three points are collinear, similar to pleftofline
   except that the result has to be exactly zero. */
#define pcollinearwithline(A, B, C)                            \
  ( tricrossproduct((A), (B), (C)) > -ROUNDERR                 \
    && tricrossproduct((A), (B), (C)) < ROUNDERR) /* == 0 */




/* Similar to pleftofline except that if they are on the same line,
   this will return 0 (so that it is not on the left). Therefore the
   name is "proper left". */
#define ppropleftofline(A, B, C)                            \
  ( tricrossproduct((A), (B), (C)) > ROUNDERR ) /* > 0   */


#define minoftwo(A, B) ( (A)<(B)+ROUNDERR ? (A) : (B) )
#define maxoftwo(A, B) ( (A)>(B)-ROUNDERR ? (A) : (B) )
#endif
