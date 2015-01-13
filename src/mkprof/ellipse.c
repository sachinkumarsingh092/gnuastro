/*********************************************************************
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
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>


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

   IMPORTANT: x_w and y_w ARE THE AXISES IN THE C STANDARD, NOT FITS.
*/
void
encloseellipse(double a, double b, double theta_rad,
	       size_t *x_w, size_t *y_w)
{
  double t_x, t_y, max_x, max_y;
  t_x=atan(b/a*tan(theta_rad));
  t_y=atan(-1*b/a/tan(theta_rad));

  max_x=a*cos(t_x)*cos(theta_rad)+b*sin(t_x)*sin(theta_rad);
  max_y=-1*a*cos(t_y)*sin(theta_rad)+b*sin(t_y)*cos(theta_rad);

  /* max_x and max_y are calculated from the
     center of the ellipse. We want the final height
     and width of the box enclosing the ellipse. So
     we have to multiply them by two, then take one
     from them (for the center).

     IMPORTANT: x_w and y_w ARE THE AXISES IN THE C STANDARD, NOT FITS.
  */
  *x_w=2*( (size_t)fabs(max_x)+1 ) - 1;
  *y_w=2*( (size_t)fabs(max_y)+1 ) - 1;
}
