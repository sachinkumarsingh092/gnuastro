/*********************************************************************
Polygon related functions and macros.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
     Sachin Kumar Singh <sachinkumarsingh092@gmail.com>
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
#ifndef __GAL_POLYGON_H__
#define __GAL_POLYGON_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */



/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */



/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */



#define GAL_POLYGON_MAX_CORNERS  50
#define GAL_POLYGON_ROUND_ERR    1e-5





/***************************************************************/
/**************     Function declarations     ******************/
/***************************************************************/
void
gal_polygon_vertices_sort_convex(double *in, size_t n, size_t *ordinds);

int
gal_polygon_is_convex(double *v, size_t n);

double
gal_polygon_area(double *v, size_t n);

int
gal_polygon_is_inside(double *v, double *p, size_t n);

int
gal_polygon_is_inside_convex(double *v, double *p, size_t n);

int
gal_polygon_ppropin(double *v, double *p, size_t n);

int
gal_polygon_is_counterclockwise(double *v, size_t n);

int
gal_polygon_to_counterclockwise(double *v, size_t n);

void
gal_polygon_clip(double *s, size_t n, double *c, size_t m,
                 double *o, size_t *numcrn);

void
gal_polygon_vertices_sort(double *in, size_t n, size_t *ordinds);

__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_POLYGON_H__ */
