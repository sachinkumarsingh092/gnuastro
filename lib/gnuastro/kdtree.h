/*********************************************************************
kdtree -- Create kd-tree and nearest neighbour searches.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Sachin Kumar Singh <sachinkumarsingh092@gmail.com>
Contributing author(s):
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Copyright (C) 2020, Free Software Foundation, Inc.

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
#ifndef __GAL_KDTREE_H__
#define __GAL_KDTREE_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>


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



gal_data_t *
gal_kdtree_create(gal_data_t *coords_raw, size_t *root);

size_t
gal_kdtree_nearest_neighbour(gal_data_t *coords_raw, gal_data_t *kdtree,
                             size_t root, double *point, double *least_dist);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_KDTREE_H__ */
