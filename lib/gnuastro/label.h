/*********************************************************************
label -- Work on labeled (positive integer valued) datasets.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
#ifndef __GAL_LABEL_H__
#define __GAL_LABEL_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>
#include <gnuastro/tile.h>

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





/* Constants for the clump over-segmentation. */
#define GAL_LABEL_INIT      -1
#define GAL_LABEL_RIVER     -2
#define GAL_LABEL_TMPCHECK  -3





/* Functions. */
gal_data_t *
gal_label_indexs(gal_data_t *labels, size_t numlabs, size_t minmapsize,
                 int quietmmap);

size_t
gal_label_watershed(gal_data_t *values, gal_data_t *indexs,
                    gal_data_t *label, size_t *topinds, int min0_max1);

void
gal_label_clump_significance(gal_data_t *values, gal_data_t *std,
                             gal_data_t *label, gal_data_t *indexs,
                             struct gal_tile_two_layer_params *tl,
                             size_t numclumps, size_t minarea, int variance,
                             int keepsmall, gal_data_t *sig,
                             gal_data_t *sigind);

void
gal_label_grow_indexs(gal_data_t *labels, gal_data_t *indexs, int withrivers,
                      int connectivity);




__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_LABEL_H__ */
