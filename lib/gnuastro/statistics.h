/*********************************************************************
Statistical functions.
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
#ifndef __GAL_STATISTICS_H__
#define __GAL_STATISTICS_H__

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



/* Maximum number of tests for sigma-clipping convergence. */
#define GAL_STATISTICS_SIG_CLIP_MAX_CONVERGE 50

/* Least acceptable mode symmetricity.*/
#define GAL_STATISTICS_MODE_GOOD_SYM         0.2f





enum bin_status
{
  GAL_STATISTICS_BINS_INVALID,           /* ==0 by C standard.  */

  GAL_STATISTICS_BINS_REGULAR,
  GAL_STATISTICS_BINS_IRREGULAR,
};


/****************************************************************
 ********               Simple statistics                 *******
 ****************************************************************/

gal_data_t *
gal_statistics_number(gal_data_t *input);

gal_data_t *
gal_statistics_minimum(gal_data_t *input);

gal_data_t *
gal_statistics_maximum(gal_data_t *input);

gal_data_t *
gal_statistics_sum(gal_data_t *input);

gal_data_t *
gal_statistics_mean(gal_data_t *input);

gal_data_t *
gal_statistics_std(gal_data_t *input);

gal_data_t *
gal_statistics_mean_std(gal_data_t *input);

gal_data_t *
gal_statistics_median(gal_data_t *input, int inplace);

size_t
gal_statistics_quantile_index(size_t size, double quantile);

gal_data_t *
gal_statistics_quantile(gal_data_t *input, double quantile, int inplace);

size_t
gal_statistics_quantile_function_index(gal_data_t *input, gal_data_t *value,
                                       int inplace);

gal_data_t *
gal_statistics_quantile_function(gal_data_t *input, gal_data_t *value,
                                 int inplace);

gal_data_t *
gal_statistics_unique(gal_data_t *input, int inplace);





/****************************************************************
 ********                     Mode                        *******
 ****************************************************************/

gal_data_t *
gal_statistics_mode(gal_data_t *input, float errorstd, int inplace);

gal_data_t *
gal_statistics_mode_mirror_plots(gal_data_t *input, gal_data_t *value,
                                 size_t numbins, int inplace,
                                 double *mirror_val);





/****************************************************************
 ********                      Sort                       *******
 ****************************************************************/

int
gal_statistics_is_sorted(gal_data_t *input, int updateflags);

void
gal_statistics_sort_increasing(gal_data_t *input);

void
gal_statistics_sort_decreasing(gal_data_t *input);

gal_data_t *
gal_statistics_no_blank_sorted(gal_data_t *input, int inplace);





/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
gal_data_t *
gal_statistics_regular_bins(gal_data_t *data, gal_data_t *range,
                            size_t numbins, double onebinstart);

gal_data_t *
gal_statistics_histogram(gal_data_t *data, gal_data_t *bins,
                         int normalize, int maxhistone);

gal_data_t *
gal_statistics_cfp(gal_data_t *data, gal_data_t *bins, int normalize);





/****************************************************************
 *****************         Outliers          ********************
 ****************************************************************/
gal_data_t *
gal_statistics_sigma_clip(gal_data_t *input, float multip, float param,
                          int inplace, int quiet);

gal_data_t *
gal_statistics_outlier_positive(gal_data_t *input, size_t window_size,
                                float sigma, float sigclip_multip,
                                float sigclip_param, int inplace, int quiet);

gal_data_t *
gal_statistics_outlier_flat_cfp(gal_data_t *input, size_t numprev,
                                float sigclip_multip, float sigclip_param,
                                float thresh, size_t numcontig, int inplace,
                                int quiet, size_t *index);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_STATISTICS_H__ */
