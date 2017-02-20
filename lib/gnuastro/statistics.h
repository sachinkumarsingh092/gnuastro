/*********************************************************************
Statistical functions.
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



/* Maximum number of tests for sigma-clipping convergence */
#define GAL_STATISTICS_MAX_SIG_CLIP_CONVERGE 50





/* Enumerators */
enum is_sorted_outputs
{
  GAL_STATISTICS_SORTED_NOT,             /* ==0 by C standard. */

  GAL_STATISTICS_SORTED_INCREASING,
  GAL_STATISTICS_SORTED_DECREASING,
};


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

gal_data_t *
gal_statistsics_quantile(gal_data_t *input, float quantile, int inplace);

size_t
gal_statistics_quantile_index(size_t size, float quant);





/****************************************************************
 ********                      Sort                       *******
 ****************************************************************/

int
gal_statistics_is_sorted(gal_data_t *data);

void
gal_statistics_sort_increasing(gal_data_t *data);

void
gal_statistics_sort_decreasing(gal_data_t *data);

gal_data_t *
gal_statistics_no_blank_sorted(gal_data_t *input, int inplace);



/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
gal_data_t *
gal_statistics_regular_bins(gal_data_t *data, gal_data_t *range,
                            size_t numbins, float onebinstart);

gal_data_t *
gal_statistics_histogram(gal_data_t *data, gal_data_t *bins,
                         int normalize, int maxhistone);

gal_data_t *
gal_statistics_cfp(gal_data_t *data, gal_data_t *bins, int normalize);





/****************************************************************
 *****************        Sigma clip         ********************
 ****************************************************************/
gal_data_t *
gal_statistics_sigma_clip(gal_data_t *input, float multip, float param,
                          int quiet);





/****************************************************************/
/*************         Identify outliers         ****************/
/****************************************************************/
void
gal_statistics_remove_outliers_flat_cdf(float *sorted, size_t *outsize);





/****************************************************************/
/*************               Mode                ****************/
/****************************************************************/
#define GAL_STATISTICS_MODE_LOW_QUANTILE  0.01f
#define GAL_STATISTICS_MODE_HIGH_QUANTILE 0.51f

#define GAL_STATISTICS_MODE_SYM_GOOD        0.2f
#define GAL_STATISTICS_MODE_LOW_QUANT_GOOD  0.02f

#define GAL_STATISTICS_MODE_SYMMETRICITY_LOW_QUANT 0.01f

struct gal_statistics_mode_params
{
  float     *sorted;   /* Sorted array to be used.                */
  size_t       size;   /* Number of elements in the sorted array. */
  size_t       lowi;   /* Lower quantile of interval.             */
  size_t       midi;   /* Middle quantile of interval.            */
  size_t       midd;   /* Middle index of inteval.                */
  size_t      highi;   /* Higher quantile of interval.            */
  float   tolerance;   /* Tolerance level to terminate search.    */
  size_t   numcheck;   /* Number of pixels after mode to check.   */
  size_t   interval;   /* Interval to check pixels.               */
  float   errorstdm;   /* Multiple of standard deviation.         */
};

void
gal_statistics_mode_mirror_plots(float *sorted, size_t size,
                                 size_t mirrorindex, float min, float max,
                                 size_t numbins, char *histsname,
                                 char *cfpsname, float mirrorplotdist);

float
gal_statistics_mode_value_from_sym(float *sorted, size_t size,
                                   size_t modeindex, float sym);

void
gal_statistics_mode_index_in_sorted(float *sorted, size_t size,
                                    float errorstdm, size_t *modeindex,
                                    float *modesym);


__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_STATISTICS_H__ */
