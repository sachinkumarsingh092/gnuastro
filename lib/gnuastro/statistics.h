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


/****************************************************************
 *****************    Mininum and Maximum    ********************
 ****************************************************************/
void
gal_statistics_long_non_blank_min(long *in, size_t size, long *min,
                                  long blank);

void
gal_statistics_long_non_blank_max(long *in, size_t size, long *max,
                                  long blank);

void
gal_statistics_float_min(float *in, size_t size, float *min);

void
gal_statistics_float_max(float *in, size_t size, float *max);

void
gal_statistics_double_min(double *in, size_t size, double *min);

double
gal_statistics_double_min_return(double *in, size_t size);

void
gal_statistics_double_max(double *in, size_t size, double *max);

double
gal_statistics_double_max_return(double *in, size_t size);

void
gal_statistics_float_max_masked(float *in, unsigned char *mask, size_t size,
                                float *max);

void
gal_statistics_float_second_max(float *in, size_t size, float *secondmax);

void
gal_statistics_float_second_min(float *in, size_t size, float *secondmin);

void
gal_statistics_f_min_max(float *in, size_t size, float *min, float *max);

void
gal_statistics_d_min_max(double *in, size_t size, double *min, double *max);

void
gal_statistics_d_max_with_index(double *in, size_t size, double *max,
                                size_t *index);

void
gal_statistics_f_max_with_index(float *in, size_t size, float *max,
                                size_t *index);

void
gal_statistics_d_min_with_index(double *in, size_t size, double *min,
                                size_t *index);

void
gal_statistics_f_min_with_index(float *in, size_t size, float *min,
                                size_t *index);





/****************************************************************
 *****************            Sum            ********************
 ****************************************************************/
float
gal_statistics_float_sum(float *in, size_t size);

float
gal_statistics_float_sum_num(float *in, size_t *size);

float
gal_statistics_float_sum_squared(float *in, size_t size);

float
gal_statistics_float_sum_mask(float *in, unsigned char *mask, size_t size,
                              size_t *nsize);

float
gal_statistics_float_sum_mask_l(float *in, long *mask, size_t size,
                                size_t *nsize);

float
gal_statistics_float_sum_squared_mask(float *in, unsigned char *mask,
                                      size_t size, size_t *nsize);

float
gal_statistics_float_sum_squared_mask_l(float *in, long *mask,
                                        size_t size, size_t *nsize);





/****************************************************************
 *****************      Average and          ********************
 ****************    Standard deviation      ********************
 ****************************************************************/
float
gal_statistics_float_average(float *in, size_t size);

double
gal_statistics_double_average(double *in, size_t size);

void
gal_statistics_f_ave(float *in, size_t size, float *ave, unsigned char *mask);

void
gal_statistics_f_ave_l(float *in, size_t size, float *ave, long *mask);

void
gal_statistics_f_ave_std(float *in, size_t size, float *ave,
                         float *std, unsigned char *mask);

void
gal_statistics_f_ave_std_l(float *in, size_t size, float *ave,
                           float *std, long *mask);





/****************************************************************
 *****************           Median            ******************
 ****************************************************************/
float
gal_statistics_median(float *array, size_t insize);

double
gal_statistics_median_double_in_place(double *array, size_t insize);





/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
void
gal_statistics_set_bins(float *sorted, size_t size, size_t numbins,
                        float min, float max, float onebinvalue,
                        float quant, float **obins);

void
gal_statistics_histogram(float *sorted, size_t size, float *bins,
                         size_t numbins, int normhist, int maxhistone);

void
gal_statistics_cumulative_fp(float *sorted, size_t size, float *bins,
                             size_t numbins, int normcfp);

void
gal_statistics_save_hist(float *sorted, size_t size, size_t numbins,
                         char *filename, char *comment);





/****************************************************************
 *****************         Quantiles         ********************
 ****************************************************************/
size_t
gal_statistics_index_from_quantile(size_t size, float quant);





/****************************************************************
 *****************        Sigma clip         ********************
 ****************************************************************/
int
gal_statistics_sigma_clip_converge(float *array, int o1_n0, size_t num_elem,
                                   float sigma_multiple, float accuracy,
                                   float *outave, float *outmed,
                                   float *outstd, int print);

int
gal_statistics_sigma_clip_certain_num(float *array, int o1_n0,
                                      size_t num_elem, float sigma_multiple,
                                      size_t numtimes, float *outave,
                                      float *outmed, float *outstd,
                                      int print);





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
