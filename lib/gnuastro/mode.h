/*********************************************************************
mode -- Find the mode of a distribution.
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
#ifndef __GAL_MODE_H__
#define __GAL_MODE_H__

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



#define GAL_MODE_LOW_QUANTILE  0.01f
#define GAL_MODE_HIGH_QUANTILE 0.51f

#define GAL_MODE_SYM_GOOD        0.2f
#define GAL_MODE_LOW_QUANT_GOOD  0.02f

#define GAL_MODE_SYMMETRICITY_LOW_QUANT 0.01f

#define GAL_MODE_GOLDEN_RATIO          1.618034f
#define GAL_MODE_TWO_TAKE_GOLDEN_RATIO 0.38197f

#define GAL_MODE_MIRROR_IS_ABOVE_RESULT    (size_t)(-1)

struct gal_mode_params
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
gal_mode_make_mirror_plots(float *sorted, size_t size, size_t mirrorindex,
                           float min, float max, size_t numbins,
                           char *histsname, char *cfpsname,
                           float mirrorplotdist);

float
gal_mode_value_from_sym(float *sorted, size_t size, size_t modeindex,
                        float sym);

void
gal_mode_index_in_sorted(float *sorted, size_t size, float errorstdm,
                         size_t *modeindex, float *modesym);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_MODE_H__ */
