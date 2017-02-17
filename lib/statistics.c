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
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/qsort.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>

#include <checkset.h>

#include "mode.h"










/****************************************************************
 ********               Simple statistics                 *******
 ****        (wrappers for functions in `arithmetic.h')      ****
 ****************************************************************/
#define SIMP_FLAGS GAL_ARITHMETIC_NUMOK

gal_data_t *
gal_statistics_number(gal_data_t *data)
{
  return gal_arithmetic(GAL_ARITHMETIC_OP_NUMVAL, SIMP_FLAGS, data);
}

gal_data_t *
gal_statistics_minimum(gal_data_t *data)
{
  return gal_arithmetic(GAL_ARITHMETIC_OP_MINVAL, SIMP_FLAGS, data);
}

gal_data_t *
gal_statistics_maximum(gal_data_t *data)
{
  return gal_arithmetic(GAL_ARITHMETIC_OP_MAXVAL, SIMP_FLAGS, data);
}

gal_data_t *
gal_statistics_sum(gal_data_t *data)
{
  return gal_arithmetic(GAL_ARITHMETIC_OP_SUMVAL, SIMP_FLAGS, data);
}

gal_data_t *
gal_statistics_mean(gal_data_t *data)
{
  return gal_arithmetic(GAL_ARITHMETIC_OP_MEANVAL, SIMP_FLAGS, data);
}

gal_data_t *
gal_statistics_std(gal_data_t *data)
{
  return gal_arithmetic(GAL_ARITHMETIC_OP_STDVAL, SIMP_FLAGS, data);
}

gal_data_t *
gal_statistics_median(gal_data_t *data)
{
  return gal_arithmetic(GAL_ARITHMETIC_OP_MEDIANVAL, SIMP_FLAGS, data);
}






















/****************************************************************
 ********                      Sort                       *******
 ****************************************************************/
/* Check if the given dataset is sorted. Output values are:

     - 0: Dataset is not sorted.
     - 1: Dataset is sorted and increasing or equal.
     - 2: dataset is sorted and decreasing.                  */

#define IS_SORTED(IT) {                                                 \
    IT *aa=data->array, *a=data->array, *af=a+data->size-1;             \
    if(a[1]>=a[0]) do if( *(a+1) < *a ) break; while(++a<af);           \
    else           do if( *(a+1) > *a ) break; while(++a<af);           \
    return ( a==af                                                      \
             ? ( aa[1]>=aa[0]                                           \
                 ? GAL_STATISTICS_SORTED_INCREASING                     \
                 : GAL_STATISTICS_SORTED_DECREASING )                   \
             : GAL_STATISTICS_SORTED_NOT );                             \
  }

int
gal_statistics_is_sorted(gal_data_t *data)
{
  /* A one-element dataset can be considered, sorted, so we'll just return
     1 (for sorted and increasing). */
  if(data->size==1) return GAL_STATISTICS_SORTED_INCREASING;

  /* Do the check. */
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:     IS_SORTED( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      IS_SORTED( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    IS_SORTED( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     IS_SORTED( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    IS_SORTED( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     IS_SORTED( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    IS_SORTED( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     IS_SORTED( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   IS_SORTED( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   IS_SORTED( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_is_sorted'", data->type);
    }

  /* Control shouldn't reach this point. */
  error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can fix the "
        "problem. For some reason, control has reached the end of "
        "`gal_statistics_is_sorted'", PACKAGE_BUGREPORT);
  return -1;
}





/* This function is ignorant to blank values, if you want to make sure
   there is no blank values, you can call `gal_blank_remove' first. */
#define STATISTICS_SORT(QSORT_F) {                                        \
    qsort(data->array, data->size, gal_data_sizeof(data->type), QSORT_F); \
  }
void
gal_statistics_sort_increasing(gal_data_t *data)
{
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:
      STATISTICS_SORT(gal_qsort_uint8_increasing);    break;
    case GAL_DATA_TYPE_INT8:
      STATISTICS_SORT(gal_qsort_int8_increasing);     break;
    case GAL_DATA_TYPE_UINT16:
      STATISTICS_SORT(gal_qsort_uint16_increasing);   break;
    case GAL_DATA_TYPE_INT16:
      STATISTICS_SORT(gal_qsort_int16_increasing);    break;
    case GAL_DATA_TYPE_UINT32:
      STATISTICS_SORT(gal_qsort_uint32_increasing);   break;
    case GAL_DATA_TYPE_INT32:
      STATISTICS_SORT(gal_qsort_int32_increasing);    break;
    case GAL_DATA_TYPE_UINT64:
      STATISTICS_SORT(gal_qsort_uint64_increasing);   break;
    case GAL_DATA_TYPE_INT64:
      STATISTICS_SORT(gal_qsort_int64_increasing);    break;
    case GAL_DATA_TYPE_FLOAT32:
      STATISTICS_SORT(gal_qsort_float32_increasing);  break;
    case GAL_DATA_TYPE_FLOAT64:
      STATISTICS_SORT(gal_qsort_float64_increasing);  break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_sort_increasing'", data->type);
    }
}





/* See explanations above `gal_statistics_sort_increasing'. */
void
gal_statistics_sort_decreasing(gal_data_t *data)
{
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:
      STATISTICS_SORT(gal_qsort_uint8_decreasing);    break;
    case GAL_DATA_TYPE_INT8:
      STATISTICS_SORT(gal_qsort_int8_decreasing);     break;
    case GAL_DATA_TYPE_UINT16:
      STATISTICS_SORT(gal_qsort_uint16_decreasing);   break;
    case GAL_DATA_TYPE_INT16:
      STATISTICS_SORT(gal_qsort_int16_decreasing);    break;
    case GAL_DATA_TYPE_UINT32:
      STATISTICS_SORT(gal_qsort_uint32_decreasing);   break;
    case GAL_DATA_TYPE_INT32:
      STATISTICS_SORT(gal_qsort_int32_decreasing);    break;
    case GAL_DATA_TYPE_UINT64:
      STATISTICS_SORT(gal_qsort_uint64_decreasing);   break;
    case GAL_DATA_TYPE_INT64:
      STATISTICS_SORT(gal_qsort_int64_decreasing);    break;
    case GAL_DATA_TYPE_FLOAT32:
      STATISTICS_SORT(gal_qsort_float32_decreasing);  break;
    case GAL_DATA_TYPE_FLOAT64:
      STATISTICS_SORT(gal_qsort_float64_decreasing);  break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_sort_decreasing'", data->type);
    }
}




















/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
/* Generate an array of regularly spaced elements. For a 1D dataset, the
   output will be 1D, for 2D, it will be 2D.

   Input arguments:

     * The `data' set you want to apply the bins to. This is only necessary
       if the range argument is not complete, see below. If `range' has all
       the necessary information, you can pass a NULL pointer for `data'.

     * The `range' data structure keeps the desired range along each
       dimension of the input data structure, it has to be in float32
       type. Note that if

         - If you want the full range of the dataset (in any dimensions,
           then just set `range' to NULL and the range will be specified
           from the minimum and maximum value of the dataset.

         - If there is one element for each dimension in range, then it is
           viewed as a quantile (Q), and the range will be: `Q to 1-Q'.

         - If there are two elements for each dimension in range, then they
           are assumed to be your desired minimum and maximum values. When
           either of the two are NaN, the minimum and maximum will be
           calculated for it.

     * The number of bins: must be larger than 0.

     * `onebinstart' A desired value for onebinstart. Note that with this
        option, the bins won't start and end exactly on the given range
        values, it will be slightly shifted to accommodate this
        request.  */
gal_data_t *
gal_statistics_regular_bins(gal_data_t *data, gal_data_t *range,
                            size_t numbins, float onebinstart)
{
  size_t i;
  gal_data_t *bins, *tmp;
  float *b, *ra, min, max, hbw, diff, binwidth;


  /* Some sanity checks. */
  if(numbins==0)
    error(EXIT_FAILURE, 0, "`numbins' in `gal_statistics_regular_bins' "
          "cannot be given a value of 0");
  if(range && range->type!=GAL_DATA_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "gal_statistics_regular_bins currently on works "
          "on ranges of type float32");
  if(data->ndim!=1)
    error(EXIT_FAILURE, 0, "gal_statistics_regular_bins currently on works "
          "in 1D data");
  if(data->type!=GAL_DATA_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "gal_statistics_regular_bins currently on works "
          "on float32 type data");


  /* Set the minimum and maximum values. */
  if(range && range->size)
    {
      ra=range->array;
      if( (range->size)%2 )
        error(EXIT_FAILURE, 0, "Quantile ranges are not implemented in "
              "`gal_statistics_regular_bins' yet.");
      else
        {
          if( isnan(ra[0]) )
            {
              tmp=gal_statistics_minimum(data);
              min=*((float *)(tmp->array));
              gal_data_free(tmp);
            }
          else min=ra[0];
          if( isnan(ra[1]) )                       /* When the maximum    */
            {                                      /* isn't set, we'll    */
              tmp=gal_statistics_maximum(data);    /* Add a very small    */
              max=*((float *)(tmp->array)) + 1e-5; /* value so all the    */
              gal_data_free(tmp);                  /* points are counted. */
            }
          else max=ra[1];
        }
    }
  else
    {
      tmp=gal_statistics_minimum(data);
      min=*((float *)(tmp->array));
      gal_data_free(tmp);
      tmp=gal_statistics_maximum(data);
      max=*((float *)(tmp->array));
      gal_data_free(tmp);
    }


  /* Allocate the space for the bins. */
  bins=gal_data_alloc(NULL, GAL_DATA_TYPE_FLOAT32, 1, &numbins, NULL,
                      0, data->minmapsize, "bin_center", data->unit,
                      "Center value of each bin.");


  /* Set central bin values. */
  b=bins->array;
  hbw = ( binwidth=(max-min)/numbins )/2;
  for(i=0;i<numbins;++i) b[i] = min + i*binwidth + hbw;


  /* Go over all the bins and stop when the sign of the two sides
     of one bin are different. */
  if( !isnan(onebinstart) )
    {
      for(i=0;i<numbins-1;++i)
        if( (b[i]-hbw) < onebinstart && (b[i+1]-hbw) > onebinstart) break;
      if( i != numbins-1 )
        {
          diff=onebinstart-b[i];
          for(i=0;i<numbins;++i)
            b[i]+=diff;
        }
    }

  /* For a check
  printf("min: %g\n", min);
  printf("max: %g\n", max);
  printf("binwidth: %g\n", binwidth);
  for(i=0;i<numbins;++i)
    printf("%zu: %.4f\n", i, b[i]);
  */

  /* Set the status of the bins to regular and return. */
  bins->status=GAL_STATISTICS_BINS_REGULAR;
  return bins;
}





/* Make a histogram of all the elements in the given dataset with bin
   values that are defined in the `bins' structure (see
   `gal_statistics_regular_bins'). */
gal_data_t *
gal_statistics_histogram(gal_data_t *data, gal_data_t *bins,
                         int normalize, int maxhistone)
{
  size_t i, *h;
  double ref=NAN;
  gal_data_t *hist;
  float *f, *ff, min, max, binwidth;


  /* Some basic sanity checks for now (until it is generalized). */
  if(data->ndim!=1)
    error(EXIT_FAILURE, 0, "gal_statistics_histogram currently on works "
          "in 1D data");
  if(data->type!=GAL_DATA_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "gal_statistics_histogram currently on works "
          "on float32 type data");


  /* Check if the bins are regular or not. For irregular bins, we can
     either use the old implementation, or GSL's histogram
     functionality. */
  if(bins->status!=GAL_STATISTICS_BINS_REGULAR)
    error(EXIT_FAILURE, 0, "the input bins to `gal_statistics_histogram' "
          "are not regular. Currently it is only implemented for regular "
          "bins");


  /* Check if normalize and maxhistone are not called together. */
  if(normalize && maxhistone)
    error(EXIT_FAILURE, 0, "only one of `normalize' and `maxhistone' may "
          "be given to `gal_statistics_histogram'");

  /* Allocate the histogram (note that we are clearning it. */
  hist=gal_data_alloc(NULL, GAL_DATA_TYPE_SIZE_T, bins->ndim, bins->dsize,
                      NULL, 1, data->minmapsize, "hist_number", "counts",
                      "Number of data points within each bin.");


  /* Set the minimum and maximum values: */
  f=bins->array;
  binwidth=f[1]-f[0];
  max = f[bins->size - 1] + binwidth/2;
  min = f[0]              - binwidth/2;


  /* Go through all the elements and find out which bin they belong to. */
  h=hist->array;
  f=data->array;
  for(i=0;i<data->size;++i)
    if(f[i]>=min && f[i]<max)
      {
        ++h[ (size_t)((f[i]-min)/binwidth) ];
        /*printf("%-10.3f%zu\n", f[i], (size_t)((f[i]-min)/binwidth) );*/
      }


  /* Find the reference to correct the histogram if necessary. */
  if(normalize)
    {
      /* Set the reference. */
      ref=0.0f;
      hist=gal_data_copy_to_new_type_free(hist, GAL_DATA_TYPE_FLOAT32);
      ff=(f=hist->array)+hist->size; do ref += *f++;   while(f<ff);

      /* Correct the name, units and comments. */
      free(hist->name); free(hist->unit); free(hist->comment);
      gal_checkset_allocate_copy("hist_normalized", &hist->name);
      gal_checkset_allocate_copy("frac", &hist->unit);
      gal_checkset_allocate_copy("Normalized histogram value for this bin",
                                 &hist->comment);
    }
  if(maxhistone)
    {
      /* Calculate the reference. */
      ref=-FLT_MAX;
      hist=gal_data_copy_to_new_type_free(hist, GAL_DATA_TYPE_FLOAT32);
      ff=(f=hist->array)+hist->size;
      do ref = *f>ref ? *f : ref; while(++f<ff);

      /* Correct the name, units and comments. */
      free(hist->name); free(hist->unit); free(hist->comment);
      gal_checkset_allocate_copy("hist_maxone", &hist->name);
      gal_checkset_allocate_copy("frac", &hist->unit);
      gal_checkset_allocate_copy("Fractional histogram value for this bin "
                                 "when maximum bin value is 1.0",
                                 &hist->comment);
    }


  /* Correct the histogram if necessary. */
  if( !isnan(ref) )
    { ff=(f=hist->array)+hist->size; do *f++ /= ref;   while(f<ff); }

  /* Return the histogram. */
  return hist;
}





/* Make a cumulative frequency plot (CFP) of all the elements in the given
   dataset with bin values that are defined in the `bins' structure (see
   `gal_statistics_regular_bins').

   The CFP is built from the histogram: in each bin, the value is the sum
   of all previous bins in the histogram. Thus, if you have already
   calculated the histogram before calling this function, you can pass it
   onto this function as the data structure in `bins->next'. If
   `bin->next!=NULL', then it is assumed to be the histogram. If it is
   NULL, then the histogram will be calculated internally and freed after
   the job is finished.

   When a histogram is given and it is normalized, the CFP will also be
   normalized (even if the normalized flag is not set here): note that a
   normalized CFP's maximum value is 1.
*/
gal_data_t *
gal_statistics_cfp(gal_data_t *data, gal_data_t *bins, int normalize)
{
  double sum;
  float *f, *ff, *hf;
  gal_data_t *hist, *cfp;
  size_t *s, *sf, *hs, sums;


  /* Some basic sanity checks for now (until it is generalized). */
  if(data->ndim!=1)
    error(EXIT_FAILURE, 0, "gal_statistics_cfp currently on works "
          "in 1D data");
  if(data->type!=GAL_DATA_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "gal_statistics_cfp currently on works "
          "on float32 type data");


  /* Check if the bins are regular or not. For irregular bins, we can
     either use the old implementation, or GSL's histogram
     functionality. */
  if(bins->status!=GAL_STATISTICS_BINS_REGULAR)
    error(EXIT_FAILURE, 0, "the input bins to `gal_statistics_cfp' "
          "are not regular. Currently it is only implemented for regular "
          "bins");


  /* Prepare the histogram. */
  hist = ( bins->next
           ? bins->next
           : gal_statistics_histogram(data, bins, 0, 0) );


  /* If the histogram has float32 type it was given by the user and is
     either normalized or its maximum was set to 1. We can only use it if
     it was normalized. If its maximum is 1, then we must ignore it and
     build the histogram again.*/
  if(hist->type==GAL_DATA_TYPE_FLOAT32)
    {
      sum=0.0f;
      ff=(f=hist->array)+hist->size; do sum += *f++;   while(f<ff);
      if(sum!=1.0f)
        hist=gal_statistics_histogram(data, bins, 0, 0);
    }

  /* Allocate the cumulative frequency plot's necessary space. */
  cfp=gal_data_alloc( NULL, hist->type, bins->ndim, bins->dsize,
                      NULL, 1, data->minmapsize,
                      ( hist->type==GAL_DATA_TYPE_FLOAT32
                        ? "cfp_normalized" : "cfp_number" ),
                      ( hist->type==GAL_DATA_TYPE_FLOAT32
                        ? "frac" : "count" ),
                      ( hist->type==GAL_DATA_TYPE_FLOAT32
                        ? "Fraction of data elements from the start to this "
                        "bin (inclusive)."
                        : "Number of data elements from the start to this "
                        "bin (inclusive).") );


  /* Fill in the cumulative frequency plot. */
  switch(hist->type)
    {
    case GAL_DATA_TYPE_SIZE_T:
      sums=0; hs=hist->array; sf=(s=cfp->array)+cfp->size;
      do sums = (*s += *hs++ + sums); while(++s<sf);
      break;

    case GAL_DATA_TYPE_FLOAT32:
      sum=0.0f; hf=hist->array; ff=(f=cfp->array)+cfp->size;
      do sum = (*f += *hf++ + sum);  while(++f<ff);
      break;

    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_statistics_cfp'", cfp->type);
    }


  /* Normalize the CFP if the user asked for it and it wasn't normalized
     until now. */
  if(normalize && cfp->type==GAL_DATA_TYPE_SIZE_T)
    {
      /* Find the sum, then divide the plot by it. Note that the sum must
         come from the histogram, not the CFP!*/
      sums=0;
      cfp=gal_data_copy_to_new_type_free(cfp, GAL_DATA_TYPE_FLOAT32);
      sf=(s=hist->array)+hist->size; do sums += *s++;   while(s<sf);
      ff=(f=cfp->array)+cfp->size;   do *f++ /= sums;   while(f<ff);

      /* Correct the name, units and comments. */
      free(cfp->name); free(cfp->unit); free(cfp->comment);
      gal_checkset_allocate_copy("cfp_normalized", &cfp->name);
      gal_checkset_allocate_copy("frac", &cfp->unit);
      gal_checkset_allocate_copy("Fraction of data elements from the start "
                                 "to this bin (inclusive).", &cfp->comment);
    }

  /* If the histogram was allocated here, free it. */
  if(hist!=bins->next) gal_data_free(hist);
  return cfp;
}




















/****************************************************************
 *****************         Quantiles         ********************
 ****************************************************************/
/* Find the index corresponding to a certain quantile, considering the
   rounding that might be needed. */
size_t
gal_statistics_index_from_quantile(size_t size, float quant)
{
  float floatindex;

  if(quant>1.0f)
    error(EXIT_FAILURE, 0, "the quantile in gal_statistics_index_from_quantile "
          "(statistics.c) Should be smaller");

  /* Find the index of the quantile. */
  floatindex=(float)size*quant;

  /*
  printf("quant: %f, size: %zu, findex: %f\n", quant, size, floatindex);
  */
  /* Note that in the conversion from float to size_t, the floor
     integer value of the float will be used. */
  if( floatindex - (int)floatindex > 0.5 )
    return floatindex+1;
  else
    return floatindex;
}




















/****************************************************************
 *****************       Sigma clip          ********************
 ****************************************************************/
/* This function will repeatedly sigma clip the data and return the
   median. It is assumed that the data have been ordered.

   o1_n0: Ordered (1), not ordered (0).
*/
int
gal_statistics_sigma_clip_converge(float *array, int o1_n0, size_t num_elem,
                                   float sigma_multiple, float accuracy,
                                   float *outave, float *outmed,
                                   float *outstd, int print)
{
  printf("\n ... in gal_statistics_sigma_clip_converge ... \n");
  exit(1);
#if 0
  size_t counter=0;
  float *start, *oldstart, *dpt;
  float ave=*outave, med=*outmed, std=*outstd;
  float oldstd=0, oldmed=0, oldave=0, *orderedarray;

  if(o1_n0==0)
    {
      gal_array_float_copy(array, num_elem, &orderedarray);
      qsort(orderedarray, num_elem, sizeof*orderedarray,
            gal_qsort_float32_increasing);
    }
  else orderedarray=array;

  start=orderedarray;
  while(counter<GAL_STATISTICS_MAX_SIG_CLIP_CONVERGE)
    {
      oldstart=start;

      med=*(start+num_elem/2);
      gal_statistics_f_ave_std(start, num_elem, &ave, &std, NULL);

      if(print)
        printf("      %zu: %f  %f  %f  %zu\n",
               counter+1, med, ave, std, num_elem);

      /* It might happen that ave and std are NaN. If so, stop the
         process here (the user has not given a mask and some pixels
         have nan values!). */
      if(isnan(ave) || isnan(std))
        return 0;

      /* Normally, oldstd should be larger than std, because the
         possible outliers have been removed. If it is not, it means
         that we have clipped too much! */
      if(counter>0 && (oldstd-std)/std<accuracy)
        {
          *outstd=oldstd; *outave=oldave; *outmed=oldmed;
          return 1;
        }

      for(dpt=start; dpt<start+num_elem; ++dpt)
        if (*dpt>med-sigma_multiple*std)
          {
            start=dpt;
            break;
          }

      for(dpt=oldstart+num_elem-1;dpt>start;dpt--)
        if (*dpt<med+sigma_multiple*std)
          {
            num_elem=dpt-start+1;
            break;
          }

      oldave=ave;
      oldmed=med;
      oldstd=std;
      ++counter;
    }
#endif
  return 0;
}





/* This function will do a certain number of sigma clips and return
   the final average, median and std. o1_n0: 1: initially ordered. 2:
   initially not ordered.*/
int
gal_statistics_sigma_clip_certain_num(float *array, int o1_n0, size_t num_elem,
                                      float sigma_multiple, size_t numtimes,
                                      float *outave, float *outmed,
                                      float *outstd, int print)
{
  printf("\n ... in gal_statistics_sigma_clip_certain_num ... \n");
  exit(1);
#if 0
  size_t counter=0;
  float ave=*outave, med=*outmed, std=*outstd;
  float *start, *oldstart, *dpt, *orderedarray;

  if(o1_n0==0)
    {
      gal_array_float_copy(array, num_elem, &orderedarray);
      qsort(orderedarray, num_elem, sizeof*orderedarray,
            gal_qsort_float32_increasing);
    }
  else orderedarray=array;

  start=orderedarray;
  for(counter=0;counter<numtimes;++counter)
    {
      oldstart=start;

      med=*(start+num_elem/2);
      gal_statistics_f_ave_std(start, num_elem, &ave, &std, NULL);

      if(print)
        printf("      %zu: %f  %f  %f  %zu\n",
               counter+1, med, ave, std, num_elem);

      /* It might happen that ave and std are nan. If so, stop the
         process here (the user has not given a mask and some pixels
         have nan values!). */
      if(isnan(ave) || isnan(std))
        return 0;


      for(dpt=start; dpt<start+num_elem; ++dpt)
        if (*dpt>med-sigma_multiple*std)
          {
            start=dpt;
            break;
          }


      for(dpt=oldstart+num_elem-1;dpt>start;dpt--)
        if (*dpt<med+sigma_multiple*std)
          {
            num_elem=dpt-start+1;
            break;
          }
    }

  if(o1_n0==0)
    free(orderedarray);

  *outave=ave;
  *outmed=med;
  *outstd=std;
#endif
  return 1;
}



















/****************************************************************/
/*************         Identify outliers         ****************/
/****************************************************************/
/* Using the cumulative distribution function this funciton will
   remove outliers from a dataset. */
void
gal_statistics_remove_outliers_flat_cdf(float *sorted, size_t *outsize)
{
  printf("\n ... in gal_statistics_remove_outliers_flat_cdf ... \n");
  exit(1);
#if 0
  int firstfound=0;
  size_t size=*outsize, i, maxind;
  float *slopes, minslope, maxslope;

  /* Find a slopes array, think of the cumulative frequency plot when
     you want to think about slopes. */
  errno=0; slopes=malloc(size*sizeof *slopes);
  if(slopes==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for slopes in "
          "gal_statistics_remove_outliers_flat_cdf (statistics.c)",
          size*sizeof *slopes);

  /* Calcuate the slope of the CDF and put it in the slopes array. */
  for(i=1;i<size-1;++i)
    slopes[i]=2/(sorted[i+1]-sorted[i-1]);

  /* Find the position of the maximum slope, note that around the
     distribution mode, the difference between the values varies less,
     so two neighbouring elements have the closest values, hence the
     largest slope (when their difference is in the denominator). */
  gal_statistics_f_max_with_index(slopes+1, size-2, &maxslope, &maxind);

  /* Find the minimum slope from the second element (for the first the
     slope is not defined. NOTE; maxind is one smaller than it should
     be because the input array to find it began from the second
     element. */
  gal_statistics_float_second_min(slopes+1, maxind+1, &minslope);

  /* Find the second place where the slope falls below `minslope`
     after the maximum position. When found, add it with one to
     account for error. Note that incase there are no outliers by this
     definition, then the final size will be equal to size! */
  for(i=maxind+1;i<size-1;++i)
    if(slopes[i]<minslope)
      {
        if(firstfound)
          break;
        else
          firstfound=1;
      }
  *outsize=i+1;

  /*
  for(i=0;i<size;++i)
    printf("%zu\t%.3f\t%.3f\n", i, arr[i], slopes[i]);
  printf("\n\nPlace to cut off for outliers is: %zu\n\n", *outsize);
  */

  free(slopes);
#endif
}




















/****************************************************************/
/*************          Mode calculation         ****************/
/****************************************************************/
void
gal_statistics_mode_mirror_plots(float *sorted, size_t size, size_t mirrorindex,
                                 float min, float max, size_t numbins,
                                 char *histsname, char *cfpsname,
                                 float mirrorplotdist)
{
  printf("\n... in gal_statistics_mode_mirror_plots ...\n");
  exit(1);
#if 0
  FILE *fp;
  size_t i, msize;
  float *out, maxhist=-FLT_MAX, maxcfp, d;
  int normhist=0, maxhistone=0, normcfp=0;
  float *bins, *mirror, *actual, mf, onebinvalue=0.0f;


  /* Find the index of the mirror and make the mirror array: */
  mf=sorted[mirrorindex];
  gal_array_float_copy(sorted, size, &actual);
  makemirrored(sorted, mirrorindex, &mirror, &msize);


  /* Set the best range if asked for, such that the mirror is on the
     1/3 of the image scale. */
  if(mirrorplotdist!=0.0f)
    {
      min=actual[gal_statistics_index_from_quantile(size, 0.001f)];
      max=mf+mirrorplotdist*(mf-min);
    }


  /* set the mirror pixel to have the value of zero.*/
  min-=mf;
  max-=mf;
  gal_array_fsum_const(actual, size, -1.0f*mf);
  gal_array_fsum_const(mirror, msize, -1.0f*mf);


  /* Allocate space for the 3 column array keeping the histograms:*/
  errno=0;
  out=malloc(numbins*3*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for out in "
          "gal_mode_make_mirror_plots (mode.c)", numbins*3*sizeof *out);


  /* Define the bin sides: */
  gal_statistics_set_bins(actual, size, numbins, min,
                          max, onebinvalue, 0, &bins);


  /* Find the histogram of the actual data and put it in out. Note
     that maxhistone=0, because here we want to use one value for both
     histograms so they are comparable. */
  gal_statistics_histogram(actual, size, bins, numbins,
                           normhist, maxhistone);
  for(i=0;i<numbins;++i)
    if(bins[i*2+1]>maxhist) maxhist=bins[i*2+1];
  for(i=0;i<numbins;++i)
    { out[i*3]=bins[i*2]; out[i*3+1]=bins[i*2+1]/maxhist; bins[i*2+1]=0.0f;}
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */
  d=(out[3]-out[0])/2;


  /* Find the histogram of the mirrired distribution and put it in
     out: */
  gal_statistics_histogram(mirror, msize, bins, numbins, normhist,
                           maxhistone);
  for(i=0;i<numbins;++i)
    { out[i*3+2]=bins[i*2+1]/maxhist; bins[i*2+1]=0.0f;}
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */


  /* Print out the histogram: */
  errno=0;
  fp=fopen(histsname, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "could not open file %s", histsname);
  fprintf(fp, "# Histogram of actual and mirrored distributions.\n");
  fprintf(fp, "# Column 0: Value in the middle of this bin.\n");
  fprintf(fp, "# Column 1: Input data.\n");
  fprintf(fp, "# Column 2: Mirror distribution.\n");
  for(i=0;i<numbins;++i)
    fprintf(fp, "%-25.6f%-25.6f%-25.6f\n", out[i*3]+d,
            out[i*3+1], out[i*3+2]);
  fclose(fp);




  /* Find the cumulative frequency plots of the two distributions: */
  gal_statistics_cumulative_fp(actual, size, bins, numbins, normcfp);
  for(i=0;i<numbins;++i)
    { out[i*3+1]=bins[i*2+1]; bins[i*2+1]=0.0f; }
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */
  gal_statistics_cumulative_fp(mirror, msize, bins, numbins, normcfp);
  for(i=0;i<numbins;++i)
    { out[i*3+2]=bins[i*2+1]; bins[i*2+1]=0.0f;}
  bins[i*2+1]=0.0f; /* bins[] actually has numbins+1 elements. */


  /* Since the two cumultiave frequency plots have to be on scale, see
     which one is larger and divided both CFPs by the size of the
     larger one. Then print the CFPs. */
  if(size>msize) maxcfp=size;
  else maxcfp=msize;
  errno=0;
  fp=fopen(cfpsname, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "could not open file %s", cfpsname);
  fprintf(fp, "# Cumulative frequency plot (average index in bin) of\n"
          "# Actual and mirrored distributions.\n");
  fprintf(fp, "# Column 0: Value in the middle of this bin.\n");
  fprintf(fp, "# Column 1: Actual data.\n");
  fprintf(fp, "# Column 2: Mirror distribution.\n");
  for(i=0;i<numbins;++i)
    fprintf(fp, "%-25.6f%-25.6f%-25.6f\n", out[i*3],
            out[i*3+1]/maxcfp, out[i*3+2]/maxcfp);
  fclose(fp);


  /* Clean up. */
  free(out);
  free(bins);
  free(mirror);
  free(actual);
#endif
}





/* It happens that you have the symmetricity and you want the flux
   value at that point, this function will do that job. In practice,
   it just finds bf from the equation to calculate symmetricity in
   modesymmetricity. */
float
gal_statistics_mode_value_from_sym(float *sorted, size_t size,
                                   size_t modeindex, float sym)
{
  float mf=sorted[modeindex];
  float af=
    sorted[gal_statistics_index_from_quantile(2*modeindex+1,
                           GAL_STATISTICS_MODE_SYMMETRICITY_LOW_QUANT)];
  return sym*(mf-af)+mf;
}





/* Find the quantile of the mode of a sorted distribution. The return
   value is either 0 (not accurate) or 1 (accurate). Accuracy is
   defined based on the difference between the maximum and minimum
   maxdiffs that were found during the golden section search.

   A good mode will have:

   modequant=(float)(modeindex)/(float)size;
   modesym>GAL_MODE_SYM_GOOD && modequant>GAL_MODE_LOW_QUANT_GOOD
*/
void
gal_statistics_mode_index_in_sorted(float *sorted, size_t size, float errorstdm,
                                    size_t *modeindex, float *modesym)
{
  struct gal_statistics_mode_params mp;

  /* Initialize the gal_mode_params structure: */
  mp.size=size;
  mp.sorted=sorted;
  mp.tolerance=0.01;
  mp.numcheck=size/2;
  mp.errorstdm=errorstdm;
  if(mp.numcheck>1000)
    mp.interval=mp.numcheck/1000;
  else mp.interval=1;
  mp.lowi  = gal_statistics_index_from_quantile(size,
                                 GAL_STATISTICS_MODE_LOW_QUANTILE);
  mp.highi = gal_statistics_index_from_quantile(size,
                                 GAL_STATISTICS_MODE_HIGH_QUANTILE);
  mp.midi  = (((float)mp.highi
               + MODE_GOLDEN_RATIO*(float)mp.lowi)
              /(1+MODE_GOLDEN_RATIO));
  mp.midd  = mirrormaxdiff(mp.sorted, mp.size, mp.midi, mp.numcheck,
                           mp.interval, mp.errorstdm);

  /* Do the golden section search and find the resulting
     symmetricity. */
  *modeindex=modegoldenselection(&mp);
  modesymmetricity(sorted, size, *modeindex, errorstdm, modesym);
}
