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
#include <gnuastro/statistics.h>

#include "mode.h"






/****************************************************************
 ********                      Sort                       *******
 ****************************************************************/
void
gal_statistics_sort()
{

}



/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
/* Set the bin lower values for all the bins. If the minimum and
   maximum are equal, then use the quantile.

   If the value of onebinvalue is NaN, then nothing will happen,
   however, if it is not a NaN, all the bins will be shifted such that
   the lower values of one of the bins is placed on this value (if it
   is in the range of the data).
*/
void
gal_statistics_set_bins(float *sorted, size_t size, size_t numbins,
                        float min, float max, float onebinvalue,
                        float quant, float **obins)
{
  size_t i;
  float diff, *bins, binwidth;

  /* Allocate space for the array. The extra bin is only for internal
     purposes (so the loops for the histogram and CFP can see the end
     of the last bin). It will never be seen by the user. */
  errno=0;
  bins=*obins=calloc((numbins+1)*2, sizeof *bins);
  if(bins==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for bins in gal_statistics_set_bins "
          "(statistics.c)", (numbins+1)*2*sizeof *bins);

  /* If the range is not defined, find it and set the bin width. */
  if(min==max)
    {
      if(quant!=0.0f)
        {
          min=sorted[ gal_statistics_index_from_quantile(size, quant)   ];
          max=sorted[ gal_statistics_index_from_quantile(size, 1-quant) ];
        }
      else
        {
          min=sorted[0];
          max=sorted[size-1];
        }
    }
  binwidth=(max-min)/numbins;

  /* Set all the bin smaller sides: */
  for(i=0;i<numbins+1;++i)
    bins[i*2]=min+i*binwidth;

  /* Go over all the bins and stop when the sign of the two sides
     of one bin are different. */
  if(isnan(onebinvalue)==0)
    {
      for(i=0;i<numbins;++i)
        if(bins[i*2]<onebinvalue && bins[(i+1)*2]>onebinvalue) break;
      if(i!=numbins)
        {
          diff=onebinvalue-bins[i*2];
          for(i=0;i<numbins+1;++i)
            bins[i*2]+=diff;
        }
    }

  /* For a check
  for(i=0;i<numbins;++i)
    printf("%zu: %.4f\n", i+1, bins[i*2]);
  */
}





void
gal_statistics_histogram(float *sorted, size_t size, float *bins,
                         size_t numbins, int normhist, int maxhistone)
{
  float max=-FLT_MAX;
  size_t histrow=0, i;

  if((long)numbins<=0)
    error(EXIT_FAILURE, 0, "the number of bins in gal_statistics_histogram "
          "(statistics.h) must be >0.  You have given asked for %ld",
          (long)numbins);

  /* Fill the histogram. */
  for(i=0;i<size;++i)
    {
      /* For a check:
      printf("%f:\n  histrow: %zu, numbins: %zu\n",
             sorted[i], histrow, numbins);
      */
      /* Data has not reached bins yet: */
      if(sorted[i]<bins[histrow*2]) continue;

      /* Jump bins until sorted[i] is smaller than the larger value of
         one bin. If we are on the last bin (where
         histrow==numbins-1), then you don't have to increase histrow
         any more, as soon as sorted[i] becomes larger than the lagest
         bin, quit the search. The 1e-6f is to account for floating
         point error. NOTE: the last interval is closed on both
         sides.*/
      if(histrow==numbins-1)
        { if(sorted[i]>bins[numbins*2]+1e-6f) break; }
      else
        while(sorted[i]>=bins[(histrow+1)*2])
          if(++histrow==numbins-1) break;

      /* For a check:
      printf("  histrow: %zu\n", histrow);
      */
      ++bins[histrow*2+1];
    }

  /* In case a normalized histogram is desired: */
  if(normhist)
    for(i=0;i<numbins;++i)
      bins[i*2+1]/=size;

  /* In case the maximum value is to become one. */
  if(maxhistone)
    {
      for(i=0;i<numbins;++i)
        if(bins[i*2+1]>max)
          max=bins[i*2+1];
      for(i=0;i<numbins;++i)
        bins[i*2+1]/=max;
    }

  /* In case you want to see the histogram:
  for(i=0;i<numbins;++i)
    printf("%zu: %.4f %.4F\n", i+1, bins[i*2], bins[i*2+1]);
  */
}





void
gal_statistics_cumulative_fp(float *sorted, size_t size, float *bins,
                             size_t numbins, int normcfp)
{
  float prevind=0;
  size_t cfprow=0, i, numinds=0;


  /* Fill the Cumulative frequency plot. The steps are just like the
     histogram above so we won't explain similar concepts here
     again. */
  for(i=0;i<size;++i)
    {
      if(sorted[i]<bins[cfprow*2]) continue;

      if(cfprow==numbins-1)
        { if(sorted[i]>bins[numbins*2]+1e-6f) break; }
      else
        while(sorted[i]>=bins[(cfprow+1)*2])
          {
            /* We have not yet reached the last bin. But we have
               reached the sorted[i] that is larger than the current
               bin and we want to switch to the next bin. So we have
               to finalize the current bin value.  bins[cfprow*2+1]
               has kept the sum of all the indexs that lie within this
               bin. Using numinds, we also kept a count of how many
               indexs there was.

               In case there were no indexs in this bin, then we have
               to pass on the previous value (since this bin did not
               contribute any data). We set numinds=0 and start
               working on the next bin. */
            if(numinds>0)
              prevind=bins[cfprow*2+1]/=numinds; /* Divide by num indexs */
            else
              bins[cfprow*2+1]=prevind;
            numinds=0;
            if(++cfprow==numbins-1) break;
          }

      /* We can't plot every index, also in every bin, there might be
         a sharp gradient. So in order to make things smooth, the
         value given to each bin will be the average of all the indexs
         that like over that bin. Here we keep the sum and number of
         indexs so in the end we can find the average. */
      bins[cfprow*2+1]+=i;
      ++numinds;
    }


  /* The last bin was not finalized. So we will do that here: */
  bins[cfprow*2+1] = numinds>0 ? bins[cfprow*2+1]/numinds : prevind;


  /* For a normalized CFP: */
  if(normcfp)
    for(i=0;i<numbins;++i)
      bins[i*2+1]/=size;


  /* In a CFP, all bins that are possibly left behind (there was no
     data to fill them) should get the same value as the lastly filled
     value. They should not be zero. Note that this has to be done
     after the possible normalization. */
  for(i=cfprow;i<numbins;++i)
    bins[i*2+1]=bins[(cfprow-1)*2+1];


  /* In case you want to see the CFP:
  for(i=0;i<numbins;++i)
    printf("%zu: %.4f %.4F\n", i+1, bins[i*2], bins[i*2+1]);
  */
}





void
gal_statistics_save_hist(float *sorted, size_t size, size_t numbins,
                         char *filename, char *comment)
{
  FILE *fp;
  size_t i;
  float onebinvalue=NAN;
  int normhist=0, maxhistone=0;
  float d, *bins, min=0.0f, max=0.0f, quant=0.0f;

  /* Set the bin sides: */
  gal_statistics_set_bins(sorted, size, numbins, min, max,
                          onebinvalue, quant, &bins);

  /* Set the size of half a bin width:*/
  d=(bins[2]-bins[0])/2;

  /* Fill the histogram: */
  gal_statistics_histogram(sorted, size, bins, numbins, normhist, maxhistone);

  /* Open the file for writing and save the histogram: */
  errno=0;
  fp=fopen(filename, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "couldn't open file %s", filename);
  fprintf(fp, "%s\n", comment);
  fprintf(fp, "# The input %zu points binned in %zu bins\n#\n",
          size, numbins);
  fprintf(fp, "# Column 0: Value in the middle of this bin.\n");
  fprintf(fp, "# Column 1: Number of points in this bin.\n");
  for(i=0;i<numbins;++i)
    fprintf(fp, "%-15.6f%.0f\n", bins[i*2]+d, bins[i*2+1]);
  fclose(fp);
  free(bins);
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
