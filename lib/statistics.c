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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <stdlib.h>

#include "forqsort.h"
#include "statistics.h"
#include "arraymanip.h"









/****************************************************************
 *****************    Mininum and Maximum    ********************
 ****************************************************************/
void
floatmin(float *in, size_t size, float *min)
{
  float tmin=FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(!isnan(*in) && *in<tmin) tmin=*in;
  while(++in<fpt);
  *min=tmin;
}





void
floatmax(float *in, size_t size, float *max)
{
  float tmax=-FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(!isnan(*in) && *in>tmax) tmax=*in;
  while(++in<fpt);
  *max=tmax;
}





void
doublemin(double *in, size_t size, double *min)
{
  double tmin=FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(!isnan(*in) && *in<tmin) tmin=*in;
  while(++in<fpt);
  *min=tmin;
}





void
doublemax(double *in, size_t size, double *max)
{
  double tmax=-FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(!isnan(*in) && *in>tmax) tmax=*in;
  while(++in<fpt);
  *max=tmax;
}




void
floatmaxmasked(float *in, unsigned char *mask, size_t size, float *max)
{
  float tmax=-FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(*mask++==0 && *in>tmax)
      tmax=*in;
  while(++in<fpt);
  *max=tmax;
}





void
floatsecondmax(float *in, size_t size, float *secondmax)
{
  float smax=-FLT_MAX, max=-FLT_MAX, *fpt;
  fpt=in+size;
  do
    {
      if(!isnan(*in) && *in>max)
	{
	  smax=max;
	  max=*in;
	}
      else if(!isnan(*in) && *in>smax) smax=*in;
    }
  while(++in<fpt);
  *secondmax=smax;
}





void
floatsecondmin(float *in, size_t size, float *secondmin)
{
  float smin=FLT_MAX, min=FLT_MAX, *fpt;
  fpt=in+size;
  do
    {
      if(!isnan(*in) && *in<min)
	{
	  smin=min;
	  min=*in;
	}
      else if(!isnan(*in) && *in<smin) smin=*in;
    }
  while(++in<fpt);
  *secondmin=smin;
}





void
fminmax(float *in, size_t size, float *min, float *max)
{
  float tmin=FLT_MAX, tmax=-FLT_MAX, *f, *fpt;

  fpt=(f=in)+size;
  do
    {
      if (*f>tmax) tmax=*f;
      if (*f<tmin) tmin=*f;
    }
  while(++f<fpt);


  /* In case there was at least one NAN value */
  if(isnan(tmin) || isnan(tmax))
    {
      tmin=FLT_MAX;
      tmax=-FLT_MAX;
      fpt=(f=in)+size;
      do
        if(isnan(*f)==0)
          {
            if (*f>tmax) tmax=*f;
            if (*f<tmin) tmin=*f;
          }
      while(++f<fpt);
    }

  if(tmin==FLT_MAX || tmax==-FLT_MAX)
    *min=*max=NAN;
  else
    {
      *max=tmax;
      *min=tmin;
    }

  *max=tmax;
  *min=tmin;
}





void
dminmax(double *in, size_t size, double *min, double *max)
{
  double tmin=FLT_MAX, tmax=-FLT_MAX, *d, *dpt;

  dpt=(d=in)+size;
  do
    {
      if (*d>tmax) tmax=*d;
      if (*d<tmin) tmin=*d;
    }
  while(++d<dpt);


  /* In case there was at least one NAN value */
  if(isnan(tmin) || isnan(tmax))
    {
      tmin=FLT_MAX;
      tmax=-FLT_MAX;
      dpt=(d=in)+size;
      do
        if(isnan(*d)==0)
          {
            if (*d>tmax) tmax=*d;
            if (*d<tmin) tmin=*d;
          }
      while(++d<dpt);
    }

  if(tmin==FLT_MAX || tmax==-FLT_MAX)
    *min=*max=NAN;
  else
    {
      *max=tmax;
      *min=tmin;
    }

  *max=tmax;
  *min=tmin;
}




void
dmax_withindex(double *in, size_t size, double *max, size_t *index)
{
  size_t tindex=0;
  double *fpt, *pt=in, tmax=-FLT_MAX;

  fpt=pt+size;
  do
    if(!isnan(*pt) && *pt>tmax)
      {
	tmax=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *max=tmax;
}





void
fmax_withindex(float *in, size_t size,
	       float *max, size_t *index)
{
  size_t tindex=0;
  float *pt=in, *fpt, tmax=-FLT_MAX;

  fpt=pt+size;
  do
    if(!isnan(*pt) && *pt>tmax)
      {
	tmax=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *max=tmax;
}





void
dmin_withindex(double *in, size_t size,
	       double *min, size_t *index)
{
  size_t tindex=0;
  double *pt=in, *fpt, tmin=FLT_MAX;

  fpt=pt+size;
  do
    if(!isnan(*pt) && *pt<tmin)
      {
	tmin=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *min=tmin;
}





void
fmin_withindex(float *in, size_t size,
	       float *min, size_t *index)
{
  size_t tindex=0;
  float *pt=in, *fpt, tmin=FLT_MAX;

  fpt=pt+size;
  do
    if(!isnan(*pt) && *pt<tmin)
      {
	tmin=*pt;
	tindex=pt-in;
      }
  while(++pt<fpt);
  *index=tindex;
  *min=tmin;
}




















/****************************************************************
 *****************            Sum            ********************
 ****************************************************************/
float
floatsum(float *in, size_t size)
{
  float *fpt;
  double sum=0;
  fpt=in+size;
  do
    if(!isnan(*in))
      sum+=*in;
  while(++in<fpt);
  return sum;
}





float
floatsumsquared(float *in, size_t size)
{
  float *fpt;
  double sum=0;
  fpt=in+size;
  do
    if(!isnan(*in))
      sum+=*in * *in;
  while(++in<fpt);
  return sum;
}





/* Sum over all elements of the array that are not covered by a
   mask. Any non-zero masked pixel is considered to be a masked
   pixel. */
float
floatsummask(float *in, unsigned char *mask,
	     size_t size, size_t *nsize)
{
  double sum=0;
  size_t counter=0;
  unsigned char *pt=mask, *fpt;

  fpt=pt+size;
  do
    if(*pt==0)
      {
	sum+=in[pt-mask];
	++counter;
      }
  while(++pt<fpt);

  *nsize=counter;
  return sum;
}





float
floatsummaskl(float *in, long *mask,
              size_t size, size_t *nsize)
{
  double sum=0;
  size_t counter=0;
  long *pt=mask, *fpt;

  fpt=pt+size;
  do
    if(*pt==0)
      {
	sum+=in[pt-mask];
	++counter;
      }
  while(++pt<fpt);

  *nsize=counter;
  return sum;
}





float
floatsumsquaredmask(float *in, unsigned char *mask,
                    size_t size, size_t *nsize)
{
  double sum=0;
  size_t counter=0;
  unsigned char *pt=mask, *fpt;

  fpt=pt+size;
  do
    if(*pt==0)
      {
	sum+=in[pt-mask] * in[pt-mask];
	++counter;
      }
  while(++pt<fpt);

  *nsize=counter;
  return sum;
}






float
floatsumsquaredmaskl(float *in, long *mask,
                     size_t size, size_t *nsize)
{
  double sum=0;
  size_t counter=0;
  long *pt=mask, *fpt;

  fpt=pt+size;
  do
    if(*pt==0)
      {
	sum+=in[pt-mask] * in[pt-mask];
	++counter;
      }
  while(++pt<fpt);

  *nsize=counter;
  return sum;
}





















/****************************************************************
 *****************      Average and          ********************
 ****************    Standard deviation      ********************
 ****************************************************************/
void
fave(float *in, size_t size, float *ave, unsigned char *mask)
{
  float sum;
  size_t nsize;
  if(mask==NULL)
    sum=floatsum(in, size);
  else
    {
      sum=floatsummask(in, mask, size, &nsize);
      size=nsize;
    }
  *ave=sum/size;
}





void
favel(float *in, size_t size, float *ave, long *mask)
{
  float sum;
  size_t nsize;
  if(mask==NULL)
    sum=floatsum(in, size);
  else
    {
      sum=floatsummaskl(in, mask, size, &nsize);
      size=nsize;
    }
  *ave=sum/size;
}





/* Find the average and standard deviation of an array, assuming that
   there is a mask array. Any mask array pixel that is not zero will
   not be included in the average and standard deviations.  Here the
   mask is assumed to be unsigned char.  */
void
favestd(float *in, size_t size, float *ave, float *std, unsigned char *mask)
{
  size_t nsize1, nsize2;
  float sum, sum2;
  if(mask)
    {
      sum=floatsummask(in, mask, size, &nsize1);
      sum2=floatsumsquaredmask(in, mask, size, &nsize2);
      if(nsize1!=nsize2)
	error(EXIT_FAILURE, 0, "A bug in favestd (lib/statistics.h). "
	      "Somehow the number of masked pixels is measured "
	      "differently. Please contact us so we can find the cause.");
      size=nsize1;
    }
  else
    {
      sum=floatsum(in, size);
      sum2=floatsumsquared(in, size);
    }
  *ave=sum/size;
  *std=sqrt( (sum2-sum*sum/size)/size );
}





/* Similar to favestd, but when the mask is assumed to be a long
   array.  */
void
favestdl(float *in, size_t size, float *ave, float *std, long *mask)
{
  size_t nsize1, nsize2;
  float sum, sum2;
  if(mask==NULL)
    {
      sum=floatsum(in, size);
      sum2=floatsumsquared(in, size);
    }
  else
    {
      sum=floatsummaskl(in, mask, size, &nsize1);
      sum2=floatsumsquaredmaskl(in, mask, size, &nsize2);
      if(nsize1!=nsize2)
	error(EXIT_FAILURE, 0, "A bug in favestl (lib/statistics.h). "
	      "Somehow the number of masked pixels is measured "
	      "differently. Please contact us so we can find the cause.");
      size=nsize1;
    }
  *ave=sum/size;
  *std=sqrt( (sum2-sum*sum/size)/size );
}





/* Find the average and standard deviation of all pixels on a
   region. A region in a larger image is defined by its starting pixel
   (`start`), its height (s0) and width (s1). This function will find
   the sum and sumsquared of all the nonmasked (==0 in mask[]) and non
   marked (==0 in byt) pixels in the region. */
void
floatavestdmaskbyt0inregion(float *in, unsigned char *byt,
			    unsigned char *mask, size_t startind,
			    size_t s0, size_t s1, size_t is1,
			    float *ave, float *std)
{
  float *i;
  size_t r, size=0;
  double sum=0, sumsq=0;
  unsigned char *m, *b, *fb;

  for(r=0;r<s0;++r)
    {
      i=in+startind;
      m=mask+startind;
      fb=(b=byt+startind)+s1;
      do
	{
	if(*m++==0 && *b==0)	/* `m` will definitely be checked and */
	  {			/* incremented, while `b` might not.  */
	    size++;
	    sum   += *i;
	    sumsq += *i * *i;
	  }
	++i;
	}
      while(++b<fb);
      startind+=is1;
    }

  *ave = sum/size;
  *std = sqrt( (sumsq-sum*sum/size)/size );
}




















/****************************************************************
 ********     Histogram and Cumulative Frequency Plot     *******
 ****************************************************************/
/* Set the bin lower values for all the bins. If the minimum and
   maximum are equal, then use the quantile. */
#define CHECKBINS 0
void
setbins(float *sorted, size_t size, size_t numbins, float min,
	float max, int binonzero, float quant, float **obins)
{
  size_t i;
  float tosubtract, *bins, binwidth;

  /* Allocate space for the array. The extra bin is only for internal
     purposes (so the loops for the histogram and CFP can see the end
     of the last bin). It will never be seen by the user. */
  errno=0;
  bins=*obins=calloc((numbins+1)*2, sizeof *bins);
  if(bins==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for bins in setbins (statistics.c)",
          (numbins+1)*2*sizeof *bins);

  /* If the range is not defined, find it and set the bin width. */
  if(min==max)
    {
      if(quant!=0.0f)
	{
	  min=sorted[ indexfromquantile(size, quant)   ];
	  max=sorted[ indexfromquantile(size, 1-quant) ];
	}
      else
	{
	  min=sorted[0];
	  max=sorted[size-1];
	  /* The number of elements in each bin is counted by those
	     equal or smaller to the smaller bin side and smaller than
	     the right of the bin. Therefore, if the maximum bin side
	     equals the maximum element value, it will not be
	     counted. So we slightly increase the maximum before
	     calculating the bin width */
	  max+=(max-min)/10000;
	}
    }
  binwidth=(max-min)/numbins;

  /* Set all the bin smaller sides: */
  for(i=0;i<numbins+1;++i)
    bins[i*2]=min+i*binwidth;

  /* If one bin is to be placed on zero. */
  if(binonzero)
    {
      for(i=0;i<numbins+1;++i)
	if(bins[i*2]>=0.0f) break;
      tosubtract=bins[i*2];
      for(i=0;i<numbins+1;++i)
	bins[i*2]-=tosubtract;
    }

  /* In case you want to check the bins: */
  if(CHECKBINS)
    for(i=0;i<numbins;++i)
      printf("%lu: %.4f\n", i+1, bins[i*2]);
}





#define CHECKHIST 0
void
histogram(float *sorted, size_t size, float *bins, size_t numbins,
	  int normhist, int maxhistone)
{
  float max=-FLT_MAX;
  size_t histrow, i;

  /* Fill the histogram. */
  histrow=0;
  for(i=0;i<size;++i)
    {
      if(sorted[i]<bins[histrow*2]) continue;
      while (sorted[i]>=bins[(histrow+1)*2])
	if(++histrow>=numbins) break;
      if(histrow>=numbins) break;
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

  /* In case you want to see the histogram: */
  if(CHECKHIST)
    for(i=0;i<numbins;++i)
      printf("%lu: %.4f %.4F\n", i+1, bins[i*2], bins[i*2+1]);
}





void
cumulativefp(float *sorted, size_t size, float *bins, size_t numbins,
	     int normcfp)
{
  float prevind=0;
  size_t cfprow=0, i, numinds=0;

  /* Fill the Cumulative frequency plot: */
  for(i=0;i<size;++i)
    {
      if(sorted[i]<bins[cfprow*2]) continue;
      while (sorted[i]>=bins[(cfprow+1)*2])
	{
	  if(numinds>0)
	    prevind=bins[cfprow*2+1]/=numinds; /* Divide by num indexs */
	  else
	    bins[cfprow*2+1]=prevind;
	  numinds=0;
	  if(++cfprow>=numbins)
	    break;
	}
      if(cfprow>=numbins) break;
      bins[cfprow*2+1]+=i;	/* Sum of indexs (see above for average) */
      ++numinds;
    }

  /* For a normalized CFP: */
  if(normcfp)
    for(i=0;i<numbins;++i)
      bins[i*2+1]/=size;

  /* In case you want to see the CFP:
  for(i=0;i<numbins;++i)
    printf("%lu: %.4f %.4F\n", i+1, bins[i*2], bins[i*2+1]);
  */
}





void
savehist(float *sorted, size_t size, size_t numbins,
	 char *filename, char *histname, size_t id)
{
  FILE *fp;
  size_t i;
  int binonzero=0, normhist=0, maxhistone=0;
  float d, *bins, min=0.0f, max=0.0f, quant=0.0f;

  /* Set the bin sides: */
  setbins(sorted, size, numbins, min, max, binonzero, quant, &bins);

  /* Set the size of half a bin width:*/
  d=(bins[2]-bins[0])/2;

  /* Fill the histogram: */
  histogram(sorted, size, bins, numbins, normhist, maxhistone);

  /* Open the file for writing and save the histogram: */
  errno=0;
  fp=fopen(filename, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "Couldn't open file %s", filename);
  fprintf(fp, "# %s: S/N histogram of %s in mesh %lu.\n",
	  PACKAGE_STRING, histname, id);
  fprintf(fp, "# %lu points in %lu bins\n", size, numbins);
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
indexfromquantile(size_t size, float quant)
{
  float floatindex;

  if(quant>1.0f)
    error(EXIT_FAILURE, 0, "The quantile in indexfromquantile (statistics.c) "
          "Should be smaller.");

  /* Find the index of the quantile. */
  floatindex=(float)size*quant;
  /*
  printf("quant: %f, size: %lu, findex: %f\n", quant, size, floatindex);
  */
  /* Note that in the conversion from float to size_t, the floor
     integer value of the float will be used. */
  if( floatindex - (int)floatindex > 0.5 )
    return floatindex+1;
  else
    return floatindex;
}




















/****************************************************************
 *****************        Sigma clip         ********************
 ****************************************************************/
/* This function will repeatedly sigma clip the data and return the
   median. It is assumed that the data have been ordered.

   o1_n0: Ordered (1), not ordered (0).
*/
int
sigmaclip_converge(float *array, int o1_n0, size_t num_elem,
		   float sigma_multiple, float accuracy,
		   float *outave, float *outmed, float *outstd)
{
  size_t counter=0;
  float *start, *oldstart, *dpt;
  float ave=*outave, med=*outmed, std=*outstd;
  float oldstd=0, oldmed=0, oldave=0, *orderedarray;

  if(o1_n0==0)
    {
      floatcopy(array, num_elem, &orderedarray);
      qsort(orderedarray, num_elem, sizeof*orderedarray,
	    floatincreasing);
    }
  else orderedarray=array;

  start=orderedarray;
  while(1)
    {
      oldstart=start;

      med=*(start+num_elem/2);
      favestd(start, num_elem, &ave, &std, NULL);

      /*
      printf("%lu: mean: %.2f, med: %.2f, std= %.2f, num inside: %lu\n",
	     counter+1, ave, med, std, num_elem);
      */

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
  return 0;
}





/* This function will do a certain number of sigma clips and return
   the final average, median and std. o1_n0: 1: initially ordered. 2:
   initially not ordered.*/
int
sigmaclip_certainnum(float *array, int o1_n0, size_t num_elem,
		     float sigma_multiple, size_t numtimes,
		     float *outave, float *outmed, float *outstd)
{
  size_t counter=0;
  float ave=*outave, med=*outmed, std=*outstd;
  float *start, *oldstart, *dpt, *orderedarray;

  if(o1_n0==0)
    {
      floatcopy(array, num_elem, &orderedarray);
      qsort(orderedarray, num_elem, sizeof*orderedarray,
	    floatincreasing);
    }
  else orderedarray=array;

  start=orderedarray;
  for(counter=0;counter<numtimes+1;++counter)
    {
      oldstart=start;

      med=*(start+num_elem/2);
      favestd(start, num_elem, &ave, &std, NULL);

      /*
      printf("%lu: mean: %.2f, med: %.2f, std= %.2f, num inside: %lu\n",
             counter+1, ave, med, std, num_elem);
      */

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
  return 1;
}
