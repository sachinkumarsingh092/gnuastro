/*********************************************************************
Functions for simple statistical analysis.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "stats.h"
#include "config.h"
#include "attaavv.h"
#include "forqsort.h"
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
    if(*in<tmin) tmin=*in;
  while(++in<fpt);
  *min=tmin;
}





void
floatmax(float *in, size_t size, float *max)
{
  float tmax=MINFD, *fpt;
  fpt=in+size;
  do
    if(*in>tmax) tmax=*in;
  while(++in<fpt);
  *max=tmax;
}





void
floatmaxmasked(float *in, unsigned char *mask, size_t size, float *max)
{
  float tmax=MINFD, *fpt;
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
  float smax=MINFD, max=MINFD, *fpt;
  fpt=in+size;
  do
    {
      if(*in>max)
	{
	  smax=max;
	  max=*in;
	}
      else if(*in>smax) smax=*in;
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
      if(*in<min)
	{
	  smin=min;
	  min=*in;
	}
      else if(*in<smin) smin=*in;
    }
  while(++in<fpt);
  *secondmin=smin;
}





void
fminmax(float *in, size_t size, float *min, float *max)
{
  float tmin=FLT_MAX, tmax=MINFD, *fpt;
  fpt=in+size;
  do
    {
      if     (*in>tmax) tmax=*in;
      else if(*in<tmin) tmin=*in;
    }
  while(++in<fpt);
  *max=tmax;
  *min=tmin;
}





void
dmax_withindex(double *in, size_t size,
        double *max, size_t *index)
{
  size_t tindex=0;
  double *fpt, *pt=in, tmax=MINFD;

  fpt=pt+size;
  do
    if(*pt>tmax)
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
  float *pt=in, *fpt, tmax=MINFD;

  fpt=pt+size;
  do
    if(*pt>tmax)
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
  double *pt=in, *fpt, tmin=MAXFD;

  fpt=pt+size;
  do
    if(*pt<tmin)
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
  float *pt=in, *fpt, tmin=MAXFD;

  fpt=pt+size;
  do
    if(*pt<tmin)
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
  float *pt, *fpt;
  double sum=0;
  fpt=in+size;
  pt=in;
  do
    sum+=*pt;
  while(++pt<fpt);
  return sum;
}





float
floatsumsquared(float *in, size_t size)
{
  float *fpt;
  double sum=0;
  fpt=in+size;
  do
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
favestd(float *in, size_t size, float *ave, float *std,
	unsigned char *mask)
{
  size_t nsize1, nsize2;
  float sum, sum2;
  if(mask)
    {
      sum=floatsummask(in, mask, size, &nsize1);
      sum2=floatsumsquaredmask(in, mask, size, &nsize2);
      assert(nsize1==nsize2);
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
favestdl(float *in, size_t size, float *ave, float *std,
    long *mask)
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
      assert(nsize1==nsize2);
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





void
floatavestdmaskbyt0inregionsclip(float *in, unsigned char *byt,
				 unsigned char *mask, size_t startind,
				 size_t s0, size_t s1, size_t is1,
				 size_t numback, float *ave, float *std)
{
  size_t r, size=0;
  float *i, *pixs, med;
  unsigned char *m, *b, *fb;

  assert( (pixs=malloc(numback*sizeof *pixs))!=NULL );

  for(r=0;r<s0;++r)
    {
      i=in+startind;
      m=mask+startind;
      fb=(b=byt+startind)+s1;
      do
	{
	  if(*m++==0 && *b==0) /* `m` will definitely be checked and */
	    pixs[size++]=*i;   /* incremented, while `b` might not.  */
	  ++i;
	}
      while(++b<fb);
      startind+=is1;
    }
  assert(size==numback);

  /* Sort the array for sigma clipping. */
  qsort(pixs, size, sizeof *pixs, floatincreasing);
  sigmaclip_converge(pixs, 1, size, 3, 0.2, ave, &med, std);

  free(pixs);
}





















/****************************************************************
 *****************         Histogram        *********************
 ****************************************************************/
/* Set the bin lower values for all the bins: */
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
  assert( (bins=*obins=calloc((numbins+1)*2, sizeof *bins))!=NULL );

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
  float max=MINFD;
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





#define CHECKCFP 0
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
	    prevind=bins[cfprow*2+1]/=numinds;
	  else
	    bins[cfprow*2+1]=prevind;
	  numinds=0;
	  if(++cfprow>=numbins)
	    break;
	}
      if(cfprow>=numbins) break;
      bins[cfprow*2+1]+=i;	/* Sum of indexs */
      ++numinds;
    }

  /* For a normalized CFP: */
  if(normcfp)
    for(i=0;i<numbins;++i)
      bins[i*2+1]/=size;

  /* In case you want to see the CFP: */
  if(CHECKCFP)
    for(i=0;i<numbins;++i)
      printf("%lu: %.4f %.4F\n", i+1, bins[i*2], bins[i*2+1]);
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
  assert( (fp=fopen(filename, "w"))!=NULL );
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
 *****************      Remove outliers      ********************
 ****************************************************************/
/* Using the cumulative distribution function this funciton will
   remove outliers from a dataset. */
void
removeoutliers_flatcdf(float *arr, size_t *outsize)
{
  int firstfound=0;
  size_t size=*outsize, i, maxind;
  float *slopes, minslope, maxslope;

  qsort(arr, size, sizeof*arr, floatincreasing);

  /* Find a slopes array, think of the cumulative frequency plot when
     you want to think about slopes. */
  assert( (slopes=malloc(size*sizeof *slopes))!=NULL );
  for(i=1;i<size-1;++i) slopes[i]=2/(arr[i+1]-arr[i-1]);

  /* Find the position of the maximum slope, note that around the
     distribution mode, the difference between the values varies less,
     so two neighbouring elements have the closest values, hence the
     largest slope (when their difference is in the denominator). */
  fmax_withindex(slopes+1, size-2, &maxslope, &maxind);

  /* Find the minimum slope from the second element (for the first the
     slope is not defined. NOTE; maxind is one smaller than it should
     be because the input array to find it began from the second
     element. */
  floatsecondmin(slopes+1, maxind+1, &minslope);

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
    printf("%lu\t%.3f\t%.3f\n", i, arr[i], slopes[i]);
  printf("\n\nPlace to cut off for outliers is: %lu\n\n", *outsize);
  */
}




















/****************************************************************
 *****************        Sigma clip         ********************
 ****************************************************************/
/* This function will repeatedly sigma clip the data and return the
   median. It is assumed that the data have been ordered.  */
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
      floatarrcpy(array, num_elem, &orderedarray);
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

      /* It might happen that ave and std are nan. If so, stop the
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
  fprintf(stderr, PACKAGE": sigmaclip_converge() could not converge\n");
  exit(0);
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
      floatarrcpy(array, num_elem, &orderedarray);
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




















/****************************************************************
 *****************         Quantiles         ********************
 ****************************************************************/
/* Find the index corresponding to a certain quantile, considering the
   rounding that might be needed. */
size_t
indexfromquantile(size_t size, float quant)
{
  float floatindex;

  assert(quant<=1.0f);

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





/* Find the value corresponding to a certain quantile.  If a mask is
   supplied (mask!=NULL) then any pixel that is masked will not be
   included in the quantile finding procedure.*/
void
valuefromquantile(float *data, size_t size, float quant,
		  float *quantflux, unsigned char *mask)
{
  int chvalue=0;
  size_t fsize;
  float *tmp, ave=0;

  assert(quant>=0 && quant<=1);

  floatarrcpymask(data, size, mask, &fsize, &tmp);

  if(quant==1)
    {
      floatmax(tmp, fsize, quantflux);
      free(tmp);
      return;
    }

  /* This is to increase the speed for a symmetric
     or positively skewed distribution, which is the case
     I have so far. */
  if(quant<=0.5)
    {
      chvalue=1;
      fave(tmp, fsize, &ave, mask);
      floatsetabovetomax(tmp, fsize, ave);
    }
  else if(quant>0.8)
    {
      chvalue=1;
      fave(tmp, fsize, &ave, mask);
      floatsetbelowtomin(tmp, fsize, ave);
    }

  qsort(tmp, fsize, sizeof*tmp, floatincreasing);

  *quantflux=tmp[indexfromquantile(fsize, quant)];

  if(chvalue==1 && *quantflux==ave)
    {
      printf("\n\n\tWarning from valuefromquantile()");
      printf("in stats.c.\n\tThe approximation did not work. ");
      printf("Speed will be slower!\n\n");
      free(tmp);
      floatarrcpymask(data, size, mask, &fsize, &tmp);
      qsort(tmp, fsize, sizeof(float), floatincreasing);
      *quantflux=tmp[ indexfromquantile(fsize, quant) ];
    }

  free(tmp);
}









void
multivaluefromquantile(float *data, size_t size, float *quants,
		       float *quantfluxs, size_t numquants,
		       unsigned char *mask)
{
  float *tmp;
  size_t i, fsize;

  for(i=0;i<numquants;++i)
    assert(quants[i]>=0 && quants[i]<=1);

  floatarrcpymask(data, size, mask, &fsize, &tmp);

  qsort(tmp, fsize, sizeof*tmp, floatincreasing);

  for(i=0;i<numquants;++i)
    quantfluxs[i]=tmp[ indexfromquantile(fsize, quants[i]) ];

  free(tmp);
}





/* Similar to valuefromquantile(), but the input array will not be
   copied. Since this function will sort the array, the result is that
   the input array will be sorted after this function is run.*/
void
valuefromquantile_nocopy(float *data, size_t size,
        float quant, float *quantflux, unsigned char *mask)
{
  size_t fsize;

  assert(quant>=0 && quant<=1);

  if(quant==1)
    {
      floatmax(data, size, quantflux);
      return;
    }

  if(mask!=NULL)
    {
      removemasked(&data, size, mask, &fsize);
      size=fsize;
    }

  qsort(data, size, sizeof(float), floatincreasing);
  *quantflux=data[ indexfromquantile(size, quant) ];
}





/* Find quantile corresponding to a certain value.  If a mask is
   supplied (mask!=NULL) then any pixel that is masked will not be
   included in the quantile finding procedure.  */
void
quantilefromvalue(float *data, size_t size,
        float *quant, float quantflux, unsigned char *mask)
{
  float *tmp;
  size_t i, fsize;

  floatarrcpymask(data, size, mask, &fsize, &tmp);

  qsort(tmp, fsize, sizeof(float), floatincreasing);

  for(i=0;i<fsize;++i)
    if(tmp[i]>quantflux)
      {
	*quant=(float)i/(float)fsize;
	break;
      }
  if(i==fsize)
    {
      printf("Input flux (%f)>largest flux (%f)\n",
	     quantflux, tmp[fsize-1]);
      *quant=1;
    }

  free(tmp);
}





/* Similar to quantilefromvalue(), but the input array will not be
   copied. Since this function will sort the array, the result is
   that the input array will be sorted after this function is run.  */
void
quantilefromvalue_nocopy(float *data, size_t size,
        float *quant, float quantflux, unsigned char *mask)
{
  size_t i, fsize;

  if(mask!=NULL)
    {
      removemasked(&data, size, mask, &fsize);
      size=fsize;
    }

  qsort(data, size, sizeof(float), floatincreasing);

  for(i=0;i<size;++i)
    if(data[i]>quantflux)
      {
	*quant=(float)i/(float)size;
	break;
      }
  if(i==size)
    {
      printf("Input flux (%f)>largest flux (%f)\n",
	     quantflux, data[size-1]);
      *quant=1;
    }
}





/* Similar to quantilefromvalue(), but the input array will be assumed
   to be already sorted. */
void
quantilefromvalue_sorted(float *sorted, size_t size,
			 float *quant, float quantflux)
{
  size_t i;

  for(i=1;i<size;++i)
    if(sorted[i]>quantflux)
      {
	if(sorted[i-1]-quantflux < quantflux-sorted[i])
	  *quant=((float)(i-1))/((float)size);
	else
	  *quant=((float)i)/((float)size);
	break;
      }
  if(i==size)
    *quant=1;
}
