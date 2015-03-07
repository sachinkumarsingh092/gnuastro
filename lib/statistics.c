/*********************************************************************
Statistical functions.
This is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
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

#include "statistics.h"










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
  float tmax=-FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(*in>tmax) tmax=*in;
  while(++in<fpt);
  *max=tmax;
}





void
doublemin(double *in, size_t size, double *min)
{
  double tmin=FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(*in<tmin) tmin=*in;
  while(++in<fpt);
  *min=tmin;
}





void
doublemax(double *in, size_t size, double *max)
{
  double tmax=-FLT_MAX, *fpt;
  fpt=in+size;
  do
    if(*in>tmax) tmax=*in;
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
dmax_withindex(double *in, size_t size,
        double *max, size_t *index)
{
  size_t tindex=0;
  double *fpt, *pt=in, tmax=-FLT_MAX;

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
  float *pt=in, *fpt, tmax=-FLT_MAX;

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
  double *pt=in, *fpt, tmin=FLT_MAX;

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
  float *pt=in, *fpt, tmin=FLT_MAX;

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
