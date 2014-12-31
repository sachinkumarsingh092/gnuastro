/*********************************************************************
Functions to manipulate arrays.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>

#include "stats.h"
#include "arraymanip.h"










/****************************************************************
 *****************   Initialize an array     ********************
 ****************************************************************/
void
longinit(long *in, size_t size, const long v)
{
  long *fpt;
  fpt=in+size;
  do
    *in=v;
  while(++in<fpt);
}





void
longinitonregion(long *in, const long v, size_t start,
		 size_t s0, size_t s1, size_t is1)
{
  size_t r;
  long *p, *fp;

  for(r=0;r<s0;++r)
    {
      fp = (p=in+start) + s1;
      do
	*p=v;
      while(++p<fp);
      start+=is1;
    }
}





void
floatinit(float *in, size_t size, const float v)
{
  float *fpt;
  fpt=in+size;
  do
    *in=v;
  while(++in<fpt);
}





void
ucharinit(unsigned char *in, size_t size,
          const unsigned char v)
{
  unsigned char *fpt;
  fpt=in+size;
  do
    *in=v;
  while(++in<fpt);
}





void
ucharinitonregion(unsigned char *in, const unsigned char v,
		  size_t start, size_t s0, size_t s1, size_t is1)
{
  size_t r;
  unsigned char *p, *fp;

  for(r=0;r<s0;++r)
    {
      fp = (p=in+start) + s1;
      do
	*p=v;
      while(++p<fp);
      start+=is1;
    }
}




















/****************************************************************
 *****************       Change values       ********************
 ****************************************************************/
void
changefvalueinarray(float *in, size_t size,
        float from, float to)
{
  float *pt=in, *fpt;
  fpt=in+size;
  for(;pt<fpt;++pt)
    if(*pt==from) *pt=to;
}





void
truncfarray(float *in, size_t size, float thresh)
{
  float *pt=in, *fpt;
  fpt=in+size;
  for(;pt<fpt;++pt)
    if(*pt>thresh) *pt=thresh;
}





void
setabovetozerofarray(float *in, size_t size, float thresh)
{
  float *pt=in, *fpt;
  fpt=in+size;
  for(;pt<fpt;++pt)
    if(*pt>thresh)
      *pt=0;
}





void
changelvalueinarray(long *in, size_t size,
                    long from, long to)
{
  long *pt=in, *fpt;
  fpt=in+size;
  for(;pt<fpt;++pt)
    if(*pt==from) *pt=to;
}





void
floatinvert(float *in, size_t size, float base)
{
    float *pt=in;
    for(;pt<in+size;++pt)
        *pt=base-*pt;
}




















/*********************************************************************
 **********************      Change type:       **********************
 *********************************************************************/
void
convertftd(float *f, size_t size, double **d)
{
  double *td, *dp;
  float *fp=f, *lfp;

  td=malloc(size*sizeof*td);
  assert(td!=NULL);

  lfp=fp+size;  dp=td;
  do
    *dp++=*fp;
  while(++fp<lfp);

  *d=td;
}





/* Convert any of the known types to float. The assertion that the
   output array be NULL before comming into this array is to make sure
   it was not allocated before this.  */
void
converttofloat(void *in, int bitpix, size_t size, float **out,
	       int byt_code, int short_code, int long_code,
	       int float_code, int double_code)
{
  float *o;
  unsigned char *u, *up;
  short         *s, *sp;
  long          *l, *lp;
  float         *f, *fp;
  double        *d, *dp;

  if(bitpix==float_code)
    {
      *out=in;
      return;
    }

  assert( *out==NULL );
  assert( (o=*out=malloc(size*sizeof *o))!=NULL  );

  if(bitpix==byt_code)
    {u=in; up=u+size; do *o++=*u; while(++u<up);}
  else if(bitpix==short_code)
    {s=in; sp=s+size; do *o++=*s; while(++s<sp);}
  else if(bitpix==long_code)
    {l=in; lp=l+size; do *o++=*l; while(++l<lp);}
  else if(bitpix==float_code)
    {f=in; fp=f+size; do *o++=*f; while(++f<fp);}
  else if(bitpix==double_code)
    {d=in; dp=d+size; do *o++=*d; while(++d<dp);}
  else
    {
      printf("\n\nbitpix=%d Not recognized in converting to float!\n",
	     bitpix);
      printf("Recognized bitpix values are:\n");
      printf("\tUnsigned char: %d\n", byt_code);
      printf("\tShort: %d\n", short_code);
      printf("\tLong: %d\n", long_code);
      printf("\tFloat: %d\n", float_code);
      printf("\tDouble: %d\n", double_code);
      exit(EXIT_FAILURE);
    }
}





/* Similar to convertofloat(), but convert to unsigne char. */
void
converttouchar(void *in, int bitpix, size_t size, unsigned char **out,
	       int byt_code, int short_code, int long_code,
	       int float_code, int double_code)
{
  unsigned char *o;
  unsigned char *u, *up;
  short         *s, *sp;
  long          *l, *lp;
  float         *f, *fp;
  double        *d, *dp;

  if(bitpix==byt_code)
    {
      *out=in;
      return;
    }

  assert( *out==NULL );
  assert( (o=*out=malloc(size*sizeof *o))!=NULL  );

  if(bitpix==byt_code)
    {u=in; up=u+size; do *o++=*u>0; while(++u<up);}
  else if(bitpix==short_code)
    {s=in; sp=s+size; do *o++=*s>0; while(++s<sp);}
  else if(bitpix==long_code)
    {l=in; lp=l+size; do *o++=*l>0; while(++l<lp);}
  else if(bitpix==float_code)
    {f=in; fp=f+size; do *o++=*f>0; while(++f<fp);}
  else if(bitpix==double_code)
    {d=in; dp=d+size; do *o++=*d>0; while(++d<dp);}
  else
    {
      printf("\n\nbitpix=%d Not recognized in converting to uchar!\n",
	     bitpix);
      printf("Recognized bitpix values are:\n");
      printf("\tUnsigned char: %d\n", byt_code);
      printf("\tShort: %d\n", short_code);
      printf("\tLong: %d\n", long_code);
      printf("\tFloat: %d\n", float_code);
      printf("\tDouble: %d\n", double_code);
      exit(EXIT_FAILURE);
    }
}





/* Similar to convertouchar(), but there is a threshold to define a
   good(=2) vs. bad(=1) uchar array. Any pixel above
   `goodbadthresh` is considered to be a bad (=1) pixel while  */
void
converttouchargoodbad(void *in, int bitpix, size_t size, unsigned char **out,
		      float goodbadthresh, int byt_code, int short_code,
		      int long_code, int float_code, int double_code)
{
  unsigned char *o;
  unsigned char *u, *up, ugb=goodbadthresh;
  short         *s, *sp, sgb=goodbadthresh;
  long          *l, *lp, lgb=goodbadthresh;
  float         *f, *fp, fgb=goodbadthresh;
  double        *d, *dp, dgb=goodbadthresh;

  if(bitpix==byt_code)
    {
      *out=in;
      return;
    }

  assert( *out==NULL );
  assert( (o=*out=malloc(size*sizeof *o))!=NULL  );

  if(bitpix==byt_code)
    {
      u=in; up=u+size;
      do {if(*u) *o=*u>ugb?1:2; else *o=0; ++o;}
      while(++u<up);
    }
  else if(bitpix==short_code)
    {
      s=in; sp=s+size;
      do {if(*s) *o=*s>sgb?1:2; else *o=0; ++o;}
      while(++s<sp);
    }
  else if(bitpix==long_code)
    {
      l=in; lp=l+size;
      do {if(*l) *o=*l>lgb?1:2; else *o=0; ++o;}
      while(++l<lp);
    }
  else if(bitpix==float_code)
    {
      f=in; fp=f+size;
      do {if(*f) *o=*f>fgb?1:2; else *o=0; ++o;}
      while(++f<fp);
    }
  else if(bitpix==double_code)
    {
      d=in; dp=d+size;
      do {if(*d) *o=*d>dgb?1:2; else *o=0; ++o;}
      while(++d<dp);
    }
  else
    {
      printf("\n\nbitpix=%d Not recognized in converting to uchar!\n",
	     bitpix);
      printf("Recognized bitpix values are:\n");
      printf("\tUnsigned char: %d\n", byt_code);
      printf("\tShort: %d\n", short_code);
      printf("\tLong: %d\n", long_code);
      printf("\tFloat: %d\n", float_code);
      printf("\tDouble: %d\n", double_code);
      exit(EXIT_FAILURE);
    }
}




















/*********************************************************************
 **********************      Remove masked      **********************
 *********************************************************************/
void
removemasked(float **in, size_t size, unsigned char *mask,
             size_t *nsize)
{
  float *tmp;
  size_t i, counter=0;

  tmp=*in;
  for(i=0;i<size;++i)
    if(mask[i]==0)
      tmp[counter++]=tmp[i];
  *nsize=counter;

  *in=(float *)realloc(*in, counter*sizeof(float));
  assert(*in!=NULL);
}




















/*********************************************************************
 **********************       Make copy         **********************
 *********************************************************************/
void
floatarrcpynomalloc(float *in, size_t size, float *out)
{
  float *fpt;
  fpt=out+size;
  do
    *out=*in++;
  while(++out<fpt);
}




/* Copy a region of a larger array into a smaller array. *in is the
   larger array and *out is the smaller array.*/
void
floatarrcpyallregionltsnomalloc(float *in, float *out, size_t start,
				size_t s0, size_t s1, size_t is1)
{
  size_t r;
  float *p, *fp;

  for(r=0;r<s0;++r)
    {
      fp = (p=in+start) + s1;
      do
	*out++=*p;
      while(++p<fp);
      start+=is1;
    }
}





/* Copy a region (of size s0xs1 starting fron index `start`) of the
   larger array `in` to the smaller array `out`. Space for out must
   already be allocated prior to calling this function. If a mask
   array is provided (mask!=NULL), then only copy those elements that
   are not masked. Note that the mask is the same size as the input
   image.

   It is upto the user to make sure that the desired box is inside the
   image. This function was primarily made to operate on a mesh
   grid. By design, all the meshes fit inside the image completely.*/
void
floatarrcpyregionnomalloc(float *in, float *out, size_t start,
			  size_t s0, size_t s1, size_t is1,
			  unsigned char *mask, size_t *outnum)
{
  size_t r, n=0;
  float *p, *fp;
  unsigned char *m;

  if(mask)
    for(r=0;r<s0;++r)
      {
	m=mask+start;
	fp = (p=in+start) + s1;
	do
	  if(*m++==0)
	    {
	      *out++=*p;
	      ++n;
	    }
	while(++p<fp);
	start+=is1;
      }
  else
    {
      for(r=0;r<s0;++r)
	{
	  fp = (p=in+start) + s1;
	  do
	    *out++=*p;
	  while(++p<fp);
	  start+=is1;
	}
      n=s0*s1;
    }
  *outnum=n;
}





void
floatarrcpy(float *in, size_t size, float **out)
{
  float *ipt, *opt, *fpt;
  ipt=in;
  *out=malloc(sizeof(float)*size);
  assert(*out!=NULL);
  fpt=*out+size;
  opt=*out;
  do
    *opt=*ipt++;
  while(++opt<fpt);
}





/* Copy the array without copying the masked elements. */
void
floatarrcpymask(float *in, size_t size,
        unsigned char *mask, size_t *nsize,
        float **out)
{
  floatarrcpy(in, size, out);

  if(mask==NULL)
    {
      *nsize=size;
      return;
    }
  else removemasked(out, size, mask, nsize);
}





void
uchararrcpy(unsigned char *in, size_t size, unsigned char **out)
{
  unsigned char *ipt, *opt, *fpt;
  ipt=in;
  *out=malloc(size*sizeof(unsigned char));
  assert(*out!=NULL);
  fpt=*out+size;
  opt=*out;
  do
    *opt=*ipt++;
  while(++opt<fpt);
}





/* Copy a region of a larger array into a smaller array. *in is the
   larger array and *out is the smaller array.*/
void
uchararrcpyallregionltsnomalloc(unsigned char *in, unsigned char *out,
				size_t start, size_t s0, size_t s1,
				size_t is1)
{
  size_t r;
  unsigned char *p, *fp;

  for(r=0;r<s0;++r)
    {
      fp = (p=in+start) + s1;
      do
	*out++=*p;
      while(++p<fp);
      start+=is1;
    }
}





/* Copy a region of a smaller array into a larger array. *in is the
   smaller array and *out is the larger array. */
void
uchararrcpyallregionstlnomalloc(unsigned char *in, unsigned char *out,
				size_t start, size_t s0, size_t s1,
				size_t is1)
{
  size_t r;
  unsigned char *p, *fp;

  for(r=0;r<s0;++r)
    {
      fp = (p=out+start) + s1;
      do
	*p=*in++;
      while(++p<fp);
      start+=is1;
    }
}





void
longarrcpy(long *in, size_t size, long **out)
{
  long *opt, *fpt;
  assert( (opt=*out=malloc(size*sizeof **out) )!=NULL );
  fpt=*out+size;
  do
    *opt=*in++;
  while(++opt<fpt);
}




















/*********************************************************************
 **********************      Mask an array      **********************
 *********************************************************************/
void
maskfarray(float *in, unsigned char *mask, size_t size,
	   unsigned char f1_b0)
{
  float min, *fpt=in, *ffpt;
  unsigned char *mpt=mask;

  floatmin(in, size, &min);

  ffpt=fpt+size;
  do
    {
      if(*mpt==f1_b0)
	*fpt=min;
      ++mpt;
    }
  while(++fpt<ffpt);
}





void
bytmaskfarray(float *in, unsigned char *byt, unsigned char *mask,
	      size_t size, unsigned char f1_b0)
{
  float min, *fpt=in, *ffpt;

  floatmin(in, size, &min);
  ffpt=fpt+size;
  do
    {
      if(*mask++ || *byt==f1_b0)
	*fpt=min;
      byt++;
    }
  while(++fpt<ffpt);
}





void
masklfarray(float *in, long *mask, size_t size, long f1_b0)
{
  long *mpt=mask;
  float min, *fpt=in, *ffpt;

  floatmin(in, size, &min);

  ffpt=fpt+size;
  if(f1_b0==0)
    do
      {
	if(*mpt==0) *fpt=min;
	++mpt;
      }
    while(++fpt<ffpt);
  else if(f1_b0==1)
    do
      {
	if(*mpt>0) *fpt=min;
	++mpt;
      }
    while(++fpt<ffpt);
}


















/*********************************************************************
 **********************   Multiply or Sum with  **********************
 *********************************************************************/
void
floatarrmwith(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in *= a;
  while(++in<fpt);
}





void
floatarrswith(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in += a;
  while(++in<fpt);
}




















/*********************************************************************
 **********************      Change values:     **********************
 *********************************************************************/
void
floatsetbelowtozero(float *in, size_t size, float min)
{
  float *fpt;
  fpt=in+size;
  do
    if(*in<min) *in=0;
  while(++in<fpt);
}





void
floatsetbelowtomin(float *in, size_t size, float min)
{
  float *fpt;
  fpt=in+size;
  do
    if(*in<min) *in=min;
  while(++in<fpt);
}





void
floatsetabovetomax(float *in, size_t size, float max)
{
  float *fpt;
  fpt=in+size;
  do
    if(*in>max) *in=max;
  while(++in<fpt);
}




















/*********************************************************************
 **********************      Shrink Array       **********************
 *********************************************************************/
/* Check to see if a box defined by the two points (x1,y1) and (x2,y2)
   is inside an array of size size1 and size2. If it doesn't overlap,
   then x1=x2 and y1=y2.*/
void
checkifinarray(int *x1, int *y1, int *x2, int *y2, int s0, int s1)
{
  int temp;
  if(*x1==*x2 && *y1==*y2) return;

  if(*x2<*x1){temp=*x1;*x1=*x2;*x2=temp;}
  if(*y2<*y1){temp=*y1;*y1=*y2;*y2=temp;}

  if(*x1<0) *x1=0;    if(*x1>s0) *x1=s0;
  if(*y1<0) *y1=0;    if(*y1>s1) *y1=s1;
  if(*x2<0) *x2=0;    if(*x2>s0) *x2=s0;
  if(*y2<0) *y2=0;    if(*y2>s1) *y2=s1;
}





/* We have a large array of size (size1*size2). We want to shrink this
    array, such that (x1,y1) comes down to point (0,0) and the new
    array now only extends to the old (x2,y2). So the size of the new
    array is: w1*w2 where w1=x2-x1 and w2=y2-y1.

    If the desired region is totally out of the array, a NULL pointer
    is returned.*/
void
floatshrinkarray(float **in, int size1, int size2,
		 int x1, int y1, int x2, int y2)
{
  float *ifpt, *ofpt, *rowstart;
  size_t i, ux1, uy1, us2, w1, w2;

  checkifinarray(&x1, &y1, &x2, &y2, size1, size2);
  if(x1==x2 || y1==y2) 		/* The required region does not */
    {				/* overlap with the array. */
      free(*in);
      *in=NULL;
      return;
    }
  /* The region covers the whole image, no need for the next step. */
  if(x1==0 && y1==0 && x2==size1 && y2==size2) return;

  w1=x2-x1;  w2=y2-y1;
  ux1=x1; uy1=y1; us2=size2;  /* The inputs are int (can be negative,
				 which is allowed: will become zero).
				 but pointers are unsigned, so to
				 faciliate the process in the loop,
				 they are converted to size_t. */
  for(i=0;i<w1;++i)
    {
      ofpt=rowstart=*in+i*w2;
      ifpt=*in+(ux1+i)*us2+uy1;
      do
	*ofpt=*ifpt++;
      while(++ofpt<rowstart+w2);

    }
  *in=(float *)realloc(*in, w1*w2*sizeof(float));
  assert(*in!=NULL);
}





/* similar to floatshrinkarray(), but the old array is kept untouched
   and a new one is created to keep the cropped region. */
void
floatshrinkarraytonew(float *in, int size1, int size2,int x1, int y1,
		      int x2, int y2, float **out)
{
  float *ifpt, *ofpt, *rowstart;
  size_t i, ux1, uy1, us2, w1, w2;

  checkifinarray(&x1, &y1, &x2, &y2, size1, size2);
  if(x1==x2 || y1==y2)
    {
      *out=NULL;
      return;
    }
  if(x1==0 && y1==0 && x2==size1 && y2==size2)
    {
      floatarrcpy(in, size1*size2, out);
      return;
    }

  w1=x2-x1;  w2=y2-y1;
  *out=malloc(w1*w2*sizeof **out);
  assert(*out!=NULL);

  ux1=x1; uy1=y1; us2=size2;
  for(i=0;i<w1;++i)
    {
      ofpt=rowstart=*out+i*w2;
      ifpt=in+(ux1+i)*us2+uy1;
      do
	*ofpt=*ifpt++;
      while(++ofpt<rowstart+w2);
    }
}





void
longshrinkarraytonew(long *in, int size1, int size2,
		      int x1, int y1, int x2, int y2,
		      long **out)
{
  long *ifpt, *ofpt, *rowstart;
  size_t i, ux1, uy1, us2, w1, w2;

  checkifinarray(&x1, &y1, &x2, &y2, size1, size2);
  if(x1==x2 || y1==y2)
    {
      *out=NULL;
      return;
    }
  if(x1==0 && y1==0 && x2==size1 && y2==size2)
    {
      longarrcpy(in, size1*size2, out);
      return;
    }

  w1=x2-x1;  w2=y2-y1;
  *out=malloc(w1*w2*sizeof **out);
  assert(*out!=NULL);

  ux1=x1; uy1=y1; us2=size2;

  for(i=0;i<w1;++i)
    {
      ofpt=rowstart=*out+i*w2;
      ifpt=in+(ux1+i)*us2+uy1;
      do
	*ofpt=*ifpt++;
      while(++ofpt<rowstart+w2);

    }
}




















/*********************************************************************
 **********************      Merge arrays       **********************
 *********************************************************************/
/* a and b are two column arrays with "numrows" rows. It is assumed
   that their first columns are equal. This function will take the
   second column of b and add it to a.  */
void
floatvmerge(float *a, float *b, size_t numrows, size_t numcols,
	    float **out)
{
  size_t i;
  float *tmp;

  assert(numcols==3);

  tmp=malloc(numcols*numrows*sizeof(float));
  assert(tmp!=NULL);

  for(i=0;i<numrows;++i)
    {
      if(b[i*2]!=a[i*2])
	{
	  printf("\n\n\tError in floatvmerge() (arraymanip.c)\n");
	  printf("\t\tThe first columns are not equal on row %lu\n\n",
		 i);
	  exit(EXIT_FAILURE);
	}
      else
	{
	  tmp[i*numcols]=a[i*2];
	  tmp[i*numcols+1]=a[i*2+1];
	  tmp[i*numcols+2]=b[i*2+1];
	}
    }
  *out=tmp;
}




















/*********************************************************************
 **********************       Print array       **********************
 *********************************************************************/
/* Print a 2D array into an ascii file. Note that s0 and s1 are one
   larger than the actual number of columns and number of rows.  */
void
printfarray(float *array, size_t s0, size_t s1,
	    char *comment, char *filename, int decimals)
{
  FILE *fp;
  double *d;
  size_t i,j;
  char fmt[20];

  sprintf(fmt, "%%-20.%df", decimals);

  convertftd(array, s0*s1, &d);	/* For printing larger  */
  fp=fopen(filename, "w"); 	/* than 6 decimals. */

  fprintf(fp, "%s", comment);
  for(i=0;i<s0;++i)
    {
      for(j=0;j<s1;++j)
	{
	  fprintf(fp, fmt, d[i*s1+j]);
	}
      fprintf(fp, "\n");
    }
  free(d);
  fclose(fp);
}
