/*********************************************************************
Functions to manipulate arrays.
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
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>

#include "arraymanip.h"




/*********************************************************************
 **********************       Initialize        **********************
 *********************************************************************/
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




void
longinit(long *in, size_t size, const long v)
{
  long *end=in+size;
  do *in++=v; while(in<end);
}





void
longinitonregion(long *in, const long v, size_t start, size_t s0,
                 size_t s1, size_t is1)
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




















/*********************************************************************
 **********************       Copy array        **********************
 *********************************************************************/
void
ucharcopy(unsigned char *in, size_t size, unsigned char **out)
{
  unsigned char *fp=in+size, *o;

  errno=0;
  o=*out=malloc(size*sizeof *out);
  if(*out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for copying", size);
  do *o++=*in; while(++in<fp);
}





void
floatcopy(float *in, size_t size, float **out)
{
  float *fp=in+size, *o;

  errno=0;
  o=*out=malloc(size*sizeof *out);
  if(*out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for copying", size);
  do *o++=*in; while(++in<fp);
}





void
floatcopyvalues(float *in, size_t size, float **out)
{
  float *fp=in+size, *o=*out;
  do *o++=*in; while(++in<fp);
}




















/*********************************************************************
 **********************         Values          **********************
 *********************************************************************/
void
fsetconst(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in = a;
  while(++in<fpt);
}





void
freplacevalue(float *in, size_t size, float from, float to)
{
  float *fpt;
  fpt=in+size;
  if(isnan(from))
    {
      do
        *in = isnan(*in) ? to : *in;
      while(++in<fpt);
    }
  else
    {
      do
        *in = *in==from ? to : *in;
      while(++in<fpt);
    }
}





/* Move all the non-NaN elements in the array to the start of the
   array so that the non-NaN alements are contiguous. This is useful
   for cases where you want to sort the data.*/
void
nonans(float *in, size_t *size)
{
  size_t outsize=0;
  float *f=in, *fp=in+*size;
  do
    if(!isnan(*in))
      {
        *f++=*in;
        ++outsize;
      }
  while(++in<fp);
  *size=outsize;
}



















/*********************************************************************
 **********************   Multiply or Sum with  **********************
 *********************************************************************/
void
fmultipconst(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in *= a;
  while(++in<fpt);
}





void
fsumconst(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in += a;
  while(++in<fpt);
}





float *
fsumarrays_return(float *in1, float *in2, size_t size)
{
  float *out, *o, *op;

  errno=0;
  o=out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for out in fsumarrays "
          "(arraymanip.c)", size*sizeof *out);

  op=o+size;
  do *o = *in1++ + *in2++; while(++o<op);

  return out;
}





void
dmultipconst(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;
  do *in *= a; while(++in<dpt);
}




void
dmultiparrays(double *in1, double *in2, size_t size)
{
  double *dpt;
  dpt=in1+size;
  do *in1 *= *in2++; while(++in1<dpt);
}





void
ddivideconst(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;
  do *in /= a; while(++in<dpt);
}





void
dconstdivide(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;
  do *in = a / *in; while(++in<dpt);
}




void
ddividearrays(double *in1, double *in2, size_t size)
{
  double *dpt;
  dpt=in1+size;
  do *in1 /= *in2++; while(++in1<dpt);
}




void
dsumconst(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;
  do *in += a; while(++in<dpt);
}





void
dsumarrays(double *in1, double *in2, size_t size)
{
  double *dpt;
  dpt=in1+size;
  do *in1 += *in2++; while(++in1<dpt);
}



void
dsubtractconst(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;
  do *in -= a; while(++in<dpt);
}





void
dconstsubtract(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;
  do *in = a - *in; while(++in<dpt);
}





void
dsubtractarrays(double *in1, double *in2, size_t size)
{
  double *dpt;
  dpt=in1+size;
  do *in1 -= *in2++; while(++in1<dpt);
}



void
dpowerconst(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;

  /* Since simply multiplying or taking the square root are more
     efficient than calling the power function. */
  if(a==2.0f)      do *in *= *in; while(++in<dpt);
  else if(a==0.5f) do *in = sqrt(*in); while(++in<dpt);
  else             do *in = pow(*in, a); while(++in<dpt);
}





void
dconstpower(double *in, size_t size, double a)
{
  double *dpt;
  dpt=in+size;
  do *in = pow(a, *in); while(++in<dpt);
}





void
dpowerarrays(double *in1, double *in2, size_t size)
{
  double *dpt;
  dpt=in1+size;
  do *in1 = pow(*in1, *in2++); while(++in1<dpt);
}





void
dlogarray(double *in, size_t size)
{
  double *dpt;
  dpt=in+size;

  do *in = log(*in); while(++in<dpt);
}





void
dlog10array(double *in, size_t size)
{
  double *dpt;
  dpt=in+size;

  do *in = log10(*in); while(++in<dpt);
}





void
dabsarray(double *in, size_t size)
{
  double *dpt;
  dpt=in+size;

  do *in = fabs(*in); while(++in<dpt);
}
