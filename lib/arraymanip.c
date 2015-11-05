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
gal_arraymanip_uchar_init_on_region(unsigned char *in, const unsigned char v,
                                    size_t start, size_t s0, size_t s1,
                                    size_t is1)
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
gal_arraymanip_long_init(long *in, size_t size, const long v)
{
  long *end=in+size;
  do *in++=v; while(in<end);
}





void
gal_arraymanip_long_init_on_region(long *in, const long v, size_t start,
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




















/*********************************************************************
 **********************       Copy array        **********************
 *********************************************************************/
void
gal_arraymanip_uchar_copy(unsigned char *in, size_t size, unsigned char **out)
{
  unsigned char *fp=in+size, *o;

  errno=0;
  o=*out=malloc(size*sizeof *out);
  if(*out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for copying", size);
  do *o++=*in; while(++in<fp);
}





void
gal_arraymanip_float_copy(float *in, size_t size, float **out)
{
  float *fp=in+size, *o;

  errno=0;
  o=*out=malloc(size*sizeof *out);
  if(*out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for copying", size);
  do *o++=*in; while(++in<fp);
}





void
gal_arraymanip_float_copy_values(float *in, size_t size, float **out)
{
  float *fp=in+size, *o=*out;
  do *o++=*in; while(++in<fp);
}




















/*********************************************************************
 **********************         Values          **********************
 *********************************************************************/
void
gal_arraymanip_fset_const(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in = a;
  while(++in<fpt);
}





void
gal_arraymanip_freplace_value(float *in, size_t size, float from, float to)
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
gal_arraymanip_no_nans(float *in, size_t *size)
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
gal_arraymanip_fmultip_const(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in *= a;
  while(++in<fpt);
}





void
gal_arraymanip_fsum_const(float *in, size_t size, float a)
{
  float *fpt;
  fpt=in+size;
  do
    *in += a;
  while(++in<fpt);
}





float *
gal_arraymanip_fsum_arrays(float *in1, float *in2, size_t size)
{
  float *out, *o, *op;

  errno=0;
  o=out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for out in gal_arraymanip_fsum_arrays"
          " (arraymanip.c)", size*sizeof *out);

  op=o+size;
  do *o = *in1++ + *in2++; while(++o<op);

  return out;
}
