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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
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
 **********************    Initialize values    **********************
 *********************************************************************/
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




















/*********************************************************************
 **********************    Initialize values    **********************
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
