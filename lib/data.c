/*********************************************************************
data -- Structure and functions to represent/work with data
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/data.h>










/*********************************************************************/
/*************      Type size and allocation       *******************/
/*********************************************************************/
size_t
gal_data_sizeof(int type)
{
  /* Allocate space for the array to keep the image. */
  switch(type)
    {
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support TBIT "
            "datatype, please get in touch with us to implement it.");

      /* The parenthesis after sizeof is not a function, it is actually a
         type cast, so we have put a space between size of and the
         parenthesis to highlight this. In C, `sizeof' is an operator, not
         a function.*/
    case GAL_DATA_TYPE_UCHAR:
      return sizeof (unsigned char);

    case GAL_DATA_TYPE_LOGICAL: case GAL_DATA_TYPE_CHAR:
      return sizeof (char);

    case GAL_DATA_TYPE_STRING:
      return sizeof (char *);

    case GAL_DATA_TYPE_USHORT:
      return sizeof (unsigned short);

    case GAL_DATA_TYPE_SHORT:
      return sizeof (short);

    case GAL_DATA_TYPE_UINT:
      return sizeof (unsigned int);

    case GAL_DATA_TYPE_INT:
      return sizeof (int);

    case GAL_DATA_TYPE_ULONG:
      return sizeof (unsigned long);

    case GAL_DATA_TYPE_LONG:
      return sizeof (long);

    case GAL_DATA_TYPE_LONGLONG:
      return sizeof (long long);

    case GAL_DATA_TYPE_FLOAT:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0, "`float` is not 32 bits on this machine");
      return sizeof (float);

    case GAL_DATA_TYPE_DOUBLE:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0, "`double` is not 64 bits on this machine");
      return sizeof (double);

    case GAL_DATA_TYPE_COMPLEX:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0, "`float` is not 32 bits on this machine");
      return sizeof (gsl_complex_float);

    case GAL_DATA_TYPE_DCOMPLEX:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0, "`double` is not 64 bits on this machine");
      return sizeof (gsl_complex);

    default:
      error(EXIT_FAILURE, 0, "type value of %d not recognized in "
            "gal_data_sizeof", type);
    }

  error(EXIT_FAILURE, 0, "Control has reached the end of `gal_data_sizeof' "
        "This is a bug! Please contact us at %s so we can find the cause of "
        "the problem.", PACKAGE_BUGREPORT);
  return -1;
}





/* Allocate an array based on the value of type. Note that the argument
   `size' is the number of elements, necessary in the array, the number of
   bytes each element needs will be determined internaly by this function
   using the datatype argument, so you don't have to worry about it. */
void *
gal_data_alloc(int type, size_t size)
{
  void *array;

  errno=0;
  array=malloc( size * gal_data_sizeof(type) );
  if(array==NULL)
    error(EXIT_FAILURE, errno, "array of %zu bytes in gal_data_alloc",
          size * gal_data_sizeof(type));

  return array;
}





void *
gal_data_alloc_blank(int type)
{
  /* Define the pointers. */
  void *allocated;
  unsigned char *b;
  char *c;
  char **str;
  unsigned short *us;
  short *s;
  unsigned int *ui;
  int *i;
  unsigned long *ul;
  long *l;
  long long *L;
  float *f;
  double *d;
  gsl_complex_float *cx;
  gsl_complex *dcx;

  /* Allocate the space for the blank value: */
  allocated=gal_data_alloc(1, type);

  /* Put the blank value into it. */
  errno=0;
  switch(type)
    {
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support blank "
            "values for `GAL_DATA_TYPE_BIT', please get in touch with "
            "us to see how we can implement it.");

    case GAL_DATA_TYPE_UCHAR:
      b=allocated;
      *b=GAL_DATA_BLANK_UCHAR;
      return b;

      /* CFITSIO says "int for keywords, char for table columns". Here we
         are only assuming table columns. So in practice this also applies
         to TSBYTE.*/
    case GAL_DATA_TYPE_CHAR: case GAL_DATA_TYPE_LOGICAL:
      c=allocated;
      *c=GAL_DATA_BLANK_CHAR;
      return c;

    case GAL_DATA_TYPE_STRING:
      str=allocated;
      *str=GAL_DATA_BLANK_STRING;
      return str;

    case GAL_DATA_TYPE_USHORT:
      us=allocated;
      *us=GAL_DATA_BLANK_USHORT;
      return us;

    case GAL_DATA_TYPE_SHORT:
      s=allocated;
      *s=GAL_DATA_BLANK_SHORT;
      return s;

    case GAL_DATA_TYPE_UINT:
      ui=allocated;
      *ui=GAL_DATA_BLANK_UINT;
      return ui;

    case GAL_DATA_TYPE_INT:
      i=allocated;
      *i=GAL_DATA_BLANK_INT;
      return i;

    case GAL_DATA_TYPE_ULONG:
      ul=allocated;
      *ul=GAL_DATA_BLANK_ULONG;
      return ul;

    case GAL_DATA_TYPE_LONG:
      l=allocated;
      *l=GAL_DATA_BLANK_LONG;
      return l;

    case GAL_DATA_TYPE_LONGLONG:
      L=allocated;
      *L=GAL_DATA_BLANK_LONGLONG;
      return L;

    case GAL_DATA_TYPE_FLOAT:
      f=allocated;
      *f=GAL_DATA_BLANK_FLOAT;
      return f;

    case GAL_DATA_TYPE_DOUBLE:
      d=allocated;
      *d=GAL_DATA_BLANK_DOUBLE;
      return d;

    case GAL_DATA_TYPE_COMPLEX:
      cx=allocated;
      GSL_SET_COMPLEX(cx, GAL_DATA_BLANK_FLOAT, GAL_DATA_BLANK_FLOAT);
      return cx;

    case GAL_DATA_TYPE_DCOMPLEX:
      dcx=allocated;
      GSL_SET_COMPLEX(dcx, GAL_DATA_BLANK_DOUBLE, GAL_DATA_BLANK_DOUBLE);
      return dcx;

    default:
      error(EXIT_FAILURE, 0, "type value of %d not recognized in "
            "`gal_data_alloc_blank'", type);
    }

  return NULL;
}
