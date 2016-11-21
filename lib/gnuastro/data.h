/*********************************************************************
data -- Structure and functions to represent/work with data
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
#ifndef __GAL_DATA_H__
#define __GAL_DATA_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <limits.h>
#include <wcslib/wcs.h>


/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */



/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */





/* Macros: */

/* The maximum dimensionality of datasets. */
#define GAL_DATA_MAXDIM    30

/* Blank values: Note that for the unsigned types or small types (like
   char), the maximum value is considered as a blank value, since the
   minimum value of an unsigned type is zero and zero is often meaningful
   in contexts were unsigned values are used. */
#define GAL_DATA_BLANK_UCHAR      UCHAR_MAX
#define GAL_DATA_BLANK_CHAR       SCHAR_MAX
#define GAL_DATA_BLANK_LOGICAL    SCHAR_MAX
#define GAL_DATA_BLANK_STRING     NULL
#define GAL_DATA_BLANK_UINT       UINT_MAX
#define GAL_DATA_BLANK_INT        INT_MIN
#define GAL_DATA_BLANK_USHORT     USHRT_MAX
#define GAL_DATA_BLANK_SHORT      INT16_MIN
#define GAL_DATA_BLANK_ULONG      ULONG_MAX
#define GAL_DATA_BLANK_LONG       INT32_MIN
#define GAL_DATA_BLANK_LONGLONG   INT64_MIN
#define GAL_DATA_BLANK_FLOAT      NAN
#define GAL_DATA_BLANK_DOUBLE     NAN





/* Macros to identify the type of data. The macros in the comment
   parenthesis is the equivalent macro in CFITSIO. */
enum gal_data_alltypes
{
  GAL_DATA_BIT,        /* Bit              (TBIT).      */
  GAL_DATA_UCHAR,      /* Unsigned char    (TBYTE).     */
  GAL_DATA_CHAR,       /* char             (TSBYTE).    */
  GAL_DATA_STRING,     /* string           (TSTRING).   */
  GAL_DATA_UINT,       /* unsigned int     (TUINT).     */
  GAL_DATA_INT,        /* int              (TINT).      */
  GAL_DATA_USHORT,     /* unsigned short   (TUSHORT).   */
  GAL_DATA_SHORT,      /* short            (TSHORT).    */
  GAL_DATA_ULONG,      /* unsigned long    (TLONG).     */
  GAL_DATA_LONG,       /* long             (TLONG).     */
  GAL_DATA_LONGLONG,   /* long long        (TLONGLONG). */
  GAL_DATA_FLOAT,      /* float            (TFLOAT).    */
  GAL_DATA_DOUBLE,     /* double           (TDOUBLE).   */
  GAL_DATA_COMPLEX,    /* Complex float    (TCOMPLEX).  */
  GAL_DATA_DCOMPLEX,   /* Complex double   (TDCOMPLEX). */
};





/* Main data structure

   If mmaped==0, it is assumed that the data is allocated (using malloc).
 */
typedef struct
{
  void    *array;      /* Array keeping data elements.             */
  int       type;      /* Type of data (from `gal_data_alltypes'). */
  size_t    ndim;      /* Number of dimensions in the array.       */
  size_t size[GAL_DATA_MAXDIM]; /* Length along each dimension.    */
  int    mmapped;      /* ==1: not in physical RAM, it is mmap'd.  */
  char *mmapname;      /* File name of the mmap.                   */
  int   anyblank;      /* ==1: has blank values.                   */
  struct wcsprm *wcs;  /* WCS information for this dataset.        */
} gal_data;

__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_BOX_H__ */
