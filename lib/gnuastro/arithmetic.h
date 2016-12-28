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
#ifndef __GAL_ARITHMETIC_H__
#define __GAL_ARITHMETIC_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>

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





enum gal_arithmetic_operators
{
  GAL_ARITHMETIC_OP_PLUS,         /*   +     */
  GAL_ARITHMETIC_OP_MINUS,        /*   -     */
  GAL_ARITHMETIC_OP_MULTIPLY,     /*   *     */
  GAL_ARITHMETIC_OP_DIVIDE,       /*   /     */
  GAL_ARITHMETIC_OP_MODULO,       /*   %     */

  GAL_ARITHMETIC_OP_LT,           /*   <     */
  GAL_ARITHMETIC_OP_LE,           /*   <=    */
  GAL_ARITHMETIC_OP_GT,           /*   >     */
  GAL_ARITHMETIC_OP_GE,           /*   >=    */
  GAL_ARITHMETIC_OP_EQ,           /*   ==    */
  GAL_ARITHMETIC_OP_NE,           /*   !=    */
  GAL_ARITHMETIC_OP_AND,          /*   &&    */
  GAL_ARITHMETIC_OP_OR,           /*   ||    */
  GAL_ARITHMETIC_OP_NOT,          /*   !     */
  GAL_ARITHMETIC_OP_ISBLANK,      /* Similar to isnan() for floats. */
  GAL_ARITHMETIC_OP_WHERE,        /*   ?:    */

  GAL_ARITHMETIC_OP_BITAND,       /*   &     */
  GAL_ARITHMETIC_OP_BITOR,        /*   |     */
  GAL_ARITHMETIC_OP_BITXOR,       /*   ^     */
  GAL_ARITHMETIC_OP_BITLSH,       /*   <<    */
  GAL_ARITHMETIC_OP_BITRSH,       /*   >>    */
  GAL_ARITHMETIC_OP_BITNOT,       /*   ~     */

  GAL_ARITHMETIC_OP_ABS,          /* abs()   */
  GAL_ARITHMETIC_OP_POW,          /* pow()   */
  GAL_ARITHMETIC_OP_SQRT,         /* sqrt()  */
  GAL_ARITHMETIC_OP_LOG,          /* log()   */
  GAL_ARITHMETIC_OP_LOG10,        /* log10() */

  GAL_ARITHMETIC_OP_MINVAL,       /* Minimum value of array.               */
  GAL_ARITHMETIC_OP_MAXVAL,       /* Maximum value of array.               */
  GAL_ARITHMETIC_OP_MIN,          /* Minimum per pixel of multiple arrays. */
  GAL_ARITHMETIC_OP_MAX,          /* Maximum per pixel of multiple arrays. */
  GAL_ARITHMETIC_OP_SUM,          /* Sum per pixel of multiple arrays.     */
  GAL_ARITHMETIC_OP_AVERAGE,      /* Average per pixel of multiple arrays. */
  GAL_ARITHMETIC_OP_MEDIAN,       /* Median per pixel of multiple arrays.  */

  GAL_ARITHMETIC_OP_TO_UCHAR,     /* Convert to unsigned char.             */
  GAL_ARITHMETIC_OP_TO_CHAR,      /* Convert to char.                      */
  GAL_ARITHMETIC_OP_TO_USHORT,    /* Convert to unsigned short.            */
  GAL_ARITHMETIC_OP_TO_SHORT,     /* Convert to short.                     */
  GAL_ARITHMETIC_OP_TO_UINT,      /* Convert to unsigned int.              */
  GAL_ARITHMETIC_OP_TO_INT,       /* Convert to int.                       */
  GAL_ARITHMETIC_OP_TO_ULONG,     /* Convert to unsigned long.             */
  GAL_ARITHMETIC_OP_TO_LONG,      /* Convert to long.                      */
  GAL_ARITHMETIC_OP_TO_LONGLONG,  /* Convert to LONGLONG.                  */
  GAL_ARITHMETIC_OP_TO_FLOAT,     /* Convert to float.                     */
  GAL_ARITHMETIC_OP_TO_DOUBLE,    /* Convert to double.                    */
};





char *
gal_arithmetic_operator_string(int operator);

gal_data_t *
gal_arithmetic_convert_to_compiled_type(gal_data_t *in, unsigned char flags);

gal_data_t *
gal_arithmetic(int operator, unsigned char flags, ...);





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_ARITHMETIC_H__ */
