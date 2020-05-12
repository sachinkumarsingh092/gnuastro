/*********************************************************************
Arithmetic operations on data structures
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <gnuastro/list.h>
#include <gnuastro/blank.h>
#include <gnuastro/units.h>
#include <gnuastro/qsort.h>
#include <gnuastro/pointer.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>
#include <gnuastro/arithmetic.h>

#include <gnuastro-internal/arithmetic-internal.h>

/* Headers for each binary operator. Since they heavily involve macros,
   their compilation can be very large if they are in a single function and
   file. So there is a separate C source and header file for each of these
   functions.*/
#include <gnuastro-internal/arithmetic-lt.h>
#include <gnuastro-internal/arithmetic-le.h>
#include <gnuastro-internal/arithmetic-gt.h>
#include <gnuastro-internal/arithmetic-ge.h>
#include <gnuastro-internal/arithmetic-eq.h>
#include <gnuastro-internal/arithmetic-ne.h>
#include <gnuastro-internal/arithmetic-or.h>
#include <gnuastro-internal/arithmetic-and.h>
#include <gnuastro-internal/arithmetic-plus.h>
#include <gnuastro-internal/arithmetic-minus.h>
#include <gnuastro-internal/arithmetic-bitor.h>
#include <gnuastro-internal/arithmetic-bitand.h>
#include <gnuastro-internal/arithmetic-bitxor.h>
#include <gnuastro-internal/arithmetic-bitlsh.h>
#include <gnuastro-internal/arithmetic-bitrsh.h>
#include <gnuastro-internal/arithmetic-modulo.h>
#include <gnuastro-internal/arithmetic-divide.h>
#include <gnuastro-internal/arithmetic-multiply.h>










/***********************************************************************/
/***************             Internal checks              **************/
/***********************************************************************/
/* Some functions are only for a floating point operand, so if the input
   isn't floating point, inform the user to change the type explicitly,
   doing it implicitly/internally puts too much responsability on the
   program. */
static void
arithmetic_check_float_input(gal_data_t *in, int operator, char *numstr)
{
  switch(in->type)
    {
    case GAL_TYPE_FLOAT32:
    case GAL_TYPE_FLOAT64:
      break;
    default:
      error(EXIT_FAILURE, 0, "the %s operator can only accept single or "
            "double precision floating point numbers as its operand. The "
            "%s operand has type %s. You can use the 'float' or 'double' "
            "operators before this operator to explicitly convert to the "
            "desired precision floating point type. If the operand was "
            "originally a typed number (string of characters), add an 'f' "
            "after it so it is directly read into the proper precision "
            "floating point number (based on the number of non-zero "
            "decimals it has)", gal_arithmetic_operator_string(operator),
            numstr, gal_type_name(in->type, 1));
    }
}




















/***********************************************************************/
/***************        Unary functions/operators         **************/
/***********************************************************************/
/* Change input data structure type. */
static gal_data_t *
arithmetic_change_type(gal_data_t *data, int operator, int flags)
{
  int type=-1;
  gal_data_t *out;

  /* Set the output type. */
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_TO_UINT8:    type=GAL_TYPE_UINT8;    break;
    case GAL_ARITHMETIC_OP_TO_INT8:     type=GAL_TYPE_INT8;     break;
    case GAL_ARITHMETIC_OP_TO_UINT16:   type=GAL_TYPE_UINT16;   break;
    case GAL_ARITHMETIC_OP_TO_INT16:    type=GAL_TYPE_INT16;    break;
    case GAL_ARITHMETIC_OP_TO_UINT32:   type=GAL_TYPE_UINT32;   break;
    case GAL_ARITHMETIC_OP_TO_INT32:    type=GAL_TYPE_INT32;    break;
    case GAL_ARITHMETIC_OP_TO_UINT64:   type=GAL_TYPE_UINT64;   break;
    case GAL_ARITHMETIC_OP_TO_INT64:    type=GAL_TYPE_INT64;    break;
    case GAL_ARITHMETIC_OP_TO_FLOAT32:  type=GAL_TYPE_FLOAT32;  break;
    case GAL_ARITHMETIC_OP_TO_FLOAT64:  type=GAL_TYPE_FLOAT64;  break;
    default:
      error(EXIT_FAILURE, 0, "%s: operator value of %d not recognized",
            __func__, operator);
    }

  /* Copy to the new type. */
  out=gal_data_copy_to_new_type(data, type);

  /* Delete the input structure if the user asked for it. */
  if(flags & GAL_ARITHMETIC_FREE)
    gal_data_free(data);

  /* Return */
  return out;
}





/* Return an array of value 1 for any zero valued element and zero for any
   non-zero valued element. */
#define TYPE_CASE_FOR_NOT(CTYPE) {                                      \
    CTYPE *a=data->array, *af=a+data->size;                             \
    do *o++ = !(*a); while(++a<af);                                     \
  }

static gal_data_t *
arithmetic_not(gal_data_t *data, int flags)
{
  uint8_t *o;
  gal_data_t *out;

  /* Allocate the output array. */
  out=gal_data_alloc(NULL, GAL_TYPE_UINT8, data->ndim, data->dsize,
                     data->wcs, 0, data->minmapsize, data->quietmmap,
                     data->name, data->unit, data->comment);
  o=out->array;


  /* Go over the pixels and set the output values. */
  switch(data->type)
    {
    case GAL_TYPE_UINT8:   TYPE_CASE_FOR_NOT( uint8_t  );   break;
    case GAL_TYPE_INT8:    TYPE_CASE_FOR_NOT( int8_t   );   break;
    case GAL_TYPE_UINT16:  TYPE_CASE_FOR_NOT( uint16_t );   break;
    case GAL_TYPE_INT16:   TYPE_CASE_FOR_NOT( int16_t  );   break;
    case GAL_TYPE_UINT32:  TYPE_CASE_FOR_NOT( uint32_t );   break;
    case GAL_TYPE_INT32:   TYPE_CASE_FOR_NOT( int32_t  );   break;
    case GAL_TYPE_UINT64:  TYPE_CASE_FOR_NOT( uint64_t );   break;
    case GAL_TYPE_INT64:   TYPE_CASE_FOR_NOT( int64_t  );   break;
    case GAL_TYPE_FLOAT32: TYPE_CASE_FOR_NOT( float    );   break;
    case GAL_TYPE_FLOAT64: TYPE_CASE_FOR_NOT( double   );   break;

    case GAL_TYPE_BIT:
      error(EXIT_FAILURE, 0, "%s: bit datatypes are not yet supported, "
            "please get in touch with us to implement it.", __func__);

    default:
      error(EXIT_FAILURE, 0, "%s: type value (%d) not recognized",
            __func__, data->type);
    }

  /* Delete the input structure if the user asked for it. */
  if(flags & GAL_ARITHMETIC_FREE)
    gal_data_free(data);

  /* Return */
  return out;
}





/* Bitwise not operator. */
static gal_data_t *
arithmetic_bitwise_not(int flags, gal_data_t *in)
{
  gal_data_t *o;
  uint8_t    *iu8  = in->array,  *iu8f  = iu8  + in->size,   *ou8;
  int8_t     *ii8  = in->array,  *ii8f  = ii8  + in->size,   *oi8;
  uint16_t   *iu16 = in->array,  *iu16f = iu16 + in->size,   *ou16;
  int16_t    *ii16 = in->array,  *ii16f = ii16 + in->size,   *oi16;
  uint32_t   *iu32 = in->array,  *iu32f = iu32 + in->size,   *ou32;
  int32_t    *ii32 = in->array,  *ii32f = ii32 + in->size,   *oi32;
  uint64_t   *iu64 = in->array,  *iu64f = iu64 + in->size,   *ou64;
  int64_t    *ii64 = in->array,  *ii64f = ii64 + in->size,   *oi64;

  /* Check the type */
  switch(in->type)
    {
    case GAL_TYPE_FLOAT32:
    case GAL_TYPE_FLOAT64:
      error(EXIT_FAILURE, 0, "%s: bitwise not (one's complement) "
            "operator can only work on integer types", __func__);
    }

  /* If we want inplace output, set the output pointer to the input
     pointer, for every pixel, the operation will be independent. */
  if(flags & GAL_ARITHMETIC_INPLACE)
    o = in;
  else
    o = gal_data_alloc(NULL, in->type, in->ndim, in->dsize, in->wcs,
                       0, in->minmapsize, in->quietmmap, NULL, NULL, NULL);

  /* Start setting the types. */
  switch(in->type)
    {
    case GAL_TYPE_UINT8:
      ou8=o->array;   do  *ou8++ = ~(*iu8++);    while(iu8<iu8f);     break;

    case GAL_TYPE_INT8:
      oi8=o->array;   do  *oi8++ = ~(*ii8++);    while(ii8<ii8f);     break;

    case GAL_TYPE_UINT16:
      ou16=o->array;  do *ou16++ = ~(*iu16++);   while(iu16<iu16f);   break;

    case GAL_TYPE_INT16:
      oi16=o->array;  do *oi16++ = ~(*ii16++);   while(ii16<ii16f);   break;

    case GAL_TYPE_UINT32:
      ou32=o->array;  do *ou32++ = ~(*iu32++);   while(iu32<iu32f);   break;

    case GAL_TYPE_INT32:
      oi32=o->array;  do *oi32++ = ~(*ii32++);   while(ii32<ii32f);   break;

    case GAL_TYPE_UINT64:
      ou64=o->array;  do *ou64++ = ~(*iu64++);   while(iu64<iu64f);   break;

    case GAL_TYPE_INT64:
      oi64=o->array;  do *oi64++ = ~(*ii64++);   while(ii64<ii64f);   break;

    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, in->type);
    }


  /* Clean up (if necessary). */
  if( (flags & GAL_ARITHMETIC_FREE) && o!=in)
    gal_data_free(in);

  /* Return */
  return o;
}





/* We don't want to use the standard function for unary functions in the
   case of the absolute operator. This is because there are multiple
   versions of this function in the C library for different types, which
   can greatly improve speed. */
#define ARITHMETIC_ABS_SGN(CTYPE, FUNC) {                          \
    CTYPE *o=out->array, *a=in->array, *af=a+in->size;             \
    do *o++ = FUNC(*a); while(++a<af);                             \
  }
static gal_data_t *
arithmetic_abs(int flags, gal_data_t *in)
{
  gal_data_t *out;

  /* Set the output array. */
  if(flags & GAL_ARITHMETIC_INPLACE)
    out=in;
  else
    out = gal_data_alloc(NULL, in->type, in->ndim, in->dsize,
                         in->wcs, 0, in->minmapsize, in->quietmmap,
                         in->name, in->unit, in->comment);

  /* Put the absolute value depending on the type. */
  switch(in->type)
    {
    /* Unsigned types are already positive, so if the input is not to be
       freed (the output must be a separate array), just copy the whole
       array. */
    case GAL_TYPE_UINT8:
    case GAL_TYPE_UINT16:
    case GAL_TYPE_UINT32:
    case GAL_TYPE_UINT64:
      if(out!=in) gal_data_copy_to_allocated(in, out);
      break;

    /* For the signed types, we actually have to go over the data and
       calculate the absolute value. There are unique functions for
       different types, so we will be using them.*/
    case GAL_TYPE_INT8:    ARITHMETIC_ABS_SGN( int8_t,  abs   );  break;
    case GAL_TYPE_INT16:   ARITHMETIC_ABS_SGN( int16_t, abs   );  break;
    case GAL_TYPE_INT32:   ARITHMETIC_ABS_SGN( int32_t, labs  );  break;
    case GAL_TYPE_INT64:   ARITHMETIC_ABS_SGN( int64_t, llabs );  break;
    case GAL_TYPE_FLOAT32: ARITHMETIC_ABS_SGN( float,   fabsf );  break;
    case GAL_TYPE_FLOAT64: ARITHMETIC_ABS_SGN( double,  fabs  );  break;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, in->type);
    }

  /* Clean up and return */
  if( (flags & GAL_ARITHMETIC_FREE) && out!=in)
    gal_data_free(in);
  return out;
}





#define UNIFUNC_RUN_FUNCTION_ON_ELEMENT(OT, IT, OP){                    \
    OT *oa=o->array;                                                    \
    IT *ia=in->array, *iaf=ia + in->size;                               \
    do *oa++ = OP(*ia++); while(ia<iaf);                                \
  }

#define UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(IT, OP)                   \
  switch(o->type)                                                       \
    {                                                                   \
    case GAL_TYPE_UINT8:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint8_t,  IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_INT8:                                                 \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int8_t,   IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_UINT16:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint16_t, IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_INT16:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int16_t,  IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_UINT32:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint32_t, IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_INT32:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int32_t,  IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_UINT64:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint64_t, IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_INT64:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int64_t,  IT, OP)                 \
        break;                                                          \
    case GAL_TYPE_FLOAT32:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(float,    IT, OP)                 \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(double,   IT, OP)                 \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",         \
            "UNIARY_FUNCTION_ON_ELEMENT", in->type);                    \
    }

#define UNIARY_FUNCTION_ON_ELEMENT_OUTPUT_STRING(OP)                    \
  switch(in->type)                                                      \
    {                                                                   \
    case GAL_TYPE_UINT8:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, uint8_t,  OP)             \
        break;                                                          \
    case GAL_TYPE_INT8:                                                 \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, int8_t,   OP)             \
        break;                                                          \
    case GAL_TYPE_UINT16:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, uint16_t, OP)             \
        break;                                                          \
    case GAL_TYPE_INT16:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, int16_t,  OP)             \
        break;                                                          \
    case GAL_TYPE_UINT32:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, uint32_t, OP)             \
        break;                                                          \
    case GAL_TYPE_INT32:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, int32_t,  OP)             \
        break;                                                          \
    case GAL_TYPE_UINT64:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, uint64_t, OP)             \
        break;                                                          \
    case GAL_TYPE_INT64:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, int64_t,  OP)             \
        break;                                                          \
    case GAL_TYPE_FLOAT32:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, float,    OP)             \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char *, double,   OP)             \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",         \
            "UNIARY_FUNCTION_ON_ELEMENT_OUTPUT_STRING", in->type);      \
    }

#define UNIARY_FUNCTION_ON_ELEMENT(OP)                                  \
  switch(in->type)                                                      \
    {                                                                   \
    case GAL_TYPE_UINT8:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(uint8_t,  OP)               \
        break;                                                          \
    case GAL_TYPE_INT8:                                                 \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(int8_t,   OP)               \
        break;                                                          \
    case GAL_TYPE_UINT16:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(uint16_t, OP)               \
        break;                                                          \
    case GAL_TYPE_INT16:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(int16_t,  OP)               \
        break;                                                          \
    case GAL_TYPE_UINT32:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(uint32_t, OP)               \
        break;                                                          \
    case GAL_TYPE_INT32:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(int32_t,  OP)               \
        break;                                                          \
    case GAL_TYPE_UINT64:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(uint64_t, OP)               \
        break;                                                          \
    case GAL_TYPE_INT64:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(int64_t,  OP)               \
        break;                                                          \
    case GAL_TYPE_FLOAT32:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(float,    OP)               \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_INSET(double,   OP)               \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",         \
            "UNIARY_FUNCTION_ON_ELEMENT", in->type);                    \
    }


#define UNIFUNC_RUN_FUNCTION_ON_ELEMENT_STRING(OT, OP){                 \
    OT *oa=o->array;                                                    \
    char **ia=in->array, **iaf=ia + in->size;                           \
    do *oa++ = OP(*ia++); while(ia<iaf);                                \
}

static gal_data_t *
arithmetic_unary_function(int operator, int flags, gal_data_t *in)
{
  uint8_t otype;
  int inplace=0;
  gal_data_t *o;

  /* See if the operation should be done in place. Note that so far, the
     output of these operators is defined in the real space (floating
     point). So even if the user requested inplace opereation, if its not a
     floating point type, its not useful.*/
  if( (flags & GAL_ARITHMETIC_INPLACE)
      && (in->type==GAL_TYPE_FLOAT32 || in->type==GAL_TYPE_FLOAT64)
      && (operator != GAL_ARITHMETIC_OP_RA_TO_DEGREE
      &&  operator != GAL_ARITHMETIC_OP_DEC_TO_DEGREE
      &&  operator != GAL_ARITHMETIC_OP_DEGREE_TO_RA
      &&  operator != GAL_ARITHMETIC_OP_DEGREE_TO_DEC ) )
    inplace=1;

  if(inplace)
    {
      o = in;
      otype=in->type;
    }
  else
    {
      otype = ( in->type==GAL_TYPE_FLOAT64
                ? GAL_TYPE_FLOAT64
                : GAL_TYPE_FLOAT32 );

      /* Check for operators which have fixed output types */
      if ( operator == GAL_ARITHMETIC_OP_RA_TO_DEGREE ||
           operator == GAL_ARITHMETIC_OP_DEC_TO_DEGREE )
        otype = GAL_TYPE_FLOAT64;

      if (operator == GAL_ARITHMETIC_OP_DEGREE_TO_RA ||
          operator == GAL_ARITHMETIC_OP_DEGREE_TO_DEC)
        otype = GAL_TYPE_STRING;

      o = gal_data_alloc(NULL, otype, in->ndim, in->dsize, in->wcs,
                         0, in->minmapsize, in->quietmmap,
                         NULL, NULL, NULL);
    }

  /* Start setting the operator and operands. */
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_SQRT:
      UNIARY_FUNCTION_ON_ELEMENT( sqrt );
      break;

    case GAL_ARITHMETIC_OP_LOG:
      UNIARY_FUNCTION_ON_ELEMENT( log );
      break;

    case GAL_ARITHMETIC_OP_LOG10:
      UNIARY_FUNCTION_ON_ELEMENT( log10 );
      break;

    case GAL_ARITHMETIC_OP_RA_TO_DEGREE:
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_STRING(double, gal_units_ra_to_degree);
      break;

    case GAL_ARITHMETIC_OP_DEC_TO_DEGREE:
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT_STRING(double, gal_units_dec_to_degree);
      break;

    case GAL_ARITHMETIC_OP_DEGREE_TO_RA:
      UNIARY_FUNCTION_ON_ELEMENT_OUTPUT_STRING(gal_units_degree_to_ra);
      break;

    case GAL_ARITHMETIC_OP_DEGREE_TO_DEC:
      UNIARY_FUNCTION_ON_ELEMENT_OUTPUT_STRING(gal_units_degree_to_dec);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",
            __func__, operator);
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, 'UNIFUNC_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, and we only have 'r' and 'l'
     allocated in any case. Alternatively, when the inputs shouldn't be
     freed, the only allocated spaces are the 'r' and 'l' arrays if their
     types weren't compiled for binary operations, we can tell this from
     the pointers: if they are different from the original pointers, they
     were allocated. */
  if( (flags & GAL_ARITHMETIC_FREE) && o!=in)
    gal_data_free(in);

  /* Return */
  return o;
}





/* Call functions in the 'gnuastro/statistics' library. */
static gal_data_t *
arithmetic_from_statistics(int operator, int flags, gal_data_t *input)
{
  gal_data_t *out=NULL;
  int ip=(flags & GAL_ARITHMETIC_INPLACE) || (flags & GAL_ARITHMETIC_FREE);

  switch(operator)
    {
    case GAL_ARITHMETIC_OP_MINVAL:   out=gal_statistics_minimum(input);break;
    case GAL_ARITHMETIC_OP_MAXVAL:   out=gal_statistics_maximum(input);break;
    case GAL_ARITHMETIC_OP_NUMBERVAL:out=gal_statistics_number(input); break;
    case GAL_ARITHMETIC_OP_SUMVAL:   out=gal_statistics_sum(input);    break;
    case GAL_ARITHMETIC_OP_MEANVAL:  out=gal_statistics_mean(input);   break;
    case GAL_ARITHMETIC_OP_STDVAL:   out=gal_statistics_std(input);    break;
    case GAL_ARITHMETIC_OP_MEDIANVAL:
      out=gal_statistics_median(input, ip); break;
    default:
      error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",
            __func__, operator);
    }

  /* If the input is to be freed, then do so and return the output. */
  if( flags & GAL_ARITHMETIC_FREE ) gal_data_free(input);
  return out;
}




















/***********************************************************************/
/***************                  Metadata                **************/
/***********************************************************************/

/* The size operator. Reports the size along a given dimension. */
static gal_data_t *
arithmetic_size(int operator, int flags, gal_data_t *in, gal_data_t *arg)
{
  size_t one=1, arg_val;
  gal_data_t *usearg=NULL, *out=NULL;

  /* Sanity checks on argument (dimension number): it should be an integer,
     and have a size of 1. */
  if(arg->type==GAL_TYPE_FLOAT32 || arg->type==GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: size operator's dimension argument"
          "must have an integer type", __func__);
  if(arg->size!=1)
    error(EXIT_FAILURE, 0, "%s: size operator's dimension argument"
          "must be a single number, but it has %zu elements", __func__,
          arg->size);


  /* Convert 'arg' to 'size_t' and read it. Note that we can only free the
     'arg' array (while changing its type), when the freeing flag has been
     set. */
  if(flags & GAL_ARITHMETIC_FREE)
    {
      arg=gal_data_copy_to_new_type_free(arg, GAL_TYPE_SIZE_T);
      arg_val=*(size_t *)(arg->array);
      gal_data_free(arg);
    }
  else
    {
      usearg=gal_data_copy_to_new_type(arg, GAL_TYPE_SIZE_T);
      arg_val=*(size_t *)(usearg->array);
      gal_data_free(usearg);
    }


  /* Sanity checks on the value of the given argument.*/
  if(arg_val>in->ndim)
    error(EXIT_FAILURE, 0, "%s: size operator's dimension argument "
          "(given %zu) cannot be larger than the dimensions of the "
          "given input (%zu)", __func__, arg_val, in->ndim);
  if(arg_val==0)
    error(EXIT_FAILURE, 0, "%s: size operator's dimension argument "
          "(given %zu) cannot be zero: dimensions are counted from 1",
          __func__, arg_val);


  /* Allocate the output array and write the desired dimension. Note that
     'dsize' is in the C order, while the output must be in FITS/Fortran
     order. Also that C order starts from 0, while the FITS order starts
     from 1. */
  out=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &one, NULL, 0,
                     in->minmapsize, 0, NULL, NULL, NULL);
  *(size_t *)(out->array)=in->dsize[in->ndim-arg_val];


  /* Clean up and return */
  if(flags & GAL_ARITHMETIC_FREE)
    gal_data_free(in);
  return out;
}




















/***********************************************************************/
/***************                  Where                   **************/
/***********************************************************************/
/* When the 'iftrue' dataset only has one element and the element is blank,
   then it will be replaced with the blank value of the type of the output
   data. */
#define DO_WHERE_OPERATION(ITT, OT) {                                   \
    ITT *it=iftrue->array;                                              \
    OT b, *o=out->array, *of=o+out->size;                               \
    gal_blank_write(&b, out->type);                                     \
    if(iftrue->size==1)                                                 \
      {                                                                 \
        if( gal_blank_is(it, iftrue->type) )                            \
          {                                                             \
            do{*o = (chb && *c==cb) ? *o : (*c ? b   : *o); ++c;      } \
            while(++o<of);                                              \
          }                                                             \
        else                                                            \
          do  {*o = (chb && *c==cb) ? *o : (*c ? *it : *o); ++c;      } \
          while(++o<of);                                                \
      }                                                                 \
    else                                                                \
      do      {*o = (chb && *c==cb) ? *o : (*c ? *it : *o); ++it; ++c;} \
      while(++o<of);                                                    \
}





#define WHERE_OUT_SET(OT)                                               \
  switch(iftrue->type)                                                  \
    {                                                                   \
    case GAL_TYPE_UINT8:    DO_WHERE_OPERATION( uint8_t,  OT);  break;  \
    case GAL_TYPE_INT8:     DO_WHERE_OPERATION( int8_t,   OT);  break;  \
    case GAL_TYPE_UINT16:   DO_WHERE_OPERATION( uint16_t, OT);  break;  \
    case GAL_TYPE_INT16:    DO_WHERE_OPERATION( int16_t,  OT);  break;  \
    case GAL_TYPE_UINT32:   DO_WHERE_OPERATION( uint32_t, OT);  break;  \
    case GAL_TYPE_INT32:    DO_WHERE_OPERATION( int32_t,  OT);  break;  \
    case GAL_TYPE_UINT64:   DO_WHERE_OPERATION( uint64_t, OT);  break;  \
    case GAL_TYPE_INT64:    DO_WHERE_OPERATION( int64_t,  OT);  break;  \
    case GAL_TYPE_FLOAT32:  DO_WHERE_OPERATION( float,    OT);  break;  \
    case GAL_TYPE_FLOAT64:  DO_WHERE_OPERATION( double,   OT);  break;  \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized for the " \
            "'iftrue' dataset", "WHERE_OUT_SET", iftrue->type);         \
    }





static void
arithmetic_where(int flags, gal_data_t *out, gal_data_t *cond,
                 gal_data_t *iftrue)
{
  int chb;    /* Read as: "Condition-Has-Blank" */
  unsigned char *c=cond->array, cb=GAL_BLANK_UINT8;

  /* The condition operator has to be unsigned char. */
  if(cond->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "%s: the condition operand must be an "
          "'uint8' type, but the given condition operand has a "
          "'%s' type", __func__, gal_type_name(cond->type, 1));

  /* The dimension and sizes of the out and condition data sets must be the
     same. */
  if( gal_dimension_is_different(out, cond) )
    error(EXIT_FAILURE, 0, "%s: the output and condition datasets "
          "must have the same size", __func__);

  /* See if the condition array has blank values. */
  chb=gal_blank_present(cond, 0);

  /* Do the operation. */
  switch(out->type)
    {
    case GAL_TYPE_UINT8:         WHERE_OUT_SET( uint8_t  );      break;
    case GAL_TYPE_INT8:          WHERE_OUT_SET( int8_t   );      break;
    case GAL_TYPE_UINT16:        WHERE_OUT_SET( uint16_t );      break;
    case GAL_TYPE_INT16:         WHERE_OUT_SET( int16_t  );      break;
    case GAL_TYPE_UINT32:        WHERE_OUT_SET( uint32_t );      break;
    case GAL_TYPE_INT32:         WHERE_OUT_SET( int32_t  );      break;
    case GAL_TYPE_UINT64:        WHERE_OUT_SET( uint64_t );      break;
    case GAL_TYPE_INT64:         WHERE_OUT_SET( int64_t  );      break;
    case GAL_TYPE_FLOAT32:       WHERE_OUT_SET( float    );      break;
    case GAL_TYPE_FLOAT64:       WHERE_OUT_SET( double   );      break;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized for the 'out'",
            __func__, out->type);
    }

  /* Clean up if necessary. */
  if(flags & GAL_ARITHMETIC_FREE)
    {
      gal_data_free(cond);
      gal_data_free(iftrue);
    }
}




















/***********************************************************************/
/***************        Multiple operand operators        **************/
/***********************************************************************/
struct multioperandparams
{
  gal_data_t      *list;        /* List of input datasets.           */
  gal_data_t       *out;        /* Output dataset.                   */
  size_t           dnum;        /* Number of input dataset.          */
  int          operator;        /* Operator to use.                  */
  uint8_t     *hasblank;        /* Array of 0s or 1s for each input. */
  float              p1;        /* Sigma-cliping parameter 1.        */
  float              p2;        /* Sigma-cliping parameter 2.        */
};





#define MULTIOPERAND_MIN(TYPE) {                                        \
    size_t n, j=0;                                                      \
    TYPE t, max, *o=p->out->array;                                      \
    gal_type_max(p->list->type, &max);                                  \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        t=max;                                                          \
        j=tprm->indexs[tind];                                           \
                                                                        \
        for(i=0;i<p->dnum;++i)  /* Loop over each array. */             \
          {   /* Only for integer types, b==b. */                       \
            if( p->hasblank[i] && b==b)                                 \
              {                                                         \
                if( a[i][j] != b )                                      \
                  { t = a[i][j] < t ? a[i][j] : t; ++n; }               \
              }                                                         \
            else { t = a[i][j] < t ? a[i][j] : t; ++n; }                \
          }                                                             \
        o[j] = n ? t : b;  /* No usable elements: set to blank. */      \
      }                                                                 \
  }





#define MULTIOPERAND_MAX(TYPE) {                                        \
    size_t n, j=0;                                                      \
    TYPE t, min, *o=p->out->array;                                      \
    gal_type_min(p->list->type, &min);                                  \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        t=min;                                                          \
        j=tprm->indexs[tind];                                           \
                                                                        \
        for(i=0;i<p->dnum;++i)  /* Loop over each array. */             \
          {   /* Only for integer types, b==b. */                       \
            if( p->hasblank[i] && b==b)                                 \
              {                                                         \
                if( a[i][j] != b )                                      \
                  { t = a[i][j] > t ? a[i][j] : t; ++n; }               \
              }                                                         \
            else { t = a[i][j] > t ? a[i][j] : t; ++n; }                \
          }                                                             \
        o[j] = n ? t : b;  /* No usable elements: set to blank. */      \
      }                                                                 \
  }





#define MULTIOPERAND_NUM {                                              \
    int use;                                                            \
    size_t n, j;                                                        \
    uint32_t *o=p->out->array;                                          \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        j=tprm->indexs[tind];                                           \
                                                                        \
        for(i=0;i<p->dnum;++i)  /* Loop over each array. */             \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(p->hasblank[i])                                          \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b       ? 1 : 0 )      /* Integer */ \
                      : ( a[i][j]==a[i][j] ? 1 : 0 ) );   /* Float   */ \
            else use=1;                                                 \
                                                                        \
            /* Increment counter if necessary. */                       \
            if(use) ++n;                                                \
          }                                                             \
        o[j] = n;                                                       \
      }                                                                 \
  }





#define MULTIOPERAND_SUM {                                              \
    int use;                                                            \
    double sum;                                                         \
    size_t n, j;                                                        \
    float *o=p->out->array;                                             \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        sum=0.0f;                                                       \
        j=tprm->indexs[tind];                                           \
                                                                        \
        for(i=0;i<p->dnum;++i)  /* Loop over each array. */             \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(p->hasblank[i])                                          \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b     ? 1 : 0 )       /* Integer */  \
                      : ( a[i][j]==*a[i] ? 1 : 0 ) );    /* Float   */  \
            else use=1;                                                 \
                                                                        \
            /* Use in sum if necessary. */                              \
            if(use) { sum += a[i][j]; ++n; }                            \
          }                                                             \
        o[j] = n ? sum : b;                                             \
      }                                                                 \
  }





#define MULTIOPERAND_MEAN {                                             \
    int use;                                                            \
    double sum;                                                         \
    size_t n, j;                                                        \
    float *o=p->out->array;                                             \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        sum=0.0f;                                                       \
        j=tprm->indexs[tind];                                           \
                                                                        \
        for(i=0;i<p->dnum;++i)  /* Loop over each array. */             \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(p->hasblank[i])                                          \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b       ? 1 : 0 )     /* Integer */  \
                      : ( a[i][j]==a[i][j] ? 1 : 0 ) );  /* Float   */  \
            else use=1;                                                 \
                                                                        \
            /* Calculate the mean if necessary. */                      \
            if(use) { sum += a[i][j]; ++n; }                            \
          }                                                             \
        o[j] = n ? sum/n : b;                                           \
      }                                                                 \
  }





#define MULTIOPERAND_STD {                                              \
    int use;                                                            \
    size_t n, j;                                                        \
    double sum, sum2;                                                   \
    float *o=p->out->array;                                             \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        sum=sum2=0.0f;                                                  \
        j=tprm->indexs[tind];                                           \
                                                                        \
        for(i=0;i<p->dnum;++i)  /* Loop over each array. */             \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(p->hasblank[i])                                          \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b       ? 1 : 0 )     /* Integer */  \
                      : ( a[i][j]==a[i][j] ? 1 : 0 ) );  /* Float   */  \
            else use=1;                                                 \
                                                                        \
            /* Calculate the necessary parameters if necessary. */      \
            if(use)                                                     \
              {                                                         \
                sum2 += a[i][j] * a[i][j];                              \
                sum  += a[i][j];                                        \
                ++n;                                                    \
              }                                                         \
          }                                                             \
        o[j] = n ? sqrt( (sum2-sum*sum/n)/n ) : b;                      \
      }                                                                 \
  }





#define MULTIOPERAND_MEDIAN(TYPE, QSORT_F) {                            \
    int use;                                                            \
    size_t n, j;                                                        \
    float *o=p->out->array;                                             \
    TYPE *pixs=gal_pointer_allocate(p->list->type, p->dnum, 0,          \
                                    __func__, "pixs");                  \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        j=tprm->indexs[tind];                                           \
                                                                        \
        /* Loop over each array: 'i' is input dataset's index. */       \
        for(i=0;i<p->dnum;++i)                                          \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(p->hasblank[i])                                          \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b       ? 1 : 0 )     /* Integer */  \
                      : ( a[i][j]==a[i][j] ? 1 : 0 ) );  /* Float   */  \
            else use=1;                                                 \
                                                                        \
            /* Put the value into the array of values. */               \
            if(use) pixs[n++]=a[i][j];                                  \
          }                                                             \
                                                                        \
        /* Sort all the values for this pixel and return the median. */ \
        if(n)                                                           \
          {                                                             \
            qsort(pixs, n, sizeof *pixs, QSORT_F);                      \
            o[j] = n%2 ? pixs[n/2] : (pixs[n/2] + pixs[n/2-1])/2 ;      \
          }                                                             \
        else                                                            \
          o[j]=b;                                                       \
      }                                                                 \
                                                                        \
    /* Clean up. */                                                     \
    free(pixs);                                                         \
  }





#define MULTIOPERAND_QUANTILE(TYPE) {                                   \
    size_t n, j;                                                        \
    gal_data_t *quantile;                                               \
    TYPE *o=p->out->array;                                              \
    TYPE *pixs=gal_pointer_allocate(p->list->type, p->dnum, 0,          \
                                    __func__, "pixs");                  \
    gal_data_t *cont=gal_data_alloc(pixs, p->list->type, 1, &p->dnum,   \
                                    NULL, 0, -1, 1, NULL, NULL, NULL);  \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        j=tprm->indexs[tind];                                           \
                                                                        \
        /* Read the necessay values from each input. */                 \
        for(i=0;i<p->dnum;++i) pixs[n++]=a[i][j];                       \
                                                                        \
        /* If there are any elements, measure the  */                   \
        if(n)                                                           \
          {                                                             \
            /* Calculate the quantile and put it in the output. */      \
            quantile=gal_statistics_quantile(cont, p->p1, 1);           \
            memcpy(&o[j], quantile->array,                              \
                   gal_type_sizeof(p->list->type));                     \
            gal_data_free(quantile);                                    \
                                                                        \
            /* Since we are doing sigma-clipping in place, the size, */ \
            /* and flags need to be reset. */                           \
            cont->flag=0;                                               \
            cont->size=cont->dsize[0]=p->dnum;                          \
          }                                                             \
        else                                                            \
          o[j]=b;                                                       \
      }                                                                 \
                                                                        \
    /* Clean up. */                                                     \
    gal_data_free(cont);                                                \
  }





#define MULTIOPERAND_SIGCLIP(TYPE) {                                    \
    size_t n, j;                                                        \
    gal_data_t *sclip;                                                  \
    uint32_t *N=p->out->array;                                          \
    float *sarr, *o=p->out->array;                                      \
    TYPE *pixs=gal_pointer_allocate(p->list->type, p->dnum, 0,          \
                                    __func__, "pixs");                  \
    gal_data_t *cont=gal_data_alloc(pixs, p->list->type, 1, &p->dnum,   \
                                    NULL, 0, -1, 1, NULL, NULL, NULL);  \
                                                                        \
    /* Go over all the pixels assigned to this thread. */               \
    for(tind=0; tprm->indexs[tind] != GAL_BLANK_SIZE_T; ++tind)         \
      {                                                                 \
        /* Initialize, 'j' is desired pixel's index. */                 \
        n=0;                                                            \
        j=tprm->indexs[tind];                                           \
                                                                        \
        /* Read the necessay values from each input. */                 \
        for(i=0;i<p->dnum;++i) pixs[n++]=a[i][j];                       \
                                                                        \
        /* If there are any elements, measure the  */                   \
        if(n)                                                           \
          {                                                             \
            /* Calculate the sigma-clip and write it in. */             \
            sclip=gal_statistics_sigma_clip(cont, p->p1, p->p2, 1, 1);  \
            sarr=sclip->array;                                          \
            switch(p->operator)                                         \
              {                                                         \
              case GAL_ARITHMETIC_OP_SIGCLIP_STD:    o[j]=sarr[3]; break;\
              case GAL_ARITHMETIC_OP_SIGCLIP_MEAN:   o[j]=sarr[2]; break;\
              case GAL_ARITHMETIC_OP_SIGCLIP_MEDIAN: o[j]=sarr[1]; break;\
              case GAL_ARITHMETIC_OP_SIGCLIP_NUMBER: N[j]=sarr[0]; break;\
              default:                                                  \
                error(EXIT_FAILURE, 0, "%s: a bug! the code %d is not " \
                      "valid for sigma-clipping results", __func__,     \
                      p->operator);                                     \
              }                                                         \
            gal_data_free(sclip);                                       \
                                                                        \
            /* Since we are doing sigma-clipping in place, the size, */ \
            /* and flags need to be reset. */                           \
            cont->flag=0;                                               \
            cont->size=cont->dsize[0]=p->dnum;                          \
          }                                                             \
        else                                                            \
          o[j]=b;                                                       \
      }                                                                 \
                                                                        \
    /* Clean up. */                                                     \
    gal_data_free(cont);                                                \
  }





#define MULTIOPERAND_TYPE_SET(TYPE, QSORT_F) {                          \
    TYPE b, **a;                                                        \
    gal_data_t *tmp;                                                    \
    size_t i=0, tind;                                                   \
                                                                        \
    /* Allocate space to keep the pointers to the arrays of each. */    \
    /* Input data structure. The operators will increment these */      \
    /* pointers while parsing them. */                                  \
    errno=0;                                                            \
    a=malloc(p->dnum*sizeof *a);                                        \
    if(a==NULL)                                                         \
      error(EXIT_FAILURE, 0, "%s: %zu bytes for 'a'",                   \
            "MULTIOPERAND_TYPE_SET", p->dnum*sizeof *a);                \
                                                                        \
    /* Fill in the array pointers and the blank value for this type. */ \
    gal_blank_write(&b, p->list->type);                                 \
    for(tmp=p->list;tmp!=NULL;tmp=tmp->next)                            \
      a[i++]=tmp->array;                                                \
                                                                        \
    /* Do the operation. */                                             \
    switch(p->operator)                                                 \
      {                                                                 \
      case GAL_ARITHMETIC_OP_MIN:                                       \
        MULTIOPERAND_MIN(TYPE);                                         \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_MAX:                                       \
        MULTIOPERAND_MAX(TYPE);                                         \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_NUMBER:                                    \
        MULTIOPERAND_NUM;                                               \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_SUM:                                       \
        MULTIOPERAND_SUM;                                               \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_MEAN:                                      \
        MULTIOPERAND_MEAN;                                              \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_STD:                                       \
        MULTIOPERAND_STD;                                               \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_MEDIAN:                                    \
        MULTIOPERAND_MEDIAN(TYPE, QSORT_F);                             \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_QUANTILE:                                  \
        MULTIOPERAND_QUANTILE(TYPE);                                    \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_SIGCLIP_STD:                               \
      case GAL_ARITHMETIC_OP_SIGCLIP_MEAN:                              \
      case GAL_ARITHMETIC_OP_SIGCLIP_MEDIAN:                            \
      case GAL_ARITHMETIC_OP_SIGCLIP_NUMBER:                            \
        MULTIOPERAND_SIGCLIP(TYPE);                                     \
        break;                                                          \
                                                                        \
      default:                                                          \
        error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",   \
              "MULTIOPERAND_TYPE_SET", p->operator);                    \
      }                                                                 \
                                                                        \
    /* Clean up. */                                                     \
    free(a);                                                            \
  }





/* Worker function on each thread. */
void *
multioperand_on_thread(void *in_prm)
{
  /* Low-level definitions to be done first. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct multioperandparams *p=(struct multioperandparams *)tprm->params;

  /* Do the operation on each thread. */
  switch(p->list->type)
    {
    case GAL_TYPE_UINT8:
      MULTIOPERAND_TYPE_SET(uint8_t,   gal_qsort_uint8_i);
      break;
    case GAL_TYPE_INT8:
      MULTIOPERAND_TYPE_SET(int8_t,    gal_qsort_int8_i);
      break;
    case GAL_TYPE_UINT16:
      MULTIOPERAND_TYPE_SET(uint16_t,  gal_qsort_uint16_i);
      break;
    case GAL_TYPE_INT16:
      MULTIOPERAND_TYPE_SET(int16_t,   gal_qsort_int16_i);
      break;
    case GAL_TYPE_UINT32:
      MULTIOPERAND_TYPE_SET(uint32_t,  gal_qsort_uint32_i);
      break;
    case GAL_TYPE_INT32:
      MULTIOPERAND_TYPE_SET(int32_t,   gal_qsort_int32_i);
      break;
    case GAL_TYPE_UINT64:
      MULTIOPERAND_TYPE_SET(uint64_t,  gal_qsort_uint64_i);
      break;
    case GAL_TYPE_INT64:
      MULTIOPERAND_TYPE_SET(int64_t,   gal_qsort_int64_i);
      break;
    case GAL_TYPE_FLOAT32:
      MULTIOPERAND_TYPE_SET(float,     gal_qsort_float32_i);
      break;
    case GAL_TYPE_FLOAT64:
      MULTIOPERAND_TYPE_SET(double,    gal_qsort_float64_i);
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, p->list->type);
    }

  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* The single operator in this function is assumed to be a linked list. The
   number of operators is determined from the fact that the last node in
   the linked list must have a NULL pointer as its 'next' element. */
static gal_data_t *
arithmetic_multioperand(int operator, int flags, gal_data_t *list,
                        gal_data_t *params, size_t numthreads)
{
  size_t i=0, dnum=1;
  float p1=NAN, p2=NAN;
  uint8_t *hasblank, otype;
  struct multioperandparams p;
  gal_data_t *out, *tmp, *ttmp;


  /* For generality, 'list' can be a NULL pointer, in that case, this
     function will return a NULL pointer and avoid further processing. */
  if(list==NULL) return NULL;


  /* If any parameters are given, prepare them. */
  for(tmp=params; tmp!=NULL; tmp=tmp->next)
    {
      /* Basic sanity checks. */
      if(tmp->size>1)
        error(EXIT_FAILURE, 0, "%s: parameters must be a single number",
              __func__);
      if(tmp->type!=GAL_TYPE_FLOAT32)
        error(EXIT_FAILURE, 0, "%s: parameters must be float32 type",
              __func__);

      /* Write them */
      if(isnan(p1)) p1=((float *)(tmp->array))[0];
      else          p2=((float *)(tmp->array))[0];

      /* Operator specific, parameter sanity checks. */
      switch(operator)
        {
        case GAL_ARITHMETIC_OP_QUANTILE:
          if(p1<0 || p1>1)
            error(EXIT_FAILURE, 0, "%s: the parameter given to the 'quantile' "
                  "operator must be between (and including) 0 and 1. The "
                  "given value is: %g", __func__, p1);
          break;
        }
    }


  /* Do a simple sanity check, comparing the operand on top of the list to
     the rest of the operands within the list. */
  for(tmp=list->next;tmp!=NULL;tmp=tmp->next)
    {
      /* Increment the number of structures. */
      ++dnum;

      /* Check the types. */
      if(tmp->type!=list->type)
        error(EXIT_FAILURE, 0, "%s: the types of all operands to the %s "
              "operator must be same", __func__,
              gal_arithmetic_operator_string(operator));

      /* Check the sizes. */
      if( gal_dimension_is_different(list, tmp) )
        error(EXIT_FAILURE, 0, "%s: the sizes of all operands to the %s "
              "operator must be same", __func__,
              gal_arithmetic_operator_string(operator));
    }


  /* Set the output dataset type */
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_MIN:            otype=list->type;       break;
    case GAL_ARITHMETIC_OP_MAX:            otype=list->type;       break;
    case GAL_ARITHMETIC_OP_NUMBER:         otype=GAL_TYPE_UINT32;  break;
    case GAL_ARITHMETIC_OP_SUM:            otype=GAL_TYPE_FLOAT32; break;
    case GAL_ARITHMETIC_OP_MEAN:           otype=GAL_TYPE_FLOAT32; break;
    case GAL_ARITHMETIC_OP_STD:            otype=GAL_TYPE_FLOAT32; break;
    case GAL_ARITHMETIC_OP_MEDIAN:         otype=GAL_TYPE_FLOAT32; break;
    case GAL_ARITHMETIC_OP_QUANTILE:       otype=list->type;       break;
    case GAL_ARITHMETIC_OP_SIGCLIP_STD:    otype=GAL_TYPE_FLOAT32; break;
    case GAL_ARITHMETIC_OP_SIGCLIP_MEAN:   otype=GAL_TYPE_FLOAT32; break;
    case GAL_ARITHMETIC_OP_SIGCLIP_MEDIAN: otype=GAL_TYPE_FLOAT32; break;
    case GAL_ARITHMETIC_OP_SIGCLIP_NUMBER: otype=GAL_TYPE_UINT32;  break;
    default:
      error(EXIT_FAILURE, 0, "%s: operator code %d isn't recognized",
            __func__, operator);
    }


  /* Set the output data structure. */
  if( (flags & GAL_ARITHMETIC_INPLACE) && otype==list->type)
    out = list;                 /* The top element in the list. */
  else
    out = gal_data_alloc(NULL, otype, list->ndim, list->dsize,
                         list->wcs, 0, list->minmapsize, list->quietmmap,
                         NULL, NULL, NULL);


  /* hasblank is used to see if a blank value should be checked for each
     list element or not. */
  hasblank=gal_pointer_allocate(GAL_TYPE_UINT8, dnum, 0, __func__,
                                "hasblank");
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    hasblank[i++]=gal_blank_present(tmp, 0);


  /* Set the parameters necessary for multithreaded operation and spin them
     off to do apply the operator. */
  p.p1=p1;
  p.p2=p2;
  p.out=out;
  p.list=list;
  p.dnum=dnum;
  p.operator=operator;
  p.hasblank=hasblank;
  gal_threads_spin_off(multioperand_on_thread, &p, out->size, numthreads);


  /* Clean up and return. Note that the operation might have been done in
     place. In that case, the top most list element was used. So we need to
     check before freeing each data structure. */
  if(flags & GAL_ARITHMETIC_FREE)
    {
      tmp=list;
      while(tmp!=NULL)
        {
          ttmp=tmp->next;
          if(tmp!=out) gal_data_free(tmp);
          tmp=ttmp;
        }
      if(params) gal_list_data_free(params);
    }
  free(hasblank);
  return out;
}




















/**********************************************************************/
/****************           Binary operators          *****************/
/**********************************************************************/
/* See if we should check for blanks. When both types are floats, blanks
   don't need to be checked (the floating point standard will do the job
   for us). It is also not necessary to check blanks in bitwise operators,
   but bitwise operators have their own macro
   ('BINARY_OP_INCR_OT_RT_LT_SET') which doesn' use 'checkblanks'.*/
int
gal_arithmetic_binary_checkblank(gal_data_t *l, gal_data_t *r)
{
  return ((((l->type!=GAL_TYPE_FLOAT32    && l->type!=GAL_TYPE_FLOAT64)
            || (r->type!=GAL_TYPE_FLOAT32 && r->type!=GAL_TYPE_FLOAT64))
           && (gal_blank_present(l, 1) || gal_blank_present(r, 1)))
          ? 1 : 0 );
}





static int
arithmetic_binary_out_type(int operator, gal_data_t *l, gal_data_t *r)
{
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_PLUS:
    case GAL_ARITHMETIC_OP_MINUS:
    case GAL_ARITHMETIC_OP_MULTIPLY:
    case GAL_ARITHMETIC_OP_DIVIDE:
      return gal_type_out(l->type, r->type);

    default:
      return GAL_TYPE_UINT8;
    }
  return -1;
}





static gal_data_t *
arithmetic_binary(int operator, int flags, gal_data_t *l, gal_data_t *r)
{
  /* Read the variable arguments. 'lo' and 'ro' keep the original data, in
     case their type isn't built (based on configure options are configure
     time). */
  int32_t otype;
  gal_data_t *o=NULL;
  size_t out_size, minmapsize;
  int quietmmap=l->quietmmap && r->quietmmap;


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_ARITHMETIC_NUMOK) && (l->size==1 || r->size==1))
      && gal_dimension_is_different(l, r) )
    error(EXIT_FAILURE, 0, "%s: the non-number inputs to %s don't have the "
          "same dimension/size", __func__,
          gal_arithmetic_operator_string(operator));


  /* Set the output type. For the comparison operators, the output type is
     either 0 or 1, so we will set the output type to 'unsigned char' for
     efficient memory and CPU usage. Since the number of operators without
     a fixed output type (like the conditionals) is less, by 'default' we
     will set the output type to 'unsigned char', and if any of the other
     operatrs are given, it will be chosen based on the input types.*/
  otype=arithmetic_binary_out_type(operator, l, r);


  /* Set the output sizes. */
  minmapsize = ( l->minmapsize < r->minmapsize
                 ? l->minmapsize : r->minmapsize );
  out_size = l->size > r->size ? l->size : r->size;


  /* If we want inplace output, set the output pointer to one input. Note
     that the output type can be different from both inputs.  */
  if(flags & GAL_ARITHMETIC_INPLACE)
    {
      if     (l->type==otype && out_size==l->size)   o = l;
      else if(r->type==otype && out_size==r->size)   o = r;
    }


  /* If the output pointer was not set above for any of the possible
     reasons, allocate it. For 'mmapsize', note that since its 'size_t', it
     will always be positive. The '-1' that is recommended to give when you
     want the value in RAM is actually the largest possible memory
     location. So we just have to choose the smaller minmapsize of the two
     to decide if the output array should be in RAM or not. */
  if(o==NULL)
    o = gal_data_alloc(NULL, otype,
                       l->size>1 ? l->ndim  : r->ndim,
                       l->size>1 ? l->dsize : r->dsize,
                       l->size>1 ? l->wcs   : r->wcs,
                       0, minmapsize, quietmmap, NULL, NULL, NULL );


  /* Call the proper function for the operator. Since they heavily involve
     macros, their compilation can be very large if they are in a single
     function and file. So there is a separate C source and header file for
     each of these functions. */
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_PLUS:     arithmetic_plus(l, r, o);     break;
    case GAL_ARITHMETIC_OP_MINUS:    arithmetic_minus(l, r, o);    break;
    case GAL_ARITHMETIC_OP_MULTIPLY: arithmetic_multiply(l, r, o); break;
    case GAL_ARITHMETIC_OP_DIVIDE:   arithmetic_divide(l, r, o);   break;
    case GAL_ARITHMETIC_OP_LT:       arithmetic_lt(l, r, o);       break;
    case GAL_ARITHMETIC_OP_LE:       arithmetic_le(l, r, o);       break;
    case GAL_ARITHMETIC_OP_GT:       arithmetic_gt(l, r, o);       break;
    case GAL_ARITHMETIC_OP_GE:       arithmetic_ge(l, r, o);       break;
    case GAL_ARITHMETIC_OP_EQ:       arithmetic_eq(l, r, o);       break;
    case GAL_ARITHMETIC_OP_NE:       arithmetic_ne(l, r, o);       break;
    case GAL_ARITHMETIC_OP_AND:      arithmetic_and(l, r, o);      break;
    case GAL_ARITHMETIC_OP_OR:       arithmetic_or(l, r, o);       break;
    case GAL_ARITHMETIC_OP_BITAND:   arithmetic_bitand(l, r, o);   break;
    case GAL_ARITHMETIC_OP_BITOR:    arithmetic_bitor(l, r, o);    break;
    case GAL_ARITHMETIC_OP_BITXOR:   arithmetic_bitxor(l, r, o);   break;
    case GAL_ARITHMETIC_OP_BITLSH:   arithmetic_bitlsh(l, r, o);   break;
    case GAL_ARITHMETIC_OP_BITRSH:   arithmetic_bitrsh(l, r, o);   break;
    case GAL_ARITHMETIC_OP_MODULO:   arithmetic_modulo(l, r, o);   break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to address "
            "the problem. %d is not a valid operator code", __func__,
            PACKAGE_BUGREPORT, operator);
    }


  /* Clean up if necessary. Note that if the operation was requested to be
     in place, then the output might be one of the inputs. */
  if(flags & GAL_ARITHMETIC_FREE)
    {
      if     (o==l)       gal_data_free(r);
      else if(o==r)       gal_data_free(l);
      else              { gal_data_free(l); gal_data_free(r); }
    }


  /* Return */
  return o;
}





#define BINFUNC_RUN_FUNCTION(OT, RT, LT, OP){                           \
    LT *la=l->array;                                                    \
    RT *ra=r->array;                                                    \
    OT *oa=o->array, *of=oa + o->size;                                  \
    if(l->size==r->size) do *oa = OP(*la++, *ra++); while(++oa<of);     \
    else if(l->size==1)  do *oa = OP(*la,   *ra++); while(++oa<of);     \
    else                 do *oa = OP(*la++, *ra  ); while(++oa<of);     \
  }


#define BINFUNC_F_OPERATOR_LEFT_RIGHT_SET(RT, LT, OP)                   \
  switch(o->type)                                                       \
    {                                                                   \
    case GAL_TYPE_FLOAT32:                                              \
      BINFUNC_RUN_FUNCTION(float, RT, LT, OP);                          \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      BINFUNC_RUN_FUNCTION(double, RT, LT, OP);                         \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type %d not recognized for o->type ", \
            "BINFUNC_F_OPERATOR_LEFT_RIGHT_SET", o->type);              \
    }


#define BINFUNC_F_OPERATOR_LEFT_SET(LT, OP)                             \
  switch(r->type)                                                       \
    {                                                                   \
    case GAL_TYPE_FLOAT32:                                              \
      BINFUNC_F_OPERATOR_LEFT_RIGHT_SET(float, LT, OP);                 \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      BINFUNC_F_OPERATOR_LEFT_RIGHT_SET(double, LT, OP);                \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type %d not recognized for r->type",  \
            "BINFUNC_F_OPERATOR_LEFT_SET", r->type);                    \
    }


#define BINFUNC_F_OPERATOR_SET(OP)                                      \
  switch(l->type)                                                       \
    {                                                                   \
    case GAL_TYPE_FLOAT32:                                              \
      BINFUNC_F_OPERATOR_LEFT_SET(float, OP);                           \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      BINFUNC_F_OPERATOR_LEFT_SET(double, OP);                          \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type %d not recognized for l->type",  \
            "BINFUNC_F_OPERATOR_SET", l->type);                         \
    }


static gal_data_t *
arithmetic_binary_function_flt(int operator, int flags, gal_data_t *l,
                               gal_data_t *r)
{
  int final_otype;
  gal_data_t *o=NULL;
  size_t out_size, minmapsize;
  int quietmmap=l->quietmmap && r->quietmmap;


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_ARITHMETIC_NUMOK) && (l->size==1 || r->size==1))
      && gal_dimension_is_different(l, r) )
    error(EXIT_FAILURE, 0, "%s: the input datasets don't have the same "
          "dimension/size", __func__);

  /* Check for the types of the left and right operands. */
  arithmetic_check_float_input(l, operator, "first");
  arithmetic_check_float_input(r, operator, "second");

  /* Set the output type. */
  final_otype = gal_type_out(l->type, r->type);

  /* Set the output sizes. */
  minmapsize = ( l->minmapsize < r->minmapsize
                 ? l->minmapsize : r->minmapsize );
  out_size = l->size > r->size ? l->size : r->size;


  /* If we want inplace output, set the output pointer to one input. Note
     that the output type can be different from both inputs.  */
  if(flags & GAL_ARITHMETIC_INPLACE)
    {
      if     (l->type==final_otype && out_size==l->size)   o = l;
      else if(r->type==final_otype && out_size==r->size)   o = r;
    }


  /* If the output pointer was not set for any reason, allocate it. For
     'mmapsize', note that since its 'size_t', it will always be
     Positive. The '-1' that is recommended to give when you want the value
     in RAM is actually the largest possible memory location. So we just
     have to choose the smaller minmapsize of the two to decide if the
     output array should be in RAM or not. */
  if(o==NULL)
    o = gal_data_alloc(NULL, final_otype,
                       l->size>1 ? l->ndim  : r->ndim,
                       l->size>1 ? l->dsize : r->dsize,
                       l->size>1 ? l->wcs : r->wcs, 0, minmapsize,
                       quietmmap, NULL, NULL, NULL);


  /* Start setting the operator and operands. */
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_POW:  BINFUNC_F_OPERATOR_SET( pow  ); break;
    default:
      error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",
            __func__, operator);
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, 'BINFUNC_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, and we only have 'r' and 'l'
     allocated in any case. Alternatively, when the inputs shouldn't be
     freed, the only allocated spaces are the 'r' and 'l' arrays if their
     types weren't compiled for binary operations, we can tell this from
     the pointers: if they are different from the original pointers, they
     were allocated. */
  if(flags & GAL_ARITHMETIC_FREE)
    {
      if     (o==l)       gal_data_free(r);
      else if(o==r)       gal_data_free(l);
      else              { gal_data_free(l); gal_data_free(r); }
    }

  /* Return */
  return o;
}




















/**********************************************************************/
/****************         High-level functions        *****************/
/**********************************************************************/
/* Order is the same as in the manual. */
int
gal_arithmetic_set_operator(char *string, size_t *num_operands)
{
  int op;

  /* Simple arithmetic operators. */
  if      (!strcmp(string, "+" ))
    { op=GAL_ARITHMETIC_OP_PLUS;              *num_operands=2;  }
  else if (!strcmp(string, "-" ))
    { op=GAL_ARITHMETIC_OP_MINUS;             *num_operands=2;  }
  else if (!strcmp(string, "x" ))
    { op=GAL_ARITHMETIC_OP_MULTIPLY;          *num_operands=2;  }
  else if (!strcmp(string, "/" ))
    { op=GAL_ARITHMETIC_OP_DIVIDE;            *num_operands=2;  }
  else if (!strcmp(string, "%" ))
    { op=GAL_ARITHMETIC_OP_MODULO;            *num_operands=2;  }

  /* Mathematical Operators. */
  else if (!strcmp(string, "abs"))
    { op=GAL_ARITHMETIC_OP_ABS;               *num_operands=1;  }
  else if (!strcmp(string, "pow"))
    { op=GAL_ARITHMETIC_OP_POW;               *num_operands=2;  }
  else if (!strcmp(string, "sqrt"))
    { op=GAL_ARITHMETIC_OP_SQRT;              *num_operands=1;  }
  else if (!strcmp(string, "log"))
    { op=GAL_ARITHMETIC_OP_LOG;               *num_operands=1;  }
  else if (!strcmp(string, "log10"))
    { op=GAL_ARITHMETIC_OP_LOG10;             *num_operands=1;  }

  /* Units conversion functions */
  else if (!strcmp(string, "ra-to-degree"))
    { op=GAL_ARITHMETIC_OP_RA_TO_DEGREE;   *num_operands=1;  }
  else if (!strcmp(string, "dec-to-degree"))
    { op=GAL_ARITHMETIC_OP_DEC_TO_DEGREE;  *num_operands=1;  }
  else if (!strcmp(string, "degree-to-ra"))
    { op=GAL_ARITHMETIC_OP_DEGREE_TO_DEC;  *num_operands=1;  }
  else if (!strcmp(string, "degree-to-dec"))
    { op=GAL_ARITHMETIC_OP_DEGREE_TO_DEC;  *num_operands=1;  }

  /* Statistical/higher-level operators. */
  else if (!strcmp(string, "minvalue"))
    { op=GAL_ARITHMETIC_OP_MINVAL;            *num_operands=1;  }
  else if (!strcmp(string, "maxvalue"))
    { op=GAL_ARITHMETIC_OP_MAXVAL;            *num_operands=1;  }
  else if (!strcmp(string, "numbervalue"))
    { op=GAL_ARITHMETIC_OP_NUMBERVAL;         *num_operands=1;  }
  else if (!strcmp(string, "sumvalue"))
    { op=GAL_ARITHMETIC_OP_SUMVAL;            *num_operands=1;  }
  else if (!strcmp(string, "meanvalue"))
    { op=GAL_ARITHMETIC_OP_MEANVAL;           *num_operands=1;  }
  else if (!strcmp(string, "stdvalue"))
    { op=GAL_ARITHMETIC_OP_STDVAL;            *num_operands=1;  }
  else if (!strcmp(string, "medianvalue"))
    { op=GAL_ARITHMETIC_OP_MEDIANVAL;         *num_operands=1;  }
  else if (!strcmp(string, "min"))
    { op=GAL_ARITHMETIC_OP_MIN;               *num_operands=-1; }
  else if (!strcmp(string, "max"))
    { op=GAL_ARITHMETIC_OP_MAX;               *num_operands=-1; }
  else if (!strcmp(string, "number"))
    { op=GAL_ARITHMETIC_OP_NUMBER;            *num_operands=-1; }
  else if (!strcmp(string, "sum"))
    { op=GAL_ARITHMETIC_OP_SUM;               *num_operands=-1; }
  else if (!strcmp(string, "mean"))
    { op=GAL_ARITHMETIC_OP_MEAN;              *num_operands=-1; }
  else if (!strcmp(string, "std"))
    { op=GAL_ARITHMETIC_OP_STD;               *num_operands=-1; }
  else if (!strcmp(string, "median"))
    { op=GAL_ARITHMETIC_OP_MEDIAN;            *num_operands=-1; }
  else if (!strcmp(string, "quantile"))
    { op=GAL_ARITHMETIC_OP_QUANTILE;          *num_operands=-1; }
  else if (!strcmp(string, "sigclip-number"))
    { op=GAL_ARITHMETIC_OP_SIGCLIP_NUMBER;    *num_operands=-1; }
  else if (!strcmp(string, "sigclip-mean"))
    { op=GAL_ARITHMETIC_OP_SIGCLIP_MEAN;      *num_operands=-1; }
  else if (!strcmp(string, "sigclip-median"))
    { op=GAL_ARITHMETIC_OP_SIGCLIP_MEDIAN;    *num_operands=-1; }
  else if (!strcmp(string, "sigclip-std"))
    { op=GAL_ARITHMETIC_OP_SIGCLIP_STD;       *num_operands=-1; }

  /* The size operator */
  else if (!strcmp(string, "size"))
    { op=GAL_ARITHMETIC_OP_SIZE;              *num_operands=2;  }

  /* Conditional operators. */
  else if (!strcmp(string, "lt" ))
    { op=GAL_ARITHMETIC_OP_LT;                *num_operands=2;  }
  else if (!strcmp(string, "le"))
    { op=GAL_ARITHMETIC_OP_LE;                *num_operands=2;  }
  else if (!strcmp(string, "gt" ))
    { op=GAL_ARITHMETIC_OP_GT;                *num_operands=2;  }
  else if (!strcmp(string, "ge"))
    { op=GAL_ARITHMETIC_OP_GE;                *num_operands=2;  }
  else if (!strcmp(string, "eq"))
    { op=GAL_ARITHMETIC_OP_EQ;                *num_operands=2;  }
  else if (!strcmp(string, "ne"))
    { op=GAL_ARITHMETIC_OP_NE;                *num_operands=2;  }
  else if (!strcmp(string, "and"))
    { op=GAL_ARITHMETIC_OP_AND;               *num_operands=2;  }
  else if (!strcmp(string, "or"))
    { op=GAL_ARITHMETIC_OP_OR;                *num_operands=2;  }
  else if (!strcmp(string, "not"))
    { op=GAL_ARITHMETIC_OP_NOT;               *num_operands=1;  }
  else if (!strcmp(string, "isblank"))
    { op=GAL_ARITHMETIC_OP_ISBLANK;           *num_operands=1;  }
  else if (!strcmp(string, "where"))
    { op=GAL_ARITHMETIC_OP_WHERE;             *num_operands=3;  }

  /* Bitwise operators. */
  else if (!strcmp(string, "bitand"))
    { op=GAL_ARITHMETIC_OP_BITAND;            *num_operands=2;  }
  else if (!strcmp(string, "bitor"))
    { op=GAL_ARITHMETIC_OP_BITOR;             *num_operands=2;  }
  else if (!strcmp(string, "bitxor"))
    { op=GAL_ARITHMETIC_OP_BITXOR;            *num_operands=2;  }
  else if (!strcmp(string, "lshift"))
    { op=GAL_ARITHMETIC_OP_BITLSH;            *num_operands=2;  }
  else if (!strcmp(string, "rshift"))
    { op=GAL_ARITHMETIC_OP_BITRSH;            *num_operands=2;  }
  else if (!strcmp(string, "bitnot"))
    { op=GAL_ARITHMETIC_OP_BITNOT;            *num_operands=1;  }

  /* Type conversion. */
  else if (!strcmp(string, "uint8"))
    { op=GAL_ARITHMETIC_OP_TO_UINT8;          *num_operands=1;  }
  else if (!strcmp(string, "int8"))
    { op=GAL_ARITHMETIC_OP_TO_INT8;           *num_operands=1;  }
  else if (!strcmp(string, "uint16"))
    { op=GAL_ARITHMETIC_OP_TO_UINT16;         *num_operands=1;  }
  else if (!strcmp(string, "int16"))
    { op=GAL_ARITHMETIC_OP_TO_INT16;          *num_operands=1;  }
  else if (!strcmp(string, "uint32"))
    { op=GAL_ARITHMETIC_OP_TO_UINT32;         *num_operands=1;  }
  else if (!strcmp(string, "int32"))
    { op=GAL_ARITHMETIC_OP_TO_INT32;          *num_operands=1;  }
  else if (!strcmp(string, "uint64"))
    { op=GAL_ARITHMETIC_OP_TO_UINT64;         *num_operands=1;  }
  else if (!strcmp(string, "int64"))
    { op=GAL_ARITHMETIC_OP_TO_INT64;          *num_operands=1;  }
  else if (!strcmp(string, "float32"))
    { op=GAL_ARITHMETIC_OP_TO_FLOAT32;        *num_operands=1;  }
  else if (!strcmp(string, "float64"))
    { op=GAL_ARITHMETIC_OP_TO_FLOAT64;        *num_operands=1;  }

  /* Operator not defined. */
  else
    { op=GAL_ARITHMETIC_OP_INVALID; *num_operands=GAL_BLANK_INT; }
  return op;
}





char *
gal_arithmetic_operator_string(int operator)
{
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_PLUS:            return "+";
    case GAL_ARITHMETIC_OP_MINUS:           return "-";
    case GAL_ARITHMETIC_OP_MULTIPLY:        return "*";
    case GAL_ARITHMETIC_OP_DIVIDE:          return "/";
    case GAL_ARITHMETIC_OP_MODULO:          return "%";

    case GAL_ARITHMETIC_OP_LT:              return "<";
    case GAL_ARITHMETIC_OP_LE:              return "<=";
    case GAL_ARITHMETIC_OP_GT:              return ">";
    case GAL_ARITHMETIC_OP_GE:              return ">=";
    case GAL_ARITHMETIC_OP_EQ:              return "==";
    case GAL_ARITHMETIC_OP_NE:              return "!=";
    case GAL_ARITHMETIC_OP_AND:             return "and";
    case GAL_ARITHMETIC_OP_OR:              return "or";
    case GAL_ARITHMETIC_OP_NOT:             return "not";
    case GAL_ARITHMETIC_OP_ISBLANK:         return "isblank";
    case GAL_ARITHMETIC_OP_WHERE:           return "where";

    case GAL_ARITHMETIC_OP_BITAND:          return "bitand";
    case GAL_ARITHMETIC_OP_BITOR:           return "bitor";
    case GAL_ARITHMETIC_OP_BITXOR:          return "bitxor";
    case GAL_ARITHMETIC_OP_BITLSH:          return "lshift";
    case GAL_ARITHMETIC_OP_BITRSH:          return "rshift";
    case GAL_ARITHMETIC_OP_BITNOT:          return "bitnot";

    case GAL_ARITHMETIC_OP_ABS:             return "abs";
    case GAL_ARITHMETIC_OP_POW:             return "pow";
    case GAL_ARITHMETIC_OP_SQRT:            return "sqrt";
    case GAL_ARITHMETIC_OP_LOG:             return "log";
    case GAL_ARITHMETIC_OP_LOG10:           return "log10";

    case GAL_ARITHMETIC_OP_RA_TO_DEGREE:    return "ra-to-degree";
    case GAL_ARITHMETIC_OP_DEC_TO_DEGREE:   return "dec-to-degree";
    case GAL_ARITHMETIC_OP_DEGREE_TO_RA:    return "degree-to-ra";
    case GAL_ARITHMETIC_OP_DEGREE_TO_DEC:   return "degree-to-dec";

    case GAL_ARITHMETIC_OP_MINVAL:          return "minvalue";
    case GAL_ARITHMETIC_OP_MAXVAL:          return "maxvalue";
    case GAL_ARITHMETIC_OP_NUMBERVAL:       return "numbervalue";
    case GAL_ARITHMETIC_OP_SUMVAL:          return "sumvalue";
    case GAL_ARITHMETIC_OP_MEANVAL:         return "meanvalue";
    case GAL_ARITHMETIC_OP_STDVAL:          return "stdvalue";
    case GAL_ARITHMETIC_OP_MEDIANVAL:       return "medianvalue";

    case GAL_ARITHMETIC_OP_MIN:             return "min";
    case GAL_ARITHMETIC_OP_MAX:             return "max";
    case GAL_ARITHMETIC_OP_NUMBER:          return "number";
    case GAL_ARITHMETIC_OP_SUM:             return "sum";
    case GAL_ARITHMETIC_OP_MEAN:            return "mean";
    case GAL_ARITHMETIC_OP_STD:             return "std";
    case GAL_ARITHMETIC_OP_MEDIAN:          return "median";
    case GAL_ARITHMETIC_OP_QUANTILE:        return "quantile";
    case GAL_ARITHMETIC_OP_SIGCLIP_NUMBER:  return "sigclip-number";
    case GAL_ARITHMETIC_OP_SIGCLIP_MEDIAN:  return "sigclip-median";
    case GAL_ARITHMETIC_OP_SIGCLIP_MEAN:    return "sigclip-mean";
    case GAL_ARITHMETIC_OP_SIGCLIP_STD:     return "sigclip-number";

    case GAL_ARITHMETIC_OP_SIZE:            return "size";

    case GAL_ARITHMETIC_OP_TO_UINT8:        return "uchar";
    case GAL_ARITHMETIC_OP_TO_INT8:         return "char";
    case GAL_ARITHMETIC_OP_TO_UINT16:       return "ushort";
    case GAL_ARITHMETIC_OP_TO_INT16:        return "short";
    case GAL_ARITHMETIC_OP_TO_UINT32:       return "uint";
    case GAL_ARITHMETIC_OP_TO_INT32:        return "int";
    case GAL_ARITHMETIC_OP_TO_UINT64:       return "ulong";
    case GAL_ARITHMETIC_OP_TO_INT64:        return "long";
    case GAL_ARITHMETIC_OP_TO_FLOAT32:      return "float32";
    case GAL_ARITHMETIC_OP_TO_FLOAT64:      return "float64";

    default:                                return NULL;
    }
  return NULL;
}





gal_data_t *
gal_arithmetic(int operator, size_t numthreads, int flags, ...)
{
  va_list va;
  gal_data_t *d1, *d2, *d3, *out=NULL;

  /* Prepare the variable arguments (starting after the flags argument). */
  va_start(va, flags);

  /* Depending on the operator, do the job: */
  switch(operator)
    {

    /* Binary operators with any data type. */
    case GAL_ARITHMETIC_OP_PLUS:
    case GAL_ARITHMETIC_OP_MINUS:
    case GAL_ARITHMETIC_OP_MULTIPLY:
    case GAL_ARITHMETIC_OP_DIVIDE:
    case GAL_ARITHMETIC_OP_LT:
    case GAL_ARITHMETIC_OP_LE:
    case GAL_ARITHMETIC_OP_GT:
    case GAL_ARITHMETIC_OP_GE:
    case GAL_ARITHMETIC_OP_EQ:
    case GAL_ARITHMETIC_OP_NE:
    case GAL_ARITHMETIC_OP_AND:
    case GAL_ARITHMETIC_OP_OR:
      d1 = va_arg(va, gal_data_t *);
      d2 = va_arg(va, gal_data_t *);
      out=arithmetic_binary(operator, flags, d1, d2);
      break;

    case GAL_ARITHMETIC_OP_NOT:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_not(d1, flags);
      break;

    case GAL_ARITHMETIC_OP_ISBLANK:
      d1 = va_arg(va, gal_data_t *);
      out = gal_blank_flag(d1);
      if(flags & GAL_ARITHMETIC_FREE) gal_data_free(d1);
      break;

    case GAL_ARITHMETIC_OP_WHERE:
      d1 = va_arg(va, gal_data_t *);    /* To modify value/array.     */
      d2 = va_arg(va, gal_data_t *);    /* Condition (unsigned char). */
      d3 = va_arg(va, gal_data_t *);    /* If true value/array.       */
      arithmetic_where(flags, d1, d2, d3);
      out=d1;
      break;


    /* Unary function operators. */
    case GAL_ARITHMETIC_OP_SQRT:
    case GAL_ARITHMETIC_OP_LOG:
    case GAL_ARITHMETIC_OP_LOG10:
    case GAL_ARITHMETIC_OP_RA_TO_DEGREE:
    case GAL_ARITHMETIC_OP_DEC_TO_DEGREE:
    case GAL_ARITHMETIC_OP_DEGREE_TO_RA:
    case GAL_ARITHMETIC_OP_DEGREE_TO_DEC:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_unary_function(operator, flags, d1);
      break;

    /* Statistical operators that return one value. */
    case GAL_ARITHMETIC_OP_MINVAL:
    case GAL_ARITHMETIC_OP_MAXVAL:
    case GAL_ARITHMETIC_OP_NUMBERVAL:
    case GAL_ARITHMETIC_OP_SUMVAL:
    case GAL_ARITHMETIC_OP_MEANVAL:
    case GAL_ARITHMETIC_OP_STDVAL:
    case GAL_ARITHMETIC_OP_MEDIANVAL:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_from_statistics(operator, flags, d1);
      break;

    /* Absolute operator. */
    case GAL_ARITHMETIC_OP_ABS:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_abs(flags, d1);
      break;

    /* Multi-operand operators */
    case GAL_ARITHMETIC_OP_MIN:
    case GAL_ARITHMETIC_OP_MAX:
    case GAL_ARITHMETIC_OP_NUMBER:
    case GAL_ARITHMETIC_OP_SUM:
    case GAL_ARITHMETIC_OP_MEAN:
    case GAL_ARITHMETIC_OP_STD:
    case GAL_ARITHMETIC_OP_MEDIAN:
    case GAL_ARITHMETIC_OP_QUANTILE:
    case GAL_ARITHMETIC_OP_SIGCLIP_STD:
    case GAL_ARITHMETIC_OP_SIGCLIP_MEAN:
    case GAL_ARITHMETIC_OP_SIGCLIP_MEDIAN:
    case GAL_ARITHMETIC_OP_SIGCLIP_NUMBER:
      d1 = va_arg(va, gal_data_t *);
      d2 = va_arg(va, gal_data_t *);
      out=arithmetic_multioperand(operator, flags, d1, d2, numthreads);
      break;


    /* Binary function operators. */
    case GAL_ARITHMETIC_OP_POW:
      d1 = va_arg(va, gal_data_t *);
      d2 = va_arg(va, gal_data_t *);
      out=arithmetic_binary_function_flt(operator, flags, d1, d2);
      break;


    /* Binary operators that only work on integer types. */
    case GAL_ARITHMETIC_OP_BITAND:
    case GAL_ARITHMETIC_OP_BITOR:
    case GAL_ARITHMETIC_OP_BITXOR:
    case GAL_ARITHMETIC_OP_BITLSH:
    case GAL_ARITHMETIC_OP_BITRSH:
    case GAL_ARITHMETIC_OP_MODULO:
      d1 = va_arg(va, gal_data_t *);
      d2 = va_arg(va, gal_data_t *);
      out=arithmetic_binary(operator, flags, d1, d2);
      break;

    case GAL_ARITHMETIC_OP_BITNOT:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_bitwise_not(flags, d1);
      break;


    /* Conversion operators. */
    case GAL_ARITHMETIC_OP_TO_UINT8:
    case GAL_ARITHMETIC_OP_TO_INT8:
    case GAL_ARITHMETIC_OP_TO_UINT16:
    case GAL_ARITHMETIC_OP_TO_INT16:
    case GAL_ARITHMETIC_OP_TO_UINT32:
    case GAL_ARITHMETIC_OP_TO_INT32:
    case GAL_ARITHMETIC_OP_TO_UINT64:
    case GAL_ARITHMETIC_OP_TO_INT64:
    case GAL_ARITHMETIC_OP_TO_FLOAT32:
    case GAL_ARITHMETIC_OP_TO_FLOAT64:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_change_type(d1, operator, flags);
      break;

    /* Size operator */
    case GAL_ARITHMETIC_OP_SIZE:
      d1 = va_arg(va, gal_data_t *);
      d2 = va_arg(va, gal_data_t *);
      out=arithmetic_size(operator, flags, d1, d2);
      break;

    /* When the operator is not recognized. */
    default:
      error(EXIT_FAILURE, 0, "%s: the argument '%d' could not be "
            "interpretted as an operator", __func__, operator);
    }

  /* End the variable argument structure and return. */
  va_end(va);
  return out;
}
