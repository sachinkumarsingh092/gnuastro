/*********************************************************************
Arithmetic operations on data structures
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2018, Free Software Foundation, Inc.

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
#include <stdlib.h>
#include <stdarg.h>

#include <gnuastro/blank.h>
#include <gnuastro/qsort.h>
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
            "%s operand has type %s. You can use the `float' or `double' "
            "operators before this operator to explicitly convert to the "
            "desired precision floating point type. If the operand was "
            "originally a typed number (string of characters), add an `f' "
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
                     data->wcs, 0, data->minmapsize, data->name, data->unit,
                     data->comment);
  o=out->array;


  /* Go over the pixels and set the output values. */
  switch(data->type)
    {
    case GAL_TYPE_UINT8:   TYPE_CASE_FOR_NOT(uint8_t);   break;
    case GAL_TYPE_INT8:    TYPE_CASE_FOR_NOT(int8_t);    break;
    case GAL_TYPE_UINT16:  TYPE_CASE_FOR_NOT(uint16_t);  break;
    case GAL_TYPE_INT16:   TYPE_CASE_FOR_NOT(int16_t);   break;
    case GAL_TYPE_UINT32:  TYPE_CASE_FOR_NOT(uint32_t);  break;
    case GAL_TYPE_INT32:   TYPE_CASE_FOR_NOT(int32_t);   break;
    case GAL_TYPE_UINT64:  TYPE_CASE_FOR_NOT(uint64_t);  break;
    case GAL_TYPE_INT64:   TYPE_CASE_FOR_NOT(int64_t);   break;
    case GAL_TYPE_FLOAT32: TYPE_CASE_FOR_NOT(float);     break;
    case GAL_TYPE_FLOAT64: TYPE_CASE_FOR_NOT(double);    break;

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
                       0, in->minmapsize, NULL, NULL, NULL);

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
                         in->wcs, 0, in->minmapsize, in->name, in->unit,
                         in->comment);

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





#define UNIFUNC_RUN_FUNCTION_ON_ELEMENT(IT, OP){                        \
    IT *ia=in->array, *oa=o->array, *iaf=ia + in->size;                 \
    do *oa++ = OP(*ia++); while(ia<iaf);                                \
  }

#define UNIARY_FUNCTION_ON_ELEMENT(OP)                                  \
  switch(in->type)                                                      \
    {                                                                   \
    case GAL_TYPE_UINT8:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint8_t, OP)                      \
        break;                                                          \
    case GAL_TYPE_INT8:                                                 \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int8_t, OP)                       \
        break;                                                          \
    case GAL_TYPE_UINT16:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint16_t, OP)                     \
        break;                                                          \
    case GAL_TYPE_INT16:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int16_t, OP)                      \
        break;                                                          \
    case GAL_TYPE_UINT32:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint32_t, OP)                     \
        break;                                                          \
    case GAL_TYPE_INT32:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int32_t, OP)                      \
        break;                                                          \
    case GAL_TYPE_UINT64:                                               \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(uint64_t, OP)                     \
        break;                                                          \
    case GAL_TYPE_INT64:                                                \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int64_t, OP)                      \
        break;                                                          \
    case GAL_TYPE_FLOAT32:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(float, OP)                        \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(double, OP)                       \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",         \
            "UNIARY_FUNCTION_ON_ELEMENT", in->type);                    \
    }

static gal_data_t *
arithmetic_unary_function(int operator, int flags, gal_data_t *in)
{
  gal_data_t *o;

  /* If we want inplace output, set the output pointer to the input
     pointer, for every pixel, the operation will be independent. */
  if(flags & GAL_ARITHMETIC_INPLACE)
    o = in;
  else
    o = gal_data_alloc(NULL, in->type, in->ndim, in->dsize, in->wcs,
                       0, in->minmapsize, NULL, NULL, NULL);

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

    default:
      error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",
            __func__, operator);
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, `UNIFUNC_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, and we only have `r' and `l'
     allocated in any case. Alternatively, when the inputs shouldn't be
     freed, the only allocated spaces are the `r' and `l' arrays if their
     types weren't compiled for binary operations, we can tell this from
     the pointers: if they are different from the original pointers, they
     were allocated. */
  if( (flags & GAL_ARITHMETIC_FREE) && o!=in)
    gal_data_free(in);

  /* Return */
  return o;
}





/* Call functions in the `gnuastro/statistics' library. */
static gal_data_t *
arithmetic_from_statistics(int operator, int flags, gal_data_t *input)
{
  gal_data_t *out=NULL;
  int ip=(flags & GAL_ARITHMETIC_INPLACE) || (flags & GAL_ARITHMETIC_FREE);

  switch(operator)
    {
    case GAL_ARITHMETIC_OP_MINVAL:  out=gal_statistics_minimum(input); break;
    case GAL_ARITHMETIC_OP_MAXVAL:  out=gal_statistics_maximum(input); break;
    case GAL_ARITHMETIC_OP_NUMVAL:  out=gal_statistics_number(input);  break;
    case GAL_ARITHMETIC_OP_SUMVAL:  out=gal_statistics_sum(input);     break;
    case GAL_ARITHMETIC_OP_MEANVAL: out=gal_statistics_mean(input);    break;
    case GAL_ARITHMETIC_OP_STDVAL:  out=gal_statistics_std(input);     break;
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
/***************                  Where                   **************/
/***********************************************************************/

/* When the `iftrue' dataset only has one element and the element is blank,
   then it will be replaced with the blank value of the type of the output
   data. */
#define DO_WHERE_OPERATION(ITT, OT) {                                \
    ITT *it=iftrue->array;                                           \
    OT b, *o=out->array, *of=o+out->size;                            \
    if(iftrue->size==1)                                              \
      {                                                              \
        if( gal_blank_present(iftrue, 0) )                           \
          {                                                          \
            gal_blank_write(&b, out->type);                          \
            do { *o = *c++ ? b : *o;        } while(++o<of);         \
          }                                                          \
        else                                                         \
          do   { *o = *c++ ? *it : *o;       } while(++o<of);        \
      }                                                              \
    else                                                             \
      do       { *o = *c++ ? *it : *o; ++it; } while(++o<of);        \
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
            "`iftrue' dataset", "WHERE_OUT_SET", iftrue->type);         \
    }





static void
arithmetic_where(int flags, gal_data_t *out, gal_data_t *cond,
                 gal_data_t *iftrue)
{
  unsigned char *c=cond->array;

  /* The condition operator has to be unsigned char. */
  if(cond->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "%s: the condition operand must be an "
          "`uint8' type, but the given condition operand has a "
          "`%s' type", __func__, gal_type_name(cond->type, 1));

  /* The dimension and sizes of the out and condition data sets must be the
     same. */
  if(gal_data_dsize_is_different(out, cond))
    error(EXIT_FAILURE, 0, "%s: the output and condition data sets of the "
          "must be the same size", __func__);

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
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized for the `out'",
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
#define MULTIOPERAND_MIN(TYPE) {                                        \
    TYPE p, max;                                                        \
    size_t n, j=0;                                                      \
    gal_type_max(list->type, &max);                                     \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        p=max;                                                          \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {   /* Only for integer types, b==b. */                       \
            if( hasblank[i] && b==b)                                    \
              {                                                         \
                if( a[i][j] != b )                                      \
                  { p = a[i][j] < p ? a[i][j] : p; ++n; }               \
              }                                                         \
            else { p = a[i][j] < p ? a[i][j] : p; ++n; }                \
          }                                                             \
        *o++ = n ? p : b;  /* No usable elements: set to blank. */      \
        ++j;                                                            \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_MAX(TYPE) {                                        \
    TYPE p, min;                                                        \
    size_t n, j=0;                                                      \
    gal_type_min(list->type, &min);                                     \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        p=min;                                                          \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {   /* Only for integer types, b==b. */                       \
            if( hasblank[i] && b==b)                                    \
              {                                                         \
                if( a[i][j] != b )                                      \
                  { p = a[i][j] > p ? a[i][j] : p; ++n; }               \
              }                                                         \
            else { p = a[i][j] > p ? a[i][j] : p; ++n; }                \
          }                                                             \
        *o++ = n ? p : b;  /* No usable elements: set to blank. */      \
        ++j;                                                            \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_NUM {                                              \
    int use;                                                            \
    size_t n, j=0;                                                      \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b       ? 1 : 0 )      /* Integer */ \
                      : ( a[i][j]==a[i][j] ? 1 : 0 ) );   /* Float   */ \
            else use=1;                                                 \
                                                                        \
            /* Increment counter if necessary. */                       \
            if(use) ++n;                                                \
          }                                                             \
        *o++ = n;                                                       \
        ++j;                                                            \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_SUM {                                              \
    int use;                                                            \
    double sum;                                                         \
    size_t n, j=0;                                                      \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        sum=0.0f;                                                       \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b     ? 1 : 0 )       /* Integer */  \
                      : ( a[i][j]==*a[i] ? 1 : 0 ) );    /* Float   */  \
            else use=1;                                                 \
                                                                        \
            /* Use in sum if necessary. */                              \
            if(use) { sum += a[i][j]; ++n; }                            \
          }                                                             \
        *o++ = n ? sum : b;                                             \
        ++j;                                                            \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_MEAN {                                             \
    int use;                                                            \
    double sum;                                                         \
    size_t n, j=0;                                                      \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        sum=0.0f;                                                       \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
              use = ( b==b                                              \
                      ? ( a[i][j]!=b       ? 1 : 0 )     /* Integer */  \
                      : ( a[i][j]==a[i][j] ? 1 : 0 ) );  /* Float   */  \
            else use=1;                                                 \
                                                                        \
            /* Calculate the mean if necessary. */                      \
            if(use) { sum += a[i][j]; ++n; }                            \
          }                                                             \
        *o++ = n ? sum/n : b;                                           \
        ++j;                                                            \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_STD {                                              \
    int use;                                                            \
    size_t n, j=0;                                                      \
    double sum, sum2;                                                   \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        sum=sum2=0.0f;                                                  \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
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
        *o++ = n ? sqrt( (sum2-sum*sum/n)/n ) : b;                      \
        ++j;                                                            \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_MEDIAN(TYPE, QSORT_F) {                            \
    int use;                                                            \
    size_t n, j=0;                                                      \
    TYPE *pixs=gal_data_malloc_array(list->type, dnum, __func__, "pixs"); \
                                                                        \
    /* Loop over each pixel */                                          \
    do                                                                  \
      {                                                                 \
        /* Initialize. */                                               \
        n=0;                                                            \
                                                                        \
        /* Loop over each array. */                                     \
        for(i=0;i<dnum;++i)                                             \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
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
            *o++ = n%2 ? pixs[n/2] : (pixs[n/2] + pixs[n/2-1])/2 ;      \
          }                                                             \
        else                                                            \
          *o++=b;                                                       \
        ++j;                                                            \
      }                                                                 \
    while(o<of);                                                        \
                                                                        \
    /* Clean up. */                                                     \
    free(pixs);                                                         \
  }





#define MULTIOPERAND_TYPE_SET(TYPE, QSORT_F) {                          \
    TYPE b, **a, *o=out->array, *of=o+out->size;                        \
    size_t i=0;  /* Different from the `i' in the main function. */     \
                                                                        \
    /* Allocate space to keep the pointers to the arrays of each. */    \
    /* Input data structure. The operators will increment these */      \
    /* pointers while parsing them. */                                  \
    errno=0;                                                            \
    a=malloc(dnum*sizeof *a);                                           \
    if(a==NULL)                                                         \
      error(EXIT_FAILURE, 0, "%s: %zu bytes for `a'",                   \
            "MULTIOPERAND_TYPE_SET", dnum*sizeof *a);                   \
                                                                        \
    /* Fill in the array pointers and the blank value for this type. */ \
    gal_blank_write(&b, list->type);                                    \
    for(tmp=list;tmp!=NULL;tmp=tmp->next)                               \
      a[i++]=tmp->array;                                                \
                                                                        \
    /* Do the operation. */                                             \
    switch(operator)                                                    \
      {                                                                 \
      case GAL_ARITHMETIC_OP_MIN:                                       \
        MULTIOPERAND_MIN(TYPE);                                         \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_MAX:                                       \
        MULTIOPERAND_MAX(TYPE);                                         \
        break;                                                          \
                                                                        \
      case GAL_ARITHMETIC_OP_NUM:                                       \
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
      default:                                                          \
        error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",   \
              "MULTIOPERAND_TYPE_SET", operator);                       \
      }                                                                 \
                                                                        \
    /* Clean up. */                                                     \
    free(a);                                                            \
  }





/* The single operator in this function is assumed to be a linked list. The
   number of operators is determined from the fact that the last node in
   the linked list must have a NULL pointer as its `next' element.*/
static gal_data_t *
arithmetic_multioperand(int operator, int flags, gal_data_t *list)
{
  uint8_t *hasblank;
  size_t i=0, dnum=1;
  gal_data_t *out, *tmp, *ttmp;


  /* For generality, `list' can be a NULL pointer, in that case, this
     function will return a NULL pointer and avoid further processing. */
  if(list==NULL) return NULL;


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
      if( gal_data_dsize_is_different(list, tmp) )
        error(EXIT_FAILURE, 0, "%s: the sizes of all operands to the %s "
              "operator must be same", __func__,
              gal_arithmetic_operator_string(operator));
    }


  /* Set the output data structure. */
  if(flags & GAL_ARITHMETIC_INPLACE)
    out = list;                 /* The top element in the list. */
  else
    out = gal_data_alloc(NULL, list->type, list->ndim, list->dsize,
                         list->wcs, 0, list->minmapsize, NULL, NULL, NULL);


  /* hasblank is used to see if a blank value should be checked for each
     list element or not. */
  hasblank=gal_data_malloc_array(GAL_TYPE_UINT8, dnum, __func__, "hasblank");
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    hasblank[i++]=gal_blank_present(tmp, 0);


  /* Start the operation. */
  switch(list->type)
    {
    case GAL_TYPE_UINT8:
      MULTIOPERAND_TYPE_SET(uint8_t,   gal_qsort_uint8_increasing);
      break;
    case GAL_TYPE_INT8:
      MULTIOPERAND_TYPE_SET(int8_t,    gal_qsort_int8_increasing);
      break;
    case GAL_TYPE_UINT16:
      MULTIOPERAND_TYPE_SET(uint16_t,  gal_qsort_uint16_increasing);
      break;
    case GAL_TYPE_INT16:
      MULTIOPERAND_TYPE_SET(int16_t,   gal_qsort_int16_increasing);
      break;
    case GAL_TYPE_UINT32:
      MULTIOPERAND_TYPE_SET(uint32_t,  gal_qsort_uint32_increasing);
      break;
    case GAL_TYPE_INT32:
      MULTIOPERAND_TYPE_SET(int32_t,   gal_qsort_int32_increasing);
      break;
    case GAL_TYPE_UINT64:
      MULTIOPERAND_TYPE_SET(uint64_t,  gal_qsort_uint64_increasing);
      break;
    case GAL_TYPE_INT64:
      MULTIOPERAND_TYPE_SET(int64_t,   gal_qsort_int64_increasing);
      break;
    case GAL_TYPE_FLOAT32:
      MULTIOPERAND_TYPE_SET(float,     gal_qsort_float32_increasing);
      break;
    case GAL_TYPE_FLOAT64:
      MULTIOPERAND_TYPE_SET(double,    gal_qsort_float64_increasing);
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, list->type);
    }


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
   (`BINARY_OP_INCR_OT_RT_LT_SET') which doesn' use `checkblanks'.*/
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
  /* Read the variable arguments. `lo' and `ro' keep the original data, in
     case their type isn't built (based on configure options are configure
     time). */
  int32_t otype;
  gal_data_t *o=NULL;
  size_t out_size, minmapsize;


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_ARITHMETIC_NUMOK) && (l->size==1 || r->size==1))
      && gal_data_dsize_is_different(l, r) )
    error(EXIT_FAILURE, 0, "%s: the non-number inputs to %s don't have the "
          "same dimension/size", __func__,
          gal_arithmetic_operator_string(operator));


  /* Set the output type. For the comparison operators, the output type is
     either 0 or 1, so we will set the output type to `unsigned char' for
     efficient memory and CPU usage. Since the number of operators without
     a fixed output type (like the conditionals) is less, by `default' we
     will set the output type to `unsigned char', and if any of the other
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
     reasons, allocate it. For `mmapsize', note that since its `size_t', it
     will always be positive. The `-1' that is recommended to give when you
     want the value in RAM is actually the largest possible memory
     location. So we just have to choose the smaller minmapsize of the two
     to decide if the output array should be in RAM or not. */
  if(o==NULL)
    o = gal_data_alloc(NULL, otype,
                       l->size>1 ? l->ndim  : r->ndim,
                       l->size>1 ? l->dsize : r->dsize,
                       l->size>1 ? l->wcs   : r->wcs,
                       0, minmapsize, NULL, NULL, NULL );


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


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_ARITHMETIC_NUMOK) && (l->size==1 || r->size==1))
      && gal_data_dsize_is_different(l, r) )
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
     `mmapsize', note that since its `size_t', it will always be
     Positive. The `-1' that is recommended to give when you want the value
     in RAM is actually the largest possible memory location. So we just
     have to choose the smaller minmapsize of the two to decide if the
     output array should be in RAM or not. */
  if(o==NULL)
    o = gal_data_alloc(NULL, final_otype,
                       l->size>1 ? l->ndim  : r->ndim,
                       l->size>1 ? l->dsize : r->dsize,
                       l->size>1 ? l->wcs : r->wcs, 0, minmapsize,
                       NULL, NULL, NULL);


  /* Start setting the operator and operands. */
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_POW:  BINFUNC_F_OPERATOR_SET( pow  ); break;
    default:
      error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",
            __func__, operator);
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, `BINFUNC_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, and we only have `r' and `l'
     allocated in any case. Alternatively, when the inputs shouldn't be
     freed, the only allocated spaces are the `r' and `l' arrays if their
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
char *
gal_arithmetic_operator_string(int operator)
{
  switch(operator)
    {
    case GAL_ARITHMETIC_OP_PLUS:         return "+";
    case GAL_ARITHMETIC_OP_MINUS:        return "-";
    case GAL_ARITHMETIC_OP_MULTIPLY:     return "*";
    case GAL_ARITHMETIC_OP_DIVIDE:       return "/";
    case GAL_ARITHMETIC_OP_MODULO:       return "%";

    case GAL_ARITHMETIC_OP_LT:           return "<";
    case GAL_ARITHMETIC_OP_LE:           return "<=";
    case GAL_ARITHMETIC_OP_GT:           return ">";
    case GAL_ARITHMETIC_OP_GE:           return ">=";
    case GAL_ARITHMETIC_OP_EQ:           return "==";
    case GAL_ARITHMETIC_OP_NE:           return "!=";
    case GAL_ARITHMETIC_OP_AND:          return "and";
    case GAL_ARITHMETIC_OP_OR:           return "or";
    case GAL_ARITHMETIC_OP_NOT:          return "not";
    case GAL_ARITHMETIC_OP_ISBLANK:      return "isblank";
    case GAL_ARITHMETIC_OP_WHERE:        return "where";

    case GAL_ARITHMETIC_OP_BITAND:       return "bitand";
    case GAL_ARITHMETIC_OP_BITOR:        return "bitor";
    case GAL_ARITHMETIC_OP_BITXOR:       return "bitxor";
    case GAL_ARITHMETIC_OP_BITLSH:       return "lshift";
    case GAL_ARITHMETIC_OP_BITRSH:       return "rshift";
    case GAL_ARITHMETIC_OP_BITNOT:       return "bitnot";

    case GAL_ARITHMETIC_OP_ABS:          return "abs";
    case GAL_ARITHMETIC_OP_POW:          return "pow";
    case GAL_ARITHMETIC_OP_SQRT:         return "sqrt";
    case GAL_ARITHMETIC_OP_LOG:          return "log";
    case GAL_ARITHMETIC_OP_LOG10:        return "log10";

    case GAL_ARITHMETIC_OP_MINVAL:       return "minvalue";
    case GAL_ARITHMETIC_OP_MAXVAL:       return "maxvalue";
    case GAL_ARITHMETIC_OP_NUMVAL:       return "numvalue";
    case GAL_ARITHMETIC_OP_SUMVAL:       return "sumvalue";
    case GAL_ARITHMETIC_OP_MEANVAL:      return "meanvalue";
    case GAL_ARITHMETIC_OP_STDVAL:       return "stdvalue";
    case GAL_ARITHMETIC_OP_MEDIANVAL:    return "medianvalue";

    case GAL_ARITHMETIC_OP_MIN:          return "min";
    case GAL_ARITHMETIC_OP_MAX:          return "max";
    case GAL_ARITHMETIC_OP_NUM:          return "num";
    case GAL_ARITHMETIC_OP_SUM:          return "sum";
    case GAL_ARITHMETIC_OP_MEAN:         return "mean";
    case GAL_ARITHMETIC_OP_STD:          return "std";
    case GAL_ARITHMETIC_OP_MEDIAN:       return "median";

    case GAL_ARITHMETIC_OP_TO_UINT8:     return "uchar";
    case GAL_ARITHMETIC_OP_TO_INT8:      return "char";
    case GAL_ARITHMETIC_OP_TO_UINT16:    return "ushort";
    case GAL_ARITHMETIC_OP_TO_INT16:     return "short";
    case GAL_ARITHMETIC_OP_TO_UINT32:    return "uint";
    case GAL_ARITHMETIC_OP_TO_INT32:     return "int";
    case GAL_ARITHMETIC_OP_TO_UINT64:    return "ulong";
    case GAL_ARITHMETIC_OP_TO_INT64:     return "long";
    case GAL_ARITHMETIC_OP_TO_FLOAT32:   return "float32";
    case GAL_ARITHMETIC_OP_TO_FLOAT64:   return "float64";

    default:
      error(EXIT_FAILURE, 0, "%s: operator code %d not recognized",
            __func__, operator);
    }

  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us to fix the problem. "
        "Control has reached the end of this function. This should not have "
        "happened", __func__);
  return NULL;
}





gal_data_t *
gal_arithmetic(int operator, int flags, ...)
{
  va_list va;
  gal_data_t *d1, *d2, *d3, *out=NULL;

  /* Prepare the variable arguments (starting after the flags argument). */
  va_start(va, flags);

  /* Depending on the operator do the job: */
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
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_unary_function(operator, flags, d1);
      break;

    /* Statistical operators that return one value. */
    case GAL_ARITHMETIC_OP_MINVAL:
    case GAL_ARITHMETIC_OP_MAXVAL:
    case GAL_ARITHMETIC_OP_NUMVAL:
    case GAL_ARITHMETIC_OP_SUMVAL:
    case GAL_ARITHMETIC_OP_MEANVAL:
    case GAL_ARITHMETIC_OP_STDVAL:
    case GAL_ARITHMETIC_OP_MEDIANVAL:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_from_statistics(operator, flags, d1);
      break;

    /* Absolute operator */
    case GAL_ARITHMETIC_OP_ABS:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_abs(flags, d1);
      break;


    /* Multi-operand operators */
    case GAL_ARITHMETIC_OP_MIN:
    case GAL_ARITHMETIC_OP_MAX:
    case GAL_ARITHMETIC_OP_NUM:
    case GAL_ARITHMETIC_OP_SUM:
    case GAL_ARITHMETIC_OP_MEAN:
    case GAL_ARITHMETIC_OP_STD:
    case GAL_ARITHMETIC_OP_MEDIAN:
      d1 = va_arg(va, gal_data_t *);
      out=arithmetic_multioperand(operator, flags, d1);
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


    /* When the operator is not recognized. */
    default:
      error(EXIT_FAILURE, 0, "%s: the argument \"%d\" could not be "
            "interpretted as an operator", __func__, operator);
    }

  /* End the variable argument structure and return. */
  va_end(va);
  return out;
}
