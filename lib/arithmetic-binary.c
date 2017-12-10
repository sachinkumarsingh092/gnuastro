/*********************************************************************
Arithmetic operations on data structures.
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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/blank.h>
#include <gnuastro/arithmetic.h>

#include <gnuastro-internal/arithmetic-internal.h>





/************************************************************************/
/*************            Native type macros            *****************/
/************************************************************************/
#if GAL_CONFIG_BIN_OP_UINT8 == 1
#define BINARY_LT_IS_UINT8                                         \
  case GAL_TYPE_UINT8: BINARY_LT_SET(uint8_t);          break;
#define BINARY_LT_SET_RT_IS_UINT8(LT)                              \
  case GAL_TYPE_UINT8: BINARY_RT_LT_SET(uint8_t, LT);   break;
#else
#define BINARY_LT_IS_UINT8
#define BINARY_LT_SET_RT_IS_UINT8(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT8 == 1
#define BINARY_LT_IS_INT8                                          \
  case GAL_TYPE_INT8: BINARY_LT_SET(int8_t);            break;
#define BINARY_LT_SET_RT_IS_INT8(LT)                               \
  case GAL_TYPE_INT8: BINARY_RT_LT_SET(int8_t, LT);     break;
#else
#define BINARY_LT_IS_INT8
#define BINARY_LT_SET_RT_IS_INT8(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT16 == 1
#define BINARY_LT_IS_UINT16                                        \
  case GAL_TYPE_UINT16: BINARY_LT_SET(uint16_t);        break;
#define BINARY_LT_SET_RT_IS_UINT16(LT)                             \
  case GAL_TYPE_UINT16: BINARY_RT_LT_SET(uint16_t, LT); break;
#else
#define BINARY_LT_IS_UINT16
#define BINARY_LT_SET_RT_IS_UINT16(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT16 == 1
#define BINARY_LT_IS_INT16                                         \
  case GAL_TYPE_INT16: BINARY_LT_SET(int16_t);          break;
#define BINARY_LT_SET_RT_IS_INT16(LT)                              \
  case GAL_TYPE_INT16: BINARY_RT_LT_SET(int16_t, LT);   break;
#else
#define BINARY_LT_IS_INT16
#define BINARY_LT_SET_RT_IS_INT16(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT32 == 1
#define BINARY_LT_IS_UINT32                                        \
  case GAL_TYPE_UINT32: BINARY_LT_SET(uint32_t);        break;
#define BINARY_LT_SET_RT_IS_UINT32(LT)                             \
  case GAL_TYPE_UINT32: BINARY_RT_LT_SET(uint32_t, LT); break;
#else
#define BINARY_LT_IS_UINT32
#define BINARY_LT_SET_RT_IS_UINT32(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT32 == 1
#define BINARY_LT_IS_INT32                                         \
  case GAL_TYPE_INT32: BINARY_LT_SET(int32_t);          break;
#define BINARY_LT_SET_RT_IS_INT32(LT)                              \
  case GAL_TYPE_INT32: BINARY_RT_LT_SET(int32_t, LT);   break;
#else
#define BINARY_LT_IS_INT32
#define BINARY_LT_SET_RT_IS_INT32(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT64 == 1
#define BINARY_LT_IS_UINT64                                        \
  case GAL_TYPE_UINT64: BINARY_LT_SET(uint64_t);        break;
#define BINARY_LT_SET_RT_IS_UINT64(LT)                             \
  case GAL_TYPE_UINT64: BINARY_RT_LT_SET(uint64_t, LT); break;
#else
#define BINARY_LT_IS_UINT64
#define BINARY_LT_SET_RT_IS_UINT64(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT64 == 1
#define BINARY_LT_IS_INT64                                         \
  case GAL_TYPE_INT64: BINARY_LT_SET(int64_t);          break;
#define BINARY_LT_SET_RT_IS_INT64(LT)                              \
  case GAL_TYPE_INT64: BINARY_RT_LT_SET(int64_t, LT);   break;
#else
#define BINARY_LT_IS_INT64
#define BINARY_LT_SET_RT_IS_INT64(LT)
#endif





#if GAL_CONFIG_BIN_OP_FLOAT32 == 1
#define BINARY_LT_IS_FLOAT32                                       \
  case GAL_TYPE_FLOAT32: BINARY_LT_SET(float);          break;
#define BINARY_LT_SET_RT_IS_FLOAT32(LT)                            \
  case GAL_TYPE_FLOAT32: BINARY_RT_LT_SET(float, LT);   break;
#else
#define BINARY_LT_IS_FLOAT32
#define BINARY_LT_SET_RT_IS_FLOAT32(LT)
#endif





#if GAL_CONFIG_BIN_OP_FLOAT64 == 1
#define BINARY_LT_IS_FLOAT64                                       \
  case GAL_TYPE_FLOAT64: BINARY_LT_SET(double);         break;
#define BINARY_LT_SET_RT_IS_FLOAT64(LT)                            \
  case GAL_TYPE_FLOAT64: BINARY_RT_LT_SET(double, LT);  break;
#else
#define BINARY_LT_IS_FLOAT64
#define BINARY_LT_SET_RT_IS_FLOAT64(LT)
#endif



















/************************************************************************/
/*************              High level macros           *****************/
/************************************************************************/
/* Final step to be used by all operators and all types. */
#define BINARY_OP_OT_RT_LT_SET(OP, OT, RT, LT) {                        \
    LT lb, *la=l->array;                                                \
    RT rb, *ra=r->array;                                                \
    OT ob, *oa=o->array, *of=oa + o->size;                              \
    if(checkblank)                                                      \
      {                                                                 \
        gal_blank_write(&lb, l->type);                                  \
        gal_blank_write(&rb, r->type);                                  \
        gal_blank_write(&ob, o->type);                                  \
        do                                                              \
          {                                                             \
            if(lb==lb && rb==rb)/* Both are integers.                */ \
              *oa = (*la!=lb  && *ra!=rb)  ? *la OP *ra : ob ;          \
            else if(lb==lb)     /* Only left operand is an integer.  */ \
              *oa = (*la!=lb  && *ra==*ra) ? *la OP *ra : ob;           \
            else                /* Only right operand is an integer. */ \
              *oa = (*la==*la && *ra!=rb)  ? *la OP *ra : ob;           \
            if(l->size>1) ++la;                                         \
            if(r->size>1) ++ra;                                         \
          }                                                             \
        while(++oa<of);                                                 \
      }                                                                 \
    else                                                                \
      {                                                                 \
        if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);   \
        else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);   \
        else                 do *oa = *la++ OP *ra;   while(++oa<of);   \
      }                                                                 \
  }




/* This is for operators like `&&' and `||', where the right operator is
   not necessarily read (and thus incremented). */
#define BINARY_OP_INCR_OT_RT_LT_SET(OP, OT, RT, LT) {                   \
    LT *la=l->array;                                                    \
    RT *ra=r->array;                                                    \
    OT *oa=o->array, *of=oa + o->size;                                  \
    if(l->size==r->size) do {*oa = *la++ OP *ra; ++ra;} while(++oa<of); \
    else if(l->size==1)  do {*oa = *la   OP *ra; ++ra;} while(++oa<of); \
    else                 do  *oa = *la++ OP *ra;        while(++oa<of); \
  }





/* For operators whose type may be any of the given inputs. */
#define BINARY_OP_RT_LT_SET(OP, RT, LT)                            \
  if(o->type==l->type)                                             \
    BINARY_OP_OT_RT_LT_SET(OP, LT, RT, LT)                         \
  else                                                             \
    BINARY_OP_OT_RT_LT_SET(OP, RT, RT, LT)





/* Left and right types set, choose what to do based on operator. */
#define BINARY_RT_LT_SET(RT, LT)                                    \
  switch(operator)                                                  \
    {                                                               \
    case GAL_ARITHMETIC_OP_PLUS:                                    \
      BINARY_OP_RT_LT_SET(          +, RT, LT);                     \
      break;                                                        \
    case GAL_ARITHMETIC_OP_MINUS:                                   \
      BINARY_OP_RT_LT_SET(          -, RT, LT);                     \
      break;                                                        \
    case GAL_ARITHMETIC_OP_MULTIPLY:                                \
      BINARY_OP_RT_LT_SET(          *, RT, LT);                     \
      break;                                                        \
    case GAL_ARITHMETIC_OP_DIVIDE:                                  \
      BINARY_OP_RT_LT_SET(          /, RT, LT);                     \
      break;                                                        \
    case GAL_ARITHMETIC_OP_LT:                                      \
      BINARY_OP_OT_RT_LT_SET(       <, uint8_t, RT, LT);            \
      break;                                                        \
    case GAL_ARITHMETIC_OP_LE:                                      \
      BINARY_OP_OT_RT_LT_SET(      <=, uint8_t, RT, LT);            \
      break;                                                        \
    case GAL_ARITHMETIC_OP_GT:                                      \
      BINARY_OP_OT_RT_LT_SET(       >, uint8_t, RT, LT);            \
      break;                                                        \
    case GAL_ARITHMETIC_OP_GE:                                      \
      BINARY_OP_OT_RT_LT_SET(      >=, uint8_t, RT, LT);            \
      break;                                                        \
    case GAL_ARITHMETIC_OP_EQ:                                      \
      BINARY_OP_OT_RT_LT_SET(      ==, uint8_t, RT, LT);            \
      break;                                                        \
    case GAL_ARITHMETIC_OP_NE:                                      \
      BINARY_OP_OT_RT_LT_SET(      !=, uint8_t, RT, LT);            \
      break;                                                        \
    case GAL_ARITHMETIC_OP_AND:                                     \
      BINARY_OP_INCR_OT_RT_LT_SET( &&, uint8_t, RT, LT);            \
      break;                                                        \
    case GAL_ARITHMETIC_OP_OR:                                      \
      BINARY_OP_INCR_OT_RT_LT_SET( ||, uint8_t, RT, LT);            \
      break;                                                        \
    default:                                                        \
      error(EXIT_FAILURE, 0, "%s: operator code %d not recognized", \
            "BINARY_RT_LT_SET", operator);                          \
    }






/* Left operand type set, see what the right operand type is. */
#define BINARY_LT_SET(LT)                                          \
  switch(r->type)                                                  \
    {                                                              \
      BINARY_LT_SET_RT_IS_UINT8(LT);                               \
      BINARY_LT_SET_RT_IS_INT8(LT);                                \
      BINARY_LT_SET_RT_IS_UINT16(LT);                              \
      BINARY_LT_SET_RT_IS_INT16(LT);                               \
      BINARY_LT_SET_RT_IS_UINT32(LT);                              \
      BINARY_LT_SET_RT_IS_INT32(LT);                               \
      BINARY_LT_SET_RT_IS_UINT64(LT);                              \
      BINARY_LT_SET_RT_IS_INT64(LT);                               \
      BINARY_LT_SET_RT_IS_FLOAT32(LT);                             \
      BINARY_LT_SET_RT_IS_FLOAT64(LT);                             \
    default:                                                       \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",    \
            "BINARY_LT_SET", r->type);                             \
    }




















/************************************************************************/
/*************              Top level function          *****************/
/************************************************************************/
gal_data_t *
arithmetic_binary(int operator, int flags, gal_data_t *lo, gal_data_t *ro)
{
  /* Read the variable arguments. `lo' and `ro' keep the original data, in
     case their type isn't built (based on configure options are configure
     time). */
  int checkblank;
  int32_t otype, final_otype;
  size_t out_size, minmapsize;
  gal_data_t *l, *r, *o=NULL, *tmp_o;


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_ARITHMETIC_NUMOK) && (lo->size==1 || ro->size==1))
      && gal_data_dsize_is_different(lo, ro) )
    error(EXIT_FAILURE, 0, "%s: the non-number inputs to %s don't have the "
          "same dimension/size", __func__,
          gal_arithmetic_operator_string(operator));


  /* Set the final output type (independent of which types are
     compiled). These needs to be done before the call to
     `gal_arithmetic_convert_to_compiled_type', because that function can
     free the space of the original data structures, thus we will loose the
     original data structure information. */
  final_otype=gal_arithmetic_binary_out_type(operator, lo, ro);


  /* Make sure the input arrays have one of the compiled types. From this
     point on, until the cleaning up section of this function, we won't be
     using the `lo' and `ro' pointers. */
  l=gal_arithmetic_convert_to_compiled_type(lo, flags);
  r=gal_arithmetic_convert_to_compiled_type(ro, flags);


  /* Set the output type. For the comparison operators, the output type is
     either 0 or 1, so we will set the output type to `unsigned char' for
     efficient memory and CPU usage. Since the number of operators without
     a fixed output type (like the conditionals) is less, by `default' we
     will set the output type to `unsigned char', and if any of the other
     operatrs are given, it will be chosen based on the input types.*/
  otype=gal_arithmetic_binary_out_type(operator, l, r);


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


  /* If the output pointer was not set for any reason, allocate it. For
     `mmapsize', note that since its `size_t', it will always be
     Positive. The `-1' that is recommended to give when you want the value
     in RAM is actually the largest possible memory location. So we just
     have to choose the smaller minmapsize of the two to decide if the
     output array should be in RAM or not. */
  if(o==NULL)
    o = gal_data_alloc(NULL, otype,
                       l->size>1 ? l->ndim  : r->ndim,
                       l->size>1 ? l->dsize : r->dsize,
                       l->size>1 ? l->wcs   : r->wcs,
                       0, minmapsize, NULL, NULL, NULL );


  /* See if we should check for blanks. When both types are floats, blanks
     don't need to be checked (the floating point standard will do the job
     for us). It is also not necessary to check blanks in bitwise
     operators, but bitwise operators have their own macro
     (`BINARY_OP_INCR_OT_RT_LT_SET') which doesn' use `checkblanks'.*/
  checkblank = ((((l->type!=GAL_TYPE_FLOAT32    && l->type!=GAL_TYPE_FLOAT64)
                  || (r->type!=GAL_TYPE_FLOAT32 && r->type!=GAL_TYPE_FLOAT64))
                 && (gal_blank_present(l, 1) || gal_blank_present(r, 1)))
                ? 1 : 0 );


  /* Start setting the operator and operands. */
  switch(l->type)
    {
      BINARY_LT_IS_UINT8;
      BINARY_LT_IS_INT8;
      BINARY_LT_IS_UINT16;
      BINARY_LT_IS_INT16;
      BINARY_LT_IS_UINT32;
      BINARY_LT_IS_INT32;
      BINARY_LT_IS_UINT64;
      BINARY_LT_IS_INT64;
      BINARY_LT_IS_FLOAT32;
      BINARY_LT_IS_FLOAT64;
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, l->type);
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, `BINARY_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, so only `r' and `l' need
     freeing. Alternatively, when the inputs shouldn't be freed, the only
     allocated spaces are the `r' and `l' arrays if their types weren't
     compiled for binary operations, we can tell this from the pointers: if
     they are different from the original pointers, they were allocated. */
  if(flags & GAL_ARITHMETIC_FREE)
    {
      if     (o==l)       gal_data_free(r);
      else if(o==r)       gal_data_free(l);
      else              { gal_data_free(l); gal_data_free(r); }
    }
  else
    {
      if(l!=lo)           gal_data_free(l);
      if(r!=ro)           gal_data_free(r);
    }


  /* The type of the output dataset (`o->type') was chosen from `l' and `r'
     (copies of the orignal operands but in a compiled type, not
     necessarily the original `lo' and `ro' data structures). So we need to
     to get the final type based on the original operands and check if the
     final output needs changing.

     IMPORTANT: This has to be done after (possibly) freeing the left and
     right operands because this step can change the `o' pointer which
     they may depend on (when working in-place). */
  if( o->type != final_otype )
    {
      tmp_o=gal_data_copy_to_new_type(o, final_otype);
      gal_data_free(o);
      o=tmp_o;
    }

  /* Return */
  return o;
}
