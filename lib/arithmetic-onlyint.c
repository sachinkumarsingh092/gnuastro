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

#include <gnuastro/arithmetic.h>

#include <gnuastro-internal/arithmetic-onlyint.h>





/************************************************************************/
/*************            Native type macros            *****************/
/************************************************************************/
#if GAL_CONFIG_BIN_OP_UINT8 == 1
#define BINOIN_LT_IS_UINT8                                         \
  case GAL_TYPE_UINT8: BINOIN_LT_SET(uint8_t);          break;
#define BINOIN_LT_SET_RT_IS_UINT8(LT)                              \
  case GAL_TYPE_UINT8: BINOIN_RT_LT_SET(uint8_t, LT);   break;
#else
#define BINOIN_LT_IS_UINT8
#define BINOIN_LT_SET_RT_IS_UINT8(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT8 == 1
#define BINOIN_LT_IS_INT8                                          \
  case GAL_TYPE_INT8: BINOIN_LT_SET(int8_t);            break;
#define BINOIN_LT_SET_RT_IS_INT8(LT)                               \
  case GAL_TYPE_INT8: BINOIN_RT_LT_SET(int8_t, LT);     break;
#else
#define BINOIN_LT_IS_INT8
#define BINOIN_LT_SET_RT_IS_INT8(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT16 == 1
#define BINOIN_LT_IS_UINT16                                        \
  case GAL_TYPE_UINT16: BINOIN_LT_SET(uint16_t);        break;
#define BINOIN_LT_SET_RT_IS_UINT16(LT)                             \
  case GAL_TYPE_UINT16: BINOIN_RT_LT_SET(uint16_t, LT); break;
#else
#define BINOIN_LT_IS_UINT16
#define BINOIN_LT_SET_RT_IS_UINT16(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT16 == 1
#define BINOIN_LT_IS_INT16                                         \
  case GAL_TYPE_INT16: BINOIN_LT_SET(int16_t);          break;
#define BINOIN_LT_SET_RT_IS_INT16(LT)                              \
  case GAL_TYPE_INT16: BINOIN_RT_LT_SET(int16_t, LT);   break;
#else
#define BINOIN_LT_IS_INT16
#define BINOIN_LT_SET_RT_IS_INT16(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT32 == 1
#define BINOIN_LT_IS_UINT32                                        \
  case GAL_TYPE_UINT32: BINOIN_LT_SET(uint32_t);        break;
#define BINOIN_LT_SET_RT_IS_UINT32(LT)                             \
  case GAL_TYPE_UINT32: BINOIN_RT_LT_SET(uint32_t, LT); break;
#else
#define BINOIN_LT_IS_UINT32
#define BINOIN_LT_SET_RT_IS_UINT32(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT32 == 1
#define BINOIN_LT_IS_INT32                                         \
  case GAL_TYPE_INT32: BINOIN_LT_SET(int32_t);          break;
#define BINOIN_LT_SET_RT_IS_INT32(LT)                              \
  case GAL_TYPE_INT32: BINOIN_RT_LT_SET(int32_t, LT);   break;
#else
#define BINOIN_LT_IS_INT32
#define BINOIN_LT_SET_RT_IS_INT32(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT64 == 1
#define BINOIN_LT_IS_UINT64                                        \
  case GAL_TYPE_UINT64: BINOIN_LT_SET(uint64_t);        break;
#define BINOIN_LT_SET_RT_IS_UINT64(LT)                             \
  case GAL_TYPE_UINT64: BINOIN_RT_LT_SET(uint64_t, LT); break;
#else
#define BINOIN_LT_IS_UINT64
#define BINOIN_LT_SET_RT_IS_UINT64(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT64 == 1
#define BINOIN_LT_IS_INT64                                         \
  case GAL_TYPE_INT64: BINOIN_LT_SET(int64_t);          break;
#define BINOIN_LT_SET_RT_IS_INT64(LT)                              \
  case GAL_TYPE_INT64: BINOIN_RT_LT_SET(int64_t, LT);   break;
#else
#define BINOIN_LT_IS_INT64
#define BINOIN_LT_SET_RT_IS_INT64(LT)
#endif



















/************************************************************************/
/*************              High level macros           *****************/
/************************************************************************/
/* Final step to be used by all operators and all types. */
#define BINOIN_OP_OT_RT_LT_SET(OP, OT, RT, LT) {                   \
    LT *la=l->array;                                               \
    RT *ra=r->array;                                               \
    OT *oa=o->array, *of=oa + o->size;                             \
    if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);  \
    else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);  \
    else                 do *oa = *la++ OP *ra;   while(++oa<of);  \
  }





/* For operators whose type may be any of the given inputs. */
#define BINOIN_OP_RT_LT_SET(OP, RT, LT)                            \
  if(o->type==l->type)                                             \
    BINOIN_OP_OT_RT_LT_SET(OP, LT, RT, LT)                         \
  else                                                             \
    BINOIN_OP_OT_RT_LT_SET(OP, RT, RT, LT)





/* Left and right types set, choose what to do based on operator. */
#define BINOIN_RT_LT_SET(RT, LT)                                   \
  switch(operator)                                                 \
    {                                                              \
    case GAL_ARITHMETIC_OP_MODULO:                                 \
      BINOIN_OP_RT_LT_SET(%, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_BITAND:                                 \
      BINOIN_OP_RT_LT_SET(&, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_BITOR:                                  \
      BINOIN_OP_RT_LT_SET(|, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_BITXOR:                                 \
      BINOIN_OP_RT_LT_SET(^, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_BITLSH:                                 \
      BINOIN_OP_RT_LT_SET(<<, RT, LT);                             \
      break;                                                       \
    case GAL_ARITHMETIC_OP_BITRSH:                                 \
      BINOIN_OP_RT_LT_SET(>>, RT, LT);                             \
      break;                                                       \
    default:                                                       \
      error(EXIT_FAILURE, 0, "operator code %d not recognized in " \
            "`BINOIN_RT_LT_SET", operator);                        \
    }






/* Left operand type set, see what the right operand type is. */
#define BINOIN_LT_SET(LT)                                          \
  switch(r->type)                                                  \
    {                                                              \
      BINOIN_LT_SET_RT_IS_UINT8(LT);                               \
      BINOIN_LT_SET_RT_IS_INT8(LT);                                \
      BINOIN_LT_SET_RT_IS_UINT16(LT);                              \
      BINOIN_LT_SET_RT_IS_INT16(LT);                               \
      BINOIN_LT_SET_RT_IS_UINT32(LT);                              \
      BINOIN_LT_SET_RT_IS_INT32(LT);                               \
      BINOIN_LT_SET_RT_IS_UINT64(LT);                              \
      BINOIN_LT_SET_RT_IS_INT64(LT);                               \
    default:                                                       \
      error(EXIT_FAILURE, 0, "type code %d not recognized in "     \
            "`BINOIN_LT_SET'", r->type);                           \
    }




















/************************************************************************/
/*************              Top level function          *****************/
/************************************************************************/
gal_data_t *
arithmetic_onlyint_binary(int operator, unsigned char flags,
                          gal_data_t *lo, gal_data_t *ro)
{
  /* Read the variable arguments. `lo' and `ro' keep the original data, in
     case their type isn't built (based on configure options are configure
     time). */
  int otype, final_otype;
  size_t out_size, minmapsize;
  gal_data_t *l, *r, *o=NULL, *tmp_o;
  char *opstring=gal_arithmetic_operator_string(operator);


  /* Simple sanity check on the input sizes and types */
  if( !( (flags & GAL_ARITHMETIC_NUMOK) && (lo->size==1 || ro->size==1))
      && gal_data_dsize_is_different(lo, ro) )
    error(EXIT_FAILURE, 0, "the non-number inputs to %s don't have the "
          "same dimension/size", opstring);

  if( lo->type==GAL_TYPE_FLOAT32 || lo->type==GAL_TYPE_FLOAT64
      || ro->type==GAL_TYPE_FLOAT32 || ro->type==GAL_TYPE_FLOAT64 )
      error(EXIT_FAILURE, 0, "the %s operator can only work on integer "
            "type operands", opstring);


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


  /* Sanity check: see if the compiled type is actually an integer. */
  if( l->type>=GAL_TYPE_FLOAT32 || r->type>=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "no larger integer compiled type. The `%s' "
          "operator can only work on integer types. The left and right "
          "operands had types `%s' and `%s'.\n\nYou can use the "
          "`--enable-bin-op-XXXX' at configure time to compile a larger "
          "type (note that unsigned types are considered to be larger than "
          "signed ones). You can run the following command for more "
          "information on these options (press the `SPACE' key to go down "
          "and `q' to return to the command-line):\n\n"
          "    $ info gnuastro \"Gnuastro configure options\"\n",
          gal_arithmetic_operator_string(operator),
          gal_type_to_string(lo->type, 1), gal_type_to_string(ro->type, 1));

  /* Set the output type. */
  otype=gal_type_out(l->type, r->type);


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


  /* Start setting the operator and operands. */
  switch(l->type)
    {
      BINOIN_LT_IS_UINT8;
      BINOIN_LT_IS_INT8;
      BINOIN_LT_IS_UINT16;
      BINOIN_LT_IS_INT16;
      BINOIN_LT_IS_UINT32;
      BINOIN_LT_IS_INT32;
      BINOIN_LT_IS_UINT64;
      BINOIN_LT_IS_INT64;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`arithmetic_onlyint_binary'", l->type);
    }


  /* The type of the output dataset (`o->type') was chosen from `l' and `r'
     (copies of the orignal operands but in a compiled type, not
     necessarily the original `lo' and `ro' data structures). So we need to
     to get the final type based on the original operands and check if the
     final output needs changing. */
  if( o->type != final_otype )
    {
      tmp_o=gal_data_copy_to_new_type(o, final_otype);
      gal_data_free(o);
      o=tmp_o;
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, `BINOIN_CONVERT_TO_COMPILED_TYPE'
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

  /* Return */
  return o;
}





gal_data_t *
arithmetic_onlyint_bitwise_not(unsigned char flags, gal_data_t *in)
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
      error(EXIT_FAILURE, 0, "the bitwise not (one's complement) "
            "operator can only work on integer types");
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
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "data_arithmetic_bitwise_not", in->type);
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
