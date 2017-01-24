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






/************************************************************************/
/*************            Native type macros            *****************/
/************************************************************************/
#if GAL_CONFIG_BIN_OP_UCHAR == 1
#define BINARY_LT_IS_UCHAR                                         \
  case GAL_DATA_TYPE_UCHAR:                                        \
    BINARY_LT_SET(unsigned char);                                  \
    break;
#define BINARY_LT_SET_RT_IS_UCHAR(LT)                              \
  case GAL_DATA_TYPE_UCHAR:                                        \
    BINARY_RT_LT_SET(unsigned char, LT);                           \
    break;
#else
#define BINARY_LT_IS_UCHAR
#define BINARY_LT_SET_RT_IS_UCHAR(LT)
#endif





#if GAL_CONFIG_BIN_OP_CHAR == 1
#define BINARY_LT_IS_CHAR                                          \
  case GAL_DATA_TYPE_CHAR:                                         \
    BINARY_LT_SET(char);                                           \
    break;
#define BINARY_LT_SET_RT_IS_CHAR(LT)                               \
  case GAL_DATA_TYPE_CHAR:                                         \
    BINARY_RT_LT_SET(char, LT);                                    \
    break;
#else
#define BINARY_LT_IS_CHAR
#define BINARY_LT_SET_RT_IS_CHAR(LT)
#endif





#if GAL_CONFIG_BIN_OP_USHORT == 1
#define BINARY_LT_IS_USHORT                                        \
  case GAL_DATA_TYPE_USHORT:                                       \
    BINARY_LT_SET(unsigned short);                                 \
    break;
#define BINARY_LT_SET_RT_IS_USHORT(LT)                             \
  case GAL_DATA_TYPE_USHORT:                                       \
    BINARY_RT_LT_SET(unsigned short, LT);                          \
    break;
#else
#define BINARY_LT_IS_USHORT
#define BINARY_LT_SET_RT_IS_USHORT(LT)
#endif





#if GAL_CONFIG_BIN_OP_SHORT == 1
#define BINARY_LT_IS_SHORT                                         \
  case GAL_DATA_TYPE_SHORT:                                        \
    BINARY_LT_SET(short);                                          \
    break;
#define BINARY_LT_SET_RT_IS_SHORT(LT)                              \
  case GAL_DATA_TYPE_SHORT:                                        \
    BINARY_RT_LT_SET(short, LT);                                   \
    break;
#else
#define BINARY_LT_IS_SHORT
#define BINARY_LT_SET_RT_IS_SHORT(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT == 1
#define BINARY_LT_IS_UINT                                          \
  case GAL_DATA_TYPE_UINT:                                         \
    BINARY_LT_SET(unsigned int);                                   \
    break;
#define BINARY_LT_SET_RT_IS_UINT(LT)                               \
  case GAL_DATA_TYPE_UINT:                                         \
    BINARY_RT_LT_SET(unsigned int, LT);                            \
    break;
#else
#define BINARY_LT_IS_UINT
#define BINARY_LT_SET_RT_IS_UINT(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT == 1
#define BINARY_LT_IS_INT                                           \
  case GAL_DATA_TYPE_INT:                                          \
    BINARY_LT_SET(int);                                            \
    break;
#define BINARY_LT_SET_RT_IS_INT(LT)                                \
  case GAL_DATA_TYPE_INT:                                          \
    BINARY_RT_LT_SET(int, LT);                                     \
    break;
#else
#define BINARY_LT_IS_INT
#define BINARY_LT_SET_RT_IS_INT(LT)
#endif





#if GAL_CONFIG_BIN_OP_ULONG == 1
#define BINARY_LT_IS_ULONG                                         \
  case GAL_DATA_TYPE_ULONG:                                        \
    BINARY_LT_SET(unsigned long);                                  \
    break;
#define BINARY_LT_SET_RT_IS_ULONG(LT)                              \
  case GAL_DATA_TYPE_ULONG:                                        \
    BINARY_RT_LT_SET(unsigned long, LT);                           \
    break;
#else
#define BINARY_LT_IS_ULONG
#define BINARY_LT_SET_RT_IS_ULONG(LT)
#endif





#if GAL_CONFIG_BIN_OP_LONG == 1
#define BINARY_LT_IS_LONG                                          \
  case GAL_DATA_TYPE_LONG:                                         \
    BINARY_LT_SET(long);                                           \
    break;
#define BINARY_LT_SET_RT_IS_LONG(LT)                               \
  case GAL_DATA_TYPE_LONG:                                         \
    BINARY_RT_LT_SET(long, LT);                                    \
    break;
#else
#define BINARY_LT_IS_LONG
#define BINARY_LT_SET_RT_IS_LONG(LT)
#endif





#if GAL_CONFIG_BIN_OP_LONGLONG == 1
#define BINARY_LT_IS_LONGLONG                                      \
  case GAL_DATA_TYPE_LONGLONG:                                     \
    BINARY_LT_SET(LONGLONG);                                       \
    break;
#define BINARY_LT_SET_RT_IS_LONGLONG(LT)                           \
  case GAL_DATA_TYPE_LONGLONG:                                     \
    BINARY_RT_LT_SET(LONGLONG, LT);                                \
    break;
#else
#define BINARY_LT_IS_LONGLONG
#define BINARY_LT_SET_RT_IS_LONGLONG(LT)
#endif





#if GAL_CONFIG_BIN_OP_FLOAT == 1
#define BINARY_LT_IS_FLOAT                                         \
  case GAL_DATA_TYPE_FLOAT:                                        \
    BINARY_LT_SET(float);                                          \
    break;
#define BINARY_LT_SET_RT_IS_FLOAT(LT)                              \
  case GAL_DATA_TYPE_FLOAT:                                        \
    BINARY_RT_LT_SET(float, LT);                                   \
    break;
#else
#define BINARY_LT_IS_FLOAT
#define BINARY_LT_SET_RT_IS_FLOAT(LT)
#endif





#if GAL_CONFIG_BIN_OP_DOUBLE == 1
#define BINARY_LT_IS_DOUBLE                                        \
  case GAL_DATA_TYPE_DOUBLE:                                       \
    BINARY_LT_SET(double);                                         \
    break;
#define BINARY_LT_SET_RT_IS_DOUBLE(LT)                             \
  case GAL_DATA_TYPE_DOUBLE:                                       \
    BINARY_RT_LT_SET(double, LT);                                  \
    break;
#else
#define BINARY_LT_IS_DOUBLE
#define BINARY_LT_SET_RT_IS_DOUBLE(LT)
#endif



















/************************************************************************/
/*************              High level macros           *****************/
/************************************************************************/
/* Final step to be used by all operators and all types. */
#define BINARY_OP_OT_RT_LT_SET(OP, OT, RT, LT) {                   \
    LT *la=l->array;                                               \
    RT *ra=r->array;                                               \
    OT *oa=o->array, *of=oa + o->size;                             \
    if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);  \
    else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);  \
    else                 do *oa = *la++ OP *ra;   while(++oa<of);  \
  }




/* This is for operators like `&&' and `||', where the right operator is
   not necessarily read (and thus incremented) incremented. */
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
#define BINARY_RT_LT_SET(RT, LT)                                   \
  switch(operator)                                                 \
    {                                                              \
    case GAL_ARITHMETIC_OP_PLUS:                                   \
      BINARY_OP_RT_LT_SET(+, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_MINUS:                                  \
      BINARY_OP_RT_LT_SET(-, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_MULTIPLY:                               \
      BINARY_OP_RT_LT_SET(*, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_DIVIDE:                                 \
      BINARY_OP_RT_LT_SET(/, RT, LT);                              \
      break;                                                       \
    case GAL_ARITHMETIC_OP_LT:                                     \
      BINARY_OP_OT_RT_LT_SET(<, unsigned char, RT, LT);            \
      break;                                                       \
    case GAL_ARITHMETIC_OP_LE:                                     \
      BINARY_OP_OT_RT_LT_SET(<=, unsigned char, RT, LT);           \
      break;                                                       \
    case GAL_ARITHMETIC_OP_GT:                                     \
      BINARY_OP_OT_RT_LT_SET(>, unsigned char, RT, LT);            \
      break;                                                       \
    case GAL_ARITHMETIC_OP_GE:                                     \
      BINARY_OP_OT_RT_LT_SET(>=, unsigned char, RT, LT);           \
      break;                                                       \
    case GAL_ARITHMETIC_OP_EQ:                                     \
      BINARY_OP_OT_RT_LT_SET(==, unsigned char, RT, LT);           \
      break;                                                       \
    case GAL_ARITHMETIC_OP_NE:                                     \
      BINARY_OP_OT_RT_LT_SET(!=, unsigned char, RT, LT);           \
      break;                                                       \
    case GAL_ARITHMETIC_OP_AND:                                    \
      BINARY_OP_INCR_OT_RT_LT_SET(&&, unsigned char, RT, LT);      \
      break;                                                       \
    case GAL_ARITHMETIC_OP_OR:                                     \
      BINARY_OP_INCR_OT_RT_LT_SET(||, unsigned char, RT, LT);      \
      break;                                                       \
    default:                                                       \
      error(EXIT_FAILURE, 0, "operator code %d not recognized in " \
            "`BINARY_RT_LT_SET", operator);                        \
    }






/* Left operand type set, see what the right operand type is. */
#define BINARY_LT_SET(LT)                                          \
  switch(r->type)                                                  \
    {                                                              \
      BINARY_LT_SET_RT_IS_UCHAR(LT);                               \
      BINARY_LT_SET_RT_IS_CHAR(LT);                                \
      BINARY_LT_SET_RT_IS_USHORT(LT);                              \
      BINARY_LT_SET_RT_IS_SHORT(LT);                               \
      BINARY_LT_SET_RT_IS_UINT(LT);                                \
      BINARY_LT_SET_RT_IS_INT(LT);                                 \
      BINARY_LT_SET_RT_IS_ULONG(LT);                               \
      BINARY_LT_SET_RT_IS_LONG(LT);                                \
      BINARY_LT_SET_RT_IS_LONGLONG(LT);                            \
      BINARY_LT_SET_RT_IS_FLOAT(LT);                               \
      BINARY_LT_SET_RT_IS_DOUBLE(LT);                              \
    default:                                                       \
      error(EXIT_FAILURE, 0, "type code %d not recognized in "     \
            "`BINARY_LT_SET'", r->type);                           \
    }




















/************************************************************************/
/*************              Top level function          *****************/
/************************************************************************/
gal_data_t *
arithmetic_binary(int operator, unsigned char flags, gal_data_t *lo,
                  gal_data_t *ro)
{
  /* Read the variable arguments. `lo' and `ro' keep the original data, in
     case their type isn't built (based on configure options are configure
     time). */
  int otype, final_otype;
  size_t out_size, minmapsize;
  gal_data_t *l, *r, *o=NULL, *tmp_o;


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_ARITHMETIC_NUMOK) && (lo->size==1 || ro->size==1))
      && gal_data_dsize_is_different(lo, ro) )
    error(EXIT_FAILURE, 0, "the non-number inputs to %s don't have the "
          "same dimension/size", gal_arithmetic_operator_string(operator));


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


  /* Start setting the operator and operands. */
  switch(l->type)
    {
      BINARY_LT_IS_UCHAR;
      BINARY_LT_IS_CHAR;
      BINARY_LT_IS_USHORT;
      BINARY_LT_IS_SHORT;
      BINARY_LT_IS_UINT;
      BINARY_LT_IS_INT;
      BINARY_LT_IS_ULONG;
      BINARY_LT_IS_LONG;
      BINARY_LT_IS_LONGLONG;
      BINARY_LT_IS_FLOAT;
      BINARY_LT_IS_DOUBLE;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`data_arithmetic_binary'", l->type);
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

  /* Return */
  return o;
}
