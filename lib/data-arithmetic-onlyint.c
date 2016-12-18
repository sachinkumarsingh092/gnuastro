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

#include <gnuastro/data.h>
#include <data-arithmetic-onlyint.h>





/************************************************************************/
/*************            Native type macros            *****************/
/************************************************************************/
#if GAL_CONFIG_BIN_OP_UCHAR == 1
#define BINOIN_LT_IS_UCHAR                                         \
  case GAL_DATA_TYPE_UCHAR:                                        \
    BINOIN_LT_SET(unsigned char);                                  \
    break;
#define BINOIN_LT_SET_RT_IS_UCHAR(LT)                              \
  case GAL_DATA_TYPE_UCHAR:                                        \
    BINOIN_RT_LT_SET(unsigned char, LT);                           \
    break;
#else
#define BINOIN_LT_IS_UCHAR
#define BINOIN_LT_SET_RT_IS_UCHAR(LT)
#endif





#if GAL_CONFIG_BIN_OP_CHAR == 1
#define BINOIN_LT_IS_CHAR                                          \
  case GAL_DATA_TYPE_CHAR:                                         \
    BINOIN_LT_SET(char);                                           \
    break;
#define BINOIN_LT_SET_RT_IS_CHAR(LT)                               \
  case GAL_DATA_TYPE_CHAR:                                         \
    BINOIN_RT_LT_SET(char, LT);                                    \
    break;
#else
#define BINOIN_LT_IS_CHAR
#define BINOIN_LT_SET_RT_IS_CHAR(LT)
#endif





#if GAL_CONFIG_BIN_OP_USHORT == 1
#define BINOIN_LT_IS_USHORT                                        \
  case GAL_DATA_TYPE_USHORT:                                       \
    BINOIN_LT_SET(unsigned short);                                 \
    break;
#define BINOIN_LT_SET_RT_IS_USHORT(LT)                             \
  case GAL_DATA_TYPE_USHORT:                                       \
    BINOIN_RT_LT_SET(unsigned short, LT);                          \
    break;
#else
#define BINOIN_LT_IS_USHORT
#define BINOIN_LT_SET_RT_IS_USHORT(LT)
#endif





#if GAL_CONFIG_BIN_OP_SHORT == 1
#define BINOIN_LT_IS_SHORT                                         \
  case GAL_DATA_TYPE_SHORT:                                        \
    BINOIN_LT_SET(short);                                          \
    break;
#define BINOIN_LT_SET_RT_IS_SHORT(LT)                              \
  case GAL_DATA_TYPE_SHORT:                                        \
    BINOIN_RT_LT_SET(short, LT);                                   \
    break;
#else
#define BINOIN_LT_IS_SHORT
#define BINOIN_LT_SET_RT_IS_SHORT(LT)
#endif





#if GAL_CONFIG_BIN_OP_UINT == 1
#define BINOIN_LT_IS_UINT                                          \
  case GAL_DATA_TYPE_UINT:                                         \
    BINOIN_LT_SET(unsigned int);                                   \
    break;
#define BINOIN_LT_SET_RT_IS_UINT(LT)                               \
  case GAL_DATA_TYPE_UINT:                                         \
    BINOIN_RT_LT_SET(unsigned int, LT);                            \
    break;
#else
#define BINOIN_LT_IS_UINT
#define BINOIN_LT_SET_RT_IS_UINT(LT)
#endif





#if GAL_CONFIG_BIN_OP_INT == 1
#define BINOIN_LT_IS_INT                                           \
  case GAL_DATA_TYPE_INT:                                          \
    BINOIN_LT_SET(int);                                            \
    break;
#define BINOIN_LT_SET_RT_IS_INT(LT)                                \
  case GAL_DATA_TYPE_INT:                                          \
    BINOIN_RT_LT_SET(int, LT);                                     \
    break;
#else
#define BINOIN_LT_IS_INT
#define BINOIN_LT_SET_RT_IS_INT(LT)
#endif





#if GAL_CONFIG_BIN_OP_ULONG == 1
#define BINOIN_LT_IS_ULONG                                         \
  case GAL_DATA_TYPE_ULONG:                                        \
    BINOIN_LT_SET(unsigned long);                                  \
    break;
#define BINOIN_LT_SET_RT_IS_ULONG(LT)                              \
  case GAL_DATA_TYPE_ULONG:                                        \
    BINOIN_RT_LT_SET(unsigned long, LT);                           \
    break;
#else
#define BINOIN_LT_IS_ULONG
#define BINOIN_LT_SET_RT_IS_ULONG(LT)
#endif





#if GAL_CONFIG_BIN_OP_LONG == 1
#define BINOIN_LT_IS_LONG                                          \
  case GAL_DATA_TYPE_LONG:                                         \
    BINOIN_LT_SET(long);                                           \
    break;
#define BINOIN_LT_SET_RT_IS_LONG(LT)                               \
  case GAL_DATA_TYPE_LONG:                                         \
    BINOIN_RT_LT_SET(long, LT);                                    \
    break;
#else
#define BINOIN_LT_IS_LONG
#define BINOIN_LT_SET_RT_IS_LONG(LT)
#endif





#if GAL_CONFIG_BIN_OP_LONGLONG == 1
#define BINOIN_LT_IS_LONGLONG                                      \
  case GAL_DATA_TYPE_LONGLONG:                                     \
    BINOIN_LT_SET(LONGLONG);                                       \
    break;
#define BINOIN_LT_SET_RT_IS_LONGLONG(LT)                           \
  case GAL_DATA_TYPE_LONGLONG:                                     \
    BINOIN_RT_LT_SET(LONGLONG, LT);                                \
    break;
#else
#define BINOIN_LT_IS_LONGLONG
#define BINOIN_LT_SET_RT_IS_LONGLONG(LT)
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
    case GAL_DATA_OPERATOR_MODULO:                                 \
      BINOIN_OP_RT_LT_SET(%, RT, LT);                              \
      break;                                                       \
    case GAL_DATA_OPERATOR_BITAND:                                 \
      BINOIN_OP_RT_LT_SET(&, RT, LT);                              \
      break;                                                       \
    case GAL_DATA_OPERATOR_BITOR:                                  \
      BINOIN_OP_RT_LT_SET(|, RT, LT);                              \
      break;                                                       \
    case GAL_DATA_OPERATOR_BITXOR:                                 \
      BINOIN_OP_RT_LT_SET(^, RT, LT);                              \
      break;                                                       \
    case GAL_DATA_OPERATOR_BITLSH:                                 \
      BINOIN_OP_RT_LT_SET(<<, RT, LT);                             \
      break;                                                       \
    case GAL_DATA_OPERATOR_BITRSH:                                 \
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
      BINOIN_LT_SET_RT_IS_UCHAR(LT);                               \
      BINOIN_LT_SET_RT_IS_CHAR(LT);                                \
      BINOIN_LT_SET_RT_IS_USHORT(LT);                              \
      BINOIN_LT_SET_RT_IS_SHORT(LT);                               \
      BINOIN_LT_SET_RT_IS_UINT(LT);                                \
      BINOIN_LT_SET_RT_IS_INT(LT);                                 \
      BINOIN_LT_SET_RT_IS_ULONG(LT);                               \
      BINOIN_LT_SET_RT_IS_LONG(LT);                                \
      BINOIN_LT_SET_RT_IS_LONGLONG(LT);                            \
    default:                                                       \
      error(EXIT_FAILURE, 0, "type code %d not recognized in "     \
            "`BINOIN_LT_SET'", r->type);                           \
    }




















/************************************************************************/
/*************              Top level function          *****************/
/************************************************************************/
gal_data_t *
data_arithmetic_onlyint_binary(int operator, unsigned char flags,
                               gal_data_t *lo, gal_data_t *ro)
{
  /* Read the variable arguments. `lo' and `ro' keep the original data, in
     case their type isn't built (based on configure options are configure
     time). */
  int otype;
  size_t out_size, minmapsize;
  gal_data_t *l, *r, *o=NULL, *tmp_o;
  char *opstring=gal_data_operator_string(operator);


  /* Simple sanity check on the input sizes and types */
  if( !( (flags & GAL_DATA_ARITH_NUMOK) && (lo->size==1 || ro->size==1))
      && gal_data_dsize_is_different(lo, ro) )
    error(EXIT_FAILURE, 0, "the non-number inputs to %s don't have the "
          "same dimension/size", opstring);

  if( lo->type==GAL_DATA_TYPE_FLOAT || lo->type==GAL_DATA_TYPE_DOUBLE
      || ro->type==GAL_DATA_TYPE_FLOAT || ro->type==GAL_DATA_TYPE_DOUBLE )
      error(EXIT_FAILURE, 0, "the %s operator can only work on integer "
            "type operands", opstring);


  /* Make sure the input arrays have one of the compiled types. From this
     point on, until the cleaning up section of this function, we won't be
     using the `lo' and `ro' pointers. */
  l=data_arithmetic_convert_to_compiled_type(lo, flags);
  r=data_arithmetic_convert_to_compiled_type(ro, flags);


  /* Set the output type. */
  otype=gal_data_out_type(l, r);


  /* Set the output sizes. */
  minmapsize = ( l->minmapsize < r->minmapsize
                 ? l->minmapsize : r->minmapsize );
  out_size = l->size > r->size ? l->size : r->size;


  /* If we want inplace output, set the output pointer to one input. Note
     that the output type can be different from both inputs.  */
  if(flags & GAL_DATA_ARITH_INPLACE)
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
      BINOIN_LT_IS_UCHAR;
      BINOIN_LT_IS_CHAR;
      BINOIN_LT_IS_USHORT;
      BINOIN_LT_IS_SHORT;
      BINOIN_LT_IS_UINT;
      BINOIN_LT_IS_INT;
      BINOIN_LT_IS_ULONG;
      BINOIN_LT_IS_LONG;
      BINOIN_LT_IS_LONGLONG;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`data_arithmetic_binary'", l->type);
    }


  /* The type of the output dataset (`o->type') was chosen from `l' and `r'
     (copies of the orignal operands but in a compiled type, not
     necessarily the original `lo' and `ro' data structures). So we need to
     to get the final type based on the original operands and check if the
     final output needs changing. */
  otype=gal_data_out_type(lo, ro);
  if( o->type != otype )
    {
      tmp_o=gal_data_copy_to_new_type(o, otype);
      gal_data_free(o, 0);
      o=tmp_o;
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, `BINOIN_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, so only `r' and `l' need
     freeing. Alternatively, when the inputs shouldn't be freed, the only
     allocated spaces are the `r' and `l' arrays if their types weren't
     compiled for binary operations, we can tell this from the pointers: if
     they are different from the original pointers, they were allocated. */
  if(flags & GAL_DATA_ARITH_FREE)
    {
      if     (o==l)       gal_data_free(r, 0);
      else if(o==r)       gal_data_free(l, 0);
      else              { gal_data_free(l, 0); gal_data_free(r, 0); }
    }
  else
    {
      if(l!=lo)           gal_data_free(l, 0);
      if(r!=ro)           gal_data_free(r, 0);
    }

  /* Return */
  return o;
}





gal_data_t *
data_arithmetic_bitwise_not(unsigned char flags, gal_data_t *in)
{
  gal_data_t *o;
  unsigned char    *iuc=in->array, *iucf=in->array+in->size, *ouc;
  char              *ic=in->array,  *icf=in->array+in->size,  *oc;
  unsigned short   *ius=in->array, *iusf=in->array+in->size, *ous;
  short             *is=in->array,  *isf=in->array+in->size,  *os;
  unsigned int     *iui=in->array, *iuif=in->array+in->size, *oui;
  int               *ii=in->array,  *iif=in->array+in->size,  *oi;
  unsigned long    *iul=in->array, *iulf=in->array+in->size, *oul;
  long              *il=in->array,  *ilf=in->array+in->size,  *ol;
  LONGLONG          *iL=in->array,  *iLf=in->array+in->size,  *oL;

  /* Check the type */
  switch(in->type)
    {
    case GAL_DATA_TYPE_FLOAT:
    case GAL_DATA_TYPE_DOUBLE:
      error(EXIT_FAILURE, 0, "the bitwise not (one's complement) "
            "operator can only work on integer types");
    }

  /* If we want inplace output, set the output pointer to the input
     pointer, for every pixel, the operation will be independent. */
  if(flags & GAL_DATA_ARITH_INPLACE)
    o = in;
  else
    o = gal_data_alloc(NULL, in->type, in->ndim, in->dsize, in->wcs,
                       0, in->minmapsize, NULL, NULL, NULL);

  /* Start setting the types. */
  switch(in->type)
    {
    case GAL_DATA_TYPE_UCHAR:
      ouc=o->array;   do *ouc++ = ~(*iuc++);  while(iuc<iucf);
    case GAL_DATA_TYPE_CHAR:
      oc=o->array;    do  *oc++ = ~(*ic++);   while(ic<icf);
    case GAL_DATA_TYPE_USHORT:
      ous=o->array;   do *ous++ = ~(*ius++);  while(ius<iusf);
    case GAL_DATA_TYPE_SHORT:
      os=o->array;    do  *os++ = ~(*is++);   while(is<isf);
    case GAL_DATA_TYPE_UINT:
      oui=o->array;   do *oui++ = ~(*iui++);  while(iui<iuif);
    case GAL_DATA_TYPE_INT:
      oi=o->array;    do  *oi++ = ~(*ii++);   while(ii<iif);
    case GAL_DATA_TYPE_ULONG:
      oul=o->array;   do *oul++ = ~(*iul++);  while(iul<iulf);
    case GAL_DATA_TYPE_LONG:
      ol=o->array;    do  *ol++ = ~(*il++);   while(il<ilf);
    case GAL_DATA_TYPE_LONGLONG:
      oL=o->array;    do  *oL++ = ~(*iL++);   while(iL<iLf);
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
  if( (flags & GAL_DATA_ARITH_FREE) && o!=in)
    gal_data_free(in, 0);

  /* Return */
  return o;
}
