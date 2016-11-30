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
#ifndef __GAL_ARITHMETIC_H__
#define __GAL_ARITHMETIC_H__




#define BINARY_OPERATOR_FOR_TYPE(LT, RT, OT, OP){                       \
    LT *la=l->array;                                                    \
    RT *ra=r->array;                                                    \
    OT *oa=o->array, *of=oa + o->size;                                  \
    if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);       \
    else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);       \
    else                 do *oa = *la++ OP *ra;   while(++oa<of);       \
  }




#if GAL_CONFIG_ARITH_CHAR == 1
#define BINARY_LEFT_RIGHT_DONE_CHAR(LT, RT, OP)                         \
  case GAL_DATA_TYPE_CHAR:                                              \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, char, OP);                         \
    break;
#define BINARY_LEFT_DONE_CHAR(LT, OP)                                   \
    case GAL_DATA_TYPE_CHAR:                                            \
      BINARY_LEFT_RIGHT_DONE(LT, char, OP);                             \
      break;
#define BINARY_MULTISWITCH_CHAR(OP)                                     \
    case GAL_DATA_TYPE_CHAR:                                            \
      BINARY_LEFT_DONE(char, OP);                                       \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_CHAR(LT, RT, OP)
#define BINARY_LEFT_DONE_CHAR(LT, OP)
#define BINARY_MULTISWITCH_CHAR(OP)
#endif



#if GAL_CONFIG_ARITH_USHORT == 1
#define BINARY_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)                         \
  case GAL_DATA_TYPE_USHORT:                                              \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, unsigned short, OP);                 \
    break;
#define BINARY_LEFT_DONE_USHORT(LT, OP)                                   \
    case GAL_DATA_TYPE_USHORT:                                            \
      BINARY_LEFT_RIGHT_DONE(LT, unsigned short, OP);                     \
      break;
#define BINARY_MULTISWITCH_USHORT(OP)                                     \
    case GAL_DATA_TYPE_USHORT:                                            \
      BINARY_LEFT_DONE(unsigned short, OP);                               \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)
#define BINARY_LEFT_DONE_USHORT(LT, OP)
#define BINARY_MULTISWITCH_USHORT(OP)
#endif



#if GAL_CONFIG_ARITH_SHORT == 1
#define BINARY_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)                         \
  case GAL_DATA_TYPE_SHORT:                                              \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, short, OP);                         \
    break;
#define BINARY_LEFT_DONE_SHORT(LT, OP)                                   \
    case GAL_DATA_TYPE_SHORT:                                            \
      BINARY_LEFT_RIGHT_DONE(LT, short, OP);                             \
      break;
#define BINARY_MULTISWITCH_SHORT(OP)                                     \
    case GAL_DATA_TYPE_SHORT:                                            \
      BINARY_LEFT_DONE(short, OP);                                       \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)
#define BINARY_LEFT_DONE_SHORT(LT, OP)
#define BINARY_MULTISWITCH_SHORT(OP)
#endif



#if GAL_CONFIG_ARITH_UINT == 1
#define BINARY_LEFT_RIGHT_DONE_UINT(LT, RT, OP)                         \
  case GAL_DATA_TYPE_UINT:                                              \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, unsigned int, OP);                 \
    break;
#define BINARY_LEFT_DONE_UINT(LT, OP)                                   \
    case GAL_DATA_TYPE_UINT:                                            \
      BINARY_LEFT_RIGHT_DONE(LT, unsigned int, OP);                     \
      break;
#define BINARY_MULTISWITCH_UINT(OP)                                     \
    case GAL_DATA_TYPE_UINT:                                            \
      BINARY_LEFT_DONE(unsigned int, OP);                               \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_UINT(LT, RT, OP)
#define BINARY_LEFT_DONE_UINT(LT, OP)
#define BINARY_MULTISWITCH_UINT(OP)
#endif



#if GAL_CONFIG_ARITH_INT == 1
#define BINARY_LEFT_RIGHT_DONE_INT(LT, RT, OP)                         \
  case GAL_DATA_TYPE_INT:                                              \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, int, OP);                         \
    break;
#define BINARY_LEFT_DONE_INT(LT, OP)                                   \
    case GAL_DATA_TYPE_INT:                                            \
      BINARY_LEFT_RIGHT_DONE(LT, int, OP);                             \
      break;
#define BINARY_MULTISWITCH_INT(OP)                                     \
    case GAL_DATA_TYPE_INT:                                            \
      BINARY_LEFT_DONE(int, OP);                                       \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_INT(LT, RT, OP)
#define BINARY_LEFT_DONE_INT(LT, OP)
#define BINARY_MULTISWITCH_INT(OP)
#endif



#if GAL_CONFIG_ARITH_ULONG == 1
#define BINARY_LEFT_RIGHT_DONE_ULONG(LT, RT, OP)                        \
    case GAL_DATA_TYPE_ULONG:                                           \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, unsigned long, OP);              \
      break;
#define BINARY_LEFT_DONE_ULONG(LT, OP)                                  \
    case GAL_DATA_TYPE_ULONG:                                           \
      BINARY_LEFT_RIGHT_DONE(LT, unsigned long, OP);                    \
      break;
#define BINARY_MULTISWITCH_ULONG(OP)                                    \
    case GAL_DATA_TYPE_ULONG:                                           \
      BINARY_LEFT_DONE(unsigned long, OP);                              \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_ULONG(LT, RT, OP)
#define BINARY_LEFT_DONE_ULONG(LT, OP)
#define BINARY_MULTISWITCH_ULONG(OP)
#endif



#if GAL_CONFIG_ARITH_LONGLONG == 1
#define BINARY_LEFT_RIGHT_DONE_LONGLONG(LT, RT, OP)                     \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, LONGLONG, OP);                   \
      break;
#define BINARY_LEFT_DONE_LONGLONG(LT, OP)                               \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      BINARY_LEFT_RIGHT_DONE(LT, long long, OP);                        \
      break;
#define BINARY_MULTISWITCH_LONGLONG(OP)                                 \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      BINARY_LEFT_DONE(LONGLONG, OP);                                   \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_LONGLONG(LT, RT, OP)
#define BINARY_LEFT_DONE_LONGLONG(LT, OP)
#define BINARY_MULTISWITCH_LONGLONG(OP)
#endif





#define BINARY_LEFT_RIGHT_DONE(LT, RT, OP)                              \
  switch(o->type)                                                       \
    {                                                                   \
                                                                        \
    case GAL_DATA_TYPE_UCHAR:                                           \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, unsigned char, OP);              \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_LONG:                                            \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, long, OP);                       \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, float, OP);                      \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, double, OP);                     \
      break;                                                            \
                                                                        \
    BINARY_LEFT_RIGHT_DONE_CHAR(LT, RT, OP)                             \
    BINARY_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)                            \
    BINARY_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)                           \
    BINARY_LEFT_RIGHT_DONE_INT(LT, RT, OP)                              \
    BINARY_LEFT_RIGHT_DONE_UINT(LT, RT, OP)                             \
    BINARY_LEFT_RIGHT_DONE_ULONG(LT, RT, OP)                            \
    BINARY_LEFT_RIGHT_DONE_LONGLONG(LT, RT, OP)                         \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for o->type in BINARY_LEFT_RIGHT_DONE", o->type);          \
    }





#define BINARY_LEFT_DONE(LT, OP)                                        \
  switch(r->type)                                                       \
    {                                                                   \
    case GAL_DATA_TYPE_UCHAR:                                           \
      BINARY_LEFT_RIGHT_DONE(LT, unsigned char, OP);                    \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_LONG:                                            \
      BINARY_LEFT_RIGHT_DONE(LT, long, OP);                             \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINARY_LEFT_RIGHT_DONE(LT, float, OP);                            \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINARY_LEFT_RIGHT_DONE(LT, double, OP);                           \
      break;                                                            \
                                                                        \
    BINARY_LEFT_DONE_CHAR(LT, OP)                                       \
    BINARY_LEFT_DONE_USHORT(LT, OP)                                     \
    BINARY_LEFT_DONE_SHORT(LT, OP)                                      \
    BINARY_LEFT_DONE_UINT(LT, OP)                                       \
    BINARY_LEFT_DONE_INT(LT, OP)                                        \
    BINARY_LEFT_DONE_ULONG(LT, OP)                                      \
    BINARY_LEFT_DONE_LONGLONG(LT, OP)                                   \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for r->type in BINARY_LEFT_DONE", r->type);                \
    }





/* Prepare the inputs and output for binary operations and do the job.*/
#define BINARY_INTERNAL(OP, OUT_TYPE) {                                 \
                                                                        \
  /* Read the variable arguments. */                                    \
  gal_data_t *l, *r;                                                    \
  l = va_arg(va, gal_data_t *);                                         \
  r = va_arg(va, gal_data_t *);                                         \
                                                                        \
                                                                        \
  /* Simple sanity check on the input sizes */                          \
  if( !( (flags & GAL_DATA_ARITH_NUMOK) && (l->size==1 || r->size==1))  \
      && gal_data_dsize_is_different(l, r) )                            \
    error(EXIT_FAILURE, 0, "The datasets don't have the same "          \
          "dimension/size");                                            \
                                                                        \
                                                                        \
  /* Set the output type and size. */                                   \
  out_type = OUT_TYPE ? OUT_TYPE : gal_data_out_type(l, r);             \
  out_size = l->size > r->size ? l->size : r->size;                     \
                                                                        \
                                                                        \
  /* If we want inplace output, set the output pointer to one input. */ \
  /* Note that the output type can be different from both inputs.    */ \
  if(flags & GAL_DATA_ARITH_INPLACE)                                    \
    {                                                                   \
      if(l->type==out_type && out_size==l->size)        o = l;          \
      else if(r->type==out_type && out_size==r->size)   o = r;          \
    }                                                                   \
                                                                        \
                                                                        \
  /* If the output pointer was not set for any reason, allocate it. */  \
  if(o==NULL)                                                           \
    o = gal_data_alloc(NULL, out_type,                                  \
                       l->size>1 ? l->ndim  : r->ndim,                  \
                       l->size>1 ? l->dsize : r->dsize,                 \
                       0, l->mmapped || r->mmapped);                    \
                                                                        \
                                                                        \
  /* Do the operations based on the different types. */                 \
  switch(l->type)                                                       \
    {                                                                   \
    case GAL_DATA_TYPE_UCHAR:                                           \
      BINARY_LEFT_DONE(unsigned char, OP);                              \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_LONG:                                            \
      BINARY_LEFT_DONE(long, OP);                                       \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINARY_LEFT_DONE(float, OP);                                      \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINARY_LEFT_DONE(double, OP);                                     \
      break;                                                            \
                                                                        \
    BINARY_MULTISWITCH_CHAR(OP)                                         \
    BINARY_MULTISWITCH_USHORT(OP)                                       \
    BINARY_MULTISWITCH_SHORT(OP)                                        \
    BINARY_MULTISWITCH_UINT(OP)                                         \
    BINARY_MULTISWITCH_INT(OP)                                          \
    BINARY_MULTISWITCH_ULONG(OP)                                        \
    BINARY_MULTISWITCH_LONGLONG(OP)                                     \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for l->type in BINARY_MULTISWITCH", l->type);              \
    }                                                                   \
                                                                        \
  /* Clean up. */                                                       \
  if(flags & GAL_DATA_ARITH_FREE)                                       \
    {                                                                   \
      if(o==l) gal_data_free(r);                                        \
      else if(o==r) gal_data_free(l);                                   \
      else {gal_data_free(l); gal_data_free(r);}                        \
    }                                                                   \
}





#endif
