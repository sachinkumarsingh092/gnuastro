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





#define BINARY_DEFINITIONS                                                \
  unsigned char  *luc = l->array, *ruc = r->array, *ouc = o->array, *ucf; \
  char           *lc  = l->array, *rc  = r->array, *oc  = o->array,  *cf; \
  unsigned short *lus = l->array, *rus = r->array, *ous = o->array, *usf; \
  short          *ls  = l->array, *rs  = r->array, *os  = o->array,  *sf; \
  unsigned int   *lui = l->array, *rui = r->array, *oui = o->array, *uif; \
  int            *li  = l->array, *ri  = r->array, *oi  = o->array, *iif; \
  unsigned long  *lul = l->array, *rul = r->array, *oul = o->array, *ulf; \
  long           *ll  = l->array, *rl  = r->array, *ol  = o->array,  *lf; \
  LONGLONG       *lL  = l->array, *rL  = r->array, *oL  = o->array,  *Lf; \
  float          *lff = l->array, *rf  = r->array, *of  = o->array,  *ff; \
  double         *ld  = l->array, *rd  = r->array, *od  = o->array,  *df;




#define OP_FUNC(OP,L, R)
#define BINARY_OPERATOR_FOR_TYPE(Tf, oT, lT, rT, OP)                    \
  Tf = oT + o->size;                                                    \
  if(l->size==r->size) do *oT = *lT++ OP *rT++; while(++oT<Tf);         \
  else if(l->size==1)  do *oT = *lT   OP *rT++; while(++oT<Tf);         \
  else                 do *oT = *lT++ OP *rT;   while(++oT<Tf);         \
  break;





/* Preparations for binary operators.

   Note that after `gal_data_to_same_type', we are only concerned with the
   `l' and `s' pointers. */
#define BINARY_INTERNAL(OP, OUT_TYPE)                                   \
                                                                        \
  gal_data_t *l, *r, *left, *right;                                     \
                                                                        \
                                                                        \
  /* Read the variable arguments. */                                    \
  left = va_arg(va, gal_data_t *);                                      \
  right = va_arg(va, gal_data_t *);                                     \
                                                                        \
                                                                        \
  /* Simple sanity check, then choose the common type. */               \
  if( (flags & GAL_DATA_ARITH_NUMOK) && (left->size==1 || right->size==1 ) ) \
    type=0;      /* Everything is ok, just a place-holder. */           \
  else if (gal_data_dsize_is_different(left, right))                    \
    error(EXIT_FAILURE, 0, "The datasets don't have the same "          \
          "dimension/size");                                            \
                                                                        \
                                                                        \
  /* Set the two datasets to the same type if they were different. */   \
  type=gal_data_out_type(left, right);                                  \
  gal_data_to_same_type(left, right, &l, &r, type,                      \
                        flags & GAL_DATA_ARITH_FREE );                  \
                                                                        \
                                                                        \
  /* Output can point to any one of the two arrays. */                  \
  if(flags & GAL_DATA_ARITH_INPLACE)                                    \
    o = l->size>1 ? l : r;                                              \
  else                                                                  \
    o = gal_data_alloc(NULL, type,                                      \
                       l->size>1 ? l->ndim  : r->ndim,                  \
                       l->size>1 ? l->dsize : r->dsize,                 \
                       0, l->mmapped || r->mmapped);                    \
                                                                        \
                                                                        \
  /* In block to allow definitions. */                                  \
  {                                                                     \
    BINARY_DEFINITIONS;                                                 \
    switch(l->type)                                                     \
      {                                                                 \
      case GAL_DATA_TYPE_UCHAR:                                         \
        BINARY_OPERATOR_FOR_TYPE(ucf, ouc, luc, ruc, OP);               \
                                                                        \
      case GAL_DATA_TYPE_CHAR:                                          \
        BINARY_OPERATOR_FOR_TYPE(cf,  oc,  lc,  rc,  OP);               \
                                                                        \
      case GAL_DATA_TYPE_USHORT:                                        \
        BINARY_OPERATOR_FOR_TYPE(usf, ous, lus, rus, OP);               \
                                                                        \
      case GAL_DATA_TYPE_SHORT:                                         \
        BINARY_OPERATOR_FOR_TYPE(sf,  os,  ls,  rs,  OP);               \
                                                                        \
      case GAL_DATA_TYPE_UINT:                                          \
        BINARY_OPERATOR_FOR_TYPE(uif, oui, lui, rui, OP);               \
                                                                        \
      case GAL_DATA_TYPE_INT:                                           \
        BINARY_OPERATOR_FOR_TYPE(iif, oi,  li,  ri,  OP);               \
                                                                        \
      case GAL_DATA_TYPE_ULONG:                                         \
        BINARY_OPERATOR_FOR_TYPE(ulf, oul, lul, rul, OP);               \
                                                                        \
      case GAL_DATA_TYPE_LONG:                                          \
        BINARY_OPERATOR_FOR_TYPE(lf,  ol,  ll,  rl,  OP);               \
                                                                        \
      case GAL_DATA_TYPE_LONGLONG:                                      \
        BINARY_OPERATOR_FOR_TYPE(Lf,  oL,  lL,  rL,  OP);               \
                                                                        \
      case GAL_DATA_TYPE_FLOAT:                                         \
        BINARY_OPERATOR_FOR_TYPE(ff, of,  lff,  rf,  OP);               \
                                                                        \
      case GAL_DATA_TYPE_DOUBLE:                                        \
        BINARY_OPERATOR_FOR_TYPE(df,  od,  ld,  rd,  OP);               \
                                                                        \
      default:                                                          \
        error(EXIT_FAILURE, 0, "type %d not recognized in "             \
              "data-arithmetic", type);                                 \
      }                                                                 \
  }                                                                     \
                                                                        \
    /* Clean up. */                                                     \
    if(flags & GAL_DATA_ARITH_FREE)                                     \
    {                                                                   \
      if(o==l) gal_data_free(r);                                        \
      else if(o==r) gal_data_free(l);                                   \
      else {gal_data_free(l); gal_data_free(r);}                        \
    }





#endif
