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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/qsort.h>









/***********************************************************************/
/***************        Unary functions/operators         **************/
/***********************************************************************/
/* Change input data structure type. */
gal_data_t *
data_arithmetic_change_type(gal_data_t *data, int operator,
                            unsigned char flags)
{
  int type=-1;
  gal_data_t *out;

  /* Set the output type. */
  switch(operator)
    {
    case GAL_DATA_OPERATOR_TO_UCHAR:    type=GAL_DATA_TYPE_UCHAR;    break;
    case GAL_DATA_OPERATOR_TO_CHAR:     type=GAL_DATA_TYPE_UCHAR;    break;
    case GAL_DATA_OPERATOR_TO_USHORT:   type=GAL_DATA_TYPE_USHORT;   break;
    case GAL_DATA_OPERATOR_TO_SHORT:    type=GAL_DATA_TYPE_SHORT;    break;
    case GAL_DATA_OPERATOR_TO_UINT:     type=GAL_DATA_TYPE_UINT;     break;
    case GAL_DATA_OPERATOR_TO_INT:      type=GAL_DATA_TYPE_INT;      break;
    case GAL_DATA_OPERATOR_TO_ULONG:    type=GAL_DATA_TYPE_ULONG;    break;
    case GAL_DATA_OPERATOR_TO_LONG:     type=GAL_DATA_TYPE_LONG;     break;
    case GAL_DATA_OPERATOR_TO_LONGLONG: type=GAL_DATA_TYPE_LONGLONG; break;
    case GAL_DATA_OPERATOR_TO_FLOAT:    type=GAL_DATA_TYPE_FLOAT;    break;
    case GAL_DATA_OPERATOR_TO_DOUBLE:   type=GAL_DATA_TYPE_DOUBLE;   break;
    default:
      error(EXIT_FAILURE, 0, "operator value of %d not recognized in "
            "`data_arithmetic_change_type'", operator);
    }

  /* Copy to the new type. */
  out=gal_data_copy_to_new_type(data, type);

  /* Delete the input structure if the user asked for it. */
  if(flags & GAL_DATA_ARITH_FREE)
    gal_data_free(data, 0);

  /* Return */
  return out;
}





/* Return an array of value 1 for any zero valued element and zero for any
   non-zero valued element. */
#define TYPE_CASE_FOR_NOT(TYPE, IN, IN_FINISH) {                        \
    case TYPE:                                                          \
      do *o++ = !*IN; while(++IN<IN_FINISH);                            \
      break;                                                            \
  }

gal_data_t *
data_arithmetic_not(gal_data_t *data, unsigned char flags)
{
  gal_data_t *out;

  /* 'value' will only be read from one of these based on the
     datatype. Which the caller assigned. If there is any problem, it is
     their responsability, not this function's.*/
  unsigned char     *uc = data->array,   *ucf = data->array + data->size, *o;
  char               *c = data->array,    *cf = data->array + data->size;
  unsigned short    *us = data->array,   *usf = data->array + data->size;
  short              *s = data->array,    *sf = data->array + data->size;
  unsigned int      *ui = data->array,   *uif = data->array + data->size;
  int               *in = data->array,   *inf = data->array + data->size;
  unsigned long     *ul = data->array,   *ulf = data->array + data->size;
  long               *l = data->array,    *lf = data->array + data->size;
  LONGLONG           *L = data->array,    *Lf = data->array + data->size;
  float              *f = data->array,    *ff = data->array + data->size;
  double             *d = data->array,    *df = data->array + data->size;


  /* Allocate the output array. */
  out=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, data->ndim, data->dsize,
                     data->wcs, 0, data->minmapsize, data->name, data->unit,
                     data->comment);
  o=out->array;


  /* Go over the pixels and set the output values. */
  switch(data->type)
    {

    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_UCHAR,    uc,  ucf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_CHAR,     c,   cf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_LOGICAL,  c,   cf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_USHORT,   us,  usf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_SHORT,    s,   sf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_UINT,     ui,  uif)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_INT,      in,  inf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_ULONG,    ul,  ulf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_LONG,     l,   lf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_LONGLONG, L,   Lf)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_FLOAT,    f,   ff)
    TYPE_CASE_FOR_NOT(GAL_DATA_TYPE_DOUBLE,   d,   df)

    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support bit "
            "datatype, please get in touch with us to implement it.");

    default:
      error(EXIT_FAILURE, 0, "type value (%d) not recognized "
            "in `data_arithmetic_not'", data->type);
    }

  /* Delete the input structure if the user asked for it. */
  if(flags & GAL_DATA_ARITH_FREE)
    gal_data_free(data, 0);

  /* Return */
  return out;
}





/* We don't want to use the standard function for unary functions in the
   case of the absolute operator. This is because there are multiple
   versions of this function in the C library for different types, which
   can greatly improve speed. */
gal_data_t *
data_arithmetic_abs(unsigned char flags, gal_data_t *in)
{
  gal_data_t *out;

  unsigned char  *ouc,  *uc = in->array,   *ucf = in->array + in->size;
  char            *oc,   *c = in->array,    *cf = in->array + in->size;
  unsigned short *ous,  *us = in->array,   *usf = in->array + in->size;
  short           *os,   *s = in->array,    *sf = in->array + in->size;
  unsigned int   *oui,  *ui = in->array,   *uif = in->array + in->size;
  int             *oi,  *ii = in->array,   *iif = in->array + in->size;
  unsigned long  *oul,  *ul = in->array,   *ulf = in->array + in->size;
  long            *ol,   *l = in->array,    *lf = in->array + in->size;
  LONGLONG        *oL,   *L = in->array,    *Lf = in->array + in->size;
  float           *of,   *f = in->array,    *ff = in->array + in->size;
  double          *od,   *d = in->array,    *df = in->array + in->size;

  /* Set the output array. */
  if(flags & GAL_DATA_ARITH_INPLACE)
    out=in;
  else
    out = gal_data_alloc(NULL, in->type, in->ndim, in->dsize,
                         in->wcs, 0, in->minmapsize, in->name, in->unit,
                         in->comment);

  /* Put the absolute value depending on the type. Note that the unsigned
     types are already positive, so if the input is not to be freed (the
     output must be a separate array), just copy the values.*/
  switch(in->type)
    {
    case GAL_DATA_TYPE_UCHAR:
      if(out!=in) { ouc=out->array; do *ouc++ = *uc++; while(uc<ucf); }
      break;
    case GAL_DATA_TYPE_CHAR:
      oc=out->array; do *oc++ = abs(*c++); while(c<cf);
      break;
    case GAL_DATA_TYPE_USHORT:
      if(out!=in) { ous=out->array; do *ous++ = *us++; while(us<usf); }
      break;
    case GAL_DATA_TYPE_SHORT:
      os=out->array; do *os++ = abs(*s++); while(s<sf);
      break;
    case GAL_DATA_TYPE_UINT:
      if(out!=in) { oui=out->array; do *oui++ = *ui++; while(ui<uif); }
      break;
    case GAL_DATA_TYPE_INT:
      oi=out->array; do *oi++ = abs(*ii++); while(ii<iif);
      break;
    case GAL_DATA_TYPE_ULONG:
      if(out!=in) { oul=out->array; do *oul++ = *ul++; while(ul<ulf); }
      break;
    case GAL_DATA_TYPE_LONG:
      ol=out->array; do *ol++ = labs(*l++); while(l<lf);
      break;
    case GAL_DATA_TYPE_LONGLONG:
      oL=out->array; do *oL++ = llabs(*L++); while(L<Lf);
      break;
    case GAL_DATA_TYPE_FLOAT:
      of=out->array; do *of++ = fabsf(*f++); while(f<ff);
      break;
    case GAL_DATA_TYPE_DOUBLE:
      od=out->array; do *od++ = fabs(*d++); while(d<df);
      break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`data_arithmetic_abs'", in->type);
    }

  /* Clean up and return */
  if( (flags & GAL_DATA_ARITH_FREE) && out!=in)
    gal_data_free(in, 0);
  return out;
}




















/***********************************************************************/
/***************          Checking functions              **************/
/***********************************************************************/

/* Some functions are only for a floating point operand, so if the input
   isn't floating point, inform the user to change the type explicitly,
   doing it implicitly/internally puts too much responsability on the
   program. */
static void
check_float_input(gal_data_t *in, int operator, char *numstr)
{
  switch(in->type)
    {
    case GAL_DATA_TYPE_FLOAT:
    case GAL_DATA_TYPE_DOUBLE:
      break;
    default:
      error(EXIT_FAILURE, 0, "the %s operator can only accept single or "
            "double precision floating point numbers as its operand. The "
            "%s operand has type %s. You can use the `float' or `double' "
            "operators before this operator to explicity convert to the "
            "desired precision floating point type. If the operand was "
            "originally a typed number (string of characters), add an `f' "
            "after it so it is directly read into the proper precision "
            "floating point number (based on the number of non-zero "
            "decimals it has)", gal_data_operator_string(operator), numstr,
            gal_data_type_string(in->type, 1));
    }
}




















/***********************************************************************/
/***************             Unary functions              **************/
/***********************************************************************/

/* Note that for floating point types *b!=*b (by definition of NaN). And in
   such cases, even if there are blank values, the smaller and larger
   condition checked will fail, therefore the final result will be what we
   want (to ignore the blank values). */
#define UNIFUNC_MINVALUE {                                              \
    int blankeq = (*b==*b && gal_data_has_blank(in)) ? 1 : 0;           \
    if(blankeq)                                                         \
      do if(*ia!=*b) *oa = *ia < *oa ? *ia : *oa; while(++ia<iaf);      \
    else                                                                \
      do *oa = *ia < *oa ? *ia : *oa; while(++ia<iaf);                  \
  }

#define UNIFUNC_MAXVALUE {                                              \
    int blankeq = (*b==*b && gal_data_has_blank(in)) ? 1 : 0;           \
    if(blankeq)                                                         \
      do if(*ia!=*b) *oa = *ia > *oa ? *ia : *oa; while(++ia<iaf);      \
    else                                                                \
      do *oa = *ia > *oa ? *ia : *oa; while(++ia<iaf);                  \
  }





#define UNIFUNC_RUN_FUNCTION_ON_ELEMENT(IT, OP){                        \
    IT *ia=in->array, *oa=o->array, *iaf=ia + in->size;                 \
    do *oa++ = OP(*ia++); while(ia<iaf);                                \
  }

#define UNIFUNC_RUN_FUNCTION_ON_ARRAY(IT){                              \
    IT *b = gal_data_alloc_blank(in->type);                             \
    IT *ia=in->array, *oa=o->array, *iaf=ia + in->size;                 \
    switch(operator)                                                    \
      {                                                                 \
      case GAL_DATA_OPERATOR_MINVAL:                                    \
        UNIFUNC_MINVALUE;                                               \
        break;                                                          \
      case GAL_DATA_OPERATOR_MAXVAL:                                    \
        UNIFUNC_MAXVALUE;                                               \
        break;                                                          \
      default:                                                          \
        error(EXIT_FAILURE, 0, "the operator code %d is not "           \
              "recognized in UNIFUNC_RUN_FUNCTION_ON_ARRAY", operator); \
      }                                                                 \
    free(b);                                                            \
  }





#define UNIARY_FUNCTION_ON_ELEMENT(OP)                                  \
  switch(in->type)                                                      \
    {                                                                   \
    case GAL_DATA_TYPE_UCHAR:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(unsigned char, OP)                \
      break;                                                            \
    case GAL_DATA_TYPE_CHAR:                                            \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(char, OP)                         \
      break;                                                            \
    case GAL_DATA_TYPE_USHORT:                                          \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(unsigned short, OP)               \
      break;                                                            \
    case GAL_DATA_TYPE_SHORT:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(short, OP)                        \
      break;                                                            \
    case GAL_DATA_TYPE_UINT:                                            \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(unsigned int, OP)                 \
      break;                                                            \
    case GAL_DATA_TYPE_INT:                                             \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(int, OP)                          \
      break;                                                            \
    case GAL_DATA_TYPE_ULONG:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(unsigned long, OP)                \
      break;                                                            \
    case GAL_DATA_TYPE_LONG:                                            \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(long, OP)                         \
      break;                                                            \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(LONGLONG, OP)                     \
      break;                                                            \
    case GAL_DATA_TYPE_FLOAT:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(float, OP)                        \
      break;                                                            \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      UNIFUNC_RUN_FUNCTION_ON_ELEMENT(double, OP)                       \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type code %d not recognized in "          \
            "`UNIFUNC_PER_ELEMENT'", in->type);                         \
    }





#define UNIARY_FUNCTION_ON_ARRAY                                        \
  switch(in->type)                                                      \
    {                                                                   \
    case GAL_DATA_TYPE_UCHAR:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(unsigned char)                      \
      break;                                                            \
    case GAL_DATA_TYPE_CHAR:                                            \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(char)                               \
      break;                                                            \
    case GAL_DATA_TYPE_USHORT:                                          \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(unsigned short)                     \
      break;                                                            \
    case GAL_DATA_TYPE_SHORT:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(short)                              \
        break;                                                          \
    case GAL_DATA_TYPE_UINT:                                            \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(unsigned int)                       \
        break;                                                          \
    case GAL_DATA_TYPE_INT:                                             \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(int)                                \
        break;                                                          \
    case GAL_DATA_TYPE_ULONG:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(unsigned long)                      \
        break;                                                          \
    case GAL_DATA_TYPE_LONG:                                            \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(long)                               \
        break;                                                          \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(LONGLONG)                           \
      break;                                                            \
    case GAL_DATA_TYPE_FLOAT:                                           \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(float)                              \
      break;                                                            \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      UNIFUNC_RUN_FUNCTION_ON_ARRAY(double)                             \
        break;                                                          \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type code %d not recognized in "          \
            "`UNIFUNC_PER_ELEMENT'", in->type);                         \
    }





gal_data_t *
data_arithmetic_unary_function(int operator, unsigned char flags,
                               gal_data_t *in)
{
  long dsize=1;
  gal_data_t *o;

  /* If we want inplace output, set the output pointer to the input
     pointer, for every pixel, the operation will be independent. */
  switch(operator)
    {

    /* Operators with only one value as output. */
    case GAL_DATA_OPERATOR_MINVAL:
      o = gal_data_alloc(NULL, in->type, 1, &dsize, NULL, 0, -1,
                         NULL, NULL, NULL);
      gal_data_type_max(o->type, o->array);
      break;
    case GAL_DATA_OPERATOR_MAXVAL:
      o = gal_data_alloc(NULL, in->type, 1, &dsize, NULL, 0, -1,
                         NULL, NULL, NULL);
      gal_data_type_min(o->type, o->array);
      break;

    /* The other operators  */
    default:
      if(flags & GAL_DATA_ARITH_INPLACE)
        o = in;
      else
        o = gal_data_alloc(NULL, in->type, in->ndim, in->dsize, in->wcs,
                           0, in->minmapsize, NULL, NULL, NULL);
    }

  /* Start setting the operator and operands. */
  switch(operator)
    {
    case GAL_DATA_OPERATOR_SQRT:
      UNIARY_FUNCTION_ON_ELEMENT( sqrt );
      break;

    case GAL_DATA_OPERATOR_LOG:
      UNIARY_FUNCTION_ON_ELEMENT( log );
      break;

    case GAL_DATA_OPERATOR_LOG10:
      UNIARY_FUNCTION_ON_ELEMENT( log10 );
      break;

    case GAL_DATA_OPERATOR_MINVAL:
    case GAL_DATA_OPERATOR_MAXVAL:
      UNIARY_FUNCTION_ON_ARRAY;
      break;

    default:
      error(EXIT_FAILURE, 0, "operator code %d not recognized in "
            "data_arithmetic_unary_function", operator);
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




















/***********************************************************************/
/***************            Binary functions              **************/
/***********************************************************************/


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
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINFUNC_RUN_FUNCTION(float, RT, LT, OP);                          \
      break;                                                            \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINFUNC_RUN_FUNCTION(double, RT, LT, OP);                         \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for o->type in BINFUNC_F_OPERATOR_LEFT_RIGHT_SET",         \
            o->type);                                                   \
    }





#define BINFUNC_F_OPERATOR_LEFT_SET(LT, OP)                             \
  switch(r->type)                                                       \
    {                                                                   \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINFUNC_F_OPERATOR_LEFT_RIGHT_SET(float, LT, OP);                 \
      break;                                                            \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINFUNC_F_OPERATOR_LEFT_RIGHT_SET(double, LT, OP);                \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for r->type in BINFUNC_F_OPERATOR_LEFT_SET", r->type);    \
    }





#define BINFUNC_F_OPERATOR_SET(OP)                                      \
  switch(l->type)                                                       \
    {                                                                   \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINFUNC_F_OPERATOR_LEFT_SET(float, OP);                           \
      break;                                                            \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINFUNC_F_OPERATOR_LEFT_SET(double, OP);                          \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for l->type in BINFUNC_F_OPERATOR_SET", l->type);          \
    }





gal_data_t *
data_arithmetic_binary_function_flt(int operator, unsigned char flags,
                                    gal_data_t *l, gal_data_t *r)
{
  int final_otype;
  gal_data_t *o=NULL;
  size_t out_size, minmapsize;


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_DATA_ARITH_NUMOK) && (l->size==1 || r->size==1))
      && gal_data_dsize_is_different(l, r) )
    error(EXIT_FAILURE, 0, "the input datasets don't have the same "
          "dimension/size in data_arithmetic_binary_function");

  /* Check for the types of the left and right operands. */
  check_float_input(l, operator, "first");
  check_float_input(r, operator, "second");

  /* Set the output type. */
  final_otype = gal_data_out_type(l, r);

  /* Set the output sizes. */
  minmapsize = ( l->minmapsize < r->minmapsize
                 ? l->minmapsize : r->minmapsize );
  out_size = l->size > r->size ? l->size : r->size;


  /* If we want inplace output, set the output pointer to one input. Note
     that the output type can be different from both inputs.  */
  if(flags & GAL_DATA_ARITH_INPLACE)
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
    case GAL_DATA_OPERATOR_POW:  BINFUNC_F_OPERATOR_SET( pow  ); break;
    default:
      error(EXIT_FAILURE, 0, "Operator code %d not recognized in "
            "data_arithmetic_binary_function", operator);
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, `BINFUNC_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, and we only have `r' and `l'
     allocated in any case. Alternatively, when the inputs shouldn't be
     freed, the only allocated spaces are the `r' and `l' arrays if their
     types weren't compiled for binary operations, we can tell this from
     the pointers: if they are different from the original pointers, they
     were allocated. */
  if(flags & GAL_DATA_ARITH_FREE)
    {
      if     (o==l)       gal_data_free(r, 0);
      else if(o==r)       gal_data_free(l, 0);
      else              { gal_data_free(l, 0); gal_data_free(r, 0); }
    }

  /* Return */
  return o;
}




















/***********************************************************************/
/***************                  Where                   **************/
/***********************************************************************/
#define DO_WHERE_OPERATION(ITT, OT) {                                \
    ITT *it=iftrue->array;                                           \
    OT *o=out->array, *of=out->array+out->size;                      \
    if(iftrue->size==1)                                              \
      do   *o = *c++ ? *it : *o;         while(++o<of);              \
    else                                                             \
      do { *o = *c++ ? *it : *o; ++it; } while(++o<of);              \
}





#define WHERE_OUT_SET(OT)                                            \
  switch(iftrue->type)                                               \
    {                                                                \
    case GAL_DATA_TYPE_UCHAR:                                        \
      DO_WHERE_OPERATION(unsigned char, OT);                         \
      break;                                                         \
    case GAL_DATA_TYPE_CHAR:                                         \
      DO_WHERE_OPERATION(char, OT);                                  \
      break;                                                         \
    case GAL_DATA_TYPE_USHORT:                                       \
      DO_WHERE_OPERATION(unsigned short, OT);                        \
      break;                                                         \
    case GAL_DATA_TYPE_SHORT:                                        \
      DO_WHERE_OPERATION(short, OT);                                 \
      break;                                                         \
    case GAL_DATA_TYPE_UINT:                                         \
      DO_WHERE_OPERATION(unsigned int, OT);                          \
      break;                                                         \
    case GAL_DATA_TYPE_INT:                                          \
      DO_WHERE_OPERATION(int, OT);                                   \
      break;                                                         \
    case GAL_DATA_TYPE_ULONG:                                        \
      DO_WHERE_OPERATION(unsigned long, OT);                         \
      break;                                                         \
    case GAL_DATA_TYPE_LONG:                                         \
      DO_WHERE_OPERATION(long, OT);                                  \
      break;                                                         \
    case GAL_DATA_TYPE_LONGLONG:                                     \
      DO_WHERE_OPERATION(LONGLONG, OT);                              \
      break;                                                         \
    case GAL_DATA_TYPE_FLOAT:                                        \
      DO_WHERE_OPERATION(float, OT);                                 \
      break;                                                         \
    case GAL_DATA_TYPE_DOUBLE:                                       \
      DO_WHERE_OPERATION(double, OT);                                \
      break;                                                         \
    default:                                                         \
      error(EXIT_FAILURE, 0, "type code %d not recognized for the "  \
            "`iftrue' dataset of `WHERE_OUT_SET'", iftrue->type);    \
    }





void
data_arithmetic_where(unsigned char flags, gal_data_t *out,
                      gal_data_t *cond, gal_data_t *iftrue)
{
  unsigned char *c=cond->array;

  /* The condition operator has to be unsigned char. */
  if(cond->type!=GAL_DATA_TYPE_UCHAR)
    error(EXIT_FAILURE, 0, "the condition operand to "
          "`data_arithmetic_where' must be an `unsigned char' type, but "
          "the given condition operator has a `%s' type",
          gal_data_type_string(cond->type, 1));

  /* The dimension and sizes of the out and condition data sets must be the
     same. */
  if(gal_data_dsize_is_different(out, cond))
    error(EXIT_FAILURE, 0, "the output and condition data sets of the "
          "`where' operator must be the same size");

  /* Do the operation. */
  switch(out->type)
    {
    case GAL_DATA_TYPE_UCHAR:
      WHERE_OUT_SET(unsigned char);
      break;
    case GAL_DATA_TYPE_CHAR:
      WHERE_OUT_SET(char);
      break;
    case GAL_DATA_TYPE_USHORT:
      WHERE_OUT_SET(unsigned short);
      break;
    case GAL_DATA_TYPE_SHORT:
      WHERE_OUT_SET(short);
      break;
    case GAL_DATA_TYPE_UINT:
      WHERE_OUT_SET(unsigned int);
      break;
    case GAL_DATA_TYPE_INT:
      WHERE_OUT_SET(int);
      break;
    case GAL_DATA_TYPE_ULONG:
      WHERE_OUT_SET(unsigned long);
      break;
    case GAL_DATA_TYPE_LONG:
      WHERE_OUT_SET(long);
      break;
    case GAL_DATA_TYPE_LONGLONG:
      WHERE_OUT_SET(LONGLONG);
      break;
    case GAL_DATA_TYPE_FLOAT:
      WHERE_OUT_SET(float);
      break;
    case GAL_DATA_TYPE_DOUBLE:
      WHERE_OUT_SET(double);
      break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized for the `out' "
            "dataset of `data_arithmetic_where'", out->type);
    }

  /* Clean up if necessary. */
  if(flags & GAL_DATA_ARITH_FREE)
    {
      gal_data_free(cond, 0);
      gal_data_free(iftrue, 0);
    }
}




















/***********************************************************************/
/***************        Multiple operand operators        **************/
/***********************************************************************/
#define MULTIOPERAND_MIN(TYPE) {                                        \
    TYPE p, max;                                                        \
    gal_data_type_max(list->type, &max);                                \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        p=max;                                                          \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {  /* Only for integer types, does *b==*b. */                 \
            if(hasblank[i] && *b==*b)                                   \
              { if( *a[i] != *b ) p = *a[i] < p ? *a[i] : p;            \
                else              p = *a[i] < p ? *a[i] : p; }          \
            ++a[i];                                                     \
          }                                                             \
        *o++=p;                                                         \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_MAX(TYPE) {                                        \
    TYPE p, min;                                                        \
    gal_data_type_min(list->type, &min);                                \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        p=min;                                                          \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {  /* Only for integer types, does *b==*b. */                 \
            if(hasblank[i] && *b==*b)                                   \
              { if( *a[i] != *b ) p = *a[i] > p ? *a[i] : p;            \
                else              p = *a[i] > p ? *a[i] : p; }          \
            ++a[i];                                                     \
          }                                                             \
        *o++=p;                                                         \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_SUM {                                              \
    double sum;                                                         \
    int n, use;                                                         \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        use=1;                                                          \
        sum=0.0f;                                                       \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
              use = ( *b==*b                                            \
                      ? ( *a[i]!=*b    ? 1 : 0 )          /* Integer */ \
                      : ( *a[i]==*a[i] ? 1 : 0 ) );       /* Float   */ \
                                                                        \
            /* a[i] must be incremented to next pixel in any case. */   \
            if(use) { sum += *a[i]++; ++n; } else ++a[i];               \
          }                                                             \
        *o++ = n ? sum : *b;                                            \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_AVERAGE {                                          \
    double sum;                                                         \
    int n, use;                                                         \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        use=1;                                                          \
        sum=0.0f;                                                       \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
              use = ( *b==*b                                            \
                      ? ( *a[i]!=*b    ? 1 : 0 )          /* Integer */ \
                      : ( *a[i]==*a[i] ? 1 : 0 ) );       /* Float   */ \
                                                                        \
            /* a[i] must be incremented to next pixel in any case. */   \
            if(use) { sum += *a[i]++; ++n; } else ++a[i];               \
          }                                                             \
        *o++ = n ? sum/n : *b;                                          \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_MEDIAN(TYPE, QSORT_F) {                            \
    TYPE *pixs;                                                         \
    int n, use;                                                         \
                                                                        \
    errno=0;                                                            \
    pixs=malloc(dnum*sizeof *pixs);                                     \
    if(pixs==NULL)                                                      \
      error(EXIT_FAILURE, 0, "%zu bytes in MULTIOPERAND_MEDIAN",        \
            dnum*sizeof *pixs);                                         \
                                                                        \
    do    /* Loop over each pixel */                                    \
      {                                                                 \
        n=0;                                                            \
        use=1;                                                          \
        for(i=0;i<dnum;++i)  /* Loop over each array. */                \
          {                                                             \
            /* Only integers and non-NaN floats: v==v is 1. */          \
            if(hasblank[i])                                             \
              use = ( *b==*b                                            \
                      ? ( *a[i]!=*b    ? 1 : 0 )          /* Integer */ \
                      : ( *a[i]==*a[i] ? 1 : 0 ) );       /* Float   */ \
                                                                        \
            /* a[i] must be incremented to next pixel in any case. */   \
            if(use) pixs[n++]=*a[i]++; else ++a[i];                     \
          }                                                             \
                                                                        \
        /* Sort all the values for this pixel and return the median. */ \
        if(n)                                                           \
          {                                                             \
            qsort(pixs, n, sizeof *pixs, QSORT_F);                      \
            *o++ = n%2 ? pixs[n/2] : (pixs[n/2] + pixs[n/2-1])/2 ;      \
          }                                                             \
        else                                                            \
          *o++=*b;                                                      \
      }                                                                 \
    while(o<of);                                                        \
  }





#define MULTIOPERAND_TYPE_SET(TYPE, QSORT_F) {                          \
    TYPE *o=out->array, *of=out->array+out->size;                       \
    TYPE **a, *b=gal_data_alloc_blank(list->type);                      \
                                                                        \
    /* Allocate space to keep the pointers to the arrays of each. */    \
    /* Input data structure. The operators will increment these */      \
    /* pointers while parsing them. */                                  \
    errno=0;                                                            \
    a=malloc(dnum*sizeof *a);                                           \
    if(a==NULL)                                                         \
      error(EXIT_FAILURE, 0, "%zu bytes for `arrays' in "               \
            "MULTIOPERAND_TYPE_SET", dnum*sizeof *a);                   \
                                                                        \
    /* Fill in the array pointers. */                                   \
    for(tmp=list;tmp!=NULL;tmp=tmp->next)                               \
      a[i++]=tmp->array;                                                \
                                                                        \
    /* Do the operation. */                                             \
    switch(operator)                                                    \
      {                                                                 \
      case GAL_DATA_OPERATOR_MIN:                                       \
        MULTIOPERAND_MIN(TYPE);                                         \
        break;                                                          \
                                                                        \
      case GAL_DATA_OPERATOR_MAX:                                       \
        MULTIOPERAND_MAX(TYPE);                                         \
        break;                                                          \
                                                                        \
      case GAL_DATA_OPERATOR_SUM:                                       \
        MULTIOPERAND_SUM;                                               \
        break;                                                          \
                                                                        \
      case GAL_DATA_OPERATOR_AVERAGE:                                   \
        MULTIOPERAND_AVERAGE;                                           \
        break;                                                          \
                                                                        \
      case GAL_DATA_OPERATOR_MEDIAN:                                    \
        MULTIOPERAND_MEDIAN(TYPE, QSORT_F);                             \
        break;                                                          \
                                                                        \
      default:                                                          \
        error(EXIT_FAILURE, 0, "the operator code %d not recognized "   \
              "in MULTIOPERAND_TYPE_SET", operator);                    \
      }                                                                 \
                                                                        \
    /* Clean up. */                                                     \
    free(b);                                                            \
    free(a);                                                            \
  }





/* The single operator in this function is assumed to be a linked list. The
   number of operators is determined from the fact that the last node in
   the linked list must have a NULL pointer as its `next' element.*/
gal_data_t *
data_arithmetic_multioperand(int operator, unsigned char flags,
                             gal_data_t *list)
{
  int *hasblank;
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
        error(EXIT_FAILURE, 0, "the types of all operands to the %s "
              "operator must be same",
              gal_data_operator_string(operator));

      /* Check the sizes. */
      if( gal_data_dsize_is_different(list, tmp) )
        error(EXIT_FAILURE, 0, "the sizes of all operands to the %s "
              "operator must be same",
              gal_data_operator_string(operator));
    }


  /* Set the output data structure. */
  if(flags & GAL_DATA_ARITH_INPLACE)
    out = list;                 /* The top element in the list. */
  else
    out = gal_data_alloc(NULL, list->type, list->ndim, list->dsize,
                         list->wcs, 0, list->minmapsize, NULL, NULL, NULL);


  /* hasblank is used to see if a blank value should be checked or not. */
  errno=0;
  hasblank=malloc(dnum*sizeof *hasblank);
  if(hasblank==NULL)
    error(EXIT_FAILURE, 0, "%zu bytes for hasblank in "
          "`data_arithmetic_multioperand", dnum*sizeof *hasblank);


  /* Fill in the hasblank array. */
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    hasblank[i++]=gal_data_has_blank(tmp);


  /* Start the operation. */
  switch(list->type)
    {
    case GAL_DATA_TYPE_UCHAR:
      MULTIOPERAND_TYPE_SET(unsigned char,  gal_qsort_uchar_increasing);
    case GAL_DATA_TYPE_CHAR:
      MULTIOPERAND_TYPE_SET(char,           gal_qsort_char_increasing);
    case GAL_DATA_TYPE_USHORT:
      MULTIOPERAND_TYPE_SET(unsigned short, gal_qsort_ushort_increasing);
    case GAL_DATA_TYPE_SHORT:
      MULTIOPERAND_TYPE_SET(short,          gal_qsort_short_increasing);
    case GAL_DATA_TYPE_UINT:
      MULTIOPERAND_TYPE_SET(unsigned int,   gal_qsort_uint_increasing);
    case GAL_DATA_TYPE_INT:
      MULTIOPERAND_TYPE_SET(int,            gal_qsort_int_increasing);
    case GAL_DATA_TYPE_ULONG:
      MULTIOPERAND_TYPE_SET(unsigned long,  gal_qsort_ulong_increasing);
    case GAL_DATA_TYPE_LONG:
      MULTIOPERAND_TYPE_SET(long,           gal_qsort_long_increasing);
    case GAL_DATA_TYPE_LONGLONG:
      MULTIOPERAND_TYPE_SET(LONGLONG,       gal_qsort_longlong_increasing);
    case GAL_DATA_TYPE_FLOAT:
      MULTIOPERAND_TYPE_SET(float,          gal_qsort_float_increasing);
    case GAL_DATA_TYPE_DOUBLE:
      MULTIOPERAND_TYPE_SET(double,         gal_qsort_double_increasing);
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`data_arithmetic_multioperand'", list->type);
    }


  /* Clean up and return. Note that the operation might have been done in
     place. In that case, the top most list element was used. So we need to
     check before freeing each data structure. */
  if(flags & GAL_DATA_ARITH_FREE)
    {
      tmp=list;
      while(tmp!=NULL)
        {
          ttmp=tmp->next;
          if(tmp!=out) gal_data_free(tmp, 0);
          tmp=ttmp;
        }
    }
  free(hasblank);
  return out;
}
