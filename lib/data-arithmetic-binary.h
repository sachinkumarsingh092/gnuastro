/*********************************************************************
Arithmetic operations on data structures.
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






/************************************************************************/
/*************      Possibly set binary types to convert    *************/
/************************************************************************/
#if GAL_CONFIG_BIN_OP_UCHAR == 1
#define BINARY_LEFT_RIGHT_DONE_UCHAR(LT, RT, OP)                        \
    case GAL_DATA_TYPE_UCHAR:                                           \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, unsigned char, OP);              \
      break;
#define BINARY_LEFT_DONE_UCHAR(LT, OP)                                  \
    case GAL_DATA_TYPE_UCHAR:                                           \
      BINARY_LEFT_RIGHT_DONE(LT, unsigned char, OP);                    \
      break;
#define BINARY_MULTISWITCH_UCHAR(OP)                                    \
    case GAL_DATA_TYPE_UCHAR:                                           \
      BINARY_LEFT_DONE(unsigned char, OP);                              \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_UCHAR(LT, RT, OP)
#define BINARY_LEFT_DONE_UCHAR(LT, OP)
#define BINARY_MULTISWITCH_UCHAR(OP)
#endif





#if GAL_CONFIG_BIN_OP_CHAR == 1
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





#if GAL_CONFIG_BIN_OP_USHORT == 1
#define BINARY_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)                       \
  case GAL_DATA_TYPE_USHORT:                                            \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, unsigned short, OP);               \
    break;
#define BINARY_LEFT_DONE_USHORT(LT, OP)                                 \
    case GAL_DATA_TYPE_USHORT:                                          \
      BINARY_LEFT_RIGHT_DONE(LT, unsigned short, OP);                   \
      break;
#define BINARY_MULTISWITCH_USHORT(OP)                                   \
    case GAL_DATA_TYPE_USHORT:                                          \
      BINARY_LEFT_DONE(unsigned short, OP);                             \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)
#define BINARY_LEFT_DONE_USHORT(LT, OP)
#define BINARY_MULTISWITCH_USHORT(OP)
#endif





#if GAL_CONFIG_BIN_OP_SHORT == 1
#define BINARY_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)                        \
  case GAL_DATA_TYPE_SHORT:                                             \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, short, OP);                        \
    break;
#define BINARY_LEFT_DONE_SHORT(LT, OP)                                  \
    case GAL_DATA_TYPE_SHORT:                                           \
      BINARY_LEFT_RIGHT_DONE(LT, short, OP);                            \
      break;
#define BINARY_MULTISWITCH_SHORT(OP)                                    \
    case GAL_DATA_TYPE_SHORT:                                           \
      BINARY_LEFT_DONE(short, OP);                                      \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)
#define BINARY_LEFT_DONE_SHORT(LT, OP)
#define BINARY_MULTISWITCH_SHORT(OP)
#endif





#if GAL_CONFIG_BIN_OP_UINT == 1
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





#if GAL_CONFIG_BIN_OP_INT == 1
#define BINARY_LEFT_RIGHT_DONE_INT(LT, RT, OP)                          \
  case GAL_DATA_TYPE_INT:                                               \
    BINARY_OPERATOR_FOR_TYPE(LT, RT, int, OP);                          \
    break;
#define BINARY_LEFT_DONE_INT(LT, OP)                                    \
    case GAL_DATA_TYPE_INT:                                             \
      BINARY_LEFT_RIGHT_DONE(LT, int, OP);                              \
      break;
#define BINARY_MULTISWITCH_INT(OP)                                      \
    case GAL_DATA_TYPE_INT:                                             \
      BINARY_LEFT_DONE(int, OP);                                        \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_INT(LT, RT, OP)
#define BINARY_LEFT_DONE_INT(LT, OP)
#define BINARY_MULTISWITCH_INT(OP)
#endif





#if GAL_CONFIG_BIN_OP_ULONG == 1
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





#if GAL_CONFIG_BIN_OP_LONG == 1
#define BINARY_LEFT_RIGHT_DONE_LONG(LT, RT, OP)                         \
    case GAL_DATA_TYPE_LONG:                                            \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, long, OP);                       \
      break;
#define BINARY_LEFT_DONE_LONG(LT, OP)                                   \
    case GAL_DATA_TYPE_LONG:                                            \
      BINARY_LEFT_RIGHT_DONE(LT, long, OP);                             \
      break;
#define BINARY_MULTISWITCH_LONG(OP)                                     \
    case GAL_DATA_TYPE_LONG:                                            \
      BINARY_LEFT_DONE(long, OP);                                       \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_LONG(LT, RT, OP)
#define BINARY_LEFT_DONE_LONG(LT, OP)
#define BINARY_MULTISWITCH_LONG(OP)
#endif





#if GAL_CONFIG_BIN_OP_LONGLONG == 1
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





#if GAL_CONFIG_BIN_OP_FLOAT == 1
#define BINARY_LEFT_RIGHT_DONE_FLOAT(LT, RT, OP)                        \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, float, OP);                      \
      break;
#define BINARY_LEFT_DONE_FLOAT(LT, OP)                                  \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINARY_LEFT_RIGHT_DONE(LT, float, OP);                            \
      break;
#define BINARY_MULTISWITCH_FLOAT(OP)                                    \
    case GAL_DATA_TYPE_FLOAT:                                           \
      BINARY_LEFT_DONE(float, OP);                                      \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_FLOAT(LT, RT, OP)
#define BINARY_LEFT_DONE_FLOAT(LT, OP)
#define BINARY_MULTISWITCH_FLOAT(OP)
#endif





#if GAL_CONFIG_BIN_OP_DOUBLE == 1
#define BINARY_LEFT_RIGHT_DONE_DOUBLE(LT, RT, OP)                       \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINARY_OPERATOR_FOR_TYPE(LT, RT, double, OP);                     \
      break;
#define BINARY_LEFT_DONE_DOUBLE(LT, OP)                                 \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINARY_LEFT_RIGHT_DONE(LT, double, OP);                           \
      break;
#define BINARY_MULTISWITCH_DOUBLE(OP)                                   \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      BINARY_LEFT_DONE(double, OP);                                     \
      break;
#else
#define BINARY_LEFT_RIGHT_DONE_DOUBLE(LT, RT, OP)
#define BINARY_LEFT_DONE_DOUBLE(LT, OP)
#define BINARY_MULTISWITCH_DOUBLE(OP)
#endif
























/************************************************************************/
/*************       Macros for specifying the type     *****************/
/************************************************************************/

#define BINARY_TYPE_FOR_CONVERT_TO_COMPILED_TYPE(intype)                \
  ntype=0;                                                              \
  switch(intype)                                                        \
    {                                                                   \
    case GAL_DATA_TYPE_UCHAR:                                           \
      if(GAL_CONFIG_BIN_OP_UCHAR) ntype=GAL_DATA_TYPE_UCHAR;            \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_USHORT)   ntype=GAL_DATA_TYPE_USHORT;  \
          else if(GAL_CONFIG_BIN_OP_SHORT)    ntype=GAL_DATA_TYPE_SHORT;   \
          else if(GAL_CONFIG_BIN_OP_UINT)     ntype=GAL_DATA_TYPE_UINT;    \
          else if(GAL_CONFIG_BIN_OP_INT)      ntype=GAL_DATA_TYPE_INT;     \
          else if(GAL_CONFIG_BIN_OP_ULONG)    ntype=GAL_DATA_TYPE_ULONG;   \
          else if(GAL_CONFIG_BIN_OP_LONG)     ntype=GAL_DATA_TYPE_LONG;    \
          else if(GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_CHAR:                                            \
      if(GAL_CONFIG_BIN_OP_CHAR) ntype=GAL_DATA_TYPE_CHAR;              \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_SHORT)    ntype=GAL_DATA_TYPE_SHORT;   \
          else if(GAL_CONFIG_BIN_OP_INT)      ntype=GAL_DATA_TYPE_INT;     \
          else if(GAL_CONFIG_BIN_OP_LONG)     ntype=GAL_DATA_TYPE_LONG;    \
          else if(GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_USHORT:                                          \
      if(GAL_CONFIG_BIN_OP_USHORT) ntype=GAL_DATA_TYPE_USHORT;          \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_UINT)     ntype=GAL_DATA_TYPE_UINT;    \
          else if(GAL_CONFIG_BIN_OP_INT)      ntype=GAL_DATA_TYPE_INT;     \
          else if(GAL_CONFIG_BIN_OP_ULONG)    ntype=GAL_DATA_TYPE_ULONG;   \
          else if(GAL_CONFIG_BIN_OP_LONG)     ntype=GAL_DATA_TYPE_LONG;    \
          else if(GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_SHORT:                                           \
      if(GAL_CONFIG_BIN_OP_SHORT) ntype=GAL_DATA_TYPE_SHORT;            \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_INT)      ntype=GAL_DATA_TYPE_INT;     \
          else if(GAL_CONFIG_BIN_OP_LONG)     ntype=GAL_DATA_TYPE_LONG;    \
          else if(GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_UINT:                                            \
      if(GAL_CONFIG_BIN_OP_UINT) ntype=GAL_DATA_TYPE_UINT;              \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_ULONG)    ntype=GAL_DATA_TYPE_ULONG;   \
          else if(GAL_CONFIG_BIN_OP_LONG)     ntype=GAL_DATA_TYPE_LONG;    \
          else if(GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_INT:                                             \
      if(GAL_CONFIG_BIN_OP_INT) ntype=GAL_DATA_TYPE_INT;                \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_LONG)     ntype=GAL_DATA_TYPE_LONG;    \
          else if(GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_ULONG:                                           \
      if(GAL_CONFIG_BIN_OP_ULONG) ntype=GAL_DATA_TYPE_ULONG;            \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_LONG:                                            \
      if(GAL_CONFIG_BIN_OP_LONG) ntype=GAL_DATA_TYPE_LONG;              \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;\
          else if(GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;   \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE;  \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      if(GAL_CONFIG_BIN_OP_LONGLONG) ntype=GAL_DATA_TYPE_LONGLONG;      \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_FLOAT)    ntype=GAL_DATA_TYPE_FLOAT;  \
          else if(GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE; \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_FLOAT:                                           \
      if(GAL_CONFIG_BIN_OP_FLOAT) ntype=GAL_DATA_TYPE_FLOAT;            \
      else                                                              \
        {                                                               \
          if     (GAL_CONFIG_BIN_OP_DOUBLE)   ntype=GAL_DATA_TYPE_DOUBLE; \
        }                                                               \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_DOUBLE:                                          \
      if(GAL_CONFIG_BIN_OP_DOUBLE) ntype=GAL_DATA_TYPE_DOUBLE;          \
      break;                                                            \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "BINARY_CONVERT_TO_COMPILED_TYPE", intype);                 \
    }





/* Note that for signed types, we won't be considering the unsigned types
   of the larger types. */
#define BINARY_CONVERT_TO_COMPILED_TYPE(in, out)                        \
                                                                        \
  /* Initialize the values that will be set. */                         \
  out=NULL;                                                             \
                                                                        \
  /* Set the values. */                                                 \
  BINARY_TYPE_FOR_CONVERT_TO_COMPILED_TYPE(in->type);                   \
                                                                        \
  /* If type is not compiled, then convert the dataset to the */        \
  /* first compiled larger type. */                                     \
  if(in->type==ntype)                                                   \
    out=in;                                                             \
  else                                                                  \
    {                                                                   \
      if(ntype)                                                         \
        {                                                               \
          out=gal_data_copy_to_new_type(in, ntype);                     \
          if(flags & GAL_DATA_ARITH_FREE)                               \
            { gal_data_free(in); in=NULL; }                             \
        }                                                               \
      else                                                              \
        {                                                               \
          char *typestring=gal_data_type_string(in->type);              \
          error(EXIT_FAILURE, 0, "The given %s type data given to "     \
                "binary operators is not compiled for native operation "\
                "and no larger types are compiled either.\n\nThe "      \
                "largest type (which can act as a fallback for any "    \
                "input type is double, so configure Gnuastro again "    \
                "with `--enable-bin-op-double' to not get this error "  \
                "any more. However, if you commonly deal with %s type " \
                "data, also enable %s with a similar option at "        \
                "configure time to greatly increase running time and "  \
                "avoid unnecessary RAM and CPU resources. Run"          \
                "`./configure --help' in Gnuastro's top source "        \
                "directory (after unpacking the tarball) for the full " \
                "list of options", typestring, typestring, typestring); \
        }                                                               \
    }





#define BINARY_OPERATOR_FOR_TYPE(LT, RT, OT, OP){                       \
    LT *la=l->array;                                                    \
    RT *ra=r->array;                                                    \
    OT *oa=o->array, *of=oa + o->size;                                  \
    if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);       \
    else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);       \
    else                 do *oa = *la++ OP *ra;   while(++oa<of);       \
  }





#define BINARY_LEFT_RIGHT_DONE(LT, RT, OP)                              \
  switch(o->type)                                                       \
    {                                                                   \
                                                                        \
    BINARY_LEFT_RIGHT_DONE_UCHAR(LT, RT, OP)                            \
    BINARY_LEFT_RIGHT_DONE_CHAR(LT, RT, OP)                             \
    BINARY_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)                            \
    BINARY_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)                           \
    BINARY_LEFT_RIGHT_DONE_INT(LT, RT, OP)                              \
    BINARY_LEFT_RIGHT_DONE_UINT(LT, RT, OP)                             \
    BINARY_LEFT_RIGHT_DONE_ULONG(LT, RT, OP)                            \
    BINARY_LEFT_RIGHT_DONE_LONG(LT, RT, OP)                             \
    BINARY_LEFT_RIGHT_DONE_LONGLONG(LT, RT, OP)                         \
    BINARY_LEFT_RIGHT_DONE_FLOAT(LT, RT, OP)                            \
    BINARY_LEFT_RIGHT_DONE_DOUBLE(LT, RT, OP)                           \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for o->type in BINARY_LEFT_RIGHT_DONE", o->type);          \
    }





#define BINARY_LEFT_DONE(LT, OP)                                        \
  switch(r->type)                                                       \
    {                                                                   \
                                                                        \
    BINARY_LEFT_DONE_UCHAR(LT, OP)                                      \
    BINARY_LEFT_DONE_CHAR(LT, OP)                                       \
    BINARY_LEFT_DONE_USHORT(LT, OP)                                     \
    BINARY_LEFT_DONE_SHORT(LT, OP)                                      \
    BINARY_LEFT_DONE_UINT(LT, OP)                                       \
    BINARY_LEFT_DONE_INT(LT, OP)                                        \
    BINARY_LEFT_DONE_ULONG(LT, OP)                                      \
    BINARY_LEFT_DONE_LONG(LT, OP)                                       \
    BINARY_LEFT_DONE_LONGLONG(LT, OP)                                   \
    BINARY_LEFT_DONE_FLOAT(LT, OP)                                      \
    BINARY_LEFT_DONE_DOUBLE(LT, OP)                                     \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for r->type in BINARY_LEFT_DONE", r->type);                \
    }




















/************************************************************************/
/*************              Top level macro             *****************/
/************************************************************************/
/* Prepare the inputs and output for binary operations and do the job.*/
#define BINARY_INTERNAL(OP, OUT_TYPE) {                                 \
                                                                        \
  /* Read the variable arguments. */                                    \
  /* `lo' and `ro' keep the original data, in case their type isn't */  \
  /* built (based on configure options are configure time). */          \
  size_t out_size, minmapsize;                                          \
  int ntype, otype, final_otype;                                        \
  gal_data_t *l, *r, *lo, *ro, *tmp_o;                                  \
                                                                        \
                                                                        \
  /* Prepare original data structures from the input arguments. */      \
  lo = va_arg(va, gal_data_t *);                                        \
  ro = va_arg(va, gal_data_t *);                                        \
                                                                        \
                                                                        \
  /* Simple sanity check on the input sizes */                          \
  if( !( (flags & GAL_DATA_ARITH_NUMOK) && (lo->size==1 || ro->size==1))\
      && gal_data_dsize_is_different(lo, ro) )                          \
    error(EXIT_FAILURE, 0, "ini BINARY_INTERNAL, the input datasets "   \
          "don't have the same dimension/size");                        \
                                                                        \
                                                                        \
  /* Set the output type and size. */                                   \
  minmapsize = ( lo->minmapsize < ro->minmapsize                        \
                 ? lo->minmapsize : ro->minmapsize );                   \
  out_size = lo->size > ro->size ? lo->size : ro->size;                 \
  final_otype = OUT_TYPE ? OUT_TYPE : gal_data_out_type(lo, ro);        \
                                                                        \
                                                                        \
  /* Make sure the input arrays have one of the compiled types. */      \
  BINARY_CONVERT_TO_COMPILED_TYPE(lo, l);                               \
  BINARY_CONVERT_TO_COMPILED_TYPE(ro, r);                               \
                                                                        \
                                                                        \
  /* Temporary output type (in case its type isn't compiled). */        \
  /* Note that the final output of this macro is put in `ntype'. */     \
  BINARY_TYPE_FOR_CONVERT_TO_COMPILED_TYPE(final_otype);                \
  otype=ntype;                                                          \
                                                                        \
                                                                        \
  /* If we want inplace output, set the output pointer to one input. */ \
  /* Note that the output type can be different from both inputs.    */ \
  if(flags & GAL_DATA_ARITH_INPLACE)                                    \
    {                                                                   \
      if     (l->type==otype && out_size==l->size)   o = l;             \
      else if(r->type==otype && out_size==r->size)   o = r;             \
    }                                                                   \
                                                                        \
                                                                        \
  /* If the output pointer was not set for any reason, allocate it. */  \
  /* For `mmapsize', note that since its `size_t', it will always be */ \
  /* Positive. The `-1' that is recommended to give when you want the */\
  /* value in RAM is actually the largest possible memory location. */  \
  /* So we just have to choose the smaller minmapsize of the two to */  \
  /* decide if the output array should be in RAM or not. */             \
  if(o==NULL)                                                           \
    o = gal_data_alloc(NULL, otype,                                     \
                       l->size>1 ? l->ndim  : r->ndim,                  \
                       l->size>1 ? l->dsize : r->dsize,                 \
                       l->size>1 ? l->wcs : r->wcs, 0, minmapsize );    \
                                                                        \
                                                                        \
  /* Do the operations based on the different types. */                 \
  switch(l->type)                                                       \
    {                                                                   \
                                                                        \
    BINARY_MULTISWITCH_UCHAR(OP)                                        \
    BINARY_MULTISWITCH_CHAR(OP)                                         \
    BINARY_MULTISWITCH_USHORT(OP)                                       \
    BINARY_MULTISWITCH_SHORT(OP)                                        \
    BINARY_MULTISWITCH_UINT(OP)                                         \
    BINARY_MULTISWITCH_INT(OP)                                          \
    BINARY_MULTISWITCH_ULONG(OP)                                        \
    BINARY_MULTISWITCH_LONG(OP)                                         \
    BINARY_MULTISWITCH_LONGLONG(OP)                                     \
    BINARY_MULTISWITCH_FLOAT(OP)                                        \
    BINARY_MULTISWITCH_DOUBLE(OP)                                       \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for l->type in BINARY_MULTISWITCH", l->type);              \
    }                                                                   \
                                                                        \
                                                                        \
  /* Clean up. Note that if the input arrays can be freed, and any of */\
  /* right or left arrays needed conversion, */                         \
  /*`BINARY_CONVERT_TO_COMPILED_TYPE' has already freed the input */    \
  /* arrays, and we only have `r' and `l' allocated in any case. */     \
  /* Alternatively, when the inputs shouldn't be freed, the only */     \
  /* allocated spaces are the `r' and `l' arrays if their types */      \
  /* weren't compiled for binary operations, we can tell this from */   \
  /* the pointers: if they are different from the original pointers, */ \
  /* they were allocated. */                                            \
  if(flags & GAL_DATA_ARITH_FREE)                                       \
    {                                                                   \
      if     (o==l)       gal_data_free(r);                             \
      else if(o==r)       gal_data_free(l);                             \
      else              { gal_data_free(l); gal_data_free(r); }         \
    }                                                                   \
  else                                                                  \
    {                                                                   \
      if(l!=lo)           gal_data_free(l);                             \
      if(r!=ro)           gal_data_free(r);                             \
    }                                                                   \
                                                                        \
                                                                        \
  /* In case otype and final_otype aren't equal, we need to convert */  \
  /* the output data structure to the proper type. */                   \
  if(otype!=final_otype)                                                \
    {                                                                   \
      tmp_o=gal_data_copy_to_new_type(o, final_otype);                  \
      gal_data_free(o);                                                 \
      o=tmp_o;                                                          \
    }                                                                   \
}





#endif
