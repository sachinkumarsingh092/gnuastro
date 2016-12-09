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
/*************      Possibly set onlyint types to convert    *************/
/************************************************************************/
#if GAL_CONFIG_BIN_OP_UCHAR == 1
#define ONLYINT_LEFT_RIGHT_DONE_UCHAR(LT, RT, OP)                       \
    case GAL_DATA_TYPE_UCHAR:                                           \
      ONLYINT_OPERATOR_FOR_TYPE(LT, RT, unsigned char, OP);             \
      break;
#define ONLYINT_LEFT_DONE_UCHAR(LT, OP)                                 \
    case GAL_DATA_TYPE_UCHAR:                                           \
      ONLYINT_LEFT_RIGHT_DONE(LT, unsigned char, OP);                   \
      break;
#define ONLYINT_MULTISWITCH_UCHAR(OP)                                   \
    case GAL_DATA_TYPE_UCHAR:                                           \
      ONLYINT_LEFT_DONE(unsigned char, OP);                             \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_UCHAR(LT, RT, OP)
#define ONLYINT_LEFT_DONE_UCHAR(LT, OP)
#define ONLYINT_MULTISWITCH_UCHAR(OP)
#endif





#if GAL_CONFIG_BIN_OP_CHAR == 1
#define ONLYINT_LEFT_RIGHT_DONE_CHAR(LT, RT, OP)                        \
  case GAL_DATA_TYPE_CHAR:                                              \
    ONLYINT_OPERATOR_FOR_TYPE(LT, RT, char, OP);                        \
    break;
#define ONLYINT_LEFT_DONE_CHAR(LT, OP)                                  \
    case GAL_DATA_TYPE_CHAR:                                            \
      ONLYINT_LEFT_RIGHT_DONE(LT, char, OP);                            \
      break;
#define ONLYINT_MULTISWITCH_CHAR(OP)                                    \
    case GAL_DATA_TYPE_CHAR:                                            \
      ONLYINT_LEFT_DONE(char, OP);                                      \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_CHAR(LT, RT, OP)
#define ONLYINT_LEFT_DONE_CHAR(LT, OP)
#define ONLYINT_MULTISWITCH_CHAR(OP)
#endif





#if GAL_CONFIG_BIN_OP_USHORT == 1
#define ONLYINT_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)                      \
  case GAL_DATA_TYPE_USHORT:                                            \
    ONLYINT_OPERATOR_FOR_TYPE(LT, RT, unsigned short, OP);              \
    break;
#define ONLYINT_LEFT_DONE_USHORT(LT, OP)                                \
    case GAL_DATA_TYPE_USHORT:                                          \
      ONLYINT_LEFT_RIGHT_DONE(LT, unsigned short, OP);                  \
      break;
#define ONLYINT_MULTISWITCH_USHORT(OP)                                  \
    case GAL_DATA_TYPE_USHORT:                                          \
      ONLYINT_LEFT_DONE(unsigned short, OP);                            \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_USHORT(LT, RT, OP)
#define ONLYINT_LEFT_DONE_USHORT(LT, OP)
#define ONLYINT_MULTISWITCH_USHORT(OP)
#endif





#if GAL_CONFIG_BIN_OP_SHORT == 1
#define ONLYINT_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)                       \
  case GAL_DATA_TYPE_SHORT:                                             \
    ONLYINT_OPERATOR_FOR_TYPE(LT, RT, short, OP);                       \
    break;
#define ONLYINT_LEFT_DONE_SHORT(LT, OP)                                 \
    case GAL_DATA_TYPE_SHORT:                                           \
      ONLYINT_LEFT_RIGHT_DONE(LT, short, OP);                           \
      break;
#define ONLYINT_MULTISWITCH_SHORT(OP)                                   \
    case GAL_DATA_TYPE_SHORT:                                           \
      ONLYINT_LEFT_DONE(short, OP);                                     \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_SHORT(LT, RT, OP)
#define ONLYINT_LEFT_DONE_SHORT(LT, OP)
#define ONLYINT_MULTISWITCH_SHORT(OP)
#endif





#if GAL_CONFIG_BIN_OP_UINT == 1
#define ONLYINT_LEFT_RIGHT_DONE_UINT(LT, RT, OP)                        \
  case GAL_DATA_TYPE_UINT:                                              \
    ONLYINT_OPERATOR_FOR_TYPE(LT, RT, unsigned int, OP);                \
    break;
#define ONLYINT_LEFT_DONE_UINT(LT, OP)                                  \
    case GAL_DATA_TYPE_UINT:                                            \
      ONLYINT_LEFT_RIGHT_DONE(LT, unsigned int, OP);                    \
      break;
#define ONLYINT_MULTISWITCH_UINT(OP)                                    \
    case GAL_DATA_TYPE_UINT:                                            \
      ONLYINT_LEFT_DONE(unsigned int, OP);                              \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_UINT(LT, RT, OP)
#define ONLYINT_LEFT_DONE_UINT(LT, OP)
#define ONLYINT_MULTISWITCH_UINT(OP)
#endif





#if GAL_CONFIG_BIN_OP_INT == 1
#define ONLYINT_LEFT_RIGHT_DONE_INT(LT, RT, OP)                         \
  case GAL_DATA_TYPE_INT:                                               \
    ONLYINT_OPERATOR_FOR_TYPE(LT, RT, int, OP);                         \
    break;
#define ONLYINT_LEFT_DONE_INT(LT, OP)                                   \
    case GAL_DATA_TYPE_INT:                                             \
      ONLYINT_LEFT_RIGHT_DONE(LT, int, OP);                             \
      break;
#define ONLYINT_MULTISWITCH_INT(OP)                                     \
    case GAL_DATA_TYPE_INT:                                             \
      ONLYINT_LEFT_DONE(int, OP);                                       \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_INT(LT, RT, OP)
#define ONLYINT_LEFT_DONE_INT(LT, OP)
#define ONLYINT_MULTISWITCH_INT(OP)
#endif





#if GAL_CONFIG_BIN_OP_ULONG == 1
#define ONLYINT_LEFT_RIGHT_DONE_ULONG(LT, RT, OP)                       \
    case GAL_DATA_TYPE_ULONG:                                           \
      ONLYINT_OPERATOR_FOR_TYPE(LT, RT, unsigned long, OP);             \
      break;
#define ONLYINT_LEFT_DONE_ULONG(LT, OP)                                 \
    case GAL_DATA_TYPE_ULONG:                                           \
      ONLYINT_LEFT_RIGHT_DONE(LT, unsigned long, OP);                   \
      break;
#define ONLYINT_MULTISWITCH_ULONG(OP)                                   \
    case GAL_DATA_TYPE_ULONG:                                           \
      ONLYINT_LEFT_DONE(unsigned long, OP);                             \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_ULONG(LT, RT, OP)
#define ONLYINT_LEFT_DONE_ULONG(LT, OP)
#define ONLYINT_MULTISWITCH_ULONG(OP)
#endif





#if GAL_CONFIG_BIN_OP_LONG == 1
#define ONLYINT_LEFT_RIGHT_DONE_LONG(LT, RT, OP)                        \
    case GAL_DATA_TYPE_LONG:                                            \
      ONLYINT_OPERATOR_FOR_TYPE(LT, RT, long, OP);                      \
      break;
#define ONLYINT_LEFT_DONE_LONG(LT, OP)                                  \
    case GAL_DATA_TYPE_LONG:                                            \
      ONLYINT_LEFT_RIGHT_DONE(LT, long, OP);                            \
      break;
#define ONLYINT_MULTISWITCH_LONG(OP)                                    \
    case GAL_DATA_TYPE_LONG:                                            \
      ONLYINT_LEFT_DONE(long, OP);                                      \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_LONG(LT, RT, OP)
#define ONLYINT_LEFT_DONE_LONG(LT, OP)
#define ONLYINT_MULTISWITCH_LONG(OP)
#endif





#if GAL_CONFIG_BIN_OP_LONGLONG == 1
#define ONLYINT_LEFT_RIGHT_DONE_LONGLONG(LT, RT, OP)                    \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      ONLYINT_OPERATOR_FOR_TYPE(LT, RT, LONGLONG, OP);                  \
      break;
#define ONLYINT_LEFT_DONE_LONGLONG(LT, OP)                              \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      ONLYINT_LEFT_RIGHT_DONE(LT, long long, OP);                       \
      break;
#define ONLYINT_MULTISWITCH_LONGLONG(OP)                                \
    case GAL_DATA_TYPE_LONGLONG:                                        \
      ONLYINT_LEFT_DONE(LONGLONG, OP);                                  \
      break;
#else
#define ONLYINT_LEFT_RIGHT_DONE_LONGLONG(LT, RT, OP)
#define ONLYINT_LEFT_DONE_LONGLONG(LT, OP)
#define ONLYINT_MULTISWITCH_LONGLONG(OP)
#endif




















/************************************************************************/
/*************       Macros for specifying the type     *****************/
/************************************************************************/

#define ONLYINT_OPERATOR_FOR_TYPE(LT, RT, OT, OP){                      \
    LT *la=l->array;                                                    \
    RT *ra=r->array;                                                    \
    OT *oa=o->array, *of=oa + o->size;                                  \
    if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);       \
    else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);       \
    else                 do *oa = *la++ OP *ra;   while(++oa<of);       \
  }





#define ONLYINT_LEFT_RIGHT_DONE(LT, RT, OP)                             \
  switch(o->type)                                                       \
    {                                                                   \
                                                                        \
      ONLYINT_LEFT_RIGHT_DONE_UCHAR(LT, RT, OP);                        \
      ONLYINT_LEFT_RIGHT_DONE_CHAR(LT, RT, OP);                         \
      ONLYINT_LEFT_RIGHT_DONE_SHORT(LT, RT, OP);                        \
      ONLYINT_LEFT_RIGHT_DONE_USHORT(LT, RT, OP);                       \
      ONLYINT_LEFT_RIGHT_DONE_INT(LT, RT, OP);                          \
      ONLYINT_LEFT_RIGHT_DONE_UINT(LT, RT, OP);                         \
      ONLYINT_LEFT_RIGHT_DONE_ULONG(LT, RT, OP);                        \
      ONLYINT_LEFT_RIGHT_DONE_LONG(LT, RT, OP);                         \
      ONLYINT_LEFT_RIGHT_DONE_LONGLONG(LT, RT, OP);                     \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for o->type in ONLYINT_LEFT_RIGHT_DONE", o->type);         \
    }





#define ONLYINT_LEFT_DONE(LT, OP)                                       \
  switch(r->type)                                                       \
    {                                                                   \
                                                                        \
      ONLYINT_LEFT_DONE_UCHAR(LT, OP);                                  \
      ONLYINT_LEFT_DONE_CHAR(LT, OP);                                   \
      ONLYINT_LEFT_DONE_USHORT(LT, OP);                                 \
      ONLYINT_LEFT_DONE_SHORT(LT, OP);                                  \
      ONLYINT_LEFT_DONE_UINT(LT, OP);                                   \
      ONLYINT_LEFT_DONE_INT(LT, OP);                                    \
      ONLYINT_LEFT_DONE_ULONG(LT, OP);                                  \
      ONLYINT_LEFT_DONE_LONG(LT, OP);                                   \
      ONLYINT_LEFT_DONE_LONGLONG(LT, OP);                               \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for r->type in ONLYINT_LEFT_DONE", r->type);               \
    }





#define ONLYINT_OPERATOR_DONE(OP)                                       \
  switch(l->type)                                                       \
    {                                                                   \
                                                                        \
      ONLYINT_MULTISWITCH_UCHAR(OP);                                    \
      ONLYINT_MULTISWITCH_CHAR(OP);                                     \
      ONLYINT_MULTISWITCH_USHORT(OP);                                   \
      ONLYINT_MULTISWITCH_SHORT(OP);                                    \
      ONLYINT_MULTISWITCH_UINT(OP);                                     \
      ONLYINT_MULTISWITCH_INT(OP);                                      \
      ONLYINT_MULTISWITCH_ULONG(OP);                                    \
      ONLYINT_MULTISWITCH_LONG(OP);                                     \
      ONLYINT_MULTISWITCH_LONGLONG(OP);                                 \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type %d not recognized in "               \
            "for l->type in data_arithmetic_onlyint", l->type);         \
    }




















/********************************************************************/
/****************        Intermediate functions       ***************/
/********************************************************************/
static int
onlyint_type_for_convert_to_compiled_type(int intype)
{
  switch(intype)
    {
    case GAL_DATA_TYPE_UCHAR:
      if(GAL_CONFIG_BIN_OP_UCHAR) return GAL_DATA_TYPE_UCHAR;
      else
        {
          if     (GAL_CONFIG_BIN_OP_USHORT)   return GAL_DATA_TYPE_USHORT;
          else if(GAL_CONFIG_BIN_OP_SHORT)    return GAL_DATA_TYPE_SHORT;
          else if(GAL_CONFIG_BIN_OP_UINT)     return GAL_DATA_TYPE_UINT;
          else if(GAL_CONFIG_BIN_OP_INT)      return GAL_DATA_TYPE_INT;
          else if(GAL_CONFIG_BIN_OP_ULONG)    return GAL_DATA_TYPE_ULONG;
          else if(GAL_CONFIG_BIN_OP_LONG)     return GAL_DATA_TYPE_LONG;
          else if(GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_CHAR:
      if(GAL_CONFIG_BIN_OP_CHAR) return GAL_DATA_TYPE_CHAR;
      else
        {
          if     (GAL_CONFIG_BIN_OP_SHORT)    return GAL_DATA_TYPE_SHORT;
          else if(GAL_CONFIG_BIN_OP_INT)      return GAL_DATA_TYPE_INT;
          else if(GAL_CONFIG_BIN_OP_LONG)     return GAL_DATA_TYPE_LONG;
          else if(GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_USHORT:
      if(GAL_CONFIG_BIN_OP_USHORT) return GAL_DATA_TYPE_USHORT;
      else
        {
          if     (GAL_CONFIG_BIN_OP_UINT)     return GAL_DATA_TYPE_UINT;
          else if(GAL_CONFIG_BIN_OP_INT)      return GAL_DATA_TYPE_INT;
          else if(GAL_CONFIG_BIN_OP_ULONG)    return GAL_DATA_TYPE_ULONG;
          else if(GAL_CONFIG_BIN_OP_LONG)     return GAL_DATA_TYPE_LONG;
          else if(GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_SHORT:
      if(GAL_CONFIG_BIN_OP_SHORT) return GAL_DATA_TYPE_SHORT;
      else
        {
          if     (GAL_CONFIG_BIN_OP_INT)      return GAL_DATA_TYPE_INT;
          else if(GAL_CONFIG_BIN_OP_LONG)     return GAL_DATA_TYPE_LONG;
          else if(GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_UINT:
      if(GAL_CONFIG_BIN_OP_UINT) return GAL_DATA_TYPE_UINT;
      else
        {
          if     (GAL_CONFIG_BIN_OP_ULONG)    return GAL_DATA_TYPE_ULONG;
          else if(GAL_CONFIG_BIN_OP_LONG)     return GAL_DATA_TYPE_LONG;
          else if(GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_INT:
      if(GAL_CONFIG_BIN_OP_INT) return GAL_DATA_TYPE_INT;
      else
        {
          if     (GAL_CONFIG_BIN_OP_LONG)     return GAL_DATA_TYPE_LONG;
          else if(GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_ULONG:
      if(GAL_CONFIG_BIN_OP_ULONG) return GAL_DATA_TYPE_ULONG;
      else
        {
          if     (GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_LONG:
      if(GAL_CONFIG_BIN_OP_LONG) return GAL_DATA_TYPE_LONG;
      else
        {
          if     (GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
        }
      break;

    case GAL_DATA_TYPE_LONGLONG:
      if(GAL_CONFIG_BIN_OP_LONGLONG) return GAL_DATA_TYPE_LONGLONG;
      break;

    default:
      error(EXIT_FAILURE, 0, "type %d not recognized in "
            "ONLYINT_CONVERT_TO_COMPILED_TYPE (note that onlyint "
            "operators only accept integer types)", intype);
    }

  return 0;
}





/* Note that for signed types, we won't be considering the unsigned types
   of the larger types. */
gal_data_t *
onlyint_convert_to_compiled_type(gal_data_t *in, unsigned char flags)
{
  int ntype;
  char *typestring;
  gal_data_t *out=NULL;

  /* Set the best compiled type. */
  ntype=onlyint_type_for_convert_to_compiled_type(in->type);

  /* If type is not compiled, then convert the dataset to the first
     compiled larger type. */
  if(in->type==ntype)
    out=in;
  else
    {
      if(ntype)
        {
          out=gal_data_copy_to_new_type(in, ntype);
          if(flags & GAL_DATA_ARITH_FREE)
            { gal_data_free(in); in=NULL; }
        }
      else
        {
          typestring=gal_data_type_string(in->type);
          error(EXIT_FAILURE, 0, "The given %s type data given to "
                "onlyint operators is not compiled for native operation "
                "and no larger types are compiled either.\n\nThe "
                "largest type (which can act as a fallback for any "
                "input type is double, so configure Gnuastro again "
                "with `--enable-bin-op-double' to not get this error "
                "any more. However, if you commonly deal with %s type "
                "data, also enable %s with a similar option at "
                "configure time to greatly increase running time and "
                "avoid unnecessary RAM and CPU resources. Run"
                "`./configure --help' in Gnuastro's top source "
                "directory (after unpacking the tarball) for the full "
                "list of options", typestring, typestring, typestring);
        }
    }

  /* Return the output data structure */
  if(out==NULL)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s, so we can fix "
          "the problem. For some reason, the `out' array in "
          "`onlyint_convert_to_compiled_type' is not set", PACKAGE_BUGREPORT);
  return out;
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
  size_t out_size, minmapsize;
  gal_data_t *l, *r, *o=NULL, *tmp_o;
  int out_type=0, otype, final_otype;


  /* Simple sanity check on the input sizes */
  if( !( (flags & GAL_DATA_ARITH_NUMOK) && (lo->size==1 || ro->size==1))
      && gal_data_dsize_is_different(lo, ro) )
    error(EXIT_FAILURE, 0, "ini ONLYINT_INTERNAL, the input datasets "
          "don't have the same dimension/size");


  /* Make sure the input arrays have one of the compiled types. */
  l=onlyint_convert_to_compiled_type(lo, flags);
  r=onlyint_convert_to_compiled_type(ro, flags);


  /* For the left/right shift bitwise operators, the length of the integer
     matters. So if the lengths of the inputs have changed with these two
     operators, then print a warning. */
  if(operator==GAL_DATA_OPERATOR_BITLSH || GAL_DATA_OPERATOR_BITRSH)
    if(lo->type!=l->type || ro->type!=r->type)
      error(EXIT_FAILURE, 0, "at least one of the input types to the bitwise "
            "left or right shift operators was not compiled. The result "
            "will thus not be what is expected. To configure Gnuastro with "
            "the respective type, use the following options at configure "
            "time: `--enable-bin-op-TYPE', run `./configure --help' to see "
            "the full list.");


  /* Set the output type. By default, `out_type' is initialized to zero,
     this means that the type of the output data structure will be
     determined based on the inputs. However, for the comparison operators,
     the output type is either 0 or 1, so we will set the output to
     unsigned character to save space and memory. Note that since this
     switch check, is only relevant for certain operators, we don't need a
     `default' statement.*/
  final_otype = out_type ? out_type : gal_data_out_type(lo, ro);
  otype=onlyint_type_for_convert_to_compiled_type(final_otype);


  /* Set the output sizes. */
  minmapsize = ( lo->minmapsize < ro->minmapsize
                 ? lo->minmapsize : ro->minmapsize );
  out_size = lo->size > ro->size ? lo->size : ro->size;


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
                       l->size>1 ? l->wcs : r->wcs, 0, minmapsize );


  /* Start setting the operator and operands. */
  switch(operator)
    {
    case GAL_DATA_OPERATOR_MODULO:     ONLYINT_OPERATOR_DONE( %  ); break;
    case GAL_DATA_OPERATOR_BITAND:     ONLYINT_OPERATOR_DONE( &  ); break;
    case GAL_DATA_OPERATOR_BITOR:      ONLYINT_OPERATOR_DONE( |  ); break;
    case GAL_DATA_OPERATOR_BITXOR:     ONLYINT_OPERATOR_DONE( ^  ); break;
    case GAL_DATA_OPERATOR_BITLSH:     ONLYINT_OPERATOR_DONE( << ); break;
    case GAL_DATA_OPERATOR_BITRSH:     ONLYINT_OPERATOR_DONE( >> ); break;
    default:
      error(EXIT_FAILURE, 0, "Operator code %d not recognized in "
            "data_arithmetic_onlyint when preparing for the operation",
            operator);
    }


  /* Clean up. Note that if the input arrays can be freed, and any of right
     or left arrays needed conversion, `ONLYINT_CONVERT_TO_COMPILED_TYPE'
     has already freed the input arrays, and we only have `r' and `l'
     allocated in any case. Alternatively, when the inputs shouldn't be
     freed, the only allocated spaces are the `r' and `l' arrays if their
     types weren't compiled for onlyint operations, we can tell this from
     the pointers: if they are different from the original pointers, they
     were allocated. */
  if(flags & GAL_DATA_ARITH_FREE)
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

  /* In case otype and final_otype aren't equal, we need to convert the
     output data structure to the proper type. */
  if(otype!=final_otype)
    {
      tmp_o=gal_data_copy_to_new_type(o, final_otype);
      gal_data_free(o);
      o=tmp_o;
    }

  /* Return */
  return o;
}
