/*********************************************************************
Arithmetic operations on data structures.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#ifndef __ARITHMETIC_BINARY_H__
#define __ARITHMETIC_BINARY_H__










/************************************************************************/
/*************             Low-level operators          *****************/
/************************************************************************/
/* Final step to be used by all operators and all types. */
#define BINARY_OP_OT_RT_LT_SET(OP, OT, LT, RT) {                        \
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





/* Blank values aren't defined for integer operators. */
#define BINARY_INT_OP_OT_RT_LT_SET(OP, OT, LT, RT) {               \
    LT *la=l->array;                                               \
    RT *ra=r->array;                                               \
    OT *oa=o->array, *of=oa + o->size;                             \
    if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);  \
    else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);  \
    else                 do *oa = *la++ OP *ra;   while(++oa<of);  \
  }





/* This is for operators like '&&' and '||', where the right operator is
   not necessarily read (and thus incremented). */
#define BINARY_OP_INCR_OT_RT_LT_SET(OP, OT, LT, RT) {                   \
    LT *la=l->array;                                                    \
    RT *ra=r->array;                                                    \
    OT *oa=o->array, *of=oa + o->size;                                  \
    if(l->size==r->size) do {*oa = *la++ OP *ra; ++ra;} while(++oa<of); \
    else if(l->size==1)  do {*oa = *la   OP *ra; ++ra;} while(++oa<of); \
    else                 do  *oa = *la++ OP *ra;        while(++oa<of); \
  }




















/************************************************************************/
/*************              Type specifiers             *****************/
/************************************************************************/

/* Flags for BINARY_OP_RT_LT_SET to identify the output type. */
enum arithmetic_binary_outtype_flags
{
  ARITHMETIC_BINARY_INVALID,                /* ==0 by C standard. */

  ARITHMETIC_BINARY_OUT_TYPE_LEFT,
  ARITHMETIC_BINARY_OUT_TYPE_RIGHT,
  ARITHMETIC_BINARY_OUT_TYPE_UINT8,
  ARITHMETIC_BINARY_OUT_TYPE_INCR_SEP,
};





/* For operators whose type may be any of the given inputs only for
   integers. For integer operators we have less options. */
#define BINARY_SET_OUT_INT(F, OP, LT, RT)                               \
  switch(F)                                                             \
    {                                                                   \
    case ARITHMETIC_BINARY_OUT_TYPE_LEFT:                               \
      BINARY_INT_OP_OT_RT_LT_SET(OP, LT, LT, RT);                       \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_RIGHT:                              \
      BINARY_INT_OP_OT_RT_LT_SET(OP, RT, LT, RT);                       \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to "   \
            "address the problem. %d not recognized for 'F'",           \
            "BINARY_SET_OUT_INT", PACKAGE_BUGREPORT, F);                \
    }





/* For operators whose type may be any of the given inputs. */
#define BINARY_SET_OUT(F, OP, LT, RT)                                   \
  switch(F)                                                             \
    {                                                                   \
    case ARITHMETIC_BINARY_OUT_TYPE_LEFT:                               \
      BINARY_OP_OT_RT_LT_SET(OP,      LT,      LT, RT);                 \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_RIGHT:                              \
      BINARY_OP_OT_RT_LT_SET(OP,      RT,      LT, RT);                 \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_UINT8:                              \
      BINARY_OP_OT_RT_LT_SET(OP,      uint8_t, LT, RT);                 \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_INCR_SEP:                           \
      BINARY_OP_INCR_OT_RT_LT_SET(OP, uint8_t, LT, RT);                 \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to "   \
            "address the problem. %d not recognized for 'F'",           \
            "BINARY_SET_OUT", PACKAGE_BUGREPORT, F);                    \
    }





/* Set the right operator type only for integers (integer operators can't
   take floating point types). So floating point types must not be in the
   list of possibilities. */
#define BINARY_SET_RT_INT(F, OP, LT)                                        \
  switch(r->type)                                                           \
    {                                                                       \
    case GAL_TYPE_UINT8:   BINARY_SET_OUT_INT( F, OP, LT, uint8_t  ) break; \
    case GAL_TYPE_INT8:    BINARY_SET_OUT_INT( F, OP, LT, int8_t   ) break; \
    case GAL_TYPE_UINT16:  BINARY_SET_OUT_INT( F, OP, LT, uint16_t ) break; \
    case GAL_TYPE_INT16:   BINARY_SET_OUT_INT( F, OP, LT, int16_t  ) break; \
    case GAL_TYPE_UINT32:  BINARY_SET_OUT_INT( F, OP, LT, uint32_t ) break; \
    case GAL_TYPE_INT32:   BINARY_SET_OUT_INT( F, OP, LT, int32_t  ) break; \
    case GAL_TYPE_UINT64:  BINARY_SET_OUT_INT( F, OP, LT, uint64_t ) break; \
    case GAL_TYPE_INT64:   BINARY_SET_OUT_INT( F, OP, LT, int64_t  ) break; \
    default:                                                                \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "       \
            "address the problem. %d is not a usable type code",            \
            "BINARY_SET_RT", PACKAGE_BUGREPORT, r->type);                   \
    }





/* Set the right operator type. */
#define BINARY_SET_RT(F, OP, LT)                                        \
  switch(r->type)                                                       \
    {                                                                   \
    case GAL_TYPE_UINT8:   BINARY_SET_OUT( F, OP, LT, uint8_t  ) break; \
    case GAL_TYPE_INT8:    BINARY_SET_OUT( F, OP, LT, int8_t   ) break; \
    case GAL_TYPE_UINT16:  BINARY_SET_OUT( F, OP, LT, uint16_t ) break; \
    case GAL_TYPE_INT16:   BINARY_SET_OUT( F, OP, LT, int16_t  ) break; \
    case GAL_TYPE_UINT32:  BINARY_SET_OUT( F, OP, LT, uint32_t ) break; \
    case GAL_TYPE_INT32:   BINARY_SET_OUT( F, OP, LT, int32_t  ) break; \
    case GAL_TYPE_UINT64:  BINARY_SET_OUT( F, OP, LT, uint64_t ) break; \
    case GAL_TYPE_INT64:   BINARY_SET_OUT( F, OP, LT, int64_t  ) break; \
    case GAL_TYPE_FLOAT32: BINARY_SET_OUT( F, OP, LT, float    ) break; \
    case GAL_TYPE_FLOAT64: BINARY_SET_OUT( F, OP, LT, double   ) break; \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "   \
            "address the problem. %d is not a usable type code",        \
            "BINARY_SET_RT", PACKAGE_BUGREPORT, r->type);               \
    }





/* Set the left operator type only for integers (integer operators can't
   take floating point types). So floating point types must not be in the
   list of possibilities. */
#define BINARY_SET_LT_INT(F, OP)                                        \
  switch(l->type)                                                       \
    {                                                                   \
    case GAL_TYPE_UINT8:   BINARY_SET_RT_INT( F, OP, uint8_t  ) break;  \
    case GAL_TYPE_INT8:    BINARY_SET_RT_INT( F, OP, int8_t   ) break;  \
    case GAL_TYPE_UINT16:  BINARY_SET_RT_INT( F, OP, uint16_t ) break;  \
    case GAL_TYPE_INT16:   BINARY_SET_RT_INT( F, OP, int16_t  ) break;  \
    case GAL_TYPE_UINT32:  BINARY_SET_RT_INT( F, OP, uint32_t ) break;  \
    case GAL_TYPE_INT32:   BINARY_SET_RT_INT( F, OP, int32_t  ) break;  \
    case GAL_TYPE_UINT64:  BINARY_SET_RT_INT( F, OP, uint64_t ) break;  \
    case GAL_TYPE_INT64:   BINARY_SET_RT_INT( F, OP, int64_t  ) break;  \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "   \
            "address the problem. %d is not a usable type code",        \
            "BINARY_SET_LT_INT", PACKAGE_BUGREPORT, l->type);           \
    }





/* Set the left operator type. */
#define BINARY_SET_LT(F, OP)                                            \
  switch(l->type)                                                       \
    {                                                                   \
    case GAL_TYPE_UINT8:   BINARY_SET_RT( F, OP, uint8_t  ) break;      \
    case GAL_TYPE_INT8:    BINARY_SET_RT( F, OP, int8_t   ) break;      \
    case GAL_TYPE_UINT16:  BINARY_SET_RT( F, OP, uint16_t ) break;      \
    case GAL_TYPE_INT16:   BINARY_SET_RT( F, OP, int16_t  ) break;      \
    case GAL_TYPE_UINT32:  BINARY_SET_RT( F, OP, uint32_t ) break;      \
    case GAL_TYPE_INT32:   BINARY_SET_RT( F, OP, int32_t  ) break;      \
    case GAL_TYPE_UINT64:  BINARY_SET_RT( F, OP, uint64_t ) break;      \
    case GAL_TYPE_INT64:   BINARY_SET_RT( F, OP, int64_t  ) break;      \
    case GAL_TYPE_FLOAT32: BINARY_SET_RT( F, OP, float    ) break;      \
    case GAL_TYPE_FLOAT64: BINARY_SET_RT( F, OP, double   ) break;      \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "   \
            "address the problem. %d is not a usable type code",        \
            "BINARY_SET_LT", PACKAGE_BUGREPORT, l->type);               \
    }



#endif
