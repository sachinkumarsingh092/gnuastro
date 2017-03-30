/*********************************************************************
Type -- Type information and basic operations.
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
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/type.h>




size_t
gal_type_sizeof(uint8_t type)
{
  /* Allocate space for the array to keep the image. */
  switch(type)
    {
    case GAL_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support bit "
            "types, please get in touch with us to implement it.");

      /* The parenthesis after sizeof is not a function, it is actually a
         type cast, so we have put a space between size of and the
         parenthesis to highlight this. In C, `sizeof' is an operator, not
         a function.*/
    case GAL_TYPE_UINT8:     return sizeof (uint8_t);
    case GAL_TYPE_INT8:      return sizeof (int8_t);
    case GAL_TYPE_UINT16:    return sizeof (uint16_t);
    case GAL_TYPE_INT16:     return sizeof (int16_t);
    case GAL_TYPE_UINT32:    return sizeof (uint32_t);
    case GAL_TYPE_INT32:     return sizeof (int32_t);
    case GAL_TYPE_UINT64:    return sizeof (uint64_t);
    case GAL_TYPE_INT64:     return sizeof (int64_t);

    case GAL_TYPE_FLOAT32:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0, "`float` is not 32 bits on this machine");
      return sizeof (float);

    case GAL_TYPE_FLOAT64:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0, "`double` is not 64 bits on this machine");
      return sizeof (double);

    case GAL_TYPE_COMPLEX32:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0, "`float` is not 32 bits on this machine");
      return sizeof (gsl_complex_float);

    case GAL_TYPE_COMPLEX64:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0, "`double` is not 64 bits on this machine");
      return sizeof (gsl_complex);

    case GAL_TYPE_STRING:
      return sizeof (char *);

    default:
      error(EXIT_FAILURE, 0, "type value of %d not recognized in "
            "gal_data_sizeof", type);
    }

  error(EXIT_FAILURE, 0, "Control has reached the end of `gal_data_sizeof' "
        "This is a bug! Please contact us at %s so we can find the cause "
        "of the problem.", PACKAGE_BUGREPORT);
  return -1;
}





char *
gal_type_to_string(uint8_t type, int long_name)
{
  switch(type)
    {
    case GAL_TYPE_BIT:
      if(long_name) return "bit";                 else return "b";

    case GAL_TYPE_UINT8:
      if(long_name) return "uint8";               else return "u8";

    case GAL_TYPE_INT8:
      if(long_name) return "int8";                else return "i8";

    case GAL_TYPE_UINT16:
      if(long_name) return "uint16";              else return "u16";

    case GAL_TYPE_INT16:
      if(long_name) return "int16";               else return "i16";

    case GAL_TYPE_UINT32:
      if(long_name) return "uint32";              else return "u32";

    case GAL_TYPE_INT32:
      if(long_name) return "int32";               else return "i32";

    case GAL_TYPE_UINT64:
      if(long_name) return "uint64";              else return "u64";

    case GAL_TYPE_INT64:
      if(long_name) return "int64";               else return "i64";

    case GAL_TYPE_FLOAT32:
      if(long_name) return "float32";             else return "f32";

    case GAL_TYPE_FLOAT64:
      if(long_name) return "float64";             else return "f64";

    case GAL_TYPE_COMPLEX32:
      if(long_name) return "complex32";           else return "c32";

    case GAL_TYPE_COMPLEX64:
      if(long_name) return "complex64";           else return "c64";

    case GAL_TYPE_STRING:
      if(long_name) return "string";              else return "str";

    case GAL_TYPE_STRLL:
      if(long_name) return "string linked list";  else return "strll";

    default:
      error(EXIT_FAILURE, 0, "type value of %d not recognized in "
            "`gal_data_type_as_string'", type);
    }

  /* Any of the cases above should return this function, so if control
     reaches here, there is a bug. */
  error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can address "
        "the problem. For some reason control has reached the end of "
        "the `gal_data_type_as_string' function. This must not happen",
        PACKAGE_BUGREPORT);
  return NULL;
}





uint8_t
gal_type_from_string(char *str)
{
  if(      !strcmp(str, "b")     || !strcmp(str, "bit") )
    return GAL_TYPE_BIT;

  else if( !strcmp(str, "u8")    || !strcmp(str, "uint8") )
    return GAL_TYPE_UINT8;

  else if( !strcmp(str, "i8")    || !strcmp(str, "int8") )
    return GAL_TYPE_INT8;

  else if( !strcmp(str, "u16")   || !strcmp(str, "uint16") )
    return GAL_TYPE_UINT16;

  else if( !strcmp(str, "i16")   || !strcmp(str, "int16") )
    return GAL_TYPE_INT16;

  else if( !strcmp(str, "u32")   || !strcmp(str, "uint32") )
    return GAL_TYPE_UINT32;

  else if( !strcmp(str, "i32")   || !strcmp(str, "int32") )
    return GAL_TYPE_INT32;

  else if( !strcmp(str, "u64")   || !strcmp(str, "uint64") )
    return GAL_TYPE_UINT64;

  else if( !strcmp(str, "i64")   || !strcmp(str, "int64") )
    return GAL_TYPE_INT64;

  else if( !strcmp(str, "f32")   || !strcmp(str, "float32") )
    return GAL_TYPE_FLOAT32;

  else if( !strcmp(str, "f64")   || !strcmp(str, "float64") )
    return GAL_TYPE_FLOAT64;

  else if( !strcmp(str, "c32")   || !strcmp(str, "complex32") )
    return GAL_TYPE_COMPLEX32;

  else if( !strcmp(str, "c64")   || !strcmp(str, "complex64") )
    return GAL_TYPE_COMPLEX64;

  else if( !strcmp(str, "str")   || !strcmp(str, "string") )
    return GAL_TYPE_STRING;

  else
    return GAL_TYPE_INVALID;

  /* Any of the cases above should return this function, so if control
     reaches here, there is a bug. */
  error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can address "
        "the problem. For some reason control has reached the end of "
        "the `gal_data_string_as_type' function. This must not happen",
        PACKAGE_BUGREPORT);
  return 0;
}





/* Put the minimum (or maximum for the `gal_data_type_max') value for the
   type in the space (that must already be allocated before the call to
   this function) pointed to by in.  */
void
gal_type_min(uint8_t type, void *in)
{
  switch(type)
    {
    case GAL_TYPE_UINT8:        *(uint8_t *)  in = 0;            break;
    case GAL_TYPE_INT8:         *(int8_t *)   in = INT8_MIN;     break;
    case GAL_TYPE_UINT16:       *(uint16_t *) in = 0;            break;
    case GAL_TYPE_INT16:        *(int16_t *)  in = INT16_MIN;    break;
    case GAL_TYPE_UINT32:       *(uint32_t *) in = 0;            break;
    case GAL_TYPE_INT32:        *(int32_t *)  in = INT32_MIN;    break;
    case GAL_TYPE_UINT64:       *(uint64_t *) in = 0;            break;
    case GAL_TYPE_INT64:        *(int64_t *)  in = INT64_MIN;    break;
    case GAL_TYPE_FLOAT32:      *(float *)    in = -FLT_MAX;     break;
    case GAL_TYPE_FLOAT64:      *(double *)   in = -DBL_MAX;     break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_data_type_min'", type);
    }
}





void
gal_type_max(uint8_t type, void *in)
{
  switch(type)
    {
    case GAL_TYPE_UINT8:        *(uint8_t *)  in = UINT8_MAX;    break;
    case GAL_TYPE_INT8:         *(int8_t *)   in = INT8_MAX;     break;
    case GAL_TYPE_UINT16:       *(uint16_t *) in = UINT16_MAX;   break;
    case GAL_TYPE_INT16:        *(int16_t *)  in = INT16_MAX;    break;
    case GAL_TYPE_UINT32:       *(uint32_t *) in = UINT32_MAX;   break;
    case GAL_TYPE_INT32:        *(int32_t *)  in = INT32_MAX;    break;
    case GAL_TYPE_UINT64:       *(uint64_t *) in = UINT64_MAX;   break;
    case GAL_TYPE_INT64:        *(int64_t *)  in = INT64_MAX;    break;
    case GAL_TYPE_FLOAT32:      *(float *)    in = FLT_MAX;      break;
    case GAL_TYPE_FLOAT64:      *(double *)   in = DBL_MAX;      break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_data_type_min'", type);
    }
}





/* Since linked lists need a different process than arrays, for functions
   that work on both, it is convenient to simiplify the check with this
   function. */
int
gal_type_is_linked_list(uint8_t type)
{
  return type==GAL_TYPE_STRLL;
}





int
gal_type_out(int first_type, int second_type)
{
  return first_type > second_type ? first_type : second_type;
}
