/*********************************************************************
blank -- Deal with blank values in a datasets
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/blank.h>

#include <checkset.h>




/* Write the blank value of the type into an already allocate space. Note
   that for STRINGS, pointer should actually be `char **'. */
void
gal_blank_write(void *ptr, uint8_t type)
{
  switch(type)
    {
    /* Numeric types */
    case GAL_DATA_TYPE_UINT8:   *(uint8_t  *)ptr = GAL_BLANK_UINT8;    break;
    case GAL_DATA_TYPE_INT8:    *(int8_t   *)ptr = GAL_BLANK_INT8;     break;
    case GAL_DATA_TYPE_UINT16:  *(uint16_t *)ptr = GAL_BLANK_UINT16;   break;
    case GAL_DATA_TYPE_INT16:   *(int16_t  *)ptr = GAL_BLANK_INT16;    break;
    case GAL_DATA_TYPE_UINT32:  *(uint32_t *)ptr = GAL_BLANK_UINT32;   break;
    case GAL_DATA_TYPE_INT32:   *(int32_t  *)ptr = GAL_BLANK_INT32;    break;
    case GAL_DATA_TYPE_UINT64:  *(uint64_t *)ptr = GAL_BLANK_UINT64;   break;
    case GAL_DATA_TYPE_INT64:   *(int64_t  *)ptr = GAL_BLANK_INT64;    break;
    case GAL_DATA_TYPE_FLOAT32: *(float    *)ptr = GAL_BLANK_FLOAT32;  break;
    case GAL_DATA_TYPE_FLOAT64: *(double   *)ptr = GAL_BLANK_FLOAT64;  break;

    /* String type. */
    case GAL_DATA_TYPE_STRING:
      gal_checkset_allocate_copy(GAL_BLANK_STRING, ptr);
      break;

    /* Complex types */
    case GAL_DATA_TYPE_COMPLEX32:
    case GAL_DATA_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "complex types are not yet supported in "
            "`gal_blank_write'");

    /* Unrecognized type. */
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_blank_write'", type);
    }
}





/* Allocate some space for the given type and put the blank value into
   it. */
void *
gal_blank_alloc_write(uint8_t type)
{
  void *out;

  /* Allocate the space to keep the blank value. */
  out=gal_data_malloc_array(type, 1);

  /* Put the blank value in the allcated space. */
  gal_blank_write(out, type);

  /* Return the allocated space. */
  return out;
}





/* Return 1 if the dataset has a blank value and zero if it doesn't. */
#define HAS_BLANK(IT) {                                         \
    IT b, *a=data->array, *af=a+data->size;                     \
    gal_blank_write(&b, data->type);                            \
    if(b==b) do if(*a==b)   return 1; while(++a<af);            \
    else     do if(*a!=*a)  return 1; while(++a<af);            \
  }
int
gal_blank_present(gal_data_t *data)
{
  char **str=data->array, **strf=str+data->size;

  /* If there is nothing in the array (its size is zero), then return 0 (no
     blank is present. */
  if(data->size==0) return 0;

  /* Go over the pixels and check: */
  switch(data->type)
    {
    /* Numeric types */
    case GAL_DATA_TYPE_UINT8:     HAS_BLANK( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      HAS_BLANK( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    HAS_BLANK( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     HAS_BLANK( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    HAS_BLANK( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     HAS_BLANK( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    HAS_BLANK( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     HAS_BLANK( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   HAS_BLANK( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   HAS_BLANK( double   );    break;

    /* String. */
    case GAL_DATA_TYPE_STRING:
      do if(!strcmp(*str++,GAL_BLANK_STRING)) return 1; while(str<strf);
      break;

    /* Complex types */
    case GAL_DATA_TYPE_COMPLEX32:
    case GAL_DATA_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "complex types are not yet supported in "
            "`gal_blank_write'");

    /* Bit */
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "bit type datasets are not yet supported in "
            "`gal_blank_present'");

    default:
      error(EXIT_FAILURE, 0, "a bug! type value (%d) not recognized "
            "in `gal_blank_present'", data->type);
    }

  /* If there was a blank value, then the function would have returned with
     a value of 1. So if it reaches here, then we can be sure that there
     was no blank values, hence, return 0. */
  return 0;
}







/* Create a dataset of the the same size as the input, but with an uint8_t
   type that has a value of 1 for data that are blank and 0 for those that
   aren't. */
#define FLAG_BLANK(IT) {                                                \
    IT b, *a=data->array;                                               \
    gal_blank_write(&b, data->type);                                    \
    if(b==b) /* Blank value can be checked with the equal comparison */ \
      do { *o = *a==b;  ++a; } while(++o<of);                           \
    else     /* Blank value will fail with the equal comparison */      \
      do { *o = *a!=*a; ++a; } while(++o<of);                           \
  }
gal_data_t *
gal_blank_flag(gal_data_t *data)
{
  uint8_t *o, *of;
  gal_data_t *out;
  char **str=data->array, **strf=str+data->size;

  /* Allocate the output array. */
  out=gal_data_alloc(NULL, GAL_DATA_TYPE_UINT8, data->ndim, data->dsize,
                     data->wcs, 0, data->minmapsize, data->name, data->unit,
                     data->comment);

  /* Set the pointers for easy looping. */
  of=(o=out->array)+data->size;

  /* Go over the pixels and set the output values. */
  switch(data->type)
    {
    /* Numeric types */
    case GAL_DATA_TYPE_UINT8:     FLAG_BLANK( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:      FLAG_BLANK( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:    FLAG_BLANK( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:     FLAG_BLANK( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:    FLAG_BLANK( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:     FLAG_BLANK( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:    FLAG_BLANK( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:     FLAG_BLANK( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:   FLAG_BLANK( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:   FLAG_BLANK( double   );    break;

    /* String. */
    case GAL_DATA_TYPE_STRING:
      do *o++ = !strcmp(*str,GAL_BLANK_STRING); while(++str<strf);
      break;

    /* Currently unsupported types. */
    case GAL_DATA_TYPE_BIT:
    case GAL_DATA_TYPE_COMPLEX32:
    case GAL_DATA_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "%s type not yet supported in `gal_blank_flag'",
            gal_data_type_as_string(data->type, 1));

    /* Bad input. */
    default:
      error(EXIT_FAILURE, 0, "type value (%d) not recognized "
            "in `gal_blank_flag'", data->type);
    }

  /* Return */
  return out;
}





/* Remove blank elements from a dataset, convert it to a 1D dataset and
   adjust the size properly. In practice this function doesn't `realloc'
   the input array, all it does is to shift the blank eleemnts to the end
   and adjust the size elements of the `gal_data_t'. */
#define BLANK_REMOVE(IT) {                                              \
    IT b, *a=data->array, *af=a+data->size, *o=data->array;             \
    if( gal_blank_present(data) )                                       \
      {                                                                 \
        gal_blank_write(&b, data->type);                                \
        if(b==b)       /* Blank value can be be checked with equal. */  \
          do if(*a!=b)  { ++num; *o++=*a; } while(++a<af);              \
        else           /* Blank value will fail on equal comparison. */ \
          do if(*a==*a) { ++num; *o++=*a; } while(++a<af);              \
      }                                                                 \
    else num=data->size;                                                \
  }
void
gal_blank_remove(gal_data_t *data)
{
  size_t num=0;

  /* Shift all non-blank elements to the start of the array. */
  switch(data->type)
    {
    case GAL_DATA_TYPE_UINT8:    BLANK_REMOVE( uint8_t  );    break;
    case GAL_DATA_TYPE_INT8:     BLANK_REMOVE( int8_t   );    break;
    case GAL_DATA_TYPE_UINT16:   BLANK_REMOVE( uint16_t );    break;
    case GAL_DATA_TYPE_INT16:    BLANK_REMOVE( int16_t  );    break;
    case GAL_DATA_TYPE_UINT32:   BLANK_REMOVE( uint32_t );    break;
    case GAL_DATA_TYPE_INT32:    BLANK_REMOVE( int32_t  );    break;
    case GAL_DATA_TYPE_UINT64:   BLANK_REMOVE( uint64_t );    break;
    case GAL_DATA_TYPE_INT64:    BLANK_REMOVE( int64_t  );    break;
    case GAL_DATA_TYPE_FLOAT32:  BLANK_REMOVE( float    );    break;
    case GAL_DATA_TYPE_FLOAT64:  BLANK_REMOVE( double   );    break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_blank_remove'", data->type);
    }

  /* Adjust the size elements of the dataset. */
  data->ndim=1;
  data->dsize[0]=data->size=num;
}





/* Print the blank value as a string. */
char *
gal_blank_as_string(uint8_t type, int width)
{
  char *blank;
  switch(type)
    {
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "bit types are not implemented in "
            "`gal_data_blank_as_string' yet.");
      break;

    case GAL_DATA_TYPE_STRING:
      if(width)
        asprintf(&blank, "%*s", width, GAL_BLANK_STRING);
      else
        asprintf(&blank, "%s", GAL_BLANK_STRING);
      break;

    case GAL_DATA_TYPE_UINT8:
      if(width)
        asprintf(&blank, "%*u", width, (uint8_t)GAL_BLANK_UINT8);
      else
        asprintf(&blank, "%u",         (uint8_t)GAL_BLANK_UINT8);
      break;

    case GAL_DATA_TYPE_INT8:
      if(width)
        asprintf(&blank, "%*d", width, (int8_t)GAL_BLANK_INT8);
      else
        asprintf(&blank, "%d",         (int8_t)GAL_BLANK_INT8);
      break;

    case GAL_DATA_TYPE_UINT16:
      if(width)
        asprintf(&blank, "%*u", width, (uint16_t)GAL_BLANK_UINT16);
      else
        asprintf(&blank, "%u",         (uint16_t)GAL_BLANK_UINT16);
      break;

    case GAL_DATA_TYPE_INT16:
      if(width)
        asprintf(&blank, "%*d", width, (int16_t)GAL_BLANK_INT16);
      else
        asprintf(&blank, "%d",         (int16_t)GAL_BLANK_INT16);
      break;

    case GAL_DATA_TYPE_UINT32:
      if(width)
        asprintf(&blank, "%*u", width, (uint32_t)GAL_BLANK_UINT32);
      else
        asprintf(&blank, "%u",         (uint32_t)GAL_BLANK_UINT32);
      break;

    case GAL_DATA_TYPE_INT32:
      if(width)
        asprintf(&blank, "%*d", width, (int32_t)GAL_BLANK_INT32);
      else
        asprintf(&blank, "%d",         (int32_t)GAL_BLANK_INT32);
      break;

    case GAL_DATA_TYPE_UINT64:
      if(width)
        asprintf(&blank, "%*lu", width, (uint64_t)GAL_BLANK_UINT64);
      else
        asprintf(&blank, "%lu",         (uint64_t)GAL_BLANK_UINT64);
      break;

    case GAL_DATA_TYPE_INT64:
      if(width)
        asprintf(&blank, "%*ld", width, (int64_t)GAL_BLANK_INT64);
      else
        asprintf(&blank, "%ld",         (int64_t)GAL_BLANK_INT64);
      break;

    case GAL_DATA_TYPE_FLOAT32:
      if(width)
        asprintf(&blank, "%*f", width, (float)GAL_BLANK_FLOAT32);
      else
        asprintf(&blank, "%f",         (float)GAL_BLANK_FLOAT32);
      break;

    case GAL_DATA_TYPE_FLOAT64:
      if(width)
        asprintf(&blank, "%*f", width, (double)GAL_BLANK_FLOAT64);
      else
        asprintf(&blank, "%f",         (double)GAL_BLANK_FLOAT64);
      break;

    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_blank_as_string'", type);
    }
  return blank;
}
