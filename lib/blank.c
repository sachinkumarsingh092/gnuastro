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




/* Write the blank value of the type into an already allocate space.

   Note that for strings, pointer should actually be `char **'. */
void
gal_blank_write(void *pointer, uint8_t type)
{
  switch(type)
    {
    case GAL_DATA_TYPE_STRING:
      gal_checkset_allocate_copy(GAL_BLANK_STRING, pointer);
      break;

    case GAL_DATA_TYPE_UINT8:
      *(uint8_t *)pointer       = GAL_BLANK_UINT8;
      break;

    case GAL_DATA_TYPE_INT8:
      *(int8_t *)pointer        = GAL_BLANK_INT8;
      break;

    case GAL_DATA_TYPE_UINT16:
      *(uint16_t *)pointer      = GAL_BLANK_UINT16;
      break;

    case GAL_DATA_TYPE_INT16:
      *(int16_t *)pointer       = GAL_BLANK_INT16;
      break;

    case GAL_DATA_TYPE_UINT32:
      *(uint32_t *)pointer      = GAL_BLANK_UINT32;
      break;

    case GAL_DATA_TYPE_INT32:
      *(int32_t *)pointer       = GAL_BLANK_INT32;
      break;

    case GAL_DATA_TYPE_UINT64:
      *(uint64_t *)pointer      = GAL_BLANK_UINT64;
      break;

    case GAL_DATA_TYPE_INT64:
      *(int64_t *)pointer       = GAL_BLANK_INT64;
      break;

    case GAL_DATA_TYPE_FLOAT32:
      *(float *)pointer         = GAL_BLANK_FLOAT32;
      break;

    case GAL_DATA_TYPE_FLOAT64:
      *(double *)pointer        = GAL_BLANK_FLOAT64;
      break;

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






/* For integers a simple equality is enough. */
#define HAS_BLANK_INT(CTYPE, BLANK) {                                   \
    CTYPE *a=data->array, *af=a+data->size;                             \
    do if(*a++ == BLANK) return 1; while(a<af);                         \
  }

/* Note that a NaN value is not equal to another NaN value, so we can't use
   the easy check for cases were the blank value is NaN. Also note that
   `isnan' is actually a macro, so it works for both float and double
   types.*/
#define HAS_BLANK_FLT(CTYPE, BLANK, MULTIP) {                           \
    CTYPE *a=data->array, *af=a+(MULTIP*data->size);                    \
    if(isnan(BLANK)) do if(isnan(*a++)) return 1; while(a<af);          \
    else             do if(*a++==BLANK) return 1; while(a<af);          \
  }

/* Return 1 if the dataset has a blank value and zero if it doesn't. */
int
gal_blank_present(gal_data_t *data)
{
  char **str=data->array, **strf=str+data->size;

  /* Go over the pixels and check: */
  switch(data->type)
    {
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support bit "
            "datatype, please get in touch with us to implement it.");

    case GAL_DATA_TYPE_UINT8:
      HAS_BLANK_INT(uint8_t, GAL_BLANK_UINT8);     break;

    case GAL_DATA_TYPE_INT8:
      HAS_BLANK_INT(int8_t, GAL_BLANK_INT8);       break;

    case GAL_DATA_TYPE_UINT16:
      HAS_BLANK_INT(uint16_t, GAL_BLANK_UINT16);   break;

    case GAL_DATA_TYPE_INT16:
      HAS_BLANK_INT(int16_t, GAL_BLANK_INT16);     break;

    case GAL_DATA_TYPE_UINT32:
      HAS_BLANK_INT(uint32_t, GAL_BLANK_UINT32);   break;

    case GAL_DATA_TYPE_INT32:
      HAS_BLANK_INT(int32_t, GAL_BLANK_INT32);     break;

    case GAL_DATA_TYPE_UINT64:
      HAS_BLANK_INT(uint64_t, GAL_BLANK_UINT64);   break;

    case GAL_DATA_TYPE_INT64:
      HAS_BLANK_INT(int64_t, GAL_BLANK_INT64);     break;

    case GAL_DATA_TYPE_FLOAT32:
      HAS_BLANK_FLT(float, GAL_BLANK_FLOAT32, 1);  break;

    case GAL_DATA_TYPE_FLOAT64:
      HAS_BLANK_FLT(double, GAL_BLANK_FLOAT64, 1); break;

    case GAL_DATA_TYPE_COMPLEX32:
      HAS_BLANK_FLT(float, GAL_BLANK_FLOAT32, 2);  break;

    case GAL_DATA_TYPE_COMPLEX64:
      HAS_BLANK_FLT(double, GAL_BLANK_FLOAT64, 2); break;

    case GAL_DATA_TYPE_STRING:
      do if(!strcmp(*str++,GAL_BLANK_STRING)) return 1; while(str<strf);
      break;

    default:
      error(EXIT_FAILURE, 0, "a bug! type value (%d) not recognized "
            "in `gal_data_has_blank'", data->type);
    }

  /* If there was a blank value, then the function would have returned with
     a value of 1. So if it reaches here, then we can be sure that there
     was no blank values, hence, return 0. */
  return 0;
}





/* For integers a simple equality is enough. */
#define FLAG_BLANK_INT(CTYPE, BLANK) {                                  \
    CTYPE *a=data->array; do *o = (*a==BLANK); while(++o<of);           \
  }

/* Note that a NaN value is not equal to another NaN value, so we can't use
   the easy check for cases were the blank value is NaN. Also note that
   `isnan' is actually a macro, so it works for both float and double
   types.*/
#define FLAG_BLANK_FLT(CTYPE, BLANK) {                                  \
    CTYPE *a=data->array;                                               \
    if(isnan(BLANK)) do *o = isnan(*a++);   while(++o<of);              \
    else             do *o = (*a++==BLANK); while(++o<of);              \
  }

#define FLAG_BLANK_COMPLEX(CTYPE, BLANK) {                              \
    CTYPE *a=data->array;                                               \
    if(isnan(BLANK))                                                    \
      do { *o=(isnan(*a) || isnan(*(a+1))); a+=2; } while(++o<of);      \
    else                                                                \
      do { *o=(*a==BLANK || *(a+1)==BLANK); a+=2; } while(++o<of);      \
  }

/* Output a data-set of the the same size as the input, but with an uint8_t
   type that has a value of 1 for data that are blank and 0 for those that
   aren't. */
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
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support bit "
            "datatype, please get in touch with us to implement it.");

    case GAL_DATA_TYPE_UINT8:
      FLAG_BLANK_INT(uint8_t, GAL_BLANK_UINT8);      break;

    case GAL_DATA_TYPE_INT8:
      FLAG_BLANK_INT(int8_t, GAL_BLANK_INT8);        break;

    case GAL_DATA_TYPE_UINT16:
      FLAG_BLANK_INT(uint16_t, GAL_BLANK_UINT16);    break;

    case GAL_DATA_TYPE_INT16:
      FLAG_BLANK_INT(int16_t, GAL_BLANK_INT16);      break;

    case GAL_DATA_TYPE_UINT32:
      FLAG_BLANK_INT(uint32_t, GAL_BLANK_UINT32);    break;

    case GAL_DATA_TYPE_INT32:
      FLAG_BLANK_INT(int32_t, GAL_BLANK_INT32);      break;

    case GAL_DATA_TYPE_UINT64:
      FLAG_BLANK_INT(uint64_t, GAL_BLANK_UINT64);    break;

    case GAL_DATA_TYPE_INT64:
      FLAG_BLANK_INT(int64_t, GAL_BLANK_INT64);      break;

    case GAL_DATA_TYPE_FLOAT32:
      FLAG_BLANK_FLT(float, GAL_BLANK_FLOAT32);      break;

    case GAL_DATA_TYPE_FLOAT64:
      FLAG_BLANK_FLT(double, GAL_BLANK_FLOAT64);     break;

    case GAL_DATA_TYPE_COMPLEX32:
      FLAG_BLANK_COMPLEX(float, GAL_BLANK_FLOAT32);  break;

    case GAL_DATA_TYPE_COMPLEX64:
      FLAG_BLANK_COMPLEX(double, GAL_BLANK_FLOAT64); break;

    case GAL_DATA_TYPE_STRING:
      do *o++ = !strcmp(*str,GAL_BLANK_STRING); while(++str<strf);
      break;

    default:
      error(EXIT_FAILURE, 0, "type value (%d) not recognized "
            "in `gal_blank_flag'", data->type);
    }

  /* Return */
  return out;
}





#define BLANK_REMOVE(IT) {                                              \
    IT b, *a=data->array, *af=a+data->size, *o=data->array;             \
    if( gal_blank_present(data) )                                       \
      {                                                                 \
        gal_blank_write(&b, data->type);                                \
        if(b==b)          /* Can be checked with equal: isn't NaN */    \
          do if(*a!=b)  { ++num; *o++=*a; } while(++a<af);              \
        else /* Blank value is NaN, so equals will fail on blank elements*/ \
          do if(*a==*a) { ++num; *o++=*a; } while(++a<af);              \
      }                                                                 \
    else num=data->size;                                                \
  }


/* Remove blank elements from a dataset, convert it to a 1D dataset and
   adjust the size properly. In practice this function doesn't `realloc'
   the input array, all it does is to shift the blank eleemnts to the end
   and adjust the size elements of the `gal_data_t'. */
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
