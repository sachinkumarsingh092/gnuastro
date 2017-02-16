/*********************************************************************
data -- Structure and functions to represent/work with data
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

#include <errno.h>
#include <error.h>
#include <fcntl.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/mman.h>

#include <gnuastro/data.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/linkedlist.h>

#include <checkset.h>










/*************************************************************
 **************            Type info           ***************
 *************************************************************/

char *
gal_data_type_as_string(uint8_t type, int long_name)
{
  switch(type)
    {
    case GAL_DATA_TYPE_BIT:
      if(long_name) return "bit";                 else return "b";

    case GAL_DATA_TYPE_UINT8:
      if(long_name) return "uint8";               else return "u8";

    case GAL_DATA_TYPE_INT8:
      if(long_name) return "int8";                else return "i8";

    case GAL_DATA_TYPE_UINT16:
      if(long_name) return "uint16";              else return "u16";

    case GAL_DATA_TYPE_INT16:
      if(long_name) return "int16";               else return "i16";

    case GAL_DATA_TYPE_UINT32:
      if(long_name) return "uint32";              else return "u32";

    case GAL_DATA_TYPE_INT32:
      if(long_name) return "int32";               else return "i32";

    case GAL_DATA_TYPE_UINT64:
      if(long_name) return "uint64";              else return "u64";

    case GAL_DATA_TYPE_INT64:
      if(long_name) return "int64";               else return "i64";

    case GAL_DATA_TYPE_FLOAT32:
      if(long_name) return "float32";             else return "f32";

    case GAL_DATA_TYPE_FLOAT64:
      if(long_name) return "float64";             else return "f64";

    case GAL_DATA_TYPE_COMPLEX32:
      if(long_name) return "complex32";           else return "c32";

    case GAL_DATA_TYPE_COMPLEX64:
      if(long_name) return "complex64";           else return "c64";

    case GAL_DATA_TYPE_STRING:
      if(long_name) return "string";              else return "str";

    case GAL_DATA_TYPE_STRLL:
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
gal_data_string_as_type(char *str)
{
  if(      !strcmp(str, "b")     || !strcmp(str, "bit") )
    return GAL_DATA_TYPE_BIT;

  else if( !strcmp(str, "u8")    || !strcmp(str, "uint8") )
    return GAL_DATA_TYPE_UINT8;

  else if( !strcmp(str, "i8")    || !strcmp(str, "int8") )
    return GAL_DATA_TYPE_INT8;

  else if( !strcmp(str, "u16")   || !strcmp(str, "uint16") )
    return GAL_DATA_TYPE_UINT16;

  else if( !strcmp(str, "i16")   || !strcmp(str, "int16") )
    return GAL_DATA_TYPE_INT16;

  else if( !strcmp(str, "u32")   || !strcmp(str, "uint32") )
    return GAL_DATA_TYPE_UINT32;

  else if( !strcmp(str, "i32")   || !strcmp(str, "int32") )
    return GAL_DATA_TYPE_INT32;

  else if( !strcmp(str, "u64")   || !strcmp(str, "uint64") )
    return GAL_DATA_TYPE_UINT64;

  else if( !strcmp(str, "i64")   || !strcmp(str, "int64") )
    return GAL_DATA_TYPE_INT64;

  else if( !strcmp(str, "f32")   || !strcmp(str, "float32") )
    return GAL_DATA_TYPE_FLOAT32;

  else if( !strcmp(str, "f64")   || !strcmp(str, "float64") )
    return GAL_DATA_TYPE_FLOAT64;

  else if( !strcmp(str, "c32")   || !strcmp(str, "complex32") )
    return GAL_DATA_TYPE_COMPLEX32;

  else if( !strcmp(str, "c64")   || !strcmp(str, "complex64") )
    return GAL_DATA_TYPE_COMPLEX64;

  else if( !strcmp(str, "str")   || !strcmp(str, "string") )
    return GAL_DATA_TYPE_STRING;

  else
    return GAL_DATA_TYPE_INVALID;

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
gal_data_type_min(uint8_t type, void *in)
{
  switch(type)
    {
    case GAL_DATA_TYPE_UINT8:    *(uint8_t *)  in = 0;            break;
    case GAL_DATA_TYPE_INT8:     *(int8_t *)   in = INT8_MIN;     break;
    case GAL_DATA_TYPE_UINT16:   *(uint16_t *) in = 0;            break;
    case GAL_DATA_TYPE_INT16:    *(int16_t *)  in = INT16_MIN;    break;
    case GAL_DATA_TYPE_UINT32:   *(uint32_t *) in = 0;            break;
    case GAL_DATA_TYPE_INT32:    *(int32_t *)  in = INT32_MIN;    break;
    case GAL_DATA_TYPE_UINT64:   *(uint64_t *) in = 0;            break;
    case GAL_DATA_TYPE_INT64:    *(int64_t *)  in = INT64_MIN;    break;
    case GAL_DATA_TYPE_FLOAT32:  *(float *)    in = -FLT_MAX;     break;
    case GAL_DATA_TYPE_FLOAT64:  *(double *)   in = -DBL_MAX;     break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_data_type_min'", type);
    }
}





void
gal_data_type_max(uint8_t type, void *in)
{
  switch(type)
    {
    case GAL_DATA_TYPE_UINT8:    *(uint8_t *)  in = UINT8_MAX;    break;
    case GAL_DATA_TYPE_INT8:     *(int8_t *)   in = INT8_MAX;     break;
    case GAL_DATA_TYPE_UINT16:   *(uint16_t *) in = UINT16_MAX;   break;
    case GAL_DATA_TYPE_INT16:    *(int16_t *)  in = INT16_MAX;    break;
    case GAL_DATA_TYPE_UINT32:   *(uint32_t *) in = UINT32_MAX;   break;
    case GAL_DATA_TYPE_INT32:    *(int32_t *)  in = INT32_MAX;    break;
    case GAL_DATA_TYPE_UINT64:   *(uint64_t *) in = UINT64_MAX;   break;
    case GAL_DATA_TYPE_INT64:    *(int64_t *)  in = INT64_MAX;    break;
    case GAL_DATA_TYPE_FLOAT32:  *(float *)    in = FLT_MAX;      break;
    case GAL_DATA_TYPE_FLOAT64:  *(double *)   in = DBL_MAX;      break;
    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_data_type_min'", type);
    }
}





/* Since linked lists need a different process than arrays, for functions
   that work on both, it is convenient to simiplify the check with this
   function. */
int
gal_data_is_linked_list(uint8_t type)
{
  return type==GAL_DATA_TYPE_STRLL;
}




















/*********************************************************************/
/*************          Size and allocation        *******************/
/*********************************************************************/
int
gal_data_dsize_is_different(gal_data_t *first, gal_data_t *second)
{
  size_t i;

  /* First make sure that the dimensionality is the same. */
  if(first->ndim!=second->ndim)
    return 1;

  /* Check if the sizes along all dimensions are the same: */
  for(i=0;i<first->ndim;++i)
    if( first->dsize[i] != second->dsize[i] )
      return 1;

  /* If it got to here, we know the dimensions have the same length. */
  return 0;
}





size_t
gal_data_sizeof(uint8_t type)
{
  /* Allocate space for the array to keep the image. */
  switch(type)
    {
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support bit "
            "types, please get in touch with us to implement it.");

      /* The parenthesis after sizeof is not a function, it is actually a
         type cast, so we have put a space between size of and the
         parenthesis to highlight this. In C, `sizeof' is an operator, not
         a function.*/
    case GAL_DATA_TYPE_UINT8:
      return sizeof (uint8_t);

    case GAL_DATA_TYPE_INT8:
      return sizeof (int8_t);

    case GAL_DATA_TYPE_UINT16:
      return sizeof (uint16_t);

    case GAL_DATA_TYPE_INT16:
      return sizeof (int16_t);

    case GAL_DATA_TYPE_UINT32:
      return sizeof (uint32_t);

    case GAL_DATA_TYPE_INT32:
      return sizeof (int32_t);

    case GAL_DATA_TYPE_UINT64:
      return sizeof (uint64_t);

    case GAL_DATA_TYPE_INT64:
      return sizeof (int64_t);

    case GAL_DATA_TYPE_FLOAT32:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0, "`float` is not 32 bits on this machine");
      return sizeof (float);

    case GAL_DATA_TYPE_FLOAT64:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0, "`double` is not 64 bits on this machine");
      return sizeof (double);

    case GAL_DATA_TYPE_COMPLEX32:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0, "`float` is not 32 bits on this machine");
      return sizeof (gsl_complex_float);

    case GAL_DATA_TYPE_COMPLEX64:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0, "`double` is not 64 bits on this machine");
      return sizeof (gsl_complex);

    case GAL_DATA_TYPE_STRING:
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





/* Copy the WCS structure from the input to the output structure. */
void
gal_data_copy_wcs(gal_data_t *in, gal_data_t *out)
{
  if(in->wcs)
    {
      /* Allocate the output WCS structure. */
      errno=0;
      out->wcs=malloc(sizeof *out->wcs);
      if(out->wcs==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for out->wcs in "
              "gal_data_copy_wcs", sizeof *out->wcs);

      /* Initialize the allocated WCS structure. The WCSLIB manual says "On
         the first invokation, and only the first invokation, wcsprm::flag
         must be set to -1 to initialize memory management"*/
      out->wcs->flag=-1;
      wcsini(1, out->ndim, out->wcs);

      /* Copy the input WCS to the output WSC structure. */
      wcscopy(1, in->wcs, out->wcs);
    }
  else
    out->wcs=NULL;
}





/* Allocate an array based on the value of type. Note that the argument
   `size' is the number of elements, necessary in the array, the number of
   bytes each element needs will be determined internaly by this function
   using the datatype argument, so you don't have to worry about it. */
void *
gal_data_malloc_array(uint8_t type, size_t size)
{
  void *array;

  errno=0;
  array=malloc( size * gal_data_sizeof(type) );
  if(array==NULL)
    error(EXIT_FAILURE, errno, "array of %zu bytes in gal_data_malloc_array",
          size * gal_data_sizeof(type));

  return array;
}





void *
gal_data_calloc_array(uint8_t type, size_t size)
{
  void *array;

  errno=0;
  array=calloc( size,  gal_data_sizeof(type) );
  if(array==NULL)
    error(EXIT_FAILURE, errno, "array of %zu bytes in gal_data_calloc_array",
          size * gal_data_sizeof(type));

  return array;
}





/* Allocate space for one blank value of the given type and put the value
   in it. */
void *
gal_data_alloc_number(uint8_t type, void *number)
{
  void *allocated;

  /* Allocate the space for the blank value: */
  allocated=gal_data_malloc_array(type, 1);

  /* Put the blank value into it. */
  errno=0;
  switch(type)
    {
    case GAL_DATA_TYPE_BIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support blank "
            "values for `GAL_DATA_TYPE_BIT', please get in touch with "
            "us to see how we can implement it.");

    case GAL_DATA_TYPE_UINT8:
      *(uint8_t *)allocated=*(uint8_t *)number;
      break;

    case GAL_DATA_TYPE_INT8:
      *(int8_t *)allocated=*(int8_t *)number;
      break;

    case GAL_DATA_TYPE_UINT16:
      *(uint16_t *)allocated=*(uint16_t *)number;
      break;

    case GAL_DATA_TYPE_INT16:
      *(int16_t *)allocated=*(int16_t *)number;
      break;

    case GAL_DATA_TYPE_UINT32:
      *(uint32_t *)allocated=*(uint32_t *)number;
      break;

    case GAL_DATA_TYPE_INT32:
      *(int32_t *)allocated=*(int32_t *)number;
      break;

    case GAL_DATA_TYPE_UINT64:
      *(uint64_t *)allocated=*(uint64_t *)number;
      break;

    case GAL_DATA_TYPE_INT64:
      *(int64_t *)allocated=*(int64_t *)number;
      break;

    case GAL_DATA_TYPE_FLOAT32:
      *(float *)allocated=*(float *)number;
      break;

    case GAL_DATA_TYPE_FLOAT64:
      *(double *)allocated=*(double *)number;
      break;

    case GAL_DATA_TYPE_COMPLEX32:
      GSL_COMPLEX_P_REAL(((gsl_complex_float *)allocated)) =
        GSL_COMPLEX_P_REAL(((gsl_complex_float *)number));
      GSL_COMPLEX_P_IMAG(((gsl_complex_float *)allocated)) =
        GSL_COMPLEX_P_IMAG(((gsl_complex_float *)number));
      break;

    case GAL_DATA_TYPE_COMPLEX64:
      GSL_COMPLEX_P_REAL(((gsl_complex *)allocated)) =
        GSL_COMPLEX_P_REAL(((gsl_complex *)number));
      GSL_COMPLEX_P_IMAG(((gsl_complex *)allocated)) =
        GSL_COMPLEX_P_IMAG(((gsl_complex *)number));
      break;

    default:
      error(EXIT_FAILURE, 0, "type value of %d not recognized in "
            "`gal_data_alloc_number'", type);
    }

  return allocated;
}





static void
gal_data_mmap(gal_data_t *data, int clear)
{
  int filedes;
  char *filename;
  unsigned char uc=0;
  size_t bsize=data->size*gal_data_sizeof(data->type);

  /* Check if the .gnuastro folder exists, write the file there. If it
     doesn't exist, then make the .gnuastro directory.*/
  gal_checkset_mkdir(".gnuastro");

  /* Set the filename */
  gal_checkset_allocate_copy("./.gnuastro/mmap_XXXXXX", &filename);

  /* Create a zero-sized file and keep its descriptor.  */
  errno=0;
  /*filedes=open(filename, O_RDWR | O_CREAT | O_EXCL | O_TRUNC );*/
  filedes=mkstemp(filename);
  if(filedes==-1)
    error(EXIT_FAILURE, errno, "%s couldn't be created", filename);


  /* Make enough space to keep the array data. */
  errno=0;
  if( lseek(filedes, bsize, SEEK_SET) == -1 )
    error(EXIT_FAILURE, errno, "%s: unable to change file position by "
          "%zu bytes", filename, bsize);


  /* Write to the newly set file position so the space is allocated. */
  if( write(filedes, &uc, bsize) == -1)
    error(EXIT_FAILURE, errno, "%s: unable to write one byte at the "
          "%zu-th position", filename, bsize);


  /* Map the memory. */
  data->array=mmap(NULL, bsize, PROT_READ | PROT_WRITE, MAP_SHARED,
                   filedes, 0);

  /* Close the file. */
  if( close(filedes) == -1 )
    error(EXIT_FAILURE, errno, "%s couldn't be closed", filename);

  /* Keep the filename. */
  data->mmapname=filename;

  /* If it was supposed to be cleared, then clear the memory. */
  if(clear) memset(data->array, 0, bsize);
}





/* Initialize the data structure.

   Some notes:

   - The `status' value is the only element that cannot be set by this
     function, it is initialized to zero.

   - If no `array' is given, a blank array of the given size will be
     allocated. If it is given the array pointer will be directly put here,
     so do not free it independently any more. If you want a separate copy
     of a dataset, you should use `gal_data_copy', not this function.

   - Space for the `name', `unit', and `comment' strings within the data
     structure are allocated here. So you can safely use literal strings,
     or statically allocated ones, or simply the strings from other data
     structures (and not have to worry about which one to free later).
*/
void
gal_data_initialize(gal_data_t *data, void *array, uint8_t type,
                    size_t ndim, size_t *dsize, struct wcsprm *wcs,
                    int clear, size_t minmapsize, char *name,
                    char *unit, char *comment)
{
  size_t i;
  gal_data_t in;

  /* Do the simple copying cases. For the display elements, set them all to
     impossible (negative) values so if not explicitly set by later steps,
     the default values are used if/when printing.*/
  data->wcs=NULL;
  data->status=0;
  data->next=NULL;
  data->ndim=ndim;
  data->type=type;
  data->mmapname=NULL;
  data->minmapsize=minmapsize;
  gal_checkset_allocate_copy(name, &data->name);
  gal_checkset_allocate_copy(unit, &data->unit);
  gal_checkset_allocate_copy(comment, &data->comment);
  data->disp_fmt=data->disp_width=data->disp_precision=-1;


  /* Copy the WCS structure. Note that the `in' data structure was just
     defined to keep this pointer to call `gal_data_copy_wcs'. */
  in.wcs=wcs;
  gal_data_copy_wcs(&in, data);


  /* Allocate space for the dsize array, only if the data are to have any
     dimensions. Note that in our convention, a number has a `ndim=1' and
     `dsize[0]=1', A 1D array also has `ndim=1', but `dsize[0]>1'. */
  if(ndim)
    {
      /* Allocate dsize. */
      errno=0;
      data->dsize=malloc(ndim*sizeof *data->dsize);
      if(data->dsize==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for data->dsize in "
              "`gal_data_alloc'", ndim*sizeof *data->dsize);


      /* Fill in the `dsize' array and in the meantime set `size': */
      data->size=1;
      for(i=0;i<ndim;++i)
        {
          /* Do a small sanity check. */
          if(dsize[i]<=0)
            error(EXIT_FAILURE, 0, "the size of a dimension cannot be zero "
                  "or negative. dsize[%zu] in `gal_data_alloc' has a value "
                  "of %ld", i, dsize[i]);

          /* Write this dimension's size, also correct the total number of
             elements. */
          data->size *= ( data->dsize[i] = dsize[i] );
        }

      /* Set the array pointer. If an non-NULL array pointer was given,
         then use it. */
      if(array)
        data->array=array;
      else
        {
          if( gal_data_sizeof(type)*data->size  > minmapsize )
            gal_data_mmap(data, clear);
          else
            {
              /* Allocate the space for the array. */
              if(clear)
                data->array = gal_data_calloc_array(data->type, data->size);
              else
                data->array = gal_data_malloc_array(data->type, data->size);
            }
        }
    }
  else
    {
      data->size=0;
      data->array=NULL;
      data->dsize=NULL;
    }
}





/* Allocate a data structure based on the given parameters. If you want to
   force the array into the hdd/ssd (mmap it), then set minmapsize=-1
   (largest possible size_t value), in this way, no file will be larger. */
gal_data_t *
gal_data_alloc(void *array, uint8_t type, size_t ndim, size_t *dsize,
               struct wcsprm *wcs, int clear, size_t minmapsize,
               char *name, char *unit, char *comment)
{
  gal_data_t *out;

  /* Allocate the space for the actual structure. */
  errno=0;
  out=malloc(sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for gal_data_t in gal_data_alloc",
          sizeof *out);

  /* Initialize the allocated array. */
  gal_data_initialize(out, array, type, ndim, dsize, wcs, clear, minmapsize,
                      name, unit, comment);

  /* Return the final structure. */
  return out;
}





/* Allocate an array of data structures and initialize all the values. */
gal_data_t *
gal_data_calloc_dataarray(size_t size)
{
  size_t i;
  gal_data_t *out;

  /* Allocate the array to keep the structures. */
  errno=0;
  out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `out' in "
          "`gal_data_calloc_dataarray'", size*sizeof *out);


  /* Set the pointers to NULL if they didn't exist and the non-pointers to
     impossible integers (so the caller knows the array is only
     allocated. `minmapsize' should be set when allocating the array and
     should be set when you run `gal_data_initialize'. */
  for(i=0;i<size;++i)
    {
      out[i].array      = NULL;
      out[i].type       = GAL_DATA_TYPE_INVALID;
      out[i].ndim       = 0;
      out[i].dsize      = NULL;
      out[i].nwcs       = 0;
      out[i].wcs        = NULL;
      out[i].mmapname   = NULL;
      out[i].name = out[i].unit = out[i].comment = NULL;
      out[i].disp_fmt = out[i].disp_width = out[i].disp_precision = -1;
    }

  /* Return the array pointer. */
  return out;
}





/* In some contexts, it is necessary for all the strings to have the same
   allocated space (when the `strlen' is different). This function will
   allocate new copies for all elements to have the same length as the
   maximum length and set all trailing elements to `\0' for those that are
   shorter than the length. The return value is the allocated space. If the
   dataset is not a string, the returned value will be -1 (largest number
   of `size_t'). */
size_t
gal_data_string_fixed_alloc_size(gal_data_t *data)
{
  size_t i, j, maxlen=0;
  char *tmp, **strarr=data->array;

  /* Return 0 if the dataset is not a string. */
  if(data->type!=GAL_DATA_TYPE_STRING)
    return -1;

  /* Get the maximum length. */
  for(i=0;i<data->size;++i)
    maxlen = strlen(strarr[i])>maxlen ? strlen(strarr[i]) : maxlen;

  /* For all elements, check the length and if they aren't equal to maxlen,
     then allocate a maxlen sized array and put the values in. */
  for(i=0;i<data->size;++i)
    {
      /* Allocate (and clear) the space for the new string. We want it to
         be cleared, so when the strings are smaller, the rest of the space
         is filled with '\0' (ASCII for 0) values.*/
      errno=0;
      tmp=calloc(maxlen+1, sizeof *strarr[i]);
      if(tmp==NULL)
        error(EXIT_FAILURE, 0, "%zu bytes for tmp in "
              "`gal_data_fixed_alloc_size_for_string'",
              maxlen+1*sizeof *strarr[i]);

      /* Put the old array into the newly allocated space. `tmp' was
         cleared (all values set to `\0', so we don't need to set the final
         one explicity after the copy.*/
      for(j=0;strarr[i][j]!='\0';++j)
        tmp[j]=strarr[i][j];

      /* Free the old array and put in the new one. */
      free(strarr[i]);
      strarr[i]=tmp;
    }

  /* Return the allocated space. */
  return maxlen+1;
}





/* Free the allocated contents of a data structure, not the structure
   itsself. The reason for this function begin separate from
   `gal_data_free) is that the data structure might be allocated as an
   array (statically like `gal_data_t da[20]', or dynamically like
   `gal_data_t *da; da=malloc(20*sizeof *da);'). In both cases, a loop will
   be necessary to delete the allocated contents of each element of the
   data structure array, but not the structure its self. After that loop,
   if the array of data structures was statically allocated, you don't have
   to do anything. If it was dynamically allocated, we just have to run
   `free(da)'.*/
void
gal_data_free_contents(gal_data_t *data)
{
  if(data==NULL)
    error(EXIT_FAILURE, 0, "the input data structure to "
          "`gal_data_free_contents' was a NULL pointer");

  /* Free all the possible allocations. */
  if(data->name)    free(data->name);
  if(data->unit)    free(data->unit);
  if(data->dsize)   free(data->dsize);
  if(data->wcs)     wcsfree(data->wcs);
  if(data->comment) free(data->comment);

  /* If the data type is string, then each element in the array is actually
     a pointer to the array of characters, so free them before freeing the
     actual array. */
  if(data->type==GAL_DATA_TYPE_STRING && data->array)
    {
      size_t i;
      char **strarr=data->array;
      for(i=0;i<data->size;++i) free(strarr[i]);
    }

  /* Free the array. */
  if(data->mmapname)
    {
      /* Delete the file keeping the array. */
      remove(data->mmapname);

      /* Free the file name space. */
      free(data->mmapname);
    }
  else
    if(data->array) free(data->array);
}





/* Free the contents of the data structure and the data structure
   itsself. */
void
gal_data_free(gal_data_t *data)
{
  if(data)
    {
      gal_data_free_contents(data);
      free(data);
    }
}




















/*********************************************************************/
/*************    Data structure as a linked list   ******************/
/*********************************************************************/
/* Add a new data structure to the top of an existing linked list of data
   structures. Note that if the new node is its self a list, all its nodes
   will be added to the list. */
void
gal_data_add_existing_to_ll(gal_data_t **list, gal_data_t *newnode)
{
  gal_data_t *tmp=newnode, *toadd;

  /* Check if newnode is itself a list or not. */
  if(newnode->next)
    {
      /* Go onto the last node in newnode's existing list. */
      while(tmp->next) tmp=tmp->next;

      /* Set the last node as the node to add to the list. */
      toadd=tmp;
    }
  else
    /* Its not a list, so just set it to `toadd'. */
    toadd=newnode;


  /* Set the next element of toadd and update what list points to.*/
  toadd->next=*list;
  *list=toadd;
}





void
gal_data_add_to_ll(gal_data_t **list, void *array, uint8_t type, size_t ndim,
                   size_t *dsize, struct wcsprm *wcs, int clear,
                   size_t minmapsize, char *name, char *unit, char *comment)
{
  gal_data_t *newnode;

  /* Put all the input information into a new data structure node. */
  newnode=gal_data_alloc(array, type, ndim, dsize, wcs, clear,
                         minmapsize, name, unit, comment);

  /* Add the new node to the list. */
  gal_data_add_existing_to_ll(list, newnode);
}





gal_data_t *
gal_data_pop_from_ll(gal_data_t **list)
{
  struct gal_data_t *out;

  /* Keep the top pointer. */
  out=*list;

  /* Move the list pointer to the next node. */
  *list=out->next;

  /* Set the next poitner of the out pointer to NULL so it isn't
     interpretted as a list any more. */
  out->next=NULL;
  return out;
}





void
gal_data_reverse_ll(gal_data_t **list)
{
  gal_data_t *popped, *in=*list, *reversed=NULL;

  /* Only do the job if the list is not NULL and has more than one node. */
  if( in && in->next )
    {
      while(in!=NULL)
        {
          popped=gal_data_pop_from_ll(&in);
          gal_data_add_existing_to_ll(&reversed, popped);
        }
      *list=reversed;
    }
}





size_t
gal_data_num_in_ll(gal_data_t *list)
{
  size_t num=0;
  while(list!=NULL)
    {
      ++num;
      list=list->next;
    }
  return num;
}





gal_data_t **
gal_data_ll_to_array_of_ptrs(gal_data_t *list, size_t *num)
{
  size_t i=0;
  gal_data_t **out;

  /* Find the number of nodes in the list. */
  *num=gal_data_num_in_ll(list);

  /* Allocate space for the array. */
  errno=0;
  out=malloc(*num * sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, 0, "%zu bytes for the output pointer in "
          "`gal_data_ll_to_array_of_ptrs'", *num * sizeof *out);

  /* Fill in the array with pointers to each data-structure. Note that we
     don't need the list pointer any more, so we can just increment it.*/
  while(list!=NULL) { out[i++]=list; list=list->next; }

  /* Return the allocated array. */
  return out;
}





void
gal_data_free_ll(gal_data_t *list)
{
  struct gal_data_t *tmp;
  while(list!=NULL)
    {
      tmp=list->next;
      gal_data_free(list);
      list=tmp;
    }
}




















/*************************************************************
 **************            Copying             ***************
 *************************************************************/

/* Only to be used in `data_copy_from_string'. */
static void
data_copy_to_string_not_parsed(char *string, void *to, uint8_t type)
{
  if( strcmp(string, GAL_BLANK_STRING) )
    gal_blank_write(to, type);
  else
    error(EXIT_FAILURE, 0, "`%s' couldn't be parsed as `%s' type",
          string, gal_data_type_as_string(type, 1));
}





/* The `from' array is an array of strings. We want to keep it as
   numbers. Note that the case where both input and output structures are
   string was */
static void
data_copy_from_string(gal_data_t *from, gal_data_t *to)
{
  size_t i;
  void *ptr;
  char **strarr=from->array, **outstrarr=to->array;

  /* Sanity check. */
  if(from->type!=GAL_DATA_TYPE_STRING)
    error(EXIT_FAILURE, 0, "`from' in `data_copy_from_string' must have "
          "a string type.");

  /* Do the copying. */
  for(i=0;i<from->size;++i)
    {
      /* Set the pointer. */
      switch(to->type)
        {
        case GAL_DATA_TYPE_UINT8:    ptr = (uint8_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_INT8:     ptr = (int8_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_UINT16:   ptr = (uint16_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_INT16:    ptr = (int16_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_UINT32:   ptr = (uint32_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_INT32:    ptr = (int32_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_UINT64:   ptr = (uint64_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_INT64:    ptr = (int64_t *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_FLOAT32:  ptr = (float *)(to->array) + i;
          break;
        case GAL_DATA_TYPE_FLOAT64:  ptr = (double *)(to->array) + i;
          break;

        case GAL_DATA_TYPE_BIT:
        case GAL_DATA_TYPE_STRLL:
        case GAL_DATA_TYPE_COMPLEX32:
        case GAL_DATA_TYPE_COMPLEX64:
          error(EXIT_FAILURE, 0, "`data_copy_from_string' currently doesn't "
                "support copying to %s type",
                gal_data_type_as_string(to->type, 1));
          break;

        default:
          error(EXIT_FAILURE, 0, "type %d not recognized for to->type in "
                "`data_copy_from_string'", to->type);
        }

      /* Read/put the input into the output. */
      if(to->type==GAL_DATA_TYPE_STRING)
        gal_checkset_allocate_copy(strarr[i], &outstrarr[i]);
      else
        {
          if( gal_data_string_to_type(&ptr, strarr[i], to->type) )
            data_copy_to_string_not_parsed(strarr[i], ptr, to->type);
        }
    }
}





/* Macros for copying to a string. */
#define COPY_TO_STR_INT(CTYPE, BLANK, FMT) {                            \
    CTYPE *a=from->array;                                               \
    for(i=0;i<from->size;++i)                                           \
      {                                                                 \
        if(a[i]!=BLANK) asprintf(&strarr[i], FMT, a[i]);                \
        else                                                            \
          gal_checkset_allocate_copy(GAL_BLANK_STRING, &strarr[i]); \
      }                                                                 \
  }

#define COPY_TO_STR_FLT(CTYPE, BLANK) {                                 \
    CTYPE *a=from->array;                                               \
    for(i=0;i<from->size;++i)                                           \
      {                                                                 \
        if(isnan(BLANK)) isblank = isnan(a[i]) ? 1 : 0;                 \
        else             isblank = a[i]==BLANK ? 1 : 0;                 \
        if(isblank==0) asprintf(&strarr[i], "%f", a[i]);                \
        else gal_checkset_allocate_copy(GAL_BLANK_STRING, &strarr[i]); \
      }                                                                 \
  }

/* Convert any given type into a string by printing it into the elements of
   the already allocated `to->array'. */
static void
data_copy_to_string(gal_data_t *from, gal_data_t *to)
{
  size_t i;
  int isblank;
  char **strarr=to->array, **instrarr=from->array;

  /* Sanity check */
  if(to->type!=GAL_DATA_TYPE_STRING)
    error(EXIT_FAILURE, 0, "`to' in `data_copy_to_string' must have a "
          "string type");

  /* Do the copying */
  switch(from->type)
    {
    case GAL_DATA_TYPE_UINT8:
      COPY_TO_STR_INT(uint8_t,  GAL_BLANK_UINT8, "%u");    break;

    case GAL_DATA_TYPE_INT8:
      COPY_TO_STR_INT(int8_t,   GAL_BLANK_INT8, "%d");     break;

    case GAL_DATA_TYPE_UINT16:
      COPY_TO_STR_INT(uint16_t, GAL_BLANK_UINT16, "%u");   break;

    case GAL_DATA_TYPE_INT16:
      COPY_TO_STR_INT(int16_t,  GAL_BLANK_INT16, "%d");    break;

    case GAL_DATA_TYPE_UINT32:
      COPY_TO_STR_INT(uint32_t, GAL_BLANK_UINT32, "%u");   break;

    case GAL_DATA_TYPE_INT32:
      COPY_TO_STR_INT(int32_t,  GAL_BLANK_INT32, "%d");    break;

    case GAL_DATA_TYPE_UINT64:
      COPY_TO_STR_INT(uint64_t, GAL_BLANK_UINT64, "%lu");  break;

    case GAL_DATA_TYPE_INT64:
      COPY_TO_STR_INT(int64_t,  GAL_BLANK_INT64, "%ld");   break;

    case GAL_DATA_TYPE_FLOAT32:
      COPY_TO_STR_FLT(float, GAL_BLANK_FLOAT32);           break;

    case GAL_DATA_TYPE_FLOAT64:
      COPY_TO_STR_FLT(double, GAL_BLANK_FLOAT32);          break;

    case GAL_DATA_TYPE_STRING:
      for(i=0;i<from->size;++i)
        gal_checkset_allocate_copy(instrarr[i], &strarr[i]);
      break;

    case GAL_DATA_TYPE_BIT:
    case GAL_DATA_TYPE_STRLL:
    case GAL_DATA_TYPE_COMPLEX32:
    case GAL_DATA_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "`data_copy_to_string' currently doesn't "
            "support copying to %s type",
            gal_data_type_as_string(from->type, 1));
      break;

    default:
      error(EXIT_FAILURE, 0, "type %d not recognized for `from->type' in "
            "`data_copy_to_string'", from->type);
    }
}





/* Copy to a new type for integers. */
#define COPY_OTYPE_ITYPE_SET_INT(otype, itype) {                        \
    itype *ia=in->array, iblank;                                        \
    otype *oa=out->array, *of=oa+out->size, oblank;                     \
                                                                        \
    /* Check if there are blank values in the input array and that */   \
    /* the types of the two structures are different. */                \
    if( in->type!=newtype && gal_blank_present(in) )                    \
      {                                                                 \
        /* Set the blank values */                                      \
        gal_blank_write(&iblank, in->type);                             \
        gal_blank_write(&oblank, newtype);                              \
                                                                        \
        /* Copy the input to the output. */                             \
        do { *oa = *ia==iblank ? oblank : *ia; ia++; } while(++oa<of);  \
      }                                                                 \
                                                                        \
    /* There were no blank elements in the input, or the input and */   \
    /* output have the same type. */                                    \
    else                                                                \
      do *oa=*ia++; while(++oa<of);                                     \
  }





/* Copy to a new type for floating point values. */
#define COPY_OTYPE_ITYPE_SET_FLT(otype, itype) {                        \
    itype *ia=in->array, iblank;                                        \
    otype *oa=out->array, *of=oa+out->size, oblank;                     \
                                                                        \
    /* Check if there are blank values in the input array and that */   \
    /* the types of the two structures are different. */                \
    if( in->type!=newtype && gal_blank_present(in) )                    \
      {                                                                 \
        /* Set the blank values */                                      \
        gal_blank_write(&iblank, in->type);                             \
        gal_blank_write(&oblank, newtype);                              \
                                                                        \
        /* When the blank value isn't NaN, then we should use the */    \
        /* equal operator to check for blank values. */                 \
        if( isnan(iblank) )                                             \
          {                                                             \
            do { *oa = isnan(*ia) ? oblank : *ia; ia++; }               \
            while(++oa<of);                                             \
          }                                                             \
        else                                                            \
          {                                                             \
            do { *oa = *ia==iblank ? oblank : *ia; ia++; }              \
            while(++oa<of);                                             \
          }                                                             \
      }                                                                 \
                                                                        \
    /* There were no blank elements in the input, or the input and */   \
    /* output have the same type. */                                    \
    else                                                                \
      do *oa=*ia++; while(++oa<of);                                     \
  }





/* gal_data_copy_to_new_type: Output type is set, now choose the input
   type. */
#define COPY_OTYPE_SET(otype)                                           \
  switch(in->type)                                                      \
    {                                                                   \
    case GAL_DATA_TYPE_UINT8:                                           \
      COPY_OTYPE_ITYPE_SET_INT(otype, uint8_t);                         \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_INT8:                                            \
      COPY_OTYPE_ITYPE_SET_INT(otype, int8_t);                          \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_UINT16:                                          \
      COPY_OTYPE_ITYPE_SET_INT(otype, uint16_t);                        \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_INT16:                                           \
      COPY_OTYPE_ITYPE_SET_INT(otype, int16_t);                         \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_UINT32:                                          \
      COPY_OTYPE_ITYPE_SET_INT(otype, uint32_t);                        \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_INT32:                                           \
      COPY_OTYPE_ITYPE_SET_INT(otype, int32_t);                         \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_UINT64:                                          \
      COPY_OTYPE_ITYPE_SET_INT(otype, uint64_t);                        \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_INT64:                                           \
      COPY_OTYPE_ITYPE_SET_INT(otype, int64_t);                         \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_FLOAT32:                                         \
      COPY_OTYPE_ITYPE_SET_FLT(otype, float);                           \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_FLOAT64:                                         \
      COPY_OTYPE_ITYPE_SET_FLT(otype, double);                          \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_STRING:                                          \
      data_copy_from_string(in, out);                                   \
      break;                                                            \
                                                                        \
    case GAL_DATA_TYPE_BIT:                                             \
    case GAL_DATA_TYPE_STRLL:                                           \
    case GAL_DATA_TYPE_COMPLEX32:                                       \
    case GAL_DATA_TYPE_COMPLEX64:                                       \
      error(EXIT_FAILURE, 0, "`gal_data_copy_to_new_type' currently "   \
            "doesn't support copying from %s type to a numeric (real) " \
            "type", gal_data_type_as_string(in->type, 1));              \
      break;                                                            \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "type code %d not recognized for "         \
            "`in->type' in COPY_OTYPE_SET", in->type);                  \
    }





/* Copy a given data structure to a new one (possibly with a new type). */
gal_data_t *
gal_data_copy_to_new_type(gal_data_t *in, uint8_t newtype)
{
  gal_data_t *out;

  /* Allocate space for the output type */
  out=gal_data_alloc(NULL, newtype, in->ndim, in->dsize, in->wcs,
                     0, in->minmapsize, in->name, in->unit, in->comment);

  /* For debugging.
  printf("in: %d (%s)\nout: %d (%s)\n\n", in->type,
         gal_data_type_as_string(in->type, 1), out->type,
         gal_data_type_as_string(out->type, 1));
  */

  /* Fill in the output array: */
  switch(newtype)
    {
    case GAL_DATA_TYPE_UINT8:   COPY_OTYPE_SET(uint8_t);         break;
    case GAL_DATA_TYPE_INT8:    COPY_OTYPE_SET(int8_t);          break;
    case GAL_DATA_TYPE_UINT16:  COPY_OTYPE_SET(uint16_t);        break;
    case GAL_DATA_TYPE_INT16:   COPY_OTYPE_SET(int16_t);         break;
    case GAL_DATA_TYPE_UINT32:  COPY_OTYPE_SET(uint32_t);        break;
    case GAL_DATA_TYPE_INT32:   COPY_OTYPE_SET(int32_t);         break;
    case GAL_DATA_TYPE_UINT64:  COPY_OTYPE_SET(uint64_t);        break;
    case GAL_DATA_TYPE_INT64:   COPY_OTYPE_SET(int64_t);         break;
    case GAL_DATA_TYPE_FLOAT32: COPY_OTYPE_SET(float);           break;
    case GAL_DATA_TYPE_FLOAT64: COPY_OTYPE_SET(double);          break;
    case GAL_DATA_TYPE_STRING:  data_copy_to_string(in, out);    break;

    case GAL_DATA_TYPE_BIT:
    case GAL_DATA_TYPE_STRLL:
    case GAL_DATA_TYPE_COMPLEX32:
    case GAL_DATA_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "`gal_data_copy_to_new_type' currently doesn't "
            "support copying to %s type",
            gal_data_type_as_string(newtype, 1));
      break;

    default:
      error(EXIT_FAILURE, 0, "type %d not recognized for "
            "for newtype in gal_data_copy_to_new_type", newtype);
    }

  /* Return the created array */
  return out;
}





gal_data_t *
gal_data_copy_to_new_type_free(gal_data_t *in, uint8_t type)
{
  gal_data_t *out;

  /* In a general application, it might happen that the type is equal with
     the type of the input. Since the job of this function is to free the
     input data set, so and the user just wants one dataset after this
     function finishes, we can safely just return the input. */
  if(type==in->type)
    return in;
  else
    {
      out=gal_data_copy_to_new_type(in, type);
      gal_data_free(in);
      return out;
    }
}





gal_data_t *
gal_data_copy(gal_data_t *in)
{
  return gal_data_copy_to_new_type(in, in->type);
}





int
gal_data_out_type(gal_data_t *first, gal_data_t *second)
{
  return first->type > second->type ? first->type : second->type;
}





/* The two input `f' and `s' datasets can be any type. But `of' and `os'
   will have type `type', if freeinputs is non-zero, then the input arrays
   will be freed if they needed to be changed to a new type. */
void
gal_data_to_same_type(gal_data_t *f,   gal_data_t *s,
                      gal_data_t **of, gal_data_t **os,
                      uint8_t type, int freeinputs)
{
  /* Change first dataset into the new type if necessary. */
  if( f->type != type )
    {
      *of=gal_data_copy_to_new_type(f, type);
      if(freeinputs)
        gal_data_free(f);
    }
  else
    *of=f;

  /* Change second dataset into the new type if necessary. */
  if( s->type != type )
    {
      *os=gal_data_copy_to_new_type(s, type);
      if(freeinputs)
        gal_data_free(s);
    }
  else
    *os=s;
}




















/*************************************************************
 **************              Write             ***************
 *************************************************************/
#define WRITE_TO_STRING(CTYPE, FMT) asprintf(&str, FMT, *(CTYPE *)ptr);

char *
gal_data_write_to_string(void *ptr, uint8_t type, int quote_if_str_has_space)
{
  char *c, *str;
  switch(type)
    {
    /* For a string we might need to make sure it has no white space
       characters, if it does, it can be printed it within quotation
       signs. */
    case GAL_DATA_TYPE_STRING:
      if(quote_if_str_has_space)
        {
          c=*(char **)ptr; while(*c!='\0') if(isspace(*c++)) break;
          if(*c=='\0') asprintf(&str, "%s",      *(char **)ptr);
          else         asprintf(&str, "\"%s\" ", *(char **)ptr);
        }
      else
        asprintf(&str, "%s", *(char **)ptr);
      break;

    case GAL_DATA_TYPE_UINT8:   WRITE_TO_STRING(uint8_t,   "%u");  break;
    case GAL_DATA_TYPE_INT8:    WRITE_TO_STRING(int8_t,    "%d");  break;
    case GAL_DATA_TYPE_UINT16:  WRITE_TO_STRING(uint16_t,  "%u");  break;
    case GAL_DATA_TYPE_INT16:   WRITE_TO_STRING(int16_t,   "%d");  break;
    case GAL_DATA_TYPE_UINT32:  WRITE_TO_STRING(uint32_t,  "%u");  break;
    case GAL_DATA_TYPE_INT32:   WRITE_TO_STRING(int32_t,   "%d");  break;
    case GAL_DATA_TYPE_UINT64:  WRITE_TO_STRING(uint64_t, "%lu");  break;
    case GAL_DATA_TYPE_INT64:   WRITE_TO_STRING(int64_t,  "%ld");  break;
    case GAL_DATA_TYPE_FLOAT32: WRITE_TO_STRING(float,   "%.6g");  break;
    case GAL_DATA_TYPE_FLOAT64: WRITE_TO_STRING(double, "%.10g");  break;

    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_data_write_to_string'", type);
    }
  return str;
}




















/*************************************************************
 **************              Read              ***************
 *************************************************************/
/* If the data structure was correctly created (the string was a number),
   then return its pointer. Otherwise, return NULL. */
gal_data_t *
gal_data_string_to_number(char *string)
{
  size_t dsize[1]={1};
  int fnz=-1, lnz=0;     /* `F'irst (or `L'ast) `N'on-`Z'ero. */
  void *ptr, *numarr;
  char *tailptr, *cp;
  uint8_t type, forcedfloat=0;

  /* Define initial spaces to keep the value. */
  uint8_t   u8;   int8_t   i8;      uint16_t u16;   int16_t i16;
  uint32_t u32;   int32_t i32;      uint64_t u64;   int64_t i64;
  float      f;   double    d;

  /* First see if the number is a double (the most generic). */
  d=strtod(string, &tailptr);
  if(*tailptr=='f') { if(tailptr[1]=='\0') forcedfloat=1; else return NULL; }
  else if (*tailptr!='\0')  return NULL;

  /* See if the number is actually an integer: */
  if( forcedfloat==0 && ceil(d) == d )
    {
      /* If the number is negative, put it in the signed types (based on
         its value). If its zero or positive, then put it in the unsigned
         types. */
      if( d < 0 )
        {
          if     (d>INT8_MIN)    {i8=d;  ptr=&i8;  type=GAL_DATA_TYPE_INT8;}
          else if(d>INT16_MIN)   {i16=d; ptr=&i16; type=GAL_DATA_TYPE_INT16;}
          else if(d>INT32_MIN)   {i32=d; ptr=&i32; type=GAL_DATA_TYPE_INT32;}
          else                   {i64=d; ptr=&i64; type=GAL_DATA_TYPE_INT64;}
        }
      else
        {
          if     (d<=UINT8_MAX)  {u8=d;  ptr=&u8;  type=GAL_DATA_TYPE_UINT8;}
          else if(d<=UINT16_MAX) {u16=d; ptr=&u16; type=GAL_DATA_TYPE_UINT16;}
          else if(d<=UINT32_MAX) {u32=d; ptr=&u32; type=GAL_DATA_TYPE_UINT32;}
          else                   {u64=d; ptr=&u64; type=GAL_DATA_TYPE_UINT64;}
        }
    }
  else
    {
      /* The maximum number of decimal digits to store in float or double
         precision floating point are:

         float:  23 mantissa bits + 1 hidden bit: log(224)÷log(10) = 7.22
         double: 52 mantissa bits + 1 hidden bit: log(253)÷log(10) = 15.95

         FLT_DIG (at least 6 in ISO C) keeps the number of digits (not zero
         before or after) that can be represented by a single precision
         floating point number. If there are more digits, then we should
         store the value as a double precision.

         Note that the number can have non-digit characters that we don't
         want, like: `.', `e', `E', `,'. */
      for(cp=string;*cp!='\0';++cp)
        if(isdigit(*cp) && *cp!='0' && fnz==-1)
          fnz=cp-string;

      /* In the previous loop, we went to the end of the string, so `cp'
         now points to its `\0'. We just have to iterate backwards! */
      for(;cp!=string;--cp)
        if(isdigit(*cp) && *cp!='0')
          {
            lnz=cp-string;
            break;
          }

      /* Calculate the number of decimal digits and decide if it the number
         should be a float or a double. */
      if( lnz-fnz < FLT_DIG || ( d<FLT_MAX && d>FLT_MIN ) )
        { f=d; ptr=&f; type=GAL_DATA_TYPE_FLOAT32; }
      else
        {      ptr=&d; type=GAL_DATA_TYPE_FLOAT64; }
    }

  /* Return the pointer to the data structure. */
  numarr=gal_data_alloc_number(type, ptr);
  return gal_data_alloc(numarr, type, 1, dsize, NULL, 0, -1,
                        NULL, NULL, NULL);
}





/* Read a string as a given data type and put a the pointer to it in
   *out. When the input `*out!=NULL', then it is assumed to be allocated
   and the value will be simply put there. If `*out==NULL', then space will
   be allocated for the given type and the string's value (in the given
   type) will be stored there.

   Note that when we are dealing with a string type, `*out' should be
   interpretted as `char **' (one element in an array of pointers to
   different strings). In other words, `out' should be `char ***'.

   This function can be used to fill in arrays of numbers from strings (in
   an already allocated data structure), or add nodes to a linked list. For
   an array, you have to pass the pointer to the `i'th element where you
   want the value to be stored, for example &(array[i]).

   If parsing was successful, it will return a 0. If there was a problem,
   it will return 1.  */
int
gal_data_string_to_type(void **out, char *string, uint8_t type)
{
  long l;
  double d;
  void *value;
  char *tailptr;
  int status=0, allocated=0;

  /* If the output is NULL, then allocate the necessary space if we are not
     dealing with a linked list. In a linked list, a NULL value is
     meaningful (it is the end of the list). */
  if( *out==NULL && !gal_data_is_linked_list(type) )
    {
      allocated=1;
      *out=gal_data_malloc_array(type, 1);
    }
  value=*out;

  /* Read the string depending on the type. */
  switch(type)
    {

    /* Linked lists, currently only string linked lists. */
    case GAL_DATA_TYPE_STRLL:
      gal_linkedlist_add_to_stll( (struct gal_linkedlist_stll **)out,
                                  string, 1);
      break;

    /* String, just allocate and copy the string and keep its pointer in
       the place `*out' points to (for strings, `*out' is `char **'). */
    case GAL_DATA_TYPE_STRING:
      gal_checkset_allocate_copy(string, value);
      break;

    /* Floating point: Read it as a double or long, then put it in the
       array. When the conversion can't be done (the string isn't a number
       for example), then just assume no blank value was given. */
    case GAL_DATA_TYPE_FLOAT32:
    case GAL_DATA_TYPE_FLOAT64:
      d=strtod(string, &tailptr);
      if(*tailptr!='\0')
        status=1;
      else
        {
          if(type==GAL_DATA_TYPE_FLOAT32) *(float *) value=d;
          else                            *(double *) value=d;
        }
      break;

    /* Integers. */
    default:
      l=strtol(string, &tailptr, 0);
      if(*tailptr!='\0')
        status=1;
      else
        switch(type)
          {
          /* The signed values can easily be put in. */
          case GAL_DATA_TYPE_INT8:         *(int8_t *)    value = l; break;
          case GAL_DATA_TYPE_INT16:        *(int16_t *)   value = l; break;
          case GAL_DATA_TYPE_INT32:        *(int32_t *)   value = l; break;
          case GAL_DATA_TYPE_INT64:        *(int64_t *)   value = l; break;

          /* For the unsigned types, the value has to be positive, so if
             the input was negative, then just return a status of one and
             don't store the value. */
          default:
            if(l<0)
              status=1;
            else
              switch(type)
                {
                case GAL_DATA_TYPE_UINT8:  *(uint8_t *)   value=l; break;
                case GAL_DATA_TYPE_UINT16: *(uint16_t *)  value=l; break;
                case GAL_DATA_TYPE_UINT32: *(uint32_t *)  value=l; break;
                case GAL_DATA_TYPE_UINT64: *(uint64_t *)  value=l; break;
                default:
                  error(EXIT_FAILURE, 0, "type code %d not recognized in "
                        "`gal_data_string_to_type'", type);
                }
          }
    }

  /* If reading was unsuccessful, then free the space if it was allocated,
     then return the status, don't touch the pointer. */
  if(status && allocated)
    {
      free(*out);
      *out=NULL;
    }
  return status;
}
