/*********************************************************************
data -- Structure and functions to represent/work with data
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <inttypes.h>

#include <gnuastro/wcs.h>
#include <gnuastro/data.h>
#include <gnuastro/tile.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>




















/*********************************************************************/
/*************              Allocation             *******************/
/*********************************************************************/
/* Allocate a data structure based on the given parameters. If you want to
   force the array into the hdd/ssd (mmap it), then set minmapsize=-1
   (largest possible size_t value), in this way, no file will be larger. */
gal_data_t *
gal_data_alloc(void *array, uint8_t type, size_t ndim, size_t *dsize,
               struct wcsprm *wcs, int clear, size_t minmapsize,
               int quietmmap, char *name, char *unit, char *comment)
{
  gal_data_t *out;

  /* Allocate the space for the actual structure. */
  errno=0;
  out=malloc(sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for gal_data_t",
          __func__, sizeof *out);

  /* Initialize the allocated array. */
  gal_data_initialize(out, array, type, ndim, dsize, wcs, clear, minmapsize,
                      quietmmap, name, unit, comment);

  /* Return the final structure. */
  return out;
}





/* Initialize the data structure.

   Some notes:

   - The 'status' value is the only element that cannot be set by this
     function, it is initialized to zero.

   - If no 'array' is given, a blank array of the given size will be
     allocated. If it is given the array pointer will be directly put here,
     so do not free it independently any more. If you want a separate copy
     of a dataset, you should use 'gal_data_copy', not this function.

   - Space for the 'name', 'unit', and 'comment' strings within the data
     structure are allocated here. So you can safely use literal strings,
     or statically allocated ones, or simply the strings from other data
     structures (and not have to worry about which one to free later).
*/
void
gal_data_initialize(gal_data_t *data, void *array, uint8_t type,
                    size_t ndim, size_t *dsize, struct wcsprm *wcs,
                    int clear, size_t minmapsize, int quietmmap,
                    char *name, char *unit, char *comment)
{
  size_t i;
  size_t data_size_limit = (size_t)(-1);

  /* Do the simple copying cases. For the display elements, set them all to
     impossible (negative) values so if not explicitly set by later steps,
     the default values are used if/when printing.*/
  data->flag       = 0;
  data->status     = 0;
  data->next       = NULL;
  data->ndim       = ndim;
  data->type       = type;
  data->block      = NULL;
  data->mmapname   = NULL;
  data->quietmmap  = quietmmap;
  data->minmapsize = minmapsize;
  gal_checkset_allocate_copy(name, &data->name);
  gal_checkset_allocate_copy(unit, &data->unit);
  gal_checkset_allocate_copy(comment, &data->comment);
  data->disp_fmt=data->disp_width=data->disp_precision=-1;


  /* Copy the WCS structure. */
  data->wcs=gal_wcs_copy(wcs);


  /* Allocate space for the dsize array, only if the data are to have any
     dimensions. Note that in our convention, a number has a 'ndim=1' and
     'dsize[0]=1', A 1D array also has 'ndim=1', but 'dsize[0]>1'. */
  if(ndim)
    {
      /* Allocate dsize. */
      errno=0;
      data->dsize=malloc(ndim*sizeof *data->dsize);
      if(data->dsize==NULL)
        error(EXIT_FAILURE, errno, "%s: %zu bytes for data->dsize",
              __func__, ndim*sizeof *data->dsize);


      /* Fill in the 'dsize' array and in the meantime set 'size': */
      data->size=1;
      for(i=0;i<ndim;++i)
        {
          /* Size along a dimension cannot be negative. */
          if(dsize[i] == 0)
            error(EXIT_FAILURE, 0, "%s: dsize[%zu]==0. The size of a "
                  "dimension cannot be zero", __func__, i);

          /* Check for possible overflow while multiplying. */
          if (dsize[i] >= data_size_limit / data->size)
            error(EXIT_FAILURE, 0, "%s: dimension %zu size is too "
                    "large %zu. Total is out of bounds",
                    __func__, i, dsize[i]);

          /* Print a warning if the size in this dimension is too
             large. May happen when the user (mistakenly) writes a negative
             value in this dimension.. */
          if (dsize[i] >= data_size_limit / 2)
            fprintf(stderr, "%s: WARNING: dsize[%zu] value %zu is probably "
                    "a mistake: it exceeds the limit %zu", __func__, i,
                    dsize[i], data_size_limit / 2);

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
          if(data->size)
            {
              if( gal_type_sizeof(type)*data->size > minmapsize )
                /* Allocate the space into disk (HDD/SSD). */
                data->array=gal_pointer_allocate_mmap(data->type, data->size,
                                                      clear, &data->mmapname,
                                                      quietmmap);
              else
                /* Allocate the space in RAM. */
                data->array = gal_pointer_allocate(data->type, data->size,
                                                   clear, __func__,
                                                   "data->array");
            }
          else data->array=NULL; /* The given size was zero! */
        }
    }
  else
    {
      data->size=0;
      data->array=NULL;
      data->dsize=NULL;
    }
}





/* Free the allocated contents of a data structure, not the structure
   itsself. The reason that this function is separate from 'gal_data_free'
   is that the data structure might be allocated as an array (statically
   like 'gal_data_t da[20]', or dynamically like 'gal_data_t *da;
   da=malloc(20*sizeof *da);'). In both cases, a loop will be necessary to
   delete the allocated contents of each element of the data structure
   array, but not the structure its self. After that loop, if the array of
   data structures was statically allocated, you don't have to do
   anything. If it was dynamically allocated, we just have to run
   'free(da)'.

   Since we aren't freeing the 'gal_data_t' its-self, after the allocated
   space for each pointer is freed, the pointer is set to NULL for safety
   (to avoid possible re-calls).
*/
void
gal_data_free_contents(gal_data_t *data)
{
  size_t i;
  char **strarr;

  if(data==NULL)
    error(EXIT_FAILURE, 0, "%s: the input data structure to "
          "'gal_data_free_contents' was a NULL pointer", __func__);

  /* Free all the possible allocations. */
  if(data->name)    { free(data->name);    data->name    = NULL; }
  if(data->unit)    { free(data->unit);    data->unit    = NULL; }
  if(data->dsize)   { free(data->dsize);   data->dsize   = NULL; }
  if(data->wcs)     { wcsfree(data->wcs);  data->wcs     = NULL; }
  if(data->comment) { free(data->comment); data->comment = NULL; }

  /* If the data type is string, then each element in the array is actually
     a pointer to the array of characters, so free them before freeing the
     actual array. */
  if(data->type==GAL_TYPE_STRING && data->array)
    {
      strarr=data->array;
      for(i=0;i<data->size;++i) if(strarr[i]) free(strarr[i]);
    }

  /* Free the array. */
  if(data->mmapname)
    {
      /* Delete the file keeping the array. */
      remove(data->mmapname);

      /* Inform the user. */
      if(!data->quietmmap)
        error(EXIT_SUCCESS, 0, "%s: deleted", data->mmapname);

      /* Free the file name space. */
      free(data->mmapname);

      /* Set the name pointer to NULL since it has been freed. */
      data->mmapname=NULL;
    }
  else
    if(data->array && data->block==NULL) free(data->array);
  data->array=NULL;
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
/*************        Array of data structures      ******************/
/*********************************************************************/
/* Allocate an array of data structures and initialize all the values. */
gal_data_t *
gal_data_array_calloc(size_t size)
{
  size_t i;
  gal_data_t *out;

  /* Allocate the array to keep the structures. */
  errno=0;
  out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for 'out'", __func__,
          size*sizeof *out);


  /* Set the pointers to NULL if they didn't exist and the non-pointers to
     impossible integers (so the caller knows the array is only
     allocated. 'minmapsize' should be set when allocating the array and
     should be set when you run 'gal_data_initialize'. */
  for(i=0;i<size;++i)
    {
      out[i].array      = NULL;
      out[i].type       = GAL_TYPE_INVALID;
      out[i].ndim       = 0;
      out[i].dsize      = NULL;
      out[i].size       = 0;
      out[i].mmapname   = NULL;
      out[i].minmapsize = -1;
      out[i].quietmmap  = 1;
      out[i].nwcs       = 0;
      out[i].wcs        = NULL;
      out[i].flag       = 0;
      out[i].status     = 0;
      out[i].next       = NULL;
      out[i].block      = NULL;
      out[i].name = out[i].unit = out[i].comment = NULL;
      out[i].disp_fmt = out[i].disp_width = out[i].disp_precision = -1;
    }

  /* Return the array pointer. */
  return out;
}





/* When you have an array of data structures. */
void
gal_data_array_free(gal_data_t *dataarr, size_t size, int free_array)
{
  size_t i;

  /* If its NULL, don't do anything. */
  if(dataarr==NULL) return;

  /* First free all the contents. */
  for(i=0;i<size;++i)
    {
      /* See if the array should be freed or not. */
      if(free_array==0)
        dataarr[i].array=NULL;

      /* Now clear the contents of the dataset. */
      gal_data_free_contents(&dataarr[i]);
    }

  /* Now you can free the whole array. */
  free(dataarr);
}




















/*************************************************************
 **************            Copying             ***************
 *************************************************************/
/* Only to be used in 'data_copy_from_string'. */
static void
data_copy_to_string_not_parsed(char *string, void *to, uint8_t type)
{
  if( strcmp(string, GAL_BLANK_STRING) )
    gal_blank_write(to, type);
  else
    error(EXIT_FAILURE, 0, "%s: '%s' couldn't be parsed as '%s' type",
          __func__, string, gal_type_name(type, 1));
}





/* The 'from' array is an array of strings. We want to keep it as
   numbers. Note that the case where both input and output structures are
   string was */
static void
data_copy_from_string(gal_data_t *from, gal_data_t *to)
{
  size_t i;
  void *ptr;
  char **strarr=from->array, **outstrarr=to->array;

  /* Sanity check. */
  if(from->type!=GAL_TYPE_STRING)
    error(EXIT_FAILURE, 0, "%s: 'from' must have a string type.", __func__);
  if(from->block)
    error(EXIT_FAILURE, 0, "%s: tiles not currently supported ('block' "
          "element must be NULL). Please contact us at %s so we can "
          "implement this feature", __func__, PACKAGE_BUGREPORT);

  /* Do the copying. */
  for(i=0;i<from->size;++i)
    {
      /* Set the pointer. */
      switch(to->type)
        {
        case GAL_TYPE_UINT8:    ptr = (uint8_t *)(to->array)  + i;   break;
        case GAL_TYPE_INT8:     ptr = (int8_t *)(to->array)   + i;   break;
        case GAL_TYPE_UINT16:   ptr = (uint16_t *)(to->array) + i;   break;
        case GAL_TYPE_INT16:    ptr = (int16_t *)(to->array)  + i;   break;
        case GAL_TYPE_UINT32:   ptr = (uint32_t *)(to->array) + i;   break;
        case GAL_TYPE_INT32:    ptr = (int32_t *)(to->array)  + i;   break;
        case GAL_TYPE_UINT64:   ptr = (uint64_t *)(to->array) + i;   break;
        case GAL_TYPE_INT64:    ptr = (int64_t *)(to->array)  + i;   break;
        case GAL_TYPE_FLOAT32:  ptr = (float *)(to->array)    + i;   break;
        case GAL_TYPE_FLOAT64:  ptr = (double *)(to->array)   + i;   break;
        case GAL_TYPE_BIT:
        case GAL_TYPE_STRLL:
        case GAL_TYPE_COMPLEX32:
        case GAL_TYPE_COMPLEX64:
          error(EXIT_FAILURE, 0, "%s: copying to %s type not currently "
                "supported", __func__, gal_type_name(to->type, 1));
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: type %d not recognized for to->type",
                __func__, to->type);
        }

      /* Read/put the input into the output. */
      if(to->type==GAL_TYPE_STRING)
        gal_checkset_allocate_copy(strarr[i], &outstrarr[i]);
      else
        {
          if( gal_type_from_string(&ptr, strarr[i], to->type) )
            data_copy_to_string_not_parsed(strarr[i], ptr, to->type);
        }
    }
}





/* Macros for copying to a string. */
#define COPY_TO_STR_INT(CTYPE, BLANK, FMT) {                            \
    CTYPE *a=from->array;                                               \
    for(i=0;i<from->size;++i)                                           \
      {                                                                 \
        if(a[i]!=BLANK)                                                 \
          {                                                             \
            if( asprintf(&strarr[i], FMT, a[i])<0 )                     \
              error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__); \
          }                                                             \
        else                                                            \
          gal_checkset_allocate_copy(GAL_BLANK_STRING, &strarr[i]);     \
      }                                                                 \
  }

#define COPY_TO_STR_FLT(CTYPE, BLANK) {                                 \
    CTYPE *a=from->array;                                               \
    for(i=0;i<from->size;++i)                                           \
      {                                                                 \
        if(isnan(BLANK)) isblank = isnan(a[i]) ? 1 : 0;                 \
        else             isblank = a[i]==BLANK ? 1 : 0;                 \
        if(isblank==0)                                                  \
          {                                                             \
            if( asprintf(&strarr[i], "%f", a[i])<0 )                    \
              error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__); \
          }                                                             \
        else gal_checkset_allocate_copy(GAL_BLANK_STRING, &strarr[i]);  \
      }                                                                 \
  }

/* Convert any given type into a string by printing it into the elements of
   the already allocated 'to->array'. */
static void
data_copy_to_string(gal_data_t *from, gal_data_t *to)
{
  size_t i;
  int isblank;
  char **strarr=to->array, **instrarr=from->array;

  /* Sanity check */
  if(to->type!=GAL_TYPE_STRING)
    error(EXIT_FAILURE, 0, "%s: 'to' must have a string type", __func__);
  if(from->block)
    error(EXIT_FAILURE, 0, "%s: tile inputs not currently supported "
          "('block' element must be NULL). Please contact us at %s so we "
          "can implement this feature", __func__, PACKAGE_BUGREPORT);

  /* Do the copying */
  switch(from->type)
    {
    case GAL_TYPE_UINT8:
      COPY_TO_STR_INT(uint8_t,  GAL_BLANK_UINT8,  "%"PRIu8);  break;

    case GAL_TYPE_INT8:
      COPY_TO_STR_INT(int8_t,   GAL_BLANK_INT8,   "%"PRId8);  break;

    case GAL_TYPE_UINT16:
      COPY_TO_STR_INT(uint16_t, GAL_BLANK_UINT16, "%"PRIu16); break;

    case GAL_TYPE_INT16:
      COPY_TO_STR_INT(int16_t,  GAL_BLANK_INT16,  "%"PRId16); break;

    case GAL_TYPE_UINT32:
      COPY_TO_STR_INT(uint32_t, GAL_BLANK_UINT32, "%"PRIu32); break;

    case GAL_TYPE_INT32:
      COPY_TO_STR_INT(int32_t,  GAL_BLANK_INT32,  "%"PRId32); break;

    case GAL_TYPE_UINT64:
      COPY_TO_STR_INT(uint64_t, GAL_BLANK_UINT64, "%"PRIu64); break;

    case GAL_TYPE_INT64:
      COPY_TO_STR_INT(int64_t,  GAL_BLANK_INT64,  "%"PRId64); break;

    case GAL_TYPE_FLOAT32:
      COPY_TO_STR_FLT(float, GAL_BLANK_FLOAT32);              break;

    case GAL_TYPE_FLOAT64:
      COPY_TO_STR_FLT(double, GAL_BLANK_FLOAT32);             break;

    case GAL_TYPE_STRING:
      for(i=0;i<from->size;++i)
        gal_checkset_allocate_copy(instrarr[i], &strarr[i]);
      break;

    case GAL_TYPE_BIT:
    case GAL_TYPE_STRLL:
    case GAL_TYPE_COMPLEX32:
    case GAL_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "%s: copying to %s type not currently supported",
            __func__, gal_type_name(from->type, 1));
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: type %d not recognized for 'from->type'",
            __func__, from->type);
    }
}





#define COPY_OT_IT_SET(OT, IT) {                                        \
    OT ob, *restrict o=out->array;                                      \
    size_t increment=0, num_increment=1;                                \
    size_t mclen=0, contig_len=in->dsize[in->ndim-1];                   \
    IT ib, *ist=NULL, *restrict i=in->array, *f=i+in->size;             \
    size_t s_e_ind[2]={0,iblock->size-1}; /* -1: this is INCLUSIVE */   \
                                                                        \
    /* If we are on a tile, the default values need to change. */       \
    if(in!=iblock)                                                      \
      ist=gal_tile_start_end_ind_inclusive(in, iblock, s_e_ind);        \
                                                                        \
    /* Constant preparations before the loop. */                        \
    if(iblock->type==out->type)                                         \
      mclen = in==iblock ? iblock->size : contig_len;                   \
    else                                                                \
      {                                                                 \
        gal_blank_write(&ob, out->type);                                \
        gal_blank_write(&ib, iblock->type);                             \
      }                                                                 \
                                                                        \
    /* Parse over the input and copy it. */                             \
    while( s_e_ind[0] + increment <= s_e_ind[1] )                       \
      {                                                                 \
        /* If we are on a tile, reset 'i' and  'f' for each round. */   \
        if(in!=iblock)                                                  \
          f = ( i = ist + increment ) + contig_len;                     \
                                                                        \
        /* When the types are the same just use memcopy, otherwise, */  \
        /* We'll have to read each number (and use internal         */  \
        /* conversion). */                                              \
        if(iblock->type==out->type)                                     \
          {                                                             \
            memcpy(o, i, mclen*gal_type_sizeof(iblock->type));          \
            o += mclen;                                                 \
          }                                                             \
        else                                                            \
          {                                                             \
            /* If the blank is a NaN value (only for floating point  */ \
            /* types), it will fail any comparison, so we'll exploit */ \
            /* this property in such cases. For other cases, a       */ \
            /* '*i==ib' is enough.                                   */ \
            if(ib==ib) do *o++ = *i==ib ? ob : *i; while(++i<f);        \
            else       do *o++ = *i!=*i ? ob : *i; while(++i<f);        \
          }                                                             \
                                                                        \
        /* Update the increment from the start of the input. */         \
        increment += ( in==iblock ? iblock->size                        \
                       : gal_tile_block_increment(iblock, in->dsize,    \
                                                  num_increment++,      \
                                                  NULL) );              \
      }                                                                 \
  }





/* gal_data_copy_to_new_type: Output type is set, now choose the input
   type. */
#define COPY_OT_SET(OT)                                                 \
  switch(iblock->type)                                                  \
    {                                                                   \
    case GAL_TYPE_UINT8:      COPY_OT_IT_SET(OT, uint8_t  );    break;  \
    case GAL_TYPE_INT8:       COPY_OT_IT_SET(OT, int8_t   );    break;  \
    case GAL_TYPE_UINT16:     COPY_OT_IT_SET(OT, uint16_t );    break;  \
    case GAL_TYPE_INT16:      COPY_OT_IT_SET(OT, int16_t  );    break;  \
    case GAL_TYPE_UINT32:     COPY_OT_IT_SET(OT, uint32_t );    break;  \
    case GAL_TYPE_INT32:      COPY_OT_IT_SET(OT, int32_t  );    break;  \
    case GAL_TYPE_UINT64:     COPY_OT_IT_SET(OT, uint64_t );    break;  \
    case GAL_TYPE_INT64:      COPY_OT_IT_SET(OT, int64_t  );    break;  \
    case GAL_TYPE_FLOAT32:    COPY_OT_IT_SET(OT, float    );    break;  \
    case GAL_TYPE_FLOAT64:    COPY_OT_IT_SET(OT, double   );    break;  \
    case GAL_TYPE_STRING:     data_copy_from_string(in, out);   break;  \
    case GAL_TYPE_BIT:                                                  \
    case GAL_TYPE_STRLL:                                                \
    case GAL_TYPE_COMPLEX32:                                            \
    case GAL_TYPE_COMPLEX64:                                            \
      error(EXIT_FAILURE, 0, "%s: copying from %s type to a numeric "   \
            "(real) type not supported", "COPY_OT_SET",                 \
            gal_type_name(in->type, 1));                                \
      break;                                                            \
                                                                        \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized for "     \
            "'in->type'", "COPY_OT_SET", in->type);                     \
    }





/* Wrapper for 'gal_data_copy_to_new_type', but will copy to the same type
   as the input. Recall that if the input is a tile (a part of the input,
   which is not-contiguous if it has more than one dimension), then the
   output will have only the elements that cover the tile.*/
gal_data_t *
gal_data_copy(gal_data_t *in)
{
  return gal_data_copy_to_new_type(in, gal_tile_block(in)->type);
}





/* Copy a given data structure to a new one with any type (for the
   output). The input can be a tile, in which case the output will be a
   contiguous patch of memory that has all the values within the input tile
   in the requested type. */
gal_data_t *
gal_data_copy_to_new_type(gal_data_t *in, uint8_t newtype)
{
  gal_data_t *out;

  /* Allocate the output datastructure. */
  out=gal_data_alloc(NULL, newtype, in->ndim, in->dsize, in->wcs,
                     0, in->minmapsize, in->quietmmap, in->name,
                     in->unit, in->comment);

  /* Fill in the output array: */
  gal_data_copy_to_allocated(in, out);

  /* Return the created array */
  return out;
}





/* Copy the input data structure into a new type and free the allocated
   space. */
gal_data_t *
gal_data_copy_to_new_type_free(gal_data_t *in, uint8_t newtype)
{
  gal_data_t *out, *iblock=gal_tile_block(in);

  /* In a general application, it might happen that the type is equal with
     the type of the input and the input isn't a tile. Since the job of
     this function is to free the input dataset, and the user just wants
     one dataset after this function finishes, we can safely just return
     the input. */
  if(newtype==iblock->type && in==iblock)
    return in;
  else
    {
      out=gal_data_copy_to_new_type(in, newtype);
      if(iblock==in)
        gal_data_free(in);
      else
        fprintf(stderr, "#####\nWarning from "
                "'gal_data_copy_to_new_type_free'\n#####\n The input "
                "dataset is a tile, not a contiguous (fully allocated) "
                "patch of memory. So it has not been freed. Please use "
                "'gal_data_copy_to_new_type' to avoid this warning.\n"
                "#####");
      return out;
    }
}





/* Copy a given dataset ('in') into an already allocated dataset 'out' (the
   actual dataset and its 'array' element). The meta-data of 'in' (except
   for 'block') will be fully copied into 'out' also. 'out->size' will be
   used to find the available space in the allocated space.

   When 'in->size != out->size' this function will behave as follows:

      'out->size < in->size': it won't re-allocate the necessary space, it
          will abort with an error, so please check before calling this
          function.

      'out->size > in->size': it will update 'out->size' and 'out->dsize'
          to be the same as the input. So if you want to re-use a
          pre-allocated space with varying input sizes, be sure to reset
          'out->size' before every call to this function. */
void
gal_data_copy_to_allocated(gal_data_t *in, gal_data_t *out)
{
  gal_data_t *iblock=gal_tile_block(in);

  /* Make sure the number of allocated elements (of whatever type) in the
     output is not smaller than the input. Note that the type is irrelevant
     because we will be doing type conversion if they differ.*/
  if( out->size < in->size  )
    error(EXIT_FAILURE, 0, "%s: the output dataset must be equal or larger "
          "than the input. the sizes are %zu and %zu respectively", __func__,
          out->size, in->size);
  if( out->ndim != in->ndim )
    error(EXIT_FAILURE, 0, "%s: the output dataset must have the same number "
          "of dimensions, the dimensions are %zu and %zu respectively",
          __func__, out->ndim, in->ndim);

  /* Free possibly allocated meta-data strings. */
  if(out->name)    free(out->name);
  if(out->unit)    free(out->unit);
  if(out->comment) free(out->comment);

  /* Write the basic meta-data. */
  out->flag           = in->flag;
  out->next           = in->next;
  out->status         = in->status;
  out->disp_width     = in->disp_width;
  out->disp_precision = in->disp_precision;
  gal_checkset_allocate_copy(in->name,    &out->name);
  gal_checkset_allocate_copy(in->unit,    &out->unit);
  gal_checkset_allocate_copy(in->comment, &out->comment);

  /* Do the copying. */
  switch(out->type)
    {
    case GAL_TYPE_UINT8:   COPY_OT_SET( uint8_t  );      break;
    case GAL_TYPE_INT8:    COPY_OT_SET( int8_t   );      break;
    case GAL_TYPE_UINT16:  COPY_OT_SET( uint16_t );      break;
    case GAL_TYPE_INT16:   COPY_OT_SET( int16_t  );      break;
    case GAL_TYPE_UINT32:  COPY_OT_SET( uint32_t );      break;
    case GAL_TYPE_INT32:   COPY_OT_SET( int32_t  );      break;
    case GAL_TYPE_UINT64:  COPY_OT_SET( uint64_t );      break;
    case GAL_TYPE_INT64:   COPY_OT_SET( int64_t  );      break;
    case GAL_TYPE_FLOAT32: COPY_OT_SET( float    );      break;
    case GAL_TYPE_FLOAT64: COPY_OT_SET( double   );      break;
    case GAL_TYPE_STRING:  data_copy_to_string(in, out); break;

    case GAL_TYPE_BIT:
    case GAL_TYPE_STRLL:
    case GAL_TYPE_COMPLEX32:
    case GAL_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "%s: copying to %s type not yet supported",
            __func__, gal_type_name(out->type, 1));
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: type %d not recognized for 'out->type'",
            __func__, out->type);
    }

  /* Correct the sizes of the output to be the same as the input. If it is
     equal, there is no problem, if not, the size information will be
     changed, so if you want to use this allocated space again, be sure to
     re-set the size parameters. */
  out->size=in->size;
  memcpy(out->dsize, in->dsize, in->ndim * sizeof *(in->dsize) );
}





/* Just a wrapper around 'gal_type_from_string_auto', to return a
   'gal_data_t' dataset hosting the allocated number. */
gal_data_t *
gal_data_copy_string_to_number(char *string)
{
  void *ptr;
  uint8_t type;
  size_t dsize=1;
  ptr=gal_type_string_to_number(string, &type);
  return ( ptr
           ? gal_data_alloc(ptr, type, 1, &dsize, NULL, 0, -1,
                            1, NULL, NULL, NULL)
           : NULL );
}
