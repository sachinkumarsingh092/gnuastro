/*********************************************************************
data -- Structure and functions to represent/work with data
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
#ifndef __GAL_DATA_H__
#define __GAL_DATA_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <math.h>
#include <limits.h>
#include <stdint.h>

#include <wcslib/wcs.h>
#include <gsl/gsl_complex.h>

#include <gnuastro/linkedlist.h>

/* When we are within Gnuastro's building process, `IN_GNUASTRO_BUILD' is
   defined. In the build process, installation information (in particular
   `GAL_CONFIG_SIZEOF_SIZE_T' that we need below) is kept in
   `config.h'. When building a user's programs, this information is kept in
   `gnuastro/config.h'. Note that all `.c' files in Gnuastro's source must
   start with the inclusion of `config.h' and that `gnuastro/config.h' is
   only created at installation time (not present during the building of
   Gnuastro).*/
#ifndef IN_GNUASTRO_BUILD
#include <gnuastro/config.h>
#endif

/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */



/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */





/* Macros to identify the type of data. */
enum gal_data_types
{
  GAL_DATA_TYPE_INVALID,     /* Invalid (=0 by C standard).             */

  GAL_DATA_TYPE_BIT,         /* 1 bit                                   */
  GAL_DATA_TYPE_UINT8,       /* 8-bit  unsigned integer.                */
  GAL_DATA_TYPE_INT8,        /* 8-bit  signed   integer.                */
  GAL_DATA_TYPE_UINT16,      /* 16-bit unsigned integer.                */
  GAL_DATA_TYPE_INT16,       /* 16-bit signed   integer.                */
  GAL_DATA_TYPE_UINT32,      /* 32-bit unsigned integer.                */
  GAL_DATA_TYPE_INT32,       /* 32-bit signed   integer.                */
  GAL_DATA_TYPE_UINT64,      /* 64-bit unsigned integer.                */
  GAL_DATA_TYPE_INT64,       /* 64-bit signed   integer.                */
  GAL_DATA_TYPE_FLOAT32,     /* 32-bit single precision floating point. */
  GAL_DATA_TYPE_FLOAT64,     /* 64-bit double precision floating point. */
  GAL_DATA_TYPE_COMPLEX32,   /* Complex 32-bit floating point.          */
  GAL_DATA_TYPE_COMPLEX64,   /* Complex 64-bit floating point.          */
  GAL_DATA_TYPE_STRING,      /* String of characters.                   */
  GAL_DATA_TYPE_STRLL,       /* Linked list of strings.                 */
};

/* `size_t' is 4 and 8 bytes on 32 and 64 bit systems respectively. In both
   cases, the standard defines `size_t' to be unsigned. During
   `./configure' the sizeof size_t was found and is stored in
   `GAL_CONFIG_SIZEOF_SIZE_T'. */
#if GAL_CONFIG_SIZEOF_SIZE_T == 4
#define GAL_DATA_TYPE_SIZE_T GAL_DATA_TYPE_UINT32
#else
#define GAL_DATA_TYPE_SIZE_T GAL_DATA_TYPE_UINT64
#endif





/* Main data structure.

   Notes
   -----

    - If mmapname==NULL, then the array is allocated (using malloc, in the
      RAM), otherwise its is mmap'd (is actually a file on the ssd/hdd).

    - minmapsize is stored in the data structure to allow any derivative
      data structures to follow the same number and decide if they should
      be mmap'd or allocated.

      - `minmapsize' ==  0:  array is definitely mmap'd.

      - `minmapsize' == -1: array is definitely in RAM.   */
typedef struct gal_data_t
{
  /* Basic information on array of data. */
  void             *array;  /* Array keeping data elements.                */
  uint8_t            type;  /* Type of data (from `gal_data_alltypes').    */
  size_t             ndim;  /* Number of dimensions in the array.          */
  size_t           *dsize;  /* Size of array along each dimension.         */
  size_t             size;  /* Total number of data-elements.              */
  char          *mmapname;  /* File name of the mmap.                      */
  size_t       minmapsize;  /* Minimum number of bytes to mmap the array.  */

  /* WCS information. */
  int                nwcs;  /* for WCSLIB: no. coord. representations.     */
  struct wcsprm      *wcs;  /* WCS information for this dataset.           */

  /* Content descriptions. */
  int              status;  /* Any context-specific status value.          */
  char              *name;  /* e.g., EXTNAME, or column, or keyword.       */
  char              *unit;  /* Units of the data.                          */
  char           *comment;  /* A more detailed description of the data.    */

  /* For printing */
  int            disp_fmt;  /* See `gal_table_diplay_formats'.             */
  int          disp_width;  /* Width of space to print in ASCII.           */
  int      disp_precision;  /* Precision to print in ASCII.                */

  /* As linked list. */
  struct gal_data_t *next;  /* To use it as a linked list if necessary.    */
} gal_data_t;




/*************************************************************
 **************        Type information        ***************
 *************************************************************/
char *
gal_data_type_as_string(uint8_t type, int long_name);

uint8_t
gal_data_string_as_type(char *str);

void
gal_data_type_min(uint8_t type, void *in);

void
gal_data_type_max(uint8_t type, void *in);

int
gal_data_is_linked_list(uint8_t type);



/*********************************************************************/
/*************         Size and allocation         *******************/
/*********************************************************************/
void
gal_data_copy_wcs(gal_data_t *in, gal_data_t *out);

int
gal_data_dsize_is_different(gal_data_t *first, gal_data_t *second);

size_t
gal_data_sizeof(uint8_t type);

void *
gal_data_malloc_array(uint8_t type, size_t size);

void *
gal_data_calloc_array(uint8_t type, size_t size);

void *
gal_data_alloc_number(uint8_t type, void *number);

void
gal_data_initialize(gal_data_t *data, void *array, uint8_t type, size_t ndim,
                    size_t *dsize, struct wcsprm *wcs, int clear,
                    size_t minmapsize, char *name, char *unit, char *comment);

gal_data_t *
gal_data_alloc(void *array, uint8_t type, size_t ndim, size_t *dsize,
               struct wcsprm *wcs, int clear, size_t minmapsize,
               char *name, char *unit, char *comment);

gal_data_t *
gal_data_calloc_dataarray(size_t size);

size_t
gal_data_string_fixed_alloc_size(gal_data_t *data);

void
gal_data_free_contents(gal_data_t *data);

void
gal_data_free(gal_data_t *data);



/*********************************************************************/
/*************    Data structure as a linked list   ******************/
/*********************************************************************/
void
gal_data_add_existing_to_ll(gal_data_t **list, gal_data_t *newnode);

void
gal_data_add_to_ll(gal_data_t **list, void *array, uint8_t type, size_t ndim,
                   size_t *dsize, struct wcsprm *wcs, int clear,
                   size_t minmapsize, char *name, char *unit, char *comment);

gal_data_t *
gal_data_pop_from_ll(struct gal_data_t **list);

void
gal_data_reverse_ll(gal_data_t **list);

size_t
gal_data_num_in_ll(struct gal_data_t *list);

gal_data_t **
gal_data_ll_to_array_of_ptrs(gal_data_t *list, size_t *num);

void
gal_data_free_ll(gal_data_t *list);



/*************************************************************
 **************            Copying             ***************
 *************************************************************/
gal_data_t *
gal_data_copy_to_new_type(gal_data_t *in, uint8_t newtype);

gal_data_t *
gal_data_copy_to_new_type_free(gal_data_t *in, uint8_t type);

gal_data_t *
gal_data_copy(gal_data_t *in);

int
gal_data_out_type(gal_data_t *first, gal_data_t *second);

void
gal_data_to_same_type(gal_data_t *f, gal_data_t *s, gal_data_t **of,
                      gal_data_t **os, uint8_t type, int freeinputs);

void
gal_data_copy_element_same_type(gal_data_t *input, size_t index, void *ptr);



/*************************************************************
 **************              Write             ***************
 *************************************************************/
char *
gal_data_write_to_string(void *ptr, uint8_t type, int quote_if_str_has_space);



/*************************************************************
 **************              Read              ***************
 *************************************************************/
gal_data_t *
gal_data_string_to_number(char *string);

int
gal_data_string_to_type(void **out, char *string, uint8_t type);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_DATA_H__ */
