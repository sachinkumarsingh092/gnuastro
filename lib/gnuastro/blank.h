/*********************************************************************
blank -- Deal with blank values in a datasets
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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
#ifndef __GAL_BLANK_H__
#define __GAL_BLANK_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>

/* When we are within Gnuastro's building process, 'IN_GNUASTRO_BUILD' is
   defined. In the build process, installation information (in particular
   'GAL_CONFIG_SIZEOF_SIZE_T' that we need below) is kept in
   'config.h'. When building a user's programs, this information is kept in
   'gnuastro/config.h'. Note that all '.c' files in Gnuastro's source must
   start with the inclusion of 'config.h' and that 'gnuastro/config.h' is
   only created at installation time (not present during the building of
   Gnuastro). */
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


/* Macros: */

/* Blank values: Note that for the unsigned types or small types (like
   char), the maximum value is considered as a blank value, since the
   minimum value of an unsigned type is zero and zero is often meaningful
   in contexts were unsigned values are used. */
#define GAL_BLANK_UINT8      UINT8_MAX
#define GAL_BLANK_INT8       INT8_MIN
#define GAL_BLANK_UINT16     UINT16_MAX
#define GAL_BLANK_INT16      INT16_MIN
#define GAL_BLANK_UINT32     UINT32_MAX
#define GAL_BLANK_INT32      INT32_MIN
#define GAL_BLANK_UINT64     UINT64_MAX
#define GAL_BLANK_INT64      INT64_MIN
#define GAL_BLANK_FLOAT32    NAN
#define GAL_BLANK_FLOAT64    NAN
#define GAL_BLANK_STRING     "n/a"

#if GAL_CONFIG_SIZEOF_SIZE_T == 4
#define GAL_BLANK_SIZE_T     GAL_BLANK_UINT32
#else
#define GAL_BLANK_SIZE_T     GAL_BLANK_UINT64
#endif

#if GAL_CONFIG_SIZEOF_LONG == 4
#define GAL_BLANK_LONG  GAL_BLANK_INT32
#define GAL_BLANK_ULONG GAL_BLANK_UINT32
#elif GAL_CONFIG_SIZEOF_LONG == 8
#define GAL_BLANK_LONG  GAL_BLANK_INT64
#define GAL_BLANK_ULONG GAL_BLANK_UINT64
#endif

#if GAL_CONFIG_SIZEOF_INT == 4
#define GAL_BLANK_INT  GAL_BLANK_INT32
#define GAL_BLANK_UINT GAL_BLANK_UINT32
#elif GAL_CONFIG_SIZEOF_INT == 2
#define GAL_BLANK_INT  GAL_BLANK_INT16
#define GAL_BLANK_UINT GAL_BLANK_UINT16
#endif


/* Functions. */
void
gal_blank_write(void *pointer, uint8_t type);

void *
gal_blank_alloc_write(uint8_t type);

void
gal_blank_initialize(gal_data_t *input);

void
gal_blank_initialize_array(void *array, size_t size, uint8_t type);

char *
gal_blank_as_string(uint8_t type, int width);

int
gal_blank_is(void *pointer, uint8_t type);

int
gal_blank_present(gal_data_t *input, int updateflag);

size_t
gal_blank_number(gal_data_t *input, int updateflag);

gal_data_t *
gal_blank_flag(gal_data_t *data);

void
gal_blank_flag_apply(gal_data_t *input, gal_data_t *flag);

void
gal_blank_remove(gal_data_t *data);

void
gal_blank_remove_realloc(gal_data_t *input);


__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_DATA_H__ */
