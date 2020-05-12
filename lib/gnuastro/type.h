/*********************************************************************
Type -- Type information and basic operations.
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
#ifndef __GAL_TYPE_H__
#define __GAL_TYPE_H__

#include <limits.h>
#include <stdint.h>

#include <gsl/gsl_complex.h>

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





/*************************************************************
 **************           Constants            ***************
 *************************************************************/
/* Macros to identify the type of data. */
enum gal_types
{
  GAL_TYPE_INVALID,         /* Invalid (=0 by C standard).             */

  GAL_TYPE_BIT,             /* 1 bit                                   */
  GAL_TYPE_UINT8,           /* 8-bit  unsigned integer.                */
  GAL_TYPE_INT8,            /* 8-bit  signed   integer.                */
  GAL_TYPE_UINT16,          /* 16-bit unsigned integer.                */
  GAL_TYPE_INT16,           /* 16-bit signed   integer.                */
  GAL_TYPE_UINT32,          /* 32-bit unsigned integer.                */
  GAL_TYPE_INT32,           /* 32-bit signed   integer.                */
  GAL_TYPE_UINT64,          /* 64-bit unsigned integer.                */
  GAL_TYPE_INT64,           /* 64-bit signed   integer.                */
  GAL_TYPE_FLOAT32,         /* 32-bit single precision floating point. */
  GAL_TYPE_FLOAT64,         /* 64-bit double precision floating point. */
  GAL_TYPE_COMPLEX32,       /* Complex 32-bit floating point.          */
  GAL_TYPE_COMPLEX64,       /* Complex 64-bit floating point.          */
  GAL_TYPE_STRING,          /* String of characters.                   */
  GAL_TYPE_STRLL,           /* Linked list of strings.                 */
};



/* Define system specific types. For example 'size_t' is 4 and 8 bytes on
   32 and 64 bit systems respectively. In both cases, the standard defines
   'size_t' to be unsigned. A similar case exists for 'long', but it is
   signed. During './configure' the sizeof 'size_t' and 'long' were found
   and are used to define an alias for these system specific types.

   Note: we are not using 'else'. This is done because by any chance, if
   the length of these types is not what is expected (4 or 8), then the
   aliases are not defined and the compiler will crash. */
#if GAL_CONFIG_SIZEOF_SIZE_T == 4
#define GAL_TYPE_SIZE_T GAL_TYPE_UINT32
#elif GAL_CONFIG_SIZEOF_SIZE_T == 8
#define GAL_TYPE_SIZE_T GAL_TYPE_UINT64
#endif

#if GAL_CONFIG_SIZEOF_LONG == 4
#define GAL_TYPE_LONG  GAL_TYPE_INT32
#define GAL_TYPE_ULONG GAL_TYPE_UINT32
#elif GAL_CONFIG_SIZEOF_LONG == 8
#define GAL_TYPE_LONG  GAL_TYPE_INT64
#define GAL_TYPE_ULONG GAL_TYPE_UINT64
#endif

#if GAL_CONFIG_SIZEOF_INT == 2
#define GAL_TYPE_INT  GAL_TYPE_INT16
#define GAL_TYPE_UINT GAL_TYPE_UINT16
#elif GAL_CONFIG_SIZEOF_INT == 4
#define GAL_TYPE_INT  GAL_TYPE_INT32
#define GAL_TYPE_UINT GAL_TYPE_UINT32
#endif







/*************************************************************
 **************         General info           ***************
 *************************************************************/
size_t
gal_type_sizeof(uint8_t type);

char *
gal_type_name(uint8_t type, int long_name);

uint8_t
gal_type_from_name(char *str);

void
gal_type_min(uint8_t type, void *in);

void
gal_type_max(uint8_t type, void *in);

int
gal_type_is_int(uint8_t type);

int
gal_type_is_list(uint8_t type);

int
gal_type_out(int first_type, int second_type);





/*************************************************************
 **************         To/from string         ***************
 *************************************************************/

char *
gal_type_bit_string(void *in, size_t size);

char *
gal_type_to_string(void *ptr, uint8_t type, int quote_if_str_has_space);

int
gal_type_from_string(void **out, char *string, uint8_t type);

void *
gal_type_string_to_number(char *string, uint8_t *type);





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_TYPE_H__ */
