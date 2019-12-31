/*********************************************************************
pointer -- facilitate working with pointers and allocation.
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
#ifndef __GAL_POINTER_H__
#define __GAL_POINTER_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <stdint.h>


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





void *
gal_pointer_increment(void *pointer, size_t increment, uint8_t type);

size_t
gal_pointer_num_between(void *earlier, void *later, uint8_t type);

void *
gal_pointer_allocate(uint8_t type, size_t size, int clear,
                     const char *funcname, const char *varname);

void *
gal_pointer_allocate_mmap(uint8_t type, size_t size, int clear,
                          char **filename, int quiet);





__END_C_DECLS    /* From C++ preparations */
#endif
