/*********************************************************************
Units -- Convert data from one unit to other.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Kartik Ohri <kartikohri13@gmail.com>
Contributing author(s):
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Copyright (C) 2020, Free Software Foundation, Inc.

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
#ifndef __GAL_UNITS_H__
#define __GAL_UNITS_H__

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




int
gal_units_extract_decimal(char *convert, const char *delimiter,
                          double *args, size_t n);

double
gal_units_ra_to_degree (char *convert);

double
gal_units_dec_to_degree (char *convert);

char *
gal_units_degree_to_ra (double decimal);

char *
gal_units_degree_to_dec (double decimal);

__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_UNITS_H__ */
