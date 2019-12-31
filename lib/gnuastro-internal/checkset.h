/*********************************************************************
Functions to check and set command line argument values and files.
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
#ifndef __GAL_CHECKSET_H__
#define __GAL_CHECKSET_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <math.h>
#include <fitsio.h>
#include <gsl/gsl_rng.h>
#include <gnuastro-internal/options.h>

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




/**************************************************************/
/**********               Environment              ************/
/**************************************************************/
gsl_rng *
gal_checkset_gsl_rng(uint8_t envseed_bool, const char **name,
                     unsigned long int *seed);



/**************************************************************/
/**********          Check fixed strings           ************/
/**************************************************************/
void
gal_checkset_known_types(char *optarg, int *type, char *filename,
                         size_t lineno);



/**************************************************************/
/**********          My String functions:          ************/
/**************************************************************/
int
gal_checkset_string_has_space(char *in);

char *
gal_checkset_malloc_cat(char *inname, char *toappend);

void
gal_checkset_allocate_copy(const char *arg, char **copy);

void
gal_checkset_allocate_copy_set(char *arg, char **copy, int *set);

char *
gal_checkset_dataset_name(char *filename, char *hdu);




/**************************************************************/
/********** Set file names and check if they exist ************/
/**************************************************************/
char *
gal_checkset_dir_part(char *filename);

char *
gal_checkset_not_dir_part(char *filename);

void
gal_checkset_check_file(char *filename);

int
gal_checkset_check_file_return(char *filename);

int
gal_checkset_writable_notexist(char *filename);

void    /* keep==0 && dontdelete==0: file will be deleted if exists.*/
gal_checkset_writable_remove(char *filename, int keep, int dontdelete);

int
gal_checkset_dir_0_file_1(char *name, int dontdelete);

char *
gal_checkset_automatic_output(struct gal_options_common_params *cp,
                              char *inname, char *suffix);

void
gal_checkset_check_dir_write_add_slash(char **dirname);

int
gal_checkset_mkdir(char *dirname);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_CHECKSET_H__ */
