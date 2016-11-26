/*********************************************************************
Functions to check and set command line argument values and files.
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
#ifndef __GAL_CHECKSET_H__
#define __GAL_CHECKSET_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <math.h>
#include <fitsio.h>


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
/*********                 Macros                **************/
/**************************************************************/
/* Check if the column number is within the boundaries of a catalog. */
#define GAL_CHECKSET_CHECK_COL_NUM_IN_CAT(INCOL,NAME) {                 \
    if( (INCOL) >= p->cs1 )                                             \
      error(EXIT_FAILURE, 0, "%s only has %zu columns while you "       \
            "have requested column %zu (counting from zero) for "       \
            "`--%s`", p->up.catname, p->cs1, (INCOL), (NAME));          \
  }





/* Check both the column number and that the values in the column are
   finitte numbers (not infinities or NAN). */
#define GAL_CHECKSET_CHECK_COL_IN_CAT(INCOL,NAME) {                     \
    size_t i;                                                           \
                                                                        \
    if( (INCOL) >= p->cs1 )                                             \
      error(EXIT_FAILURE, 0, "%s only has %zu columns while you "       \
            "have requested column %zu (counting from zero) for "       \
            "`--%s`", p->up.catname, p->cs1, (INCOL), (NAME));          \
                                                                        \
    for(i=0;i<p->cs0;++i)                                               \
      if( !isfinite(p->cat[i*p->cs1+(INCOL)]) )                         \
        error(EXIT_FAILURE, 0, "%s: column %zu (--%s) in row %zu "      \
              "could not be read as a number. See %s. Note that "       \
              "counting starts from zero",                              \
              p->up.catname, (INCOL), (NAME), i, GAL_TXTARRAY_LOG);     \
  }




#define GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE(name,string) {       \
    if(gal_checkset_string_has_space(string))                           \
      fprintf(fp, CONF_SHOWFMT"\"%s\"\n", name, string);                \
    else                                                                \
      fprintf(fp, CONF_SHOWFMT"%s\n", name, string);                    \
  }



/****************************************************************
 ************      Check and convert strings    *****************
 ****************************************************************/
int
gal_checkset_str_is_double(char *string, double *out);



/**************************************************************/
/********* Read arguments and check their values **************/
/**************************************************************/
void
gal_checkset_int_zero_or_one(char *optarg, int *var, char *lo, char so,
                             char* spack, char *filename, size_t lineno);

void
gal_checkset_int_4_or_8(char *optarg, int *var, char *lo, char so,
                        char *spack, char *filename, size_t lineno);

void
gal_checkset_int_el_zero(char *optarg, int *var, char *lo, char so,
                         char *spack, char *filename, size_t lineno);

void
gal_checkset_int_l_zero(char *optarg, int *var, char *lo, char so,
                        char *spack, char *filename, size_t lineno);

void
gal_checkset_int_smaller_equal_to(char *optarg, int *var, char *lo,
                                  char so, char *spack, char *filename,
                                  size_t lineno, long maxvalue);

void
gal_checkset_long_el_zero(char *optarg, long *var, char *lo, char so,
                          char *spack, char *filename, size_t lineno);

void
gal_checkset_any_long(char *optarg, long *var, char *lo, char so,
                      char *spack, char *filename, size_t lineno);

void
gal_checkset_sizet_el_zero(char *optarg, size_t *var, char *lo, char so,
                           char *spack, char *filename, size_t lineno);

void
gal_checkset_sizet_l_zero(char *optarg, size_t *var, char *lo, char so,
                          char *spack, char *filename, size_t lineno);

void
gal_checkset_sizet_p_odd(char *optarg, size_t *var, char *lo, char so,
                         char* spack, char *filename, size_t lineno);

void
gal_checkset_float_l_0(char *optarg, float *var, char *lo, char so,
                       char *spack, char *filename, size_t lineno);

void
gal_checkset_float_l_0_s_1(char *optarg, float *var, char *lo, char so,
                           char *spack, char *filename, size_t lineno);

void
gal_checkset_any_float(char *optarg, float *var, char *lo, char so,
                       char *spack, char *filename, size_t lineno);

void
gal_checkset_double_l_0(char *optarg, double *var, char *lo, char so,
                        char *spack, char *filename, size_t lineno);

void
gal_checkset_double_el_0(char *optarg, double *var, char *lo, char so,
                         char* spack, char *filename, size_t lineno);

void
gal_checkset_double_l_value(char *optarg, double *var, char *lo, char so,
                            char* spack, double value, char *filename,
                            size_t lineno);

void
gal_checkset_any_double(char *optarg, double *var, char *lo, char so,
                        char *spack, char *filename, size_t lineno);



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
gal_checkset_allocate_copy(char *arg, char **copy);

void
gal_checkset_allocate_copy_set(char *arg, char **copy, int *set);




/**************************************************************/
/********** Set file names and check if they exist ************/
/**************************************************************/
void
gal_checkset_check_file(char *filename);

int
gal_checkset_check_file_report(char *filename);

void
gal_checkset_check_remove_file(char *filename, int dontdelete);

int
gal_checkset_dir_0_file_1(char *name, int dontdelete);

void
gal_checkset_automatic_output(char *inname, char *suffix,
                              int removedirinfo, int dontdelete,
                              char **outname);

char *
gal_checkset_dir_part(char *input);

char *
gal_checkset_not_dir_part(char *input);

void
gal_checkset_check_dir_write_add_slash(char **dirname);

void
gal_checkset_mkdir(char *dirname);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_CHECKSET_H__ */
