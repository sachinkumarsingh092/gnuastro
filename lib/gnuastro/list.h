/*********************************************************************
Functions for linked lists.
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
#ifndef __GAL_LIST_H__
#define __GAL_LIST_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>


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





/****************************************************************
 *****************          String           ********************
 ****************************************************************/
typedef struct gal_list_str_t
{
  char *v;
  struct gal_list_str_t *next;
} gal_list_str_t;

void
gal_list_str_add(gal_list_str_t **list, char *value,
                 int allocate);

char *
gal_list_str_pop(gal_list_str_t **list);

size_t
gal_list_str_number(gal_list_str_t *list);

gal_list_str_t *
gal_list_str_last(gal_list_str_t *list);

void
gal_list_str_print(gal_list_str_t *list);

void
gal_list_str_reverse(gal_list_str_t **list);

void
gal_list_str_free(gal_list_str_t *list, int freevalue);





/****************************************************************
 *****************          int32            ********************
 ****************************************************************/
typedef struct gal_list_i32_t
{
  int32_t v;
  struct gal_list_i32_t *next;
} gal_list_i32_t;

void
gal_list_i32_add(gal_list_i32_t **list, int32_t value);

int32_t
gal_list_i32_pop(gal_list_i32_t **list);

size_t
gal_list_i32_number(gal_list_i32_t *list);

gal_list_i32_t *
gal_list_i32_last(gal_list_i32_t *list);

void
gal_list_i32_print(gal_list_i32_t *list);

void
gal_list_i32_reverse(gal_list_i32_t **list);

int32_t *
gal_list_i32_to_array(gal_list_i32_t *list, int reverse, size_t *num);

void
gal_list_i32_free(gal_list_i32_t *list);





/****************************************************************
 *****************           size_t          ********************
 ****************************************************************/
typedef struct gal_list_sizet_t
{
  size_t v;
  struct gal_list_sizet_t *next;
} gal_list_sizet_t;

void
gal_list_sizet_add(gal_list_sizet_t **list, size_t value);

size_t
gal_list_sizet_pop(gal_list_sizet_t **list);

size_t
gal_list_sizet_number(gal_list_sizet_t *list);

gal_list_sizet_t *
gal_list_sizet_last(gal_list_sizet_t *list);

void
gal_list_sizet_print(gal_list_sizet_t *list);

void
gal_list_sizet_reverse(gal_list_sizet_t **list);

size_t *
gal_list_sizet_to_array(gal_list_sizet_t *list, int reverse, size_t *num);

void
gal_list_sizet_free(gal_list_sizet_t *list);





/****************************************************************
 *****************           float           ********************
 ****************************************************************/
typedef struct gal_list_f32_t
{
    float v;
    struct gal_list_f32_t *next;
} gal_list_f32_t;

void
gal_list_f32_add(gal_list_f32_t **list, float value);

float
gal_list_f32_pop(gal_list_f32_t **list);

size_t
gal_list_f32_number(gal_list_f32_t *list);

gal_list_f32_t *
gal_list_f32_last(gal_list_f32_t *list);

void
gal_list_f32_reverse(gal_list_f32_t **list);

void
gal_list_f32_print(gal_list_f32_t *list);

float *
gal_list_f32_to_array(gal_list_f32_t *list, int reverse, size_t *num);

void
gal_list_f32_free(gal_list_f32_t *list);





/****************************************************************
 *****************          double           ********************
 ****************************************************************/
typedef struct gal_list_f64_t
{
    double v;
    struct gal_list_f64_t *next;
} gal_list_f64_t;

void
gal_list_f64_add(gal_list_f64_t **list, double value);

double
gal_list_f64_pop(gal_list_f64_t **list);

size_t
gal_list_f64_number(gal_list_f64_t *list);

gal_list_f64_t *
gal_list_f64_last(gal_list_f64_t *list);

void
gal_list_f64_print(gal_list_f64_t *list);

void
gal_list_f64_reverse(gal_list_f64_t **list);

double *
gal_list_f64_to_array(gal_list_f64_t *list, int reverse, size_t *num);

void
gal_list_f64_free(gal_list_f64_t *list);





/****************************************************************
 *****************          void *           ********************
 ****************************************************************/
typedef struct gal_list_void_t
{
  void *v;
  struct gal_list_void_t *next;
} gal_list_void_t;

void
gal_list_void_add(gal_list_void_t **list, void *value);

void *
gal_list_void_pop(gal_list_void_t **list);

size_t
gal_list_void_number(gal_list_void_t *list);

gal_list_void_t *
gal_list_void_last(gal_list_void_t *list);

void
gal_list_void_reverse(gal_list_void_t **list);

void
gal_list_void_free(gal_list_void_t *list, int freevalue);





/****************************************************************
 *****************       Ordered size_t      ********************
 ****************************************************************/
typedef struct gal_list_osizet_t
{
  size_t v;                       /* The actual value. */
  float s;                        /* The parameter to sort by. */
  struct gal_list_osizet_t *next;
} gal_list_osizet_t;

void
gal_list_osizet_add(gal_list_osizet_t **list,
                    size_t value, float tosort);

size_t
gal_list_osizet_pop(gal_list_osizet_t **list, float *sortvalue);

void
gal_list_osizet_to_sizet_free(gal_list_osizet_t *in,
                              gal_list_sizet_t **out);





/****************************************************************
 ***********     Doubly linked, ordered size_t     **************
 ****************************************************************/
typedef struct gal_list_dosizet_t
{
  size_t v;                       /* The actual value. */
  float s;                        /* The parameter to sort by. */
  struct gal_list_dosizet_t *prev;
  struct gal_list_dosizet_t *next;
} gal_list_dosizet_t;

void
gal_list_dosizet_add(gal_list_dosizet_t **largest,
                     gal_list_dosizet_t **smallest, size_t value, float tosort);

size_t
gal_list_dosizet_pop_smallest(gal_list_dosizet_t **lartest,
                              gal_list_dosizet_t **smallest, float *tosort);

void
gal_list_dosizet_print(gal_list_dosizet_t *largest,
                       gal_list_dosizet_t *smallest);

void
gal_list_dosizet_to_sizet(gal_list_dosizet_t *in,
                          gal_list_sizet_t **out);

void
gal_list_dosizet_free(gal_list_dosizet_t *largest);





/****************************************************************
 *****************        gal_data_t         ********************
 ****************************************************************/
void
gal_list_data_add(gal_data_t **list, gal_data_t *newnode);

void
gal_list_data_add_alloc(gal_data_t **list, void *array, uint8_t type,
                        size_t ndim, size_t *dsize, struct wcsprm *wcs,
                        int clear, size_t minmapsize, int quietmmap,
                        char *name, char *unit, char *comment);

gal_data_t *
gal_list_data_pop(gal_data_t **list);

void
gal_list_data_reverse(gal_data_t **list);

gal_data_t **
gal_list_data_to_array_ptr(gal_data_t *list, size_t *num);

size_t
gal_list_data_number(gal_data_t *list);

gal_data_t *
gal_list_data_last(gal_data_t *list);

void
gal_list_data_free(gal_data_t *list);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_LIST_H__ */
