/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2019, Free Software Foundation, Inc.

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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <gsl/gsl_heapsort.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/qsort.h>
#include <gnuastro/pointer.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>
#include <gnuastro/permutation.h>

#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "ui.h"
#include "arithmetic.h"



/**************************************************************/
/********     Selecting and ordering of columns      **********/
/**************************************************************/
static void
table_apply_permutation(gal_data_t *table, size_t *permutation,
                        size_t newsize, int inverse)
{
  gal_data_t *tmp;

  for(tmp=table;tmp!=NULL;tmp=tmp->next)
    {
      /* Apply the permutation. */
      if(inverse)
        gal_permutation_apply_inverse(tmp, permutation);
      else
        gal_permutation_apply(tmp, permutation);

      /* Correct the size. */
      tmp->size=tmp->dsize[0]=newsize;
    }
}





static void
table_range(struct tableparams *p)
{
  uint8_t *u;
  double *rarr;
  gal_data_t *mask;
  struct list_range *tmp;
  gal_data_t *ref, *perm, *range, *blmask;
  size_t i, g, b, *s, *sf, one=1, ngood=0;
  gal_data_t *min, *max, *ltmin, *gemax, *sum;

  int numok=GAL_ARITHMETIC_NUMOK;
  int inplace=GAL_ARITHMETIC_INPLACE;

  /* Allocate datasets for the necessary numbers and write them in. */
  min=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one, NULL, 0, -1, 1,
                     NULL, NULL, NULL);
  max=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one, NULL, 0, -1, 1,
                     NULL, NULL, NULL);
  perm=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, p->table->dsize, NULL, 0,
                      p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);
  mask=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, p->table->dsize, NULL, 1,
                      p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);

  /* Go over all the necessary range options. */
  range=p->range;
  for(tmp=p->rangecol;tmp!=NULL;tmp=tmp->next)
    {
      /* Set the minimum and maximum values. */
      rarr=range->array;
      ((double *)(min->array))[0] = rarr[0];
      ((double *)(max->array))[0] = rarr[1];

      /* Set the reference column to read values from. */
      ref=tmp->v;

      /* Find all the bad elements (smaller than the minimum, larger than
         the maximum or blank) so we can flag them. */
      ltmin=gal_arithmetic(GAL_ARITHMETIC_OP_LT, 1, numok, ref, min);
      gemax=gal_arithmetic(GAL_ARITHMETIC_OP_GE, 1, numok, ref, max);
      blmask = ( gal_blank_present(ref, 1)
                 ? gal_arithmetic(GAL_ARITHMETIC_OP_ISBLANK, 1, 0, ref)
                 : NULL );

      /* Merge all the flags into one array. */
      ltmin=gal_arithmetic(GAL_ARITHMETIC_OP_OR, 1, inplace, ltmin, gemax);
      if(blmask)
        ltmin=gal_arithmetic(GAL_ARITHMETIC_OP_OR, 1, inplace, ltmin, blmask);

      /* Add these flags to all previous flags. */
      mask=gal_arithmetic(GAL_ARITHMETIC_OP_OR, 1, inplace, mask, ltmin);

      /* For a check.
      {
        float *f=ref->array;
        uint8_t *m=mask->array;
        uint8_t *u=ltmin->array, *uf=u+ltmin->size;
        printf("\n\nInput column: %s\n", ref->name ? ref->name : "No Name");
        printf("Range: %g, %g\n", rarr[0], rarr[1]);
        printf("%-20s%-20s%-20s\n", "Value", "This mask",
               "Including previous");
        do printf("%-20f%-20u%-20u\n", *f++, *u++, *m++); while(u<uf);
        exit(0);
      }
      */

      /* Clean up. */
      gal_data_free(ltmin);
      gal_data_free(gemax);

      /* Increment pointers. */
      range=range->next;
    }

  /* Count the number of bad elements. */
  sum=gal_statistics_sum(mask);
  ngood = p->table->size - ((double *)(sum->array))[0];

  /* Define the permutation: elements within range remain on the top of
     the list, while the ones outside of will be placed after them
     (starting from the index after the last good one). */
  g=0;          /* Good indexs (starting from 0). */
  b=ngood;      /* Bad indexs (starting from total number of good). */
  u=mask->array;
  sf=(s=perm->array)+perm->size;
  do *s = *u++ ? b++ : g++; while(++s<sf);

  /* For a check
  {
    size_t i;
    double *v=ref->array;
    uint8_t *a=mask->array;
    printf("ref->type: %s\n", gal_type_name(ref->type, 1));
    for(i=0;i<ref->size;++i) printf("%u, %g\n", a[i], v[i]);
    gal_permutation_check(perm->array, perm->size);
  }
  */

  /* Apply the final permutation to the whole table. */
  table_apply_permutation(p->table, perm->array, ngood, 1);

  /* If the sort column is not in the table (the proper range has already
     been applied to it), and we need to sort the resulting columns
     afterwards, we should also apply the permutation on the sort
     column. */
  if(p->sortcol && p->sortin==0)
    table_apply_permutation(p->sortcol, perm->array, ngood, 1);

  /* Clean up. */
  i=0;
  for(tmp=p->rangecol;tmp!=NULL;tmp=tmp->next)
    { if(p->freerange[i]) {gal_data_free(tmp->v); tmp->v=NULL;} ++i; }
  ui_list_range_free(p->rangecol, 0);
  gal_data_free(mask);
  gal_data_free(perm);
  gal_data_free(sum);
  gal_data_free(min);
  gal_data_free(max);
  free(p->freerange);
}





static void
table_sort(struct tableparams *p)
{
  gal_data_t *perm;
  size_t c=0, *s, *sf;
  int (*qsortfn)(const void *, const void *)=NULL;

  /* In case there are no columns to sort, skip this function. */
  if(p->table->size==0) return;

  /* Allocate the permutation array and fill it. */
  perm=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, p->table->dsize, NULL, 0,
                      p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);
  sf=(s=perm->array)+perm->size; do *s=c++; while(++s<sf);

  /* Set the proper qsort function. */
  if(p->descending)
    switch(p->sortcol->type)
      {
      case GAL_TYPE_UINT8:   qsortfn=gal_qsort_index_single_uint8_d;   break;
      case GAL_TYPE_INT8:    qsortfn=gal_qsort_index_single_int8_d;    break;
      case GAL_TYPE_UINT16:  qsortfn=gal_qsort_index_single_uint16_d;  break;
      case GAL_TYPE_INT16:   qsortfn=gal_qsort_index_single_int16_d;   break;
      case GAL_TYPE_UINT32:  qsortfn=gal_qsort_index_single_uint32_d;  break;
      case GAL_TYPE_INT32:   qsortfn=gal_qsort_index_single_int32_d;   break;
      case GAL_TYPE_UINT64:  qsortfn=gal_qsort_index_single_uint64_d;  break;
      case GAL_TYPE_INT64:   qsortfn=gal_qsort_index_single_int64_d;   break;
      case GAL_TYPE_FLOAT32: qsortfn=gal_qsort_index_single_float32_d; break;
      case GAL_TYPE_FLOAT64: qsortfn=gal_qsort_index_single_float64_d; break;
      default:
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
              "the problem. The code `%u' wasn't recognized as a data type",
              __func__, PACKAGE_BUGREPORT, p->sortcol->type);
      }
  else
    switch(p->sortcol->type)
      {
      case GAL_TYPE_UINT8:   qsortfn=gal_qsort_index_single_uint8_i;   break;
      case GAL_TYPE_INT8:    qsortfn=gal_qsort_index_single_int8_i;    break;
      case GAL_TYPE_UINT16:  qsortfn=gal_qsort_index_single_uint16_i;  break;
      case GAL_TYPE_INT16:   qsortfn=gal_qsort_index_single_int16_i;   break;
      case GAL_TYPE_UINT32:  qsortfn=gal_qsort_index_single_uint32_i;  break;
      case GAL_TYPE_INT32:   qsortfn=gal_qsort_index_single_int32_i;   break;
      case GAL_TYPE_UINT64:  qsortfn=gal_qsort_index_single_uint64_i;  break;
      case GAL_TYPE_INT64:   qsortfn=gal_qsort_index_single_int64_i;   break;
      case GAL_TYPE_FLOAT32: qsortfn=gal_qsort_index_single_float32_i; break;
      case GAL_TYPE_FLOAT64: qsortfn=gal_qsort_index_single_float64_i; break;
      default:
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
              "the problem. The code `%u' wasn't recognized as a data type",
              __func__, PACKAGE_BUGREPORT, p->sortcol->type);
      }

  /* Sort the indexs from the values. */
  gal_qsort_index_single=p->sortcol->array;
  qsort(perm->array, perm->size, sizeof *s, qsortfn);

  /* For a check (only on float32 type `sortcol'):
  {
    float *f=p->sortcol->array;
    sf=(s=perm->array)+perm->size;
    do printf("%f\n", f[*s]); while(++s<sf);
    exit(0);
  }
  */

  /* Sort all the output columns with this permutation. */
  table_apply_permutation(p->table, perm->array, perm->size, 0);

  /* Clean up. */
  gal_data_free(perm);
  if(p->freesort) gal_data_free(p->sortcol);
}





static void
table_head_tail(struct tableparams *p)
{
  char **strarr;
  gal_data_t *col;
  size_t i, start, end;

  /* Go over all the columns and make the necessary corrections. */
  for(col=p->table;col!=NULL;col=col->next)
    {
      /* If we are dealing with strings, we'll need to free the strings
         that the columns that will not be used point to (outside the
         allocated array directly `gal_data_t'). We don't have to worry
         about the space for the actual pointers (they will be free'd by
         `free' in any case, since they are in the initially allocated
         array).*/
      if(col->type==GAL_TYPE_STRING)
        {
          /* Set the start and ending indexs. */
          start = p->head!=GAL_BLANK_SIZE_T ? p->head        : 0;
          end   = p->head!=GAL_BLANK_SIZE_T ? p->table->size : p->tail;

          /* Free their allocated spaces. */
          strarr=col->array;
          for(i=start; i<end; ++i) { free(strarr[i]); strarr[i]=NULL; }
        }

      /* For `--tail', we'll need to bring the last columns to the
         start. Note that we are using `memmove' because we want to be safe
         with overlap. */
      if(p->tail!=GAL_BLANK_SIZE_T)
        memmove(col->array,
                gal_pointer_increment(col->array, col->size - p->tail,
                                      col->type),
                p->tail*gal_type_sizeof(col->type));

      /* In any case (head or tail), the new number of column elements is
         the given value. */
      col->size = col->dsize[0] = ( p->head!=GAL_BLANK_SIZE_T
                                    ? p->head
                                    : p->tail );
    }
}





/**************************************************************/
/***************       Top function         *******************/
/**************************************************************/
void
table(struct tableparams *p)
{
  /* Apply a certain range (if required) to the output sample. */
  if(p->range) table_range(p);

  /* Sort it (if required). */
  if(p->sort) table_sort(p);

  /* If the output number of rows is limited, apply them. */
  if(p->head!=GAL_BLANK_SIZE_T || p->tail!=GAL_BLANK_SIZE_T)
    table_head_tail(p);

  /* If any operations are needed, do them. */
  if(p->outcols)
    arithmetic_operate(p);

  /* Write the output. */
  gal_table_write(p->table, NULL, p->cp.tableformat, p->cp.output,
                  "TABLE", p->colinfoinstdout);
}
