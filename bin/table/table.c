/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <gsl/gsl_heapsort.h>

#include <gnuastro/txt.h>
#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/qsort.h>
#include <gnuastro/pointer.h>
#include <gnuastro/polygon.h>
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





static gal_data_t *
table_selection_range(struct tableparams *p, gal_data_t *col)
{
  size_t one=1;
  double *darr;
  int numok=GAL_ARITHMETIC_NUMOK;
  int inplace=GAL_ARITHMETIC_INPLACE;
  gal_data_t *min=NULL, *max=NULL, *tmp, *ltmin, *gemax=NULL;

  /* First, make sure everything is OK. */
  if(p->range==NULL)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us to fix the "
          "problem at %s. 'p->range' should not be NULL at this point",
          __func__, PACKAGE_BUGREPORT);

  /* Allocations. */
  min=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one, NULL, 0, -1, 1,
                     NULL, NULL, NULL);
  max=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one, NULL, 0, -1, 1,
                     NULL, NULL, NULL);

  /* Read the range of values for this column. */
  darr=p->range->array;
  ((double *)(min->array))[0] = darr[0];
  ((double *)(max->array))[0] = darr[1];

  /* Move 'p->range' to the next element in the list and free the current
     one (we have already read its values and don't need it any more). */
  tmp=p->range;
  p->range=p->range->next;
  gal_data_free(tmp);

  /* Find all the elements outside this range (smaller than the minimum,
     larger than the maximum or blank) as separate binary flags.. */
  ltmin=gal_arithmetic(GAL_ARITHMETIC_OP_LT, 1, numok, col, min);
  gemax=gal_arithmetic(GAL_ARITHMETIC_OP_GE, 1, numok, col, max);

  /* Merge them both into one array. */
  ltmin=gal_arithmetic(GAL_ARITHMETIC_OP_OR, 1, inplace, ltmin, gemax);

  /* For a check.
  {
    size_t i;
    uint8_t *u=ltmin->array;
    for(i=0;i<ltmin->size;++i) printf("%zu: %u\n", i, u[i]);
    exit(0);
  }
  */

  /* Clean up and return. */
  gal_data_free(gemax);
  gal_data_free(min);
  gal_data_free(max);
  return ltmin;
}





/* Read column value of any type as a double for the polygon options. */
static double
selection_polygon_read_point(gal_data_t *col, size_t i)
{
  /* Check and assign the points to the points array. */
  switch(col->type)
    {
    case GAL_TYPE_INT8:    return (( int8_t   *)col->array)[i];
    case GAL_TYPE_UINT8:   return (( uint8_t  *)col->array)[i];
    case GAL_TYPE_UINT16:  return (( uint16_t *)col->array)[i];
    case GAL_TYPE_INT16:   return (( int16_t  *)col->array)[i];
    case GAL_TYPE_UINT32:  return (( uint32_t *)col->array)[i];
    case GAL_TYPE_INT32:   return (( int32_t  *)col->array)[i];
    case GAL_TYPE_UINT64:  return (( uint64_t *)col->array)[i];
    case GAL_TYPE_INT64:   return (( int64_t  *)col->array)[i];
    case GAL_TYPE_FLOAT32: return (( float    *)col->array)[i];
    case GAL_TYPE_FLOAT64: return (( double   *)col->array)[i];
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, col->type);
    }

  /* Control should not reach here. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NAN;
}





/* Mask the rows that are not in the given polygon. */
static gal_data_t *
table_selection_polygon(struct tableparams *p, gal_data_t *col1,
                        gal_data_t *col2, int in1out0)
{
  uint8_t *oarr;
  double point[2];
  gal_data_t *out=NULL;
  size_t i, psize=p->polygon->size/2;

  /* Allocate the output array: This array will have a '0' for the points
     which are inside the polygon and '1' for those that are outside of it
     (to be masked/removed from the input). */
  out=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, col1->dsize, NULL, 0, -1, 1,
                     NULL, NULL, NULL);
  oarr=out->array;

  /* Loop through all the rows in the given columns and check the points.*/
  for(i=0; i<col1->size; i++)
    {
      /* Read the column values as a double. */
      point[0]=selection_polygon_read_point(col1, i);
      point[1]=selection_polygon_read_point(col2, i);

      /* For '--inpolygon', if point is inside polygon, put 0, otherwise
         1. Note that we are building a mask for the rows that must be
         discarded, so we want '1' for the points we don't want. */
      oarr[i] = (in1out0
                 ? !gal_polygon_is_inside(p->polygon->array, point, psize)
                 :  gal_polygon_is_inside(p->polygon->array, point, psize));

      /* For a check
      printf("(%f,%f): %s, %u\n", point[0], point[1], oarr[i]);
      */
    }

  /* Return the output column. */
  return out;
}





/* Given a string dataset and a single string, return a 'uint8_t' array
   with the same size as the string dataset that has a '1' for all the
   elements that are equal. */
static gal_data_t *
table_selection_string_eq_ne(gal_data_t *column, char *reference, int e0n1)
{
  gal_data_t *out;
  uint8_t *oarr, comp;
  size_t i, size=column->size;
  char **strarr=column->array;

  /* Allocate the output binary dataset. */
  out=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &size, NULL, 0, -1, 1,
                     NULL, NULL, NULL);
  oarr=out->array;

  /* Parse the values and mark the outputs IN THE OPPOSITE manner (we are
     marking the ones that must be removed). */
  for(i=0;i<size;++i)
    {
      comp=strcmp(strarr[i], reference);
      oarr[i] = e0n1 ? (comp==0) : (comp!=0);
    }

  /* Return. */
  return out;
}





static gal_data_t *
table_selection_equal_or_notequal(struct tableparams *p, gal_data_t *col,
                                  int e0n1)
{
  void *varr;
  char **strarr;
  size_t i, one=1;
  int numok=GAL_ARITHMETIC_NUMOK;
  int inplace=GAL_ARITHMETIC_INPLACE;
  gal_data_t *eq, *out=NULL, *value=NULL;
  gal_data_t *arg = e0n1 ? p->notequal : p->equal;

  /* Note that this operator is used to make the "masked" array, so when
     'e0n1==0' the operator should be 'GAL_ARITHMETIC_OP_NE' and
     vice-versa.

     For the merging with other elements, when 'e0n1==0', we need the
     'GAL_ARITHMETIC_OP_AND', but for 'e0n1==1', it should be 'OR'. */
  int mergeop  = e0n1 ? GAL_ARITHMETIC_OP_OR : GAL_ARITHMETIC_OP_AND;
  int operator = e0n1 ? GAL_ARITHMETIC_OP_EQ : GAL_ARITHMETIC_OP_NE;

  /* First, make sure everything is OK. */
  if(arg==NULL)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us to fix the "
          "problem at %s. 'p->range' should not be NULL at this point",
          __func__, PACKAGE_BUGREPORT);

  /* To easily parse the given values. */
  strarr=arg->array;

  /* Go through the values given to this call of the option and flag the
     elements. */
  for(i=0;i<arg->size;++i)
    {
      /* Write the value  */
      if(col->type==GAL_TYPE_STRING)
        eq=table_selection_string_eq_ne(col, strarr[i], e0n1);
      else
        {
          /* Allocate the value dataset. */
          value=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one, NULL, 0, -1, 1,
                               NULL, NULL, NULL);
          varr=value->array;

          /* Read the stored string as a float64. */
          if( gal_type_from_string(&varr, strarr[i], GAL_TYPE_FLOAT64) )
            {
              fprintf(stderr, "%s couldn't be read as a number.\n", strarr[i]);
              exit(EXIT_FAILURE);
            }

          /* Mark the rows that are equal (irrespective of the column's
             original numerical datatype). */
          eq=gal_arithmetic(operator, 1, numok, col, value);
        }

      /* Merge the results with (possible) previous results. */
      if(out)
        {
          out=gal_arithmetic(mergeop, 1, inplace, out, eq);
          gal_data_free(eq);
        }
      else
        out=eq;
    }

  /* For a check.
  {
    uint8_t *u=out->array;
    for(i=0;i<out->size;++i) printf("%zu: %u\n", i, u[i]);
    exit(0);
  }
  */


  /* Move the main pointer to the next possible call of the given
     option. Note that 'arg' already points to 'p->equal' or 'p->notequal',
     so it will automatically be freed with the next step.*/
  if(e0n1) p->notequal=p->notequal->next;
  else     p->equal=p->equal->next;

  /* Clean up and return. */
  gal_data_free(value);
  gal_data_free(arg);
  return out;
}





static void
table_selection(struct tableparams *p)
{
  uint8_t *u;
  struct list_select *tmp;
  gal_data_t *mask, *addmask=NULL;
  gal_data_t *sum, *perm, *blmask;
  size_t i, g, b, *s, *sf, ngood=0;
  int inplace=GAL_ARITHMETIC_INPLACE;

  /* Allocate datasets for the necessary numbers and write them in. */
  perm=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, p->table->dsize, NULL, 0,
                      p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);
  mask=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, p->table->dsize, NULL, 1,
                      p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);

  /* Go over each selection criteria and remove the necessary elements. */
  for(tmp=p->selectcol;tmp!=NULL;tmp=tmp->next)
    {
      switch(tmp->type)
        {
        case SELECT_TYPE_RANGE:
          addmask=table_selection_range(p, tmp->col);
          break;

        /* '--inpolygon' and '--outpolygon' need two columns. */
        case SELECT_TYPE_INPOLYGON:
        case SELECT_TYPE_OUTPOLYGON:
          addmask=table_selection_polygon(p, tmp->col, tmp->next->col,
                                          tmp->type==SELECT_TYPE_INPOLYGON);
          tmp=tmp->next;
          break;

        case SELECT_TYPE_EQUAL:
          addmask=table_selection_equal_or_notequal(p, tmp->col, 0);
          break;

        case SELECT_TYPE_NOTEQUAL:
          addmask=table_selection_equal_or_notequal(p, tmp->col, 1);
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s "
                "to fix the problem. The code %d is not a recognized "
                "range identifier", __func__, PACKAGE_BUGREPORT,
                tmp->type);
        }

      /* Remove any blank elements. */
      if(gal_blank_present(tmp->col, 1))
        {
          blmask = gal_arithmetic(GAL_ARITHMETIC_OP_ISBLANK, 1, 0, tmp->col);
          addmask=gal_arithmetic(GAL_ARITHMETIC_OP_OR, 1, inplace,
                                 addmask, blmask);
          gal_data_free(blmask);
        }

      /* Add this mask array to the cumulative mask array (of all
         selections). */
      mask=gal_arithmetic(GAL_ARITHMETIC_OP_OR, 1, inplace, mask, addmask);

      /* For a check.
         {
           float *f=ref->array;
           uint8_t *m=mask->array;
           uint8_t *u=addmask->array, *uf=u+addmask->size;
           printf("\n\nInput column: %s\n", ref->name ? ref->name : "No Name");
           printf("Range: %g, %g\n", rarr[0], rarr[1]);
           printf("%-20s%-20s%-20s\n", "Value", "This mask",
           "Including previous");
           do printf("%-20f%-20u%-20u\n", *f++, *u++, *m++); while(u<uf);
           exit(0);
           }
        */

      /* Final clean up. */
      gal_data_free(addmask);
    }

  /* Find the final number of elements to print. */
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
  for(tmp=p->selectcol;tmp!=NULL;tmp=tmp->next)
    { if(p->freeselect[i]) {gal_data_free(tmp->col); tmp->col=NULL;} ++i; }
  ui_list_select_free(p->selectcol, 0);
  gal_data_free(mask);
  gal_data_free(perm);
  free(p->freeselect);
  gal_data_free(sum);
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

  /* For string columns, print a descriptive message. Note that some FITS
     tables were found that do actually have numbers stored in string
     types! */
  if(p->sortcol->type==GAL_TYPE_STRING)
    error(EXIT_FAILURE, 0, "sort column has a string type, but it can "
          "(currently) only work on numbers.\n\n"
          "TIP: if you know the columns contents are all numbers that are "
          "just stored as strings, you can use this program to save the "
          "table as a text file, modify the column meta-data (for example "
          "to type 'i32' or 'f32' instead of 'strN'), then use this "
          "program again to save it as a FITS table.\n\n"
          "For more on column metadata in plain text format, please run "
          "the following command (or see the 'Gnuastro text table format "
          "section of the book/manual):\n\n"
          "    $ info gnuastro \"gnuastro text table format\"");

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
              "the problem. The code '%u' wasn't recognized as a data type",
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
              "the problem. The code '%u' wasn't recognized as a data type",
              __func__, PACKAGE_BUGREPORT, p->sortcol->type);
      }

  /* Sort the indexs from the values. */
  gal_qsort_index_single=p->sortcol->array;
  qsort(perm->array, perm->size, sizeof *s, qsortfn);

  /* For a check (only on float32 type 'sortcol'):
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
         allocated array directly 'gal_data_t'). We don't have to worry
         about the space for the actual pointers (they will be free'd by
         'free' in any case, since they are in the initially allocated
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

      /* For '--tail', we'll need to bring the last columns to the
         start. Note that we are using 'memmove' because we want to be safe
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





/*This function concatenates two table column wise .
 It  attaches catcolumn table at the back of first table */
static void
table_catcolumn(struct tableparams *p)
{
  char *hdu=NULL;
  gal_data_t *tocat, *final;
  gal_list_str_t *filell, *hdull;
  struct gal_options_common_params *cp=&p->cp;

  /* Go over all the given files. */
  hdull=p->catcolhdu;
  for(filell=p->catcolumn; filell!=NULL; filell=filell->next)
    {
      /* Set the HDU (not necessary for non-FITS tables). */
      if(gal_fits_name_is_fits(filell->v))
        {
          if(hdull) { hdu=hdull->v; hdull=hdull->next; }
          else
            error(EXIT_FAILURE, 0, "not enough '--catcolhdu's. For every "
                  "FITS table given to '--catcolumn', a call to "
                  "'--catcolhdu' is necessary to identify its HDU/extension");
        }
      else hdu=NULL;

      /* Read the catcolumn table. */
      tocat=gal_table_read(filell->v, hdu, NULL, NULL, cp->searchin,
                           cp->ignorecase, cp->minmapsize, p->cp.quietmmap,
                           NULL);

      /* Check the number of rows. */
      if(tocat->dsize[0]!=p->table->dsize[0])
        error(EXIT_FAILURE, 0, "%s: incorrect number of rows. The table given "
              "to '--catcolumn' must have the same number of rows as the "
              "main argument (after all row-selections have been applied), "
              "but they have %zu and %zu rows respectively",
              gal_fits_name_save_as_string(filell->v, hdu), tocat->dsize[0],
              p->table->dsize[0]);

      /* Find the final column of the main table and add this table.*/
      final=gal_list_data_last(p->table);
      final->next=tocat;
    }
}




















/**************************************************************/
/***************       Top function         *******************/
/**************************************************************/
void
table(struct tableparams *p)
{
  /* Apply a certain range (if required) to the output sample. */
  if(p->selection) table_selection(p);

  /* Sort it (if required). */
  if(p->sort) table_sort(p);

  /* If the output number of rows is limited, apply them. */
  if(p->head!=GAL_BLANK_SIZE_T || p->tail!=GAL_BLANK_SIZE_T)
    table_head_tail(p);

  /* If any operations are needed, do them. */
  if(p->outcols)
    arithmetic_operate(p);

  /* Concatenate the columns of tables(if required)*/
  if(p->catcolumn) table_catcolumn(p);

  /* Write the output. */
  gal_table_write(p->table, NULL, p->cp.tableformat, p->cp.output,
                  "TABLE", p->colinfoinstdout);
}
