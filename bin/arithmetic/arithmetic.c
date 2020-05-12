/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/array.h>
#include <gnuastro/binary.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/interpolate.h>

#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "operands.h"
#include "arithmetic.h"














/***************************************************************/
/*************          Internal functions         *************/
/***************************************************************/
#define SET_NUM_OP(CTYPE) {                                \
    CTYPE a=*(CTYPE *)(numpop->array); if(a>0) return a;    }

static size_t
pop_number_of_operands(struct arithmeticparams *p, int op, char *token_string,
                       gal_data_t **params)
{
  char *cstring="first";
  size_t c, numparams=0;
  gal_data_t *tmp, *numpop;

  /* See if this operator needs any parameters. If so, pop them. */
  switch(op)
    {
    case GAL_ARITHMETIC_OP_QUANTILE:
      numparams=1;
      break;
    case GAL_ARITHMETIC_OP_SIGCLIP_STD:
    case GAL_ARITHMETIC_OP_SIGCLIP_MEAN:
    case GAL_ARITHMETIC_OP_SIGCLIP_MEDIAN:
    case GAL_ARITHMETIC_OP_SIGCLIP_NUMBER:
      numparams=2;
      break;
    }

  /* If any parameters should be read, read them. */
  *params=NULL;
  for(c=0;c<numparams;++c)
    {
      /* If it only has one element, save it as floating point and add it
         to the list. */
      tmp=operands_pop(p, token_string);
      if(tmp->size>1)
        error(EXIT_FAILURE, 0, "the %s popped operand of the '%s' "
              "operator must be a single number", cstring, token_string);
      tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
      gal_list_data_add(params, tmp);

      /* A small sanity check (none of the parameters for sigma-clipping,
         or quantile estimation can be negative. */
      if( ((float *)(tmp->array))[0]<=0.0 )
        error(EXIT_FAILURE, 0, "the %s popped operand of the '%s' "
              "operator cannot be negative", cstring, token_string);

      /* Increment the counter string. */
      cstring=c?"third":"second";
    }

  /* Check if its a number. */
  numpop=operands_pop(p, token_string);
  if(numpop->size>1)
    error(EXIT_FAILURE, 0, "the %s popped operand of the '%s' "
          "operator (number of input datasets) must be a number, not an "
          "array", cstring, token_string);

  /* Check its type and return the value. */
  switch(numpop->type)
    {
    /* For the integer types, if they are unsigned, then just pass their
       value, if they are signed, you have to make sure they are zero or
       positive. */
    case GAL_TYPE_UINT8:   SET_NUM_OP(uint8_t);     break;
    case GAL_TYPE_INT8:    SET_NUM_OP(int8_t);      break;
    case GAL_TYPE_UINT16:  SET_NUM_OP(uint16_t);    break;
    case GAL_TYPE_INT16:   SET_NUM_OP(int16_t);     break;
    case GAL_TYPE_UINT32:  SET_NUM_OP(uint32_t);    break;
    case GAL_TYPE_INT32:   SET_NUM_OP(int32_t);     break;
    case GAL_TYPE_UINT64:  SET_NUM_OP(uint64_t);    break;
    case GAL_TYPE_INT64:   SET_NUM_OP(int64_t);     break;

    /* Floating point numbers are not acceptable in this context. */
    case GAL_TYPE_FLOAT32:
    case GAL_TYPE_FLOAT64:
      error(EXIT_FAILURE, 0, "the %s popped operand of the '%s' "
            "operator (number of input datasets) must be an integer type",
            cstring, token_string);

    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, numpop->type);
    }

  /* If control reaches here, then the number must have been a negative
     value, so print an error. */
  error(EXIT_FAILURE, 0, "the %s popped operand of the '%s' operator "
        "cannot be zero or a negative number", cstring,
        token_string);
  return 0;
}




















/**********************************************************************/
/****************         Filtering operators         *****************/
/**********************************************************************/
#define ARITHMETIC_FILTER_DIM 10

struct arithmetic_filter_p
{
  int           operator;       /* The type of filtering.                */
  size_t          *fsize;       /* Filter size.                          */
  size_t        *hpfsize;       /* Positive Half-filter size.            */
  size_t        *hnfsize;       /* Negative Half-filter size.            */
  float     sclip_multip;       /* Sigma multiple in sigma-clipping.     */
  float      sclip_param;       /* Termination critera in sigma-cliping. */
  gal_data_t      *input;       /* Input dataset.                        */
  gal_data_t        *out;       /* Output dataset.                       */

  int           hasblank;       /* If the dataset has blank values.      */
};





/* Main filtering work function. */
static void *
arithmetic_filter(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct arithmetic_filter_p *afp=(struct arithmetic_filter_p *)tprm->params;
  gal_data_t *input=afp->input;

  size_t sind=-1;
  size_t ind, index, one=1;
  gal_data_t *sigclip, *result=NULL;
  size_t *hpfsize=afp->hpfsize, *hnfsize=afp->hnfsize;
  size_t *tsize, *dsize=input->dsize, *fsize=afp->fsize;
  size_t i, j, coord[ARITHMETIC_FILTER_DIM], ndim=input->ndim;
  size_t start[ARITHMETIC_FILTER_DIM], end[ARITHMETIC_FILTER_DIM];
  gal_data_t *tile=gal_data_alloc(NULL, input->type, ndim, afp->fsize, NULL,
                                  0, -1, 1, NULL, NULL, NULL);

  /* Prepare the tile. */
  free(tile->array);
  tsize=tile->dsize;
  tile->block=input;


  /* Go over all the pixels that were assigned to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading, put the index in 'ind'. */
      ind=tprm->indexs[i];

      /* Get the coordinate of the pixel. */
      gal_dimension_index_to_coord(ind, ndim, dsize, coord);

      /* See which dimensions need trimming. */
      tile->size=1;
      for(j=0;j<ndim;++j)
        {
          /* Estimate the coordinate of the filter's starting point. Note
             that we are dealing with size_t (unsigned int) type here, so
             there are no negatives. A negative result will produce an
             extremely large number, so instead of checking for negative,
             we can just see if the result of a subtraction is less than
             the width of the input. */
          if( (coord[j] - hnfsize[j] > dsize[j])
              || (coord[j] + hpfsize[j] >= dsize[j]) )
            {
              start[j] = ( (coord[j] - hnfsize[j] > dsize[j])
                           ? 0 : coord[j] - hnfsize[j] );
              end[j]   = ( (coord[j] + hpfsize[j] >= dsize[j])
                           ? dsize[j]
                           : coord[j] + hpfsize[j] + 1);
              tsize[j] = end[j] - start[j];
            }
          else  /* NOT on the edge (given requested filter width). */
            {
              tsize[j] = fsize[j];
              start[j] = coord[j] - hnfsize[j];
            }
          tile->size *= tsize[j];
        }

      /* For a test.
         printf("coord: %zu, %zu\n", coord[1]+1, coord[0]+1);
         printf("\tstart: %zu, %zu\n", start[1]+1, start[0]+1);
         printf("\ttsize: %zu, %zu\n", tsize[1], tsize[0]);
      */

      /* Set the tile's starting pointer. */
      index=gal_dimension_coord_to_index(ndim, dsize, start);
      tile->array=gal_pointer_increment(input->array, index, input->type);

      /* Do the necessary calculation. */
      switch(afp->operator)
        {
        case ARITHMETIC_OP_FILTER_MEDIAN:
          result=gal_statistics_median(tile, 0);
          break;


        case ARITHMETIC_OP_FILTER_MEAN:
          result=gal_statistics_mean(tile);
          break;


        case ARITHMETIC_OP_FILTER_SIGCLIP_MEAN:
        case ARITHMETIC_OP_FILTER_SIGCLIP_MEDIAN:
          /* Find the sigma-clipped results. */
          sigclip=gal_statistics_sigma_clip(tile, afp->sclip_multip,
                                            afp->sclip_param, 0, 1);

          /* Set the required index. */
          switch(afp->operator)
            {
            case ARITHMETIC_OP_FILTER_SIGCLIP_MEAN:   sind = 2; break;
            case ARITHMETIC_OP_FILTER_SIGCLIP_MEDIAN: sind = 1; break;
            default:
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at "
                    "%s to fix the problem. The 'afp->operator' value "
                    "%d is not recognized as sigma-clipped median or "
                    "mean", __func__, PACKAGE_BUGREPORT, afp->operator);
            }

          /* Allocate the output and write the value into it. */
          result=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL,
                                0, -1, 1, NULL, NULL, NULL);
          ((float *)(result->array))[0] =
            ((float *)(sigclip->array))[sind];

          /* Clean up. */
          gal_data_free(sigclip);
          break;


        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s "
                "to fix the problem. 'afp->operator' code %d is not "
                "recognized", PACKAGE_BUGREPORT, __func__,
                afp->operator);
        }

      /* Make sure the output array type and result's type are the
         same. */
      if(result->type!=afp->out->type)
        result=gal_data_copy_to_new_type_free(result, afp->out->type);


      /* Copy the result into the output array. */
      memcpy(gal_pointer_increment(afp->out->array, ind, afp->out->type),
             result->array, gal_type_sizeof(afp->out->type));

      /* Clean up for this pixel. */
      gal_data_free(result);
    }


  /* Clean up for this thread. */
  tile->array=NULL;
  tile->block=NULL;
  gal_data_free(tile);


  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





static void
wrapper_for_filter(struct arithmeticparams *p, char *token, int operator)
{
  int type=GAL_TYPE_INVALID;
  size_t i=0, ndim, nparams, one=1;
  struct arithmetic_filter_p afp={0};
  size_t fsize[ARITHMETIC_FILTER_DIM];
  gal_data_t *tmp, *tmp2, *zero, *comp, *params_list=NULL;
  size_t hnfsize[ARITHMETIC_FILTER_DIM], hpfsize[ARITHMETIC_FILTER_DIM];
  int issigclip=(operator==ARITHMETIC_OP_FILTER_SIGCLIP_MEAN
                 || operator==ARITHMETIC_OP_FILTER_SIGCLIP_MEDIAN);


  /* Get the input's number of dimensions. */
  afp.input=operands_pop(p, token);
  afp.operator=operator;
  ndim=afp.input->ndim;
  afp.hnfsize=hnfsize;
  afp.hpfsize=hpfsize;
  afp.fsize=fsize;


  /* A small sanity check. */
  if(ndim>ARITHMETIC_FILTER_DIM)
    error(EXIT_FAILURE, 0, "%s: currently only datasets with less than "
          "%d dimensions are acceptable. The input has %zu dimensions",
          __func__, ARITHMETIC_FILTER_DIM, ndim);


  /* A zero value for checking the value of input widths. */
  zero=gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &one, NULL, 1, -1, 1, NULL,
                      NULL, NULL);


  /* Based on the first popped operand's dimensions and the operator, of
     pop the necessary number of operands. */
  nparams = ndim + (issigclip ? 2 : 0 );
  for(i=0;i<nparams;++i)
    gal_list_data_add(&params_list, operands_pop(p, token));


  /* Make sure the parameters only have single values. */
  i=0;
  for(tmp=params_list; tmp!=NULL; tmp=tmp->next)
    {
      ++i;
      if(tmp->size!=1)
        error(EXIT_FAILURE, 0, "the parameters given to the filtering "
              "operators can only be numbers. Value number %zu has %zu "
              "elements, so its an array", i, tmp->size);
    }


  /* If this is a sigma-clipping filter, the top two operands are the
     sigma-clipping parameters. */
  if(issigclip)
    {
      /* Read the sigma-clipping multiple (first element in the list). */
      tmp=gal_list_data_pop(&params_list);
      tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
      afp.sclip_multip=*(float *)(tmp->array);
      gal_data_free(tmp);

      /* Read the sigma-clipping termination parameter. */
      tmp=gal_list_data_pop(&params_list);
      tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
      afp.sclip_param=*(float *)(tmp->array);
      gal_data_free(tmp);
    }


  /* If the input only has one element, filtering makes no sense, so don't
     waste time, just add the input onto the stack. */
  if(afp.input->size==1) afp.out=afp.input;
  else
    {
      /* Allocate an array for the size of the filter and fill it in. The
         values must be written in the inverse order since the user gives
         dimensions with the FITS standard. */
      i=ndim-1;
      for(tmp=params_list; tmp!=NULL; tmp=tmp->next)
        {
          /* Make sure the user has given an integer type. */
          if(tmp->type==GAL_TYPE_FLOAT32 || tmp->type==GAL_TYPE_FLOAT64)
            error(EXIT_FAILURE, 0, "lengths of filter along dimensions "
                  "must be integer values, not floats. The given length "
                  "along dimension %zu is a float", ndim-i);

          /* Make sure it isn't negative. */
          comp=gal_arithmetic(GAL_ARITHMETIC_OP_GT, 1, 0, tmp, zero);
          if( *(uint8_t *)(comp->array) == 0 )
            error(EXIT_FAILURE, 0, "lengths of filter along dimensions "
                  "must be positive. The given length in dimension %zu"
                  "is either zero or negative", ndim-i);
          gal_data_free(comp);

          /* Convert the input into size_t and put it into the array that
             keeps the filter size. */
          tmp2=gal_data_copy_to_new_type(tmp, GAL_TYPE_SIZE_T);
          fsize[ i ] = *(size_t *)(tmp2->array);
          gal_data_free(tmp2);

          /* If the width is larger than the input's size, change the width
             to the input's size. */
          if( fsize[i] > afp.input->dsize[i] )
            error(EXIT_FAILURE, 0, "%s: the filter size along dimension %zu "
                  "(%zu) is greater than the input's length in that "
                  "dimension (%zu)", __func__, i, fsize[i],
                  afp.input->dsize[i]);

          /* Go onto the previous dimension. */
          --i;
        }


      /* Set the half filter sizes. Note that when the size is an odd
         number, the number of pixels before and after the actual pixel are
         equal, but for an even number, we will look into one element more
         when looking before than the ones after. */
      for(i=0;i<ndim;++i)
        {
          if( fsize[i]%2 )
            hnfsize[i]=hpfsize[i]=fsize[i]/2;
          else
            { hnfsize[i]=fsize[i]/2; hpfsize[i]=fsize[i]/2-1; }
        }

      /* For a test.
      printf("fsize: %zu, %zu\n", fsize[0], fsize[1]);
      printf("hnfsize: %zu, %zu\n", hnfsize[0], hnfsize[1]);
      printf("hpfsize: %zu, %zu\n", hpfsize[0], hpfsize[1]);
      */

      /* See if the input has blank pixels. */
      afp.hasblank=gal_blank_present(afp.input, 1);


      /* Set the type of the output dataset. */
      switch(operator)
        {
        case ARITHMETIC_OP_FILTER_MEDIAN:
        case ARITHMETIC_OP_FILTER_SIGCLIP_MEDIAN:
          type=afp.input->type;
          break;

        case ARITHMETIC_OP_FILTER_MEAN:
        case ARITHMETIC_OP_FILTER_SIGCLIP_MEAN:
          type=GAL_TYPE_FLOAT64;
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to fix "
                "the problem. The 'operator' code %d is not recognized",
                PACKAGE_BUGREPORT, __func__, operator);
        }


      /* Allocate the output dataset. Note that filtering doesn't change
         the units of the dataset. */
      afp.out=gal_data_alloc(NULL, type, ndim, afp.input->dsize,
                             afp.input->wcs, 0, afp.input->minmapsize,
                             afp.input->quietmmap, NULL, afp.input->unit,
                             NULL);


      /* Spin off threads for each pixel. */
      gal_threads_spin_off(arithmetic_filter, &afp, afp.input->size,
                           p->cp.numthreads);
    }


  /* Add the output to the top of the stack. */
  operands_add(p, NULL, afp.out);

  /* Clean up and add the output on top of the stack */
  gal_data_free(zero);
  gal_data_free(afp.input);
  gal_list_data_free(params_list);
}




















/***************************************************************/
/*************            Other functions          *************/
/***************************************************************/
static int
arithmetic_binary_sanity_checks(gal_data_t *in, gal_data_t *conn,
                                char *operator)
{
  int conn_int;

  /* Do proper sanity checks on 'conn'. */
  if(conn->size!=1)
    error(EXIT_FAILURE, 0, "the first popped operand to '%s' must be a "
          "single number. However, it has %zu elements", operator,
          conn->size);
  if(conn->type==GAL_TYPE_FLOAT32 || conn->type==GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "the first popped operand to '%s' is the "
          "connectivity (a value between 1 and the number of dimensions) "
          "therefore, it must NOT be a floating point", operator);

  /* Convert the connectivity value to a 32-bit integer and read it in and
     make sure it is not larger than the number of dimensions. */
  conn=gal_data_copy_to_new_type_free(conn, GAL_TYPE_INT32);
  conn_int = *((int32_t *)(conn->array));
  if(conn_int>in->ndim)
    error(EXIT_FAILURE, 0, "the first popped operand of '%s' (%d) is "
          "larger than the number of dimensions in the second-popped "
          "operand (%zu)", operator, conn_int, in->ndim);

  /* Make sure the array has an unsigned 8-bit type. */
  if(in->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "the second popped operand of '%s' doesn't "
          "have an 8-bit unsigned integer type. It must be a binary "
          "dataset (only being equal to zero is checked). You can use "
          "the 'uint8' operator for type conversion", operator);

  /* Clean up and return the integer value of 'conn'. */
  gal_data_free(conn);
  return conn_int;
}





static void
arithmetic_erode_dilate(struct arithmeticparams *p, char *token, int op)
{
  int conn_int;

  /* Pop the two necessary operands. */
  gal_data_t *conn = operands_pop(p, token);
  gal_data_t *in   = operands_pop(p, token);

  /* Do the sanity checks. */
  conn_int=arithmetic_binary_sanity_checks(in, conn, token);

  /* Do the operation. */
  switch(op)
    {
    case ARITHMETIC_OP_ERODE:  gal_binary_erode(in,  1, conn_int, 1); break;
    case ARITHMETIC_OP_DILATE: gal_binary_dilate(in, 1, conn_int, 1); break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
            "problem. The operator code %d not recognized", __func__,
            PACKAGE_BUGREPORT, op);
    }

  /* Push the result onto the stack. */
  operands_add(p, NULL, in);
}





static void
arithmetic_connected_components(struct arithmeticparams *p, char *token)
{
  int conn_int;
  gal_data_t *out=NULL;

  /* Pop the two necessary operands. */
  gal_data_t *conn = operands_pop(p, token);
  gal_data_t *in   = operands_pop(p, token);

  /* Basic sanity checks. */
  conn_int=arithmetic_binary_sanity_checks(in, conn, token);

  /* Do the connected components labeling. */
  gal_binary_connected_components(in, &out, conn_int);

  /* Push the result onto the stack. */
  operands_add(p, NULL, out);

  /* Clean up ('conn' was freed in the sanity check). */
  gal_data_free(in);
}





static void
arithmetic_fill_holes(struct arithmeticparams *p, char *token)
{
  int conn_int;

  /* Pop the two necessary operands. */
  gal_data_t *conn = operands_pop(p, token);
  gal_data_t *in   = operands_pop(p, token);

  /* Basic sanity checks. */
  conn_int=arithmetic_binary_sanity_checks(in, conn, token);

  /* Fill the holes */
  gal_binary_holes_fill(in, conn_int, -1);

  /* Push the result onto the stack. */
  operands_add(p, NULL, in);
}





static void
arithmetic_invert(struct arithmeticparams *p, char *token)
{
  gal_data_t *in = operands_pop(p, token);

  uint8_t *u8  = in->array, *u8f  = u8  + in->size;
  uint8_t *u16 = in->array, *u16f = u16 + in->size;
  uint8_t *u32 = in->array, *u32f = u32 + in->size;
  uint8_t *u64 = in->array, *u64f = u64 + in->size;

  /* Do the inversion based on type. */
  switch(in->type)
    {
    case GAL_TYPE_UINT8:  do *u8  = UINT8_MAX-*u8;   while(++u8<u8f);   break;
    case GAL_TYPE_UINT16: do *u16 = UINT16_MAX-*u16; while(++u16<u16f); break;
    case GAL_TYPE_UINT32: do *u32 = UINT32_MAX-*u32; while(++u32<u32f); break;
    case GAL_TYPE_UINT64: do *u64 = UINT64_MAX-*u64; while(++u64<u64f); break;
    default:
      error(EXIT_FAILURE, 0, "'invert' operand has %s type. 'invert' can "
            "only take unsigned integer types.\n\nYou can use any of the "
            "'uint8', 'uint16', 'uint32', or 'uint64' operators to chage "
            "the type before calling 'invert'",
            gal_type_name(in->type, 1));
    }

  /* Push the result onto the stack. */
  operands_add(p, NULL, in);
}





static void
arithmetic_interpolate(struct arithmeticparams *p, char *token)
{
  int num_int;
  gal_data_t *interpolated;

  /* First pop the number of nearby neighbors.*/
  gal_data_t *num = operands_pop(p, token);

  /* Then pop the actual dataset to interpolate. */
  gal_data_t *in = operands_pop(p, token);

  /* Do proper sanity checks on 'num'. */
  if(num->size!=1)
    error(EXIT_FAILURE, 0, "the first popped operand to "
          "'interpolate-medianngb' must be a single number. However, "
          "it has %zu elements", num->size);
  if(num->type==GAL_TYPE_FLOAT32 || num->type==GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "the first popped operand to "
          "'interpolate-medianngb' is the number of nearby neighbors (a "
          "counter, an integer). It must NOT be a floating point.\n\n"
          "If its already an integer, but in a floating point container, "
          "you can use the 'int32' operator to convert it to a 32-bit "
          "integer for example");

  /* Convert the given number to a 32-bit integer and read it in. */
  num=gal_data_copy_to_new_type_free(num, GAL_TYPE_INT32);
  num_int = *((int32_t *)(num->array));

  /* Call the interpolation function. */
  interpolated=gal_interpolate_close_neighbors(in, NULL, p->cp.interpmetric,
                                               num_int, p->cp.numthreads,
                                               1, 0);

  /* Clean up and push the interpolated array onto the stack. */
  gal_data_free(in);
  gal_data_free(num);
  operands_add(p, NULL, interpolated);
}





static void
arithmetic_collapse(struct arithmeticparams *p, char *token, int operator)
{
  long dim;
  gal_data_t *collapsed=NULL;

  /* First popped operand is the dimension. */
  gal_data_t *dimension = operands_pop(p, token);

  /* The second popped operand is the desired input dataset. */
  gal_data_t *input = operands_pop(p, token);


  /* Small sanity check. */
  if( dimension->ndim!=1 || dimension->size!=1)
    error(EXIT_FAILURE, 0, "first popped operand of 'collapse-*' operators "
          "(dimension to collapse) must be a single number (single-element, "
          "one-dimensional dataset). But it has %zu dimension(s) and %zu "
          "element(s).", dimension->ndim, dimension->size);
  if(dimension->type==GAL_TYPE_FLOAT32 || dimension->type==GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "first popped operand of 'collapse-*' operators "
          "(dimension to collapse) must have an integer type, but it has "
          "a floating point type ('%s')", gal_type_name(dimension->type,1));
  dimension=gal_data_copy_to_new_type_free(dimension, GAL_TYPE_LONG);
  dim=((long *)(dimension->array))[0];
  if(dim<0 || dim==0)
    error(EXIT_FAILURE, 0, "first popped operand of 'collapse-*' operators "
          "(dimension to collapse) must be positive (larger than zero), it "
          "is %ld", dim);
  if(dim > input->ndim)
    error(EXIT_FAILURE, 0, "input dataset to '%s' has %zu dimension(s), "
          "but you have asked to collapse along dimension %zu", token,
          input->ndim, dim);


  /* If a WCS structure has been read, we'll need to pass it to
     'gal_dimension_collapse', so it modifies it respectively. */
  if(p->wcs_collapsed==0)
    {
      p->wcs_collapsed=1;
      input->wcs=p->refdata.wcs;
    }


  /* Run the relevant library function. */
  switch(operator)
    {
    case ARITHMETIC_OP_COLLAPSE_SUM:
      collapsed=gal_dimension_collapse_sum(input, input->ndim-dim, NULL);
      break;

    case ARITHMETIC_OP_COLLAPSE_MEAN:
      collapsed=gal_dimension_collapse_mean(input, input->ndim-dim, NULL);
      break;

    case ARITHMETIC_OP_COLLAPSE_NUMBER:
      collapsed=gal_dimension_collapse_number(input, input->ndim-dim);
      break;

    case ARITHMETIC_OP_COLLAPSE_MIN:
      collapsed=gal_dimension_collapse_minmax(input, input->ndim-dim, 0);
      break;

    case ARITHMETIC_OP_COLLAPSE_MAX:
      collapsed=gal_dimension_collapse_minmax(input, input->ndim-dim, 1);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
            "problem. The operator code %d is not recognized", __func__,
            PACKAGE_BUGREPORT, operator);
    }


  /* If a WCS structure existed, a modified WCS is now present in
     'collapsed->wcs'. So we'll let the freeing of 'input' free the old
     'p->refdata.wcs' structure and we'll put the new one there, then we'll
     set 'collapsed->wcs' to 'NULL', so the new one isn't freed. */
  if(collapsed->wcs)
    {
      p->refdata.wcs = collapsed->wcs;
      collapsed->wcs = NULL;
    }


  /* We'll also need to correct the size of the reference dataset if it
     hasn't been corrected yet. We'll use 'memcpy' to write the new 'dsize'
     values into the old ones. The dimensions have decreased, so we won't
     be writing outside of allocated space that 'p->refdata.dsize' points
     to. */
  if( p->refdata.ndim != collapsed->ndim )
    {
      p->refdata.ndim -= 1;
      memcpy( p->refdata.dsize, collapsed->dsize,
              p->refdata.ndim * (sizeof *p->refdata.dsize) );
    }


  /* Clean up and add the collapsed dataset to the top of the operands. */
  gal_data_free(input);
  gal_data_free(dimension);
  operands_add(p, NULL, collapsed);
}





void
arithmetic_tofile(struct arithmeticparams *p, char *token, int freeflag)
{
  /* Pop the top dataset. */
  gal_data_t *popped = operands_pop(p, token);
  char *filename = &token[ freeflag
                           ? OPERATOR_PREFIX_LENGTH_TOFILEFREE
                           : OPERATOR_PREFIX_LENGTH_TOFILE     ];

  /* Save it to a file. */
  popped->wcs=p->refdata.wcs;
  if(popped->ndim==1 && p->onedasimage==0)
    gal_table_write(popped, NULL, p->cp.tableformat, filename,
                    "ARITHMETIC", 0);
  else
    gal_fits_img_write(popped, filename, NULL, PROGRAM_NAME);
  if(!p->cp.quiet)
    printf(" - Write: %s\n", filename);

  /* Reset the WCS to NULL and put it back on the stack. */
  popped->wcs=NULL;
  if(freeflag)
    gal_data_free(popped);
  else
    operands_add(p, NULL, popped);
}





void
arithmetic_unique(struct arithmeticparams *p, char *token, int operator)
{
  /* Pass the popped operand to the statistics library. */
  gal_data_t *input = gal_statistics_unique(operands_pop(p, token), 1);
  operands_add(p, NULL, input);
}





void
arithmetic_add_dimension(struct arithmeticparams *p, char *token, int operator)
{
  gal_data_t *out=NULL;
  gal_data_t *tmp = operands_pop(p, token);
  size_t i, num, dsize[3], ndim=3, nbytes=0;

  /* Make sure the first operand is a number. */
  if(tmp->size!=1)
    error(EXIT_FAILURE, 0, "first popped operand to '%s' must be a "
          "number (specifying how many datasets to use)", token);

  /* Put the value into 'num'. */
  tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_SIZE_T);
  num=*(size_t *)(tmp->array);
  gal_data_free(tmp);

  /* Pop all the datasets and put them in a list. */
  for(i=0;i<num;++i)
    {
      /* Pop the operand. */
      tmp=operands_pop(p, token);

      /* Things that differ from the first dataset and the rest. */
      if(out) /* Not the first. */
        {
          /* Basic sanity checks. */
          if(tmp->type!=out->type)
            error(EXIT_FAILURE, 0, "the operands to '%s' have to have the "
                  "same data type (the inputs contain atleast two types: "
                  "'%s' and '%s')", token, gal_type_name(tmp->type, 1),
                  gal_type_name(out->type, 1));
          if( tmp->ndim!=out->ndim-1
              || tmp->dsize[0]!=out->dsize[1]
              || tmp->dsize[1]!=out->dsize[2] )
            error(EXIT_FAILURE, 0, "the operands to '%s' have to have the "
                  "same size", token);
        }
      else  /* First popped operand. */
        {
          /* First popped operand, do necessary basic checks here. */
          if(tmp->ndim!=2)
            error(EXIT_FAILURE, 0, "currently only 2-dimensional datasets "
                  "are acceptable for '%s', please get in touch with us at "
                  "%s so we add functionality for different dimensions",
                  token, PACKAGE_BUGREPORT);

          /* Allocate the output dataset. */
          dsize[0]=num;
          dsize[1]=tmp->dsize[0];
          dsize[2]=tmp->dsize[1];
          out = gal_data_alloc(NULL, tmp->type, ndim, dsize, NULL, 0,
                               p->cp.minmapsize, p->cp.quietmmap, NULL,
                               NULL, NULL);

          /* Get the number of bytes in each dataset. */
          nbytes=gal_type_sizeof(tmp->type)*tmp->size;
        }

      /* Copy the dataset into the higher-dimensional output. */
      memcpy(gal_pointer_increment(out->array, (num-1-i)*tmp->size,
                                   tmp->type),
             tmp->array, nbytes);

      /* Clean up. */
      gal_data_free(tmp);
    }

  /* Put the higher-dimensional output on the operands stack. */
  operands_add(p, NULL, out);
}




















/***************************************************************/
/*************      Reverse Polish algorithm       *************/
/***************************************************************/
static int
arithmetic_set_operator(char *string, size_t *num_operands)
{
  /* Use the library's main function for its own operators. */
  int op = gal_arithmetic_set_operator(string, num_operands);

  /* If its not a library operator, check if its an internal operator. */
  if(op==GAL_ARITHMETIC_OP_INVALID)
    {
      /* Non-library operators. */
      if      (!strcmp(string, "filter-mean"))
        { op=ARITHMETIC_OP_FILTER_MEAN;           *num_operands=0; }
      else if (!strcmp(string, "filter-median"))
        { op=ARITHMETIC_OP_FILTER_MEDIAN;         *num_operands=0; }
      else if (!strcmp(string, "filter-sigclip-mean"))
        { op=ARITHMETIC_OP_FILTER_SIGCLIP_MEAN;   *num_operands=0; }
      else if (!strcmp(string, "filter-sigclip-median"))
        { op=ARITHMETIC_OP_FILTER_SIGCLIP_MEDIAN; *num_operands=0; }
      else if (!strcmp(string, "erode"))
        { op=ARITHMETIC_OP_ERODE;                 *num_operands=0; }
      else if (!strcmp(string, "dilate"))
        { op=ARITHMETIC_OP_DILATE;                *num_operands=0; }
      else if (!strcmp(string, "connected-components"))
        { op=ARITHMETIC_OP_CONNECTED_COMPONENTS;  *num_operands=0; }
      else if (!strcmp(string, "fill-holes"))
        { op=ARITHMETIC_OP_FILL_HOLES;            *num_operands=0; }
      else if (!strcmp(string, "invert"))
        { op=ARITHMETIC_OP_INVERT;                *num_operands=0; }
      else if (!strcmp(string, "interpolate-medianngb"))
        { op=ARITHMETIC_OP_INTERPOLATE_MEDIANNGB; *num_operands=0; }
      else if (!strcmp(string, "collapse-sum"))
        { op=ARITHMETIC_OP_COLLAPSE_SUM;          *num_operands=0; }
      else if (!strcmp(string, "collapse-min"))
        { op=ARITHMETIC_OP_COLLAPSE_MIN;          *num_operands=0; }
      else if (!strcmp(string, "collapse-max"))
        { op=ARITHMETIC_OP_COLLAPSE_MAX;          *num_operands=0; }
      else if (!strcmp(string, "collapse-mean"))
        { op=ARITHMETIC_OP_COLLAPSE_MEAN;         *num_operands=0; }
      else if (!strcmp(string, "collapse-number"))
        { op=ARITHMETIC_OP_COLLAPSE_NUMBER;       *num_operands=0; }
      else if (!strcmp(string, "unique"))
        { op=ARITHMETIC_OP_UNIQUE;                *num_operands=0; }
      else if (!strcmp(string, "add-dimension"))
        { op=ARITHMETIC_OP_ADD_DIMENSION;         *num_operands=0; }
      else
        error(EXIT_FAILURE, 0, "the argument '%s' could not be "
              "interpretted as a file name, named dataset, number, "
              "or operator", string);
    }

  /* Return the operator code. */
  return op;
}





static void
arithmetic_operator_run(struct arithmeticparams *p, int operator,
                        char *operator_string, size_t num_operands)
{
  size_t i;
  unsigned int numop;
  gal_data_t *d1=NULL, *d2=NULL, *d3=NULL;
  int flags = ( GAL_ARITHMETIC_INPLACE | GAL_ARITHMETIC_FREE
                | GAL_ARITHMETIC_NUMOK );

  /* When 'num_operands!=0', the operator is in the library. */
  if(num_operands)
    {
      /* Pop the necessary number of operators. Note that the
         operators are poped from a linked list (which is
         last-in-first-out). So for the operators which need a
         specific order, the first poped operand is actally the
         last (right most, in in-fix notation) input operand.*/
      switch(num_operands)
        {
        case 1:
          d1=operands_pop(p, operator_string);
          break;

        case 2:
          d2=operands_pop(p, operator_string);
          d1=operands_pop(p, operator_string);
          break;

        case 3:
          d3=operands_pop(p, operator_string);
          d2=operands_pop(p, operator_string);
          d1=operands_pop(p, operator_string);
          break;

        case -1:
          /* This case is when the number of operands is itself an
             operand. So except for sigma-clipping (that has other
             parameters), the first popped operand must be an
             integer number, we will use that to construct a linked
             list of any number of operands within the single 'd1'
             pointer. */
          numop=pop_number_of_operands(p, operator, operator_string, &d2);
          for(i=0;i<numop;++i)
            gal_list_data_add(&d1, operands_pop(p, operator_string));
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
                "the problem. '%zu' is not recognized as an operand "
                "counter (with '%s')", __func__, PACKAGE_BUGREPORT,
                num_operands, operator_string);
        }

      /* Run the arithmetic operation. Note that 'gal_arithmetic'
         is a variable argument function (like printf). So the
         number of arguments it uses depend on the operator. So
         when the operator doesn't need three operands, the extra
         arguments will be ignored. */
      operands_add(p, NULL, gal_arithmetic(operator, p->cp.numthreads,
                                           flags, d1, d2, d3));
    }

  /* No need to call the arithmetic library, call the proper
     wrappers directly. */
  else
    {
      switch(operator)
        {
        case ARITHMETIC_OP_FILTER_MEAN:
        case ARITHMETIC_OP_FILTER_MEDIAN:
        case ARITHMETIC_OP_FILTER_SIGCLIP_MEAN:
        case ARITHMETIC_OP_FILTER_SIGCLIP_MEDIAN:
          wrapper_for_filter(p, operator_string, operator);
          break;

        case ARITHMETIC_OP_ERODE:
        case ARITHMETIC_OP_DILATE:
          arithmetic_erode_dilate(p, operator_string, operator);
          break;

        case ARITHMETIC_OP_CONNECTED_COMPONENTS:
          arithmetic_connected_components(p, operator_string);
          break;

        case ARITHMETIC_OP_FILL_HOLES:
          arithmetic_fill_holes(p, operator_string);
          break;

        case ARITHMETIC_OP_INVERT:
          arithmetic_invert(p, operator_string);
          break;

        case ARITHMETIC_OP_INTERPOLATE_MEDIANNGB:
          arithmetic_interpolate(p, operator_string);
          break;

        case ARITHMETIC_OP_COLLAPSE_SUM:
        case ARITHMETIC_OP_COLLAPSE_MIN:
        case ARITHMETIC_OP_COLLAPSE_MAX:
        case ARITHMETIC_OP_COLLAPSE_MEAN:
        case ARITHMETIC_OP_COLLAPSE_NUMBER:
          arithmetic_collapse(p, operator_string, operator);
          break;

        case ARITHMETIC_OP_UNIQUE:
          arithmetic_unique(p, operator_string, operator);
          break;

        case ARITHMETIC_OP_ADD_DIMENSION:
          arithmetic_add_dimension(p, operator_string, operator);
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! please contact us at "
                "%s to fix the problem. The code %d is not "
                "recognized for 'op'", __func__, PACKAGE_BUGREPORT,
                operator);
        }
    }
}





/* This function implements the reverse polish algorithm as explained
   in the Wikipedia page.

   NOTE that in ui.c, the input linked list of tokens was ordered to
   have the same order as what the user provided. */
void
reversepolish(struct arithmeticparams *p)
{
  gal_data_t *data;
  size_t num_operands=0;
  gal_list_str_t *token;
  char *hdu, *filename, *printnum;
  int operator=GAL_ARITHMETIC_OP_INVALID;


  /* Prepare the processing: */
  p->operands=NULL;
  p->popcounter=0;


  /* Go over each input token and do the work. */
  for(token=p->tokens;token!=NULL;token=token->next)
    {
      /* The 'tofile-' operator's string can end in a '.fits', similar to a
         FITS file input file. So, it needs to be checked before checking
         for a filename. If we have a name or number, then add it to the
         operands linked list. Otherwise, pull out two members and do the
         specified operation on them. */
      operator=GAL_ARITHMETIC_OP_INVALID;
      if( !strncmp(OPERATOR_PREFIX_TOFILE, token->v,
                   OPERATOR_PREFIX_LENGTH_TOFILE) )
        arithmetic_tofile(p, token->v, 0);
      else if( !strncmp(OPERATOR_PREFIX_TOFILEFREE, token->v,
                   OPERATOR_PREFIX_LENGTH_TOFILE) )
        arithmetic_tofile(p, token->v, 1);
      else if( !strncmp(token->v, OPERATOR_PREFIX_SET,
                        OPERATOR_PREFIX_LENGTH_SET) )
        operands_set_name(p, token->v);
      else if( gal_array_name_recognized(token->v)
          || operands_is_name(p, token->v) )
        operands_add(p, token->v, NULL);
      else if( (data=gal_data_copy_string_to_number(token->v)) )
        operands_add(p, NULL, data);
      /* Last option is an operator: the program will abort if the token
         isn't an operator. */
      else
        {
          operator=arithmetic_set_operator(token->v, &num_operands);
          arithmetic_operator_run(p, operator, token->v, num_operands);
        }

      /* Increment the token counter. */
      ++p->tokencounter;
    }


  /* If there aren't any more operands (a variable has been set but not
     used), then there is nothing to create. */
  if(p->operands==NULL)
    error(EXIT_FAILURE, 0, "no operands on the stack to write (as output)");


  /* If there is more than one node in the operands stack then the user has
     given too many operands which is an error. */
  if(p->operands->next!=NULL)
    error(EXIT_FAILURE, 0, "too many operands");


  /* If the final operand has a filename, but its 'data' element is NULL,
     then the file hasn't actually be read yet. In this case, we need to
     read the contents of the file and put the resulting dataset into the
     operands 'data' element. This can happen for example if no operators
     are called and there is only one filename as an argument (which can
     happen in scripts). */
  if(p->operands->data==NULL && p->operands->filename)
    {
      /* Read the desired image and report it if necessary. */
      hdu=p->operands->hdu;
      filename=p->operands->filename;
      if( gal_fits_name_is_fits(filename) )
        {
          /* Read the data, note that the WCS has already been set. */
          p->operands->data=gal_array_read_one_ch(filename, hdu, NULL,
                                                  p->cp.minmapsize,
                                                  p->cp.quietmmap);
          data=p->operands->data;
          data->ndim=gal_dimension_remove_extra(data->ndim, data->dsize, NULL);
          if(!p->cp.quiet) printf(" - %s (hdu %s) is read.\n", filename, hdu);
        }
      else
        error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to fix "
              "the problem. While 'operands->data' is NULL, the filename "
              "('%s') is not recognized as a FITS file", __func__,
              PACKAGE_BUGREPORT, filename);
    }


  /* If the final data structure has more than one element, write it as a
     FITS file. Otherwise, print it in the standard output. */
  data=p->operands->data;
  if(data->size==1)
    {
      /* Make the string to print the number. */
      printnum=gal_type_to_string(data->array, data->type, 0);
      printf("%s\n", printnum);

      /* Clean up. */
      free(printnum);
    }
  else
    {
      /* Put a copy of the WCS structure from the reference image, it
         will be freed while freeing 'data'. */
      data->wcs=p->refdata.wcs;
      if(data->ndim==1 && p->onedasimage==0)
        gal_table_write(data, NULL, p->cp.tableformat,
                        p->onedonstdout ? NULL : p->cp.output,
                        "ARITHMETIC", 0);
      else
        gal_fits_img_write(data, p->cp.output, NULL, PROGRAM_NAME);
      if(!p->cp.quiet)
        printf(" - Write (final): %s\n", p->cp.output);
    }


  /* Clean up, note that above, we copied the pointer to 'refdata->wcs'
     into 'data', so it is freed when freeing 'data'. */
  gal_data_free(data);
  free(p->refdata.dsize);
  gal_list_data_free(p->named);


  /* Clean up. Note that the tokens were taken from the command-line
     arguments, so the string within each token linked list must not be
     freed. */
  gal_list_str_free(p->tokens, 0);
  free(p->operands);
}



















/***************************************************************/
/*************             Top function            *************/
/***************************************************************/
void
arithmetic(struct arithmeticparams *p)
{
  /* Parse the arguments */
  reversepolish(p);
}
