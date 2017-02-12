/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/array.h>
#include <gnuastro/arithmetic.h>

#include <checkset.h>

#include "main.h"

#include "operands.h"
#include "arithmetic.h"














/***************************************************************/
/*************          Internal functions         *************/
/***************************************************************/
#define SET_NUM_OP(CTYPE) {                                \
    CTYPE a=*(CTYPE *)(data->array); if(a>0) return a;    }

static size_t
set_number_of_operands(struct imgarithparams *p, gal_data_t *data,
                       char *token_string)
{
  /* Check if its a number. */
  if(data->size>1)
    error(EXIT_FAILURE, 0, "the first popped operand to the \"%s\" "
          "operator must be a number, not an array", token_string);

  /* Check its type and return the value. */
  switch(data->type)
    {

    /* For the integer types, if they are unsigned, then just pass their
       value, if they are signed, you have to make sure they are zero or
       positive. */
    case GAL_DATA_TYPE_UINT8:   SET_NUM_OP(uint8_t);     break;
    case GAL_DATA_TYPE_INT8:    SET_NUM_OP(int8_t);      break;
    case GAL_DATA_TYPE_UINT16:  SET_NUM_OP(uint16_t);    break;
    case GAL_DATA_TYPE_INT16:   SET_NUM_OP(int16_t);     break;
    case GAL_DATA_TYPE_UINT32:  SET_NUM_OP(uint32_t);    break;
    case GAL_DATA_TYPE_INT32:   SET_NUM_OP(int32_t);     break;
    case GAL_DATA_TYPE_UINT64:  SET_NUM_OP(uint64_t);    break;
    case GAL_DATA_TYPE_INT64:   SET_NUM_OP(int64_t);     break;

    /* Floating point numbers are not acceptable in this context. */
    case GAL_DATA_TYPE_FLOAT32:
    case GAL_DATA_TYPE_FLOAT64:
      error(EXIT_FAILURE, 0, "the first popped operand to the \"%s\" "
            "operator must be an integer type", token_string);

    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`set_number_of_operands", data->type);
    }

  /* If control reaches here, then the number must have been a negative
     value, so print an error. */
  error(EXIT_FAILURE, 0, "the first popped operand to the \"%s\" operator "
        "cannot be zero or a negative number", token_string);
  return 0;
}




















/***************************************************************/
/*************      Reverse Polish algorithm       *************/
/***************************************************************/
/* This function implements the reverse polish algorithm as explained
   in the Wikipedia page.

   NOTE that in ui.c, the input linked list of tokens was ordered to
   have the same order as what the user provided. */
void
reversepolish(struct imgarithparams *p)
{
  int op=0, nop=0;
  unsigned int numop, i;
  struct gal_linkedlist_stll *token;
  gal_data_t *d1=NULL, *d2=NULL, *d3=NULL;
  unsigned char flags = ( GAL_ARITHMETIC_INPLACE | GAL_ARITHMETIC_FREE
                          | GAL_ARITHMETIC_NUMOK );


  /* Prepare the processing: */
  p->operands=NULL;
  p->addcounter=p->popcounter=0;


  /* Go over each input token and do the work. */
  for(token=p->tokens;token!=NULL;token=token->next)
    {
      /* If we have a name or number, then add it to the operands linked
         list. Otherwise, pull out two members and do the specified
         operation on them. */
      if(gal_fits_name_is_fits(token->v))
        add_operand(p, token->v, NULL);
      else if( (d1=gal_data_string_to_number(token->v)) )
        add_operand(p, NULL, d1);
      else
        {

          /* Order is the same as in the manual. */
          /* Simple arithmetic operators. */
          if      (!strcmp(token->v, "+" ))
            { op=GAL_ARITHMETIC_OP_PLUS;          nop=2;  }
          else if (!strcmp(token->v, "-" ))
            { op=GAL_ARITHMETIC_OP_MINUS;         nop=2;  }
          else if (!strcmp(token->v, "*" ))
            { op=GAL_ARITHMETIC_OP_MULTIPLY;      nop=2;  }
          else if (!strcmp(token->v, "/" ))
            { op=GAL_ARITHMETIC_OP_DIVIDE;        nop=2;  }
          else if (!strcmp(token->v, "%" ))
            { op=GAL_ARITHMETIC_OP_MODULO;        nop=2;  }

          /* Mathematical Operators. */
          else if (!strcmp(token->v, "abs"))
            { op=GAL_ARITHMETIC_OP_ABS;           nop=1;  }
          else if (!strcmp(token->v, "pow"))
            { op=GAL_ARITHMETIC_OP_POW;           nop=2;  }
          else if (!strcmp(token->v, "sqrt"))
            { op=GAL_ARITHMETIC_OP_SQRT;          nop=1;  }
          else if (!strcmp(token->v, "log"))
            { op=GAL_ARITHMETIC_OP_LOG;           nop=1;  }
          else if (!strcmp(token->v, "log10"))
            { op=GAL_ARITHMETIC_OP_LOG10;         nop=1;  }

          /* Statistical/higher-level operators. */
          else if (!strcmp(token->v, "minvalue"))
            { op=GAL_ARITHMETIC_OP_MINVAL;        nop=1;  }
          else if (!strcmp(token->v, "maxvalue"))
            { op=GAL_ARITHMETIC_OP_MAXVAL;        nop=1;  }
          else if (!strcmp(token->v, "min"))
            { op=GAL_ARITHMETIC_OP_MIN;           nop=-1; }
          else if (!strcmp(token->v, "max"))
            { op=GAL_ARITHMETIC_OP_MAX;           nop=-1; }
          else if (!strcmp(token->v, "sum"))
            { op=GAL_ARITHMETIC_OP_SUM;           nop=-1; }
          else if (!strcmp(token->v, "average"))
            { op=GAL_ARITHMETIC_OP_AVERAGE;       nop=-1; }
          else if (!strcmp(token->v, "median"))
            { op=GAL_ARITHMETIC_OP_MEDIAN;        nop=-1; }

          /* Conditional operators. */
          else if (!strcmp(token->v, "lt" ))
            { op=GAL_ARITHMETIC_OP_LT;            nop=2;  }
          else if (!strcmp(token->v, "le"))
            { op=GAL_ARITHMETIC_OP_LE;            nop=2;  }
          else if (!strcmp(token->v, "gt" ))
            { op=GAL_ARITHMETIC_OP_GT;            nop=2;  }
          else if (!strcmp(token->v, "ge"))
            { op=GAL_ARITHMETIC_OP_LE;            nop=2;  }
          else if (!strcmp(token->v, "eq"))
            { op=GAL_ARITHMETIC_OP_EQ;            nop=2;  }
          else if (!strcmp(token->v, "ne"))
            { op=GAL_ARITHMETIC_OP_NE;            nop=2;  }
          else if (!strcmp(token->v, "and"))
            { op=GAL_ARITHMETIC_OP_AND;           nop=2;  }
          else if (!strcmp(token->v, "or"))
            { op=GAL_ARITHMETIC_OP_OR;            nop=2;  }
          else if (!strcmp(token->v, "not"))
            { op=GAL_ARITHMETIC_OP_NOT;           nop=1;  }
          else if (!strcmp(token->v, "isblank"))
            { op=GAL_ARITHMETIC_OP_ISBLANK;       nop=1;  }
          else if (!strcmp(token->v, "where"))
            { op=GAL_ARITHMETIC_OP_WHERE;         nop=3;  }

          /* Bitwise operators. */
          else if (!strcmp(token->v, "bitand"))
            { op=GAL_ARITHMETIC_OP_BITAND;        nop=2;  }
          else if (!strcmp(token->v, "bitor"))
            { op=GAL_ARITHMETIC_OP_BITOR;         nop=2;  }
          else if (!strcmp(token->v, "bitxor"))
            { op=GAL_ARITHMETIC_OP_BITXOR;        nop=2;  }
          else if (!strcmp(token->v, "lshift"))
            { op=GAL_ARITHMETIC_OP_BITLSH;        nop=2;  }
          else if (!strcmp(token->v, "rshift"))
            { op=GAL_ARITHMETIC_OP_BITRSH;        nop=2;  }
          else if (!strcmp(token->v, "bitnot"))
            { op=GAL_ARITHMETIC_OP_BITNOT;        nop=1;  }

          /* Type conversion. */
          else if (!strcmp(token->v, "uint8"))
            { op=GAL_ARITHMETIC_OP_TO_UINT8;      nop=1;  }
          else if (!strcmp(token->v, "int8"))
            { op=GAL_ARITHMETIC_OP_TO_INT8;       nop=1;  }
          else if (!strcmp(token->v, "uint16"))
            { op=GAL_ARITHMETIC_OP_TO_UINT16;     nop=1;  }
          else if (!strcmp(token->v, "int16"))
            { op=GAL_ARITHMETIC_OP_TO_INT16;      nop=1;  }
          else if (!strcmp(token->v, "uint64"))
            { op=GAL_ARITHMETIC_OP_TO_UINT64;     nop=1;  }
          else if (!strcmp(token->v, "int64"))
            { op=GAL_ARITHMETIC_OP_TO_INT64;      nop=1;  }
          else if (!strcmp(token->v, "float32"))
            { op=GAL_ARITHMETIC_OP_TO_FLOAT32;    nop=1;  }
          else if (!strcmp(token->v, "float64"))
            { op=GAL_ARITHMETIC_OP_TO_FLOAT64;    nop=1;  }

          /* Finished checks with known operators */
          else
            error(EXIT_FAILURE, 0, "the argument \"%s\" could not be "
                  "interpretted as a FITS file, number, or operator",
                  token->v);


          /* Pop the necessary number of operators. Note that the operators
             are poped from a linked list (which is last-in-first-out). So
             for the operators which need a specific order, the first poped
             operand is actally the last (right most, in in-fix notation)
             input operand.*/
          switch(nop)
            {
            case 1:
              d1=pop_operand(p, token->v);
              break;

            case 2:
              d2=pop_operand(p, token->v);
              d1=pop_operand(p, token->v);
              break;

            case 3:
              d3=pop_operand(p, token->v);
              d2=pop_operand(p, token->v);
              d1=pop_operand(p, token->v);
              break;

            case -1:
              /* This case is when the number of operands is itself an
                 operand. So the first popped operand must be an integer
                 number, we will use that to construct a linked list of any
                 number of operands within the single `d1' pointer. */
              d1=NULL;
              numop=set_number_of_operands(p, pop_operand(p, token->v),
                                           token->v);
              for(i=0;i<numop;++i)
                gal_data_add_existing_to_ll(&d1, pop_operand(p, token->v));
              break;

            default:
              error(EXIT_FAILURE, 0, "No operators need %d operands", nop);
            }


          /* Run the arithmetic operation. Note that `gal_arithmetic' is a
             variable argument function (like printf). So the number of
             arguments it uses depend on the operator. So when the operator
             doesn't need three operands, the extra arguments will be
             ignored. */
          add_operand(p, NULL, gal_arithmetic(op, flags, d1, d2, d3));
        }
    }

  /* If there is more than one node in the operands stack then the user has
     given too many operands which is an error. */
  if(p->operands->next!=NULL)
    error(EXIT_FAILURE, 0, "too many operands");


  /* Set the output type. */
  d1=p->operands->data;


  /* If the final data structure has more than one element, write it as a
     FITS file. Otherwise, print it in the standard output. */
  if(d1->size==1)
    {
      /* To simplify the printing process, we will first change it to
         double, then use printf's `%g' to print it, so integers will be
         printed as an integer.  */
      d2=gal_data_copy_to_new_type(d1, GAL_DATA_TYPE_FLOAT32);
      printf("%g\n", *(double *)d2->array);
      gal_data_free(d2);
    }
  else
    {
      /* Put a copy of the WCS structure from the reference image, it
         will be freed while freeing d1. */
      d1->wcs=p->refdata.wcs;
      gal_fits_img_write(d1, p->cp.output, NULL, PROGRAM_STRING);
      if(!p->cp.quiet)
        printf(" - Output written to %s\n", p->cp.output);
    }


  /* Clean up, note that above, we copied the pointer to `refdata->wcs'
     into `d1', so it is freed when freeing d1. */
  gal_data_free(d1);
  free(p->refdata.dsize);

  /* Clean up. Note that the tokens were taken from the command-line
     arguments, so the string within each token linked list must not be
     freed. */
  gal_linkedlist_free_stll(p->tokens, 0);
  free(p->operands);
}



















/***************************************************************/
/*************             Top function            *************/
/***************************************************************/
void
imgarith(struct imgarithparams *p)
{
  /* Parse the arguments */
  reversepolish(p);
}
