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

#include <gnuastro/data.h>
#include <gnuastro/fits.h>
#include <gnuastro/array.h>
#include <gnuastro/statistics.h>

#include <checkset.h>

#include "main.h"

#include "operands.h"
#include "arithmetic.h"
















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
  struct gal_linkedlist_stll *token;
  gal_data_t *d1=NULL, *d2=NULL, *d3=NULL;
  unsigned char flags = ( GAL_DATA_ARITH_INPLACE | GAL_DATA_ARITH_FREE
                          | GAL_DATA_ARITH_NUMOK );


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
          if      (!strcmp(token->v, "+" ))
            { op=GAL_DATA_OPERATOR_PLUS;          nop=2;  }
          else if (!strcmp(token->v, "-" ))
            { op=GAL_DATA_OPERATOR_MINUS;         nop=2;  }
          else if (!strcmp(token->v, "*" ))
            { op=GAL_DATA_OPERATOR_MULTIPLY;      nop=2;  }
          else if (!strcmp(token->v, "/" ))
            { op=GAL_DATA_OPERATOR_DIVIDE;        nop=2;  }

          else if (!strcmp(token->v, "lt" ))
            { op=GAL_DATA_OPERATOR_LT;            nop=2;  }
          else if (!strcmp(token->v, "le"))
            { op=GAL_DATA_OPERATOR_LE;            nop=2;  }
          else if (!strcmp(token->v, "gt" ))
            { op=GAL_DATA_OPERATOR_GT;            nop=2;  }
          else if (!strcmp(token->v, "ge"))
            { op=GAL_DATA_OPERATOR_LE;            nop=2;  }
          else if (!strcmp(token->v, "eq"))
            { op=GAL_DATA_OPERATOR_EQ;            nop=2;  }
          else if (!strcmp(token->v, "ne"))
            { op=GAL_DATA_OPERATOR_NE;            nop=2;  }
          else if (!strcmp(token->v, "and"))
            { op=GAL_DATA_OPERATOR_AND;           nop=2;  }
          else if (!strcmp(token->v, "or"))
            { op=GAL_DATA_OPERATOR_OR;            nop=2;  }
          else if (!strcmp(token->v, "bitand"))
            { op=GAL_DATA_OPERATOR_BITAND;        nop=2;  }
          else if (!strcmp(token->v, "bitor"))
            { op=GAL_DATA_OPERATOR_BITOR;         nop=2;  }

          else if (!strcmp(token->v, "not"))
            { op=GAL_DATA_OPERATOR_NOT;           nop=1;  }
          else if (!strcmp(token->v, "isblank"))
            { op=GAL_DATA_OPERATOR_ISBLANK;       nop=1;  }
          else if (!strcmp(token->v, "where"))
            { op=GAL_DATA_OPERATOR_WHERE;         nop=3;  }

          else if (!strcmp(token->v, "abs"))
            { op=GAL_DATA_OPERATOR_ABS;           nop=1;  }
          else if (!strcmp(token->v, "pow"))
            { op=GAL_DATA_OPERATOR_POW;           nop=2;  }
          else if (!strcmp(token->v, "sqrt"))
            { op=GAL_DATA_OPERATOR_SQRT;          nop=1;  }
          else if (!strcmp(token->v, "log"))
            { op=GAL_DATA_OPERATOR_LOG;           nop=1;  }
          else if (!strcmp(token->v, "log10"))
            { op=GAL_DATA_OPERATOR_LOG10;         nop=1;  }

          else if (!strcmp(token->v, "minval"))
            { op=GAL_DATA_OPERATOR_MINVAL;        nop=1;  }
          else if (!strcmp(token->v, "maxval"))
            { op=GAL_DATA_OPERATOR_MAXVAL;        nop=1;  }
          else if (!strcmp(token->v, "min"))
            { op=GAL_DATA_OPERATOR_MIN;           nop=-1; }
          else if (!strcmp(token->v, "max"))
            { op=GAL_DATA_OPERATOR_MAX;           nop=-1; }
          else if (!strcmp(token->v, "average"))
            { op=GAL_DATA_OPERATOR_AVERAGE;       nop=-1; }
          else if (!strcmp(token->v, "median"))
            { op=GAL_DATA_OPERATOR_MEDIAN;        nop=-1; }
          else
            error(EXIT_FAILURE, 0, "the argument \"%s\" could not be "
                  "interpretted as a FITS file, number, or operator",
                  token->v);

          /* Pop the necessary number of operators. Note that the operators
             are poped from a linked list (which is last-in-first-out). So
             for the operators which need a specific order, the first poped
             operand is actally the last (right most, in in-fix notation)
             input operand.*/
          if(nop==1)
            d1=pop_operand(p, token->v);
          else if(nop==2)
          {
            d2=pop_operand(p, token->v);
            d1=pop_operand(p, token->v);
          }
          else if(nop==3)
          {
            d3=pop_operand(p, token->v);
            d2=pop_operand(p, token->v);
            d1=pop_operand(p, token->v);
          }
          else if(nop==-1)
            /* It might be possible to add a `next' pointer to the data
               structure so it can act like a linked list too. In that
               case, t */
            error(EXIT_FAILURE, 0, "Operators with an unknown number of "
                  "operands are not yet implemented.");
          else
            error(EXIT_FAILURE, 0, "No operators need %d operands", nop);

          /* Do the arithemtic operation. */
          add_operand(p, NULL, gal_data_arithmetic(op, flags, d1, d2, d3));
        }
    }

#if 0

          else if(!strcmp(token->v, "abs"))       takeabs(p);
          else if(!strcmp(token->v, "pow"))       topower(p, NULL);
          else if(!strcmp(token->v, "sqrt"))      takesqrt(p);
          else if(!strcmp(token->v, "log"))       takelog(p);
          else if(!strcmp(token->v, "log10"))     takelog10(p);
          else if(!strcmp(token->v, "minvalue"))  findmin(p);
          else if(!strcmp(token->v, "maxvalue"))  findmax(p);
          else if(!strcmp(token->v, "min")
                  || !strcmp(token->v, "max")
                  || !strcmp(token->v, "average")
                  || !strcmp(token->v, "median")) alloppixs(p, token->v);
          else if(!strcmp(token->v, "lt")
                  || !strcmp(token->v, "le")
                  || !strcmp(token->v, "gt")
                  || !strcmp(token->v, "ge")
                  || !strcmp(token->v, "eq")
                  || !strcmp(token->v, "neq"))    conditionals(p, token->v);
          else if(!strcmp(token->v, "and")
                  || !strcmp(token->v, "or"))     andor(p, token->v);
          else if(!strcmp(token->v, "not"))       notfunc(p);
          else if(!strcmp(token->v, "isblank"))   opisblank(p);
          else if(!strcmp(token->v, "where"))     where(p);
#endif


  /* If there is more than one node in the operands stack then the user has
     given too many operands which is an error. */
  if(p->operands->next!=NULL)
    error(EXIT_FAILURE, 0, "too many operands");


  /* Copy the WCS structure into the final dataset. */
  p->operands->data->wcs=p->refdata.wcs;


  /* Set the output type. */
  if(p->outtype==p->operands->data->type)
    d1=p->operands->data;
  else
    {
      d1=gal_data_copy_to_new_type(p->operands->data, p->outtype);
      gal_data_free(p->operands->data);
    }
  gal_fits_write_img(p->operands->data, p->cp.output, "Arithmetic", NULL,
                     SPACK_STRING);


  /* Clean up: */
  gal_data_free(d1);
  free(p->refdata.dsize);
  wcsfree(p->refdata.wcs);


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
