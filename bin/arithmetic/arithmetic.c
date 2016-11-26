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
#include <gnuastro/statistics.h>

#include <checkset.h>

#include "main.h"

#include "operands.h"
#include "operators.h"
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
  double number;
  gal_data_t *ddata;
  struct gal_linkedlist_stll *token;

  /* Prepare the processing: */
  p->operands=NULL;
  p->addcounter=p->popcounter=0;

  /* Go over each input token and do the work. */
  for(token=p->tokens;token!=NULL;token=token->next)
    {
      /* If we have a name or number, then add it to the operands
         linked list. Otherwise, pull out two members and do the
         specified operation on them. */
      if(gal_fits_name_is_fits(token->v))
        add_operand(p, token->v, NOOPTNUMBER, NOOPTDATA);
      else if(gal_checkset_str_is_double(token->v, &number))
        add_operand(p, NOOPTFILENAME, number, NOOPTDATA);
      else
        {
          if     (!strcmp(token->v, "+"))         sum(p);
          else if(!strcmp(token->v, "-"))         subtract(p);
          else if(!strcmp(token->v, "*"))         multiply(p);
          else if(!strcmp(token->v, "/"))         divide(p);
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
          else
            error(EXIT_FAILURE, 0, "the argument \"%s\" could not be "
                  "interpretted as an operator", token->v);
        }
    }

  /* If there is more than one node in the operands stack then the user has
     given too many operands which is an error. */
  if(p->operands->next!=NULL)
    error(EXIT_FAILURE, 0, "too many operands");

  /* If the remaining operand is an array then save the array as a
     FITS image, if not, simply print the floating point number. */
  if(p->operands->data)
    {
      /* Internally, all arrays were double type. But the user could set
         the the output type using the type option. So if the user has
         asked for anything other than a double, we will have to convert
         the arrays.*/
      if(p->outtype==GAL_DATA_TYPE_DOUBLE)
        ddata=p->operands->data;
      else
        {
          ddata=gal_data_copy_to_new_type(p->operands->data, p->outtype);
          gal_data_free(p->operands->data);
        }

      /* Copy the WCS structure into the final dataset and write it into
         the output file. */
      ddata->wcs=p->refdata.wcs;
      gal_fits_write_img(ddata, p->cp.output, "Arithmetic", NULL,
                         SPACK_STRING);

      /* Clean up: */
      gal_data_free(ddata);
      free(p->refdata.dsize);
      wcsfree(p->refdata.wcs);
    }
  else
    printf("%g\n", p->operands->number);

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
