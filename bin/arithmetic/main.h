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
#ifndef MAIN_H
#define MAIN_H

#include <gnuastro/fits.h>
#include <gnuastro/linkedlist.h>

#include <options.h>


/* Progarm name macros: */
#define PROGRAM_NAME "Arithmetic"      /* Program full name.       */
#define PROGRAM_EXEC "astarithmetic"   /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION






/* Constants: */
#define NEG_DASH_REPLACE 11 /* Vertical tab (ASCII=11) for negative dash */





/* In every node of the operand linked list, only one of the `filename' or
   `data' should be non-NULL. Otherwise it will be a bug and will cause
   problems. All the operands operate on this premise. */
struct operand
{
  char       *filename;    /* !=NULL if the operand is a filename. */
  char            *hdu;    /* !=NULL if the operand is a filename. */
  gal_data_t     *data;    /* !=NULL if the operand is a dataset.  */
  struct operand *next;    /* Pointer to next operand.             */
};






struct imgarithparams
{
  /* Other structures: */
  struct gal_options_common_params cp;  /* Common parameters.           */

  /* Input: */
  struct gal_linkedlist_stll   *hdus; /* List of all given HDU strings. */
  struct gal_linkedlist_stll *tokens; /* List of all arithmetic tokens. */
  size_t        addcounter;  /* The number of FITS images added.        */
  size_t        popcounter;  /* The number of FITS images popped.       */
  gal_data_t       refdata;  /* Container for information of the data.  */

  /* Operating mode: */

  /* Internal: */
  struct operand *operands;  /* The operands linked list.               */
  time_t           rawtime;  /* Starting time of the program.           */
};





#endif
