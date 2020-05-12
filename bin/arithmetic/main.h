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
#ifndef MAIN_H
#define MAIN_H

#include <gnuastro/fits.h>
#include <gnuastro/list.h>

#include <gnuastro-internal/options.h>


/* Progarm name macros: */
#define PROGRAM_NAME   "Arithmetic"     /* Program full name.       */
#define PROGRAM_EXEC   "astarithmetic"  /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





/* Constants: */
#define NEG_DASH_REPLACE 11 /* Vertical tab (ASCII=11) for negative dash */
#define OPERATOR_PREFIX_SET               "set-"
#define OPERATOR_PREFIX_TOFILE            "tofile-"
#define OPERATOR_PREFIX_TOFILEFREE        "tofilefree-"
#define OPERATOR_PREFIX_LENGTH_SET        strlen(OPERATOR_PREFIX_SET)
#define OPERATOR_PREFIX_LENGTH_TOFILE     strlen(OPERATOR_PREFIX_TOFILE)
#define OPERATOR_PREFIX_LENGTH_TOFILEFREE strlen(OPERATOR_PREFIX_TOFILEFREE)





/* In every node of the operand linked list, only one of the 'filename' or
   'data' should be non-NULL. Otherwise it will be a bug and will cause
   problems. All the operands operate on this premise. */
struct operand
{
  char       *filename;    /* !=NULL if the operand is a filename. */
  char            *hdu;    /* !=NULL if the operand is a filename. */
  gal_data_t     *data;    /* !=NULL if the operand is a dataset.  */
  struct operand *next;    /* Pointer to next operand.             */
};






struct arithmeticparams
{
  /* Other structures: */
  struct gal_options_common_params cp;  /* Common parameters.           */

  /* Input: */
  gal_list_str_t     *hdus;  /* List of all given HDU strings.          */
  gal_list_str_t   *tokens;  /* List of all arithmetic tokens.          */
  char            *wcsfile;  /* File to use for output's WCS.           */
  char             *wcshdu;  /* Extension to use for output's WCS.      */
  size_t        popcounter;  /* The number of FITS images popped.       */
  gal_data_t       refdata;  /* Container for information of the data.  */
  char          *globalhdu;  /* Single HDU for all inputs.              */
  uint8_t      onedasimage;  /* Write 1D outputs as an image not table. */
  uint8_t     onedonstdout;  /* Write 1D outputs on stdout, not table.  */
  gal_data_t        *named;  /* List containing variables.              */
  size_t      tokencounter;  /* Counter for finding place in tokens.    */

  /* Operating mode: */
  int        wcs_collapsed;  /* If the internal WCS is already collapsed.*/

  /* Internal: */
  struct operand *operands;  /* The operands linked list.               */
  time_t           rawtime;  /* Starting time of the program.           */
};





#endif
