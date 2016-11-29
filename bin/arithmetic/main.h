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

#include <commonparams.h>

/* Progarm name macros: */
#define SPACK           "astarithmetic"   /* Subpackage executable name. */
#define SPACK_NAME      "Arithmetic"      /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION


/* Do not use the main commonargs.h option reader for the --hdu
   option. */
#define NOTCOMMONHDU   1





/* Constants: */
#define NEGDASHREPLACE 11  /* A vertical tab (ASCII=11) for negative dash */





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





struct uiparams
{
  char        *maskname;   /* Name of mask image.                   */
  char            *mhdu;   /* Mask image HDU.                       */

  int           typeset;
  int       masknameset;
  int masknameallocated;
  int           mhduset;
};





struct imgarithparams
{
  /* Other structures: */
  struct uiparams         up;  /* User interface parameters.            */
  struct gal_commonparams cp;  /* Common parameters.                    */

  /* Input: */
  struct gal_linkedlist_stll *hdus; /* List of all given HDU strings.   */
  struct gal_linkedlist_stll *tokens; /* List of all arithmetic tokens. */
  size_t           numfits;  /* Total number of input FITS images.      */
  size_t        addcounter;  /* The number of FITS images added.        */
  size_t        popcounter;  /* The number of FITS images popped.       */
  gal_data_t       refdata;  /* Container for information of the data.  */

  /* Output: */
  int              outtype;  /* User's desired output bixpix value.     */

  /* Operating mode: */

  /* Internal: */
  struct operand *operands;  /* The operands linked list.               */
  time_t           rawtime;  /* Starting time of the program.           */
};





#endif
