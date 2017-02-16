/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

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

/* Include necessary headers */
#include <gnuastro/data.h>

#include <options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "Statistics"    /* Program full name.       */
#define PROGRAM_EXEC   "aststatistics" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







/* Main program parameters structure */
struct statisticsparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.             */
  struct gal_linkedlist_ill  *toprint; /* Values to print in one row.    */
  char          *inputname;  /* Input filename.                          */
  char             *column;  /* Column name or number if input is table. */
  float       greaterequal;  /* Only use values >= this value.           */
  float           lessthan;  /* Only use values <  this value.           */
  float      quantilerange;  /* Quantile (Q) range: from Q to 1-Q.       */

  uint8_t        asciihist;  /* Print an ASCII histogram.                */
  uint8_t        histogram;  /* Save histogram in output.                */
  uint8_t       cumulative;  /* Save cumulative distibution in output.   */
  char         *sigclipstr;  /* Multiple of sigma, and tolerance or num. */
  float        mirrorquant;  /* Quantile of mirror distribution to save. */
  size_t           numbins;  /* Number of bins in histogram or CFP.      */
  uint8_t         lowerbin;  /* Save interval lower limit, not center.   */
  uint8_t        normalize;  /* set the sum of all bins to 1.            */
  float        onebinstart;  /* Shift bins to start at this value.       */
  uint8_t        maxbinone;  /* Set the maximum bin to 1.                */

  /* Internal */
  gal_data_t        *input;  /* Input data structure.                    */
  int               isfits;  /* Input is a FITS file.                    */
  int             hdu_type;  /* Type of HDU (image or table).            */
  time_t           rawtime;  /* Starting time of the program.            */
};

#endif
