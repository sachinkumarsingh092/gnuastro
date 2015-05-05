/*********************************************************************
ImageStatistics - Get general statistics about the image.
ImgeStatistics is part of GNU Astronomy Utilities (Gnuastro) package.

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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MAIN_H
#define MAIN_H

#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astimgstat"       /* Subpackage executable name. */
#define SPACK_NAME      "ImageStatistics"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION






struct uiparams
{
  char             *inputname;  /* Name of input file.               */
  char              *maskname;  /* Name of mask image.               */
  char                  *mhdu;  /* Mask image HDU.                   */
  int             masknameset;
  int       masknameallocated;
  int                 mhduset;
  int          histnumbinsset;
  int              histminset;
  int              histmaxset;
  int            histquantset;
  int               cfpnumset;
  int               cfpminset;
  int               cfpmaxset;
  int             cfpquantset;
  int        sigclipmultipset;
  int     sigcliptoleranceset;
  int           sigclipnumset;

};





struct imgstatparams
{
  /* Other structures: */
  struct uiparams       up;  /* User interface parameters.            */
  struct commonparams   cp;  /* Common parameters.                    */

  /* Input: */
  float               *img;  /* Input image array.                    */
  size_t              size;  /* Number of non-blank data elements.    */
  int            ignoremin;  /* Ignore all data with minimum value.   */

  /* Output: */

  /* Histogram: */
  char           *histname;  /* Histogram file name.                  */
  int             normhist;  /* ==1: Normalize the histogram.         */
  int           maxhistone;  /* Scale such that max bin is one.       */
  int            binonzero;  /* Shift histogram, one bin starts at 0. */
  size_t       histnumbins;  /* Number of bins in the histogram.      */
  float            histmin;  /* Minimum value to use in histogram.    */
  float            histmax;  /* Maximum value to use in histogram.    */
  float          histquant;  /* Quantile range of histogram.          */

  /* Cumulative frequency plot: */
  char            *cfpname;  /* Cumultiave frequency plot file name.  */
  int              normcfp;  /* ==1: Normalize the CFP.               */
  int      maxcfpeqmaxhist;  /* ==1: CFP max equal to historgram max. */
  int           cfpsimhist;  /* ==1: CFP range equal to hist range.   */
  size_t            cfpnum;  /* The number of points to sample CFP.   */
  float             cfpmin;  /* Minimum value for CFP.                */
  float             cfpmax;  /* Maximum value for CFP.                */
  float           cfpquant;  /* Quantile range of CFP.                */

  /* Sigma clipping: */
  int              sigclip;  /* ==1: Do sigma clipping.               */
  float      sigclipmultip;  /* Multiple of STD in sigma clipping.    */
  float   sigcliptolerance;  /* Tolerance level in sigma clipping.    */
  size_t        sigclipnum;  /* Number of times to sigma clip.        */


  /* Operating mode: */

  /* Internal: */
  time_t           rawtime;  /* Starting time of the program.         */
};

#endif
