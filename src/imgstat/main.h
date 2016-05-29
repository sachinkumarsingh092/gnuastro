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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MAIN_H
#define MAIN_H

#include <gnuastro/commonparams.h>

/* Progarm name macros: */
#define SPACK           "astimgstat"       /* Subpackage executable name. */
#define SPACK_NAME      "ImageStatistics"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION



#define PRINTFLT  "%.6f\n"
#define PRINTINT  "%.0f\n"
#define STRVAL    "   -- %-45s%s\n"
#define FNAMEVAL  "   -- %-45s%f\n"
#define SNAMEVAL  "   -- %-45s%lu\n"
#define ASCIIHISTNUMBINS    60
#define ASCIIHISTHEIGHT     10
#define HISTSTRING     "Histogram"
#define CFPSTRING      "Cumulative Frequency Plot"



struct uiparams
{
  char             *inputname;  /* Name of input file.               */
  char              *maskname;  /* Name of mask image.               */
  char                  *mhdu;  /* Mask image HDU.                   */
  int             masknameset;
  int       masknameallocated;
  int                 mhduset;
  int           mirrordistset;
  int          onebinvalueset;
  int          histnumbinsset;
  int              histminset;
  int              histmaxset;
  int            histquantset;
  int               cfpnumset;
  int               cfpminset;
  int               cfpmaxset;
  int             cfpquantset;
  int       mirrorplotdistset;
  int        sigclipmultipset;
  int     sigcliptoleranceset;
  int           sigclipnumset;

};





struct imgstatparams
{
  /* Other structures: */
  struct uiparams         up; /* User interface parameters.             */
  struct gal_commonparams cp; /* Common parameters.                     */

  /* Input: */
  float               *img;  /* Input image array.                      */
  float            *sorted;  /* Sorted input data.                      */
  size_t              size;  /* Number of non-blank data elements.      */
  int            ignoremin;  /* Ignore all data with minimum value.     */
  float         mirrordist;  /* Distance to go out for mode check.      */

  /* Output: */
  int            asciihist;  /* ==1: print an ASCII histogram.          */
  float        onebinvalue;  /* Shift bins so one bin starts with this. */
  int             lowerbin;  /* Interval lower limit as column 1.       */
  float             mirror;  /* Mirror quantile.                        */
  char         *mirrorhist;  /* Name of mirror histogram.               */
  char          *mirrorcfp;  /* Name of mirror CFP.                     */
  char          *mhistname;  /* Name of mode mirror histogram.          */
  float     mirrorplotdist;  /* Dist. after mode to display on plot.    */
  char           *mcfpname;  /* Name of mode mirror CFP.                */
  int   histrangeformirror;  /* ==1: Use histogram range for mirror.    */

  /* Histogram: */
  char           *histname;  /* Histogram file name.                    */
  int             normhist;  /* ==1: Normalize the histogram.           */
  int           maxhistone;  /* Scale such that max bin is one.         */
  size_t       histnumbins;  /* Number of bins in the histogram.        */
  float            histmin;  /* Minimum value to use in histogram.      */
  float            histmax;  /* Maximum value to use in histogram.      */
  float          histquant;  /* Quantile range of histogram.            */

  /* Cumulative frequency plot: */
  char            *cfpname;  /* Cumultiave frequency plot file name.    */
  int              normcfp;  /* ==1: Normalize the CFP.                 */
  int      maxcfpeqmaxhist;  /* ==1: CFP max equal to historgram max.   */
  int           cfpsimhist;  /* ==1: CFP range equal to hist range.     */
  size_t            cfpnum;  /* The number of points to sample CFP.     */
  float             cfpmin;  /* Minimum value for CFP.                  */
  float             cfpmax;  /* Maximum value for CFP.                  */
  float           cfpquant;  /* Quantile range of CFP.                  */

  /* Sigma clipping: */
  int              sigclip;  /* ==1: Do sigma clipping.                 */
  float      sigclipmultip;  /* Multiple of STD in sigma clipping.      */
  float   sigcliptolerance;  /* Tolerance level in sigma clipping.      */
  size_t        sigclipnum;  /* Number of times to sigma clip.          */


  /* Operating mode: */

  /* Internal: */
  time_t           rawtime;  /* Starting time of the program.           */
};

#endif
