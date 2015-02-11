/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MAIN_H
#define MAIN_H

#include <stdint.h>

#include "config.h"
#include "commonparams.h"



/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astconvertt" /* Subpackage executable name. */
#define SPACK_NAME      "ConvertType" /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION





/* Format codes: */
#define FITSFORMAT 1
#define TXTFORMAT  2
#define JPEGFORMAT 3
#define EPSFORMAT  4



struct uiparams
{
  char         *h2;              /* Header for second FITS image.   */
  char         *h3;              /* Header for third FITS image.    */
  char         *h4;              /* Header for fourth FITS image.   */

  int        h2set;
  int        h3set;
  int        h4set;

  int   qualityset;
  int widthincmset;
  int   kincmykset;

  int   fluxlowset;
  int  fluxhighset;
  int   maxbyteset;
};





struct converttparams
{
  /* Before actual program: */
  struct commonparams  cp;      /* Common parameters.                  */
  struct     uiparams  up;      /* User interface parameters.          */

  /* Input: */
  struct stll *inputnames;      /* The names of input files.           */
  size_t        numinputs;      /* Number of input files.              */
  int           inputtype;      /* The type of the input file.         */

  /* Output: */
  int          outputtype;      /* The type of the output file.        */
  int             kincmyk;      /* ==1: Only input is K in CMYK.       */
  int             quality;      /* Quality of JPEG image.              */
  float         widthincm;      /* Width in centimeters.               */

  /* Flux: */
  float           fluxlow;      /* Lower flux truncation value.        */
  float          fluxhigh;      /* Higher flux truncation value.       */
  uint8_t         maxbyte;      /* Maximum byte value.                 */
  int           flminbyte;      /* fluxlow is minimum byte.            */
  int           fhmaxbyte;      /* fluxhigh is maximum byte.           */
  char           *convert;      /* The value conversion string.        */
  int           convfirst;      /* First convert, then truncate.       */
  int                 log;      /* ==1: Save input in log scale.       */
  int              invert;      /* ==1: Invert the image.              */

  /* INTERNAL PARAMETERS: */
  time_t          rawtime;      /* Starting time of the program.       */

  /* Input channels: */
  size_t            numch;      /* Current Channel.                    */
  double           *ch[4];      /* Data for each channel.              */
  size_t            s0[4];      /* First C axis size for each channel. */
  size_t            s1[4];      /* Second C axis size for each channel.*/
};


#endif
