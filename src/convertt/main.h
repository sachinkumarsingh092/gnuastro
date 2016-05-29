/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <stdint.h>

#include "commonparams.h"



/* Progarm name macros: */
#define SPACK           "astconvertt" /* Subpackage executable name. */
#define SPACK_NAME      "ConvertType" /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION





#define BLANKCHANNELNAME "blank"




/* Format codes: */
#define TXTFORMAT     -1
#define FITSFORMAT     1
#define JPEGFORMAT     2
#define EPSFORMAT      3
#define PDFFORMAT      4





/* This is used in converting certain values in the array. */
struct change
{
  double from;
  double to;
  struct change *next;
};





struct uiparams
{
  char         *hdu2;            /* Header for second FITS image.   */
  char         *hdu3;            /* Header for third FITS image.    */
  char         *hdu4;            /* Header for fourth FITS image.   */

  int        hdu2set;
  int        hdu3set;
  int        hdu4set;

  int     qualityset;
  int   widthincmset;
  int borderwidthset;

  int     fluxlowset;
  int    fluxhighset;
  int     maxbyteset;
};





struct converttparams
{
  /* Before actual program: */
  struct uiparams         up;  /* User interface parameters.           */
  struct gal_commonparams cp; /* Common parameters.                    */

  /* Input: */
  struct gal_linkedlist_stll *inputnames; /* The names of input files. */
  size_t          numinputs;  /* Number of input files.                */
  int             inputtype;  /* The type of the input file.           */

  /* Output: */
  int            outputtype;  /* The type of the output file.          */
  int               quality;  /* Quality of JPEG image.                */
  float           widthincm;  /* Width in centimeters.                 */
  int           borderwidth;  /* Width of border in PostScript points. */
  int                   hex;  /* Use hexadecimal not ASCII85 encoding. */

  /* Flux: */
  double            fluxlow;  /* Lower flux truncation value.          */
  double           fluxhigh;  /* Higher flux truncation value.         */
  uint8_t           maxbyte;  /* Maximum byte value.                   */
  int             flminbyte;  /* fluxlow is minimum byte.              */
  int             fhmaxbyte;  /* fluxhigh is maximum byte.             */
  struct change     *change;  /* The value conversion string.          */
  int      changeaftertrunc;  /* First convert, then truncate.         */
  int                   log;  /* ==1: Save input in log scale.         */
  int                invert;  /* ==1: Invert the image.                */

  /* INTERNAL PARAMETERS: */
  time_t            rawtime;  /* Starting time of the program.         */

  /* Input channels: */
  char            *names[4];  /* File names in input order.            */
  size_t              numch;  /* Current Channel.                      */
  int            isblank[4];  /* ==1: this channel is blank or fill.   */
  int            bitpixs[4];  /* Bitpix values for each channel.       */
  int             numnul[4];  /* Number of NUL in each channel.        */
  double             *ch[4];  /* Data for each channel.                */
  uint8_t           *ech[4];  /* 8-bit color channels.                 */
  size_t              s0[4];  /* First C axis size for each channel.   */
  size_t              s1[4];  /* Second C axis size for each channel.  */
};


#endif
