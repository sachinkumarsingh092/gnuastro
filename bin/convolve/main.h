/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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

/* Program names.  */
#define PROGRAM_NAME   "Convolve"      /* Program full name.       */
#define PROGRAM_EXEC   "astconvolve"   /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





/* Macros */
#define CONVFLOATINGPOINTERR 1e-10





/* Enumerators */
enum complex_to_real
{
  COMPLEX_TO_REAL_INVALID,           /* ==0 by C standard. */

  COMPLEX_TO_REAL_SPEC,
  COMPLEX_TO_REAL_PHASE,
  COMPLEX_TO_REAL_REAL,
};


enum domain_codes
{
  CONVOLVE_DOMAIN_INVALID,           /* ==0 by C standard. */

  CONVOLVE_DOMAIN_SPATIAL,
  CONVOLVE_DOMAIN_FREQUENCY,
};





/* Processing parameters structure */
struct convolveparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.              */
  char             *filename;  /* Name of input file.                     */
  char           *kernelname;  /* File name of kernel.                    */
  char                 *khdu;  /* HDU of kernel.                          */
  uint8_t       nokernelflip;  /* Do not flip the kernel.                 */
  uint8_t       nokernelnorm;  /* Do not normalize the kernel.            */
  double        minsharpspec;  /* Deconvolution: min spect. of sharp img. */
  uint8_t     checkfreqsteps;  /* View the frequency domain steps.        */
  size_t               *tile;  /* Size of tiles along each dim. (C order).*/
  size_t        *numchannels;  /* No. of tiles along each dim. (C order). */
  float        remainderfrac;  /* Frac. of remainers in each dim to cut.  */
  uint8_t         convoverch;  /* Convolve over channel borders.          */
  uint8_t         checktiles;  /* Tile IDs in an img, the size of input.  */
  char            *domainstr;  /* String value specifying domain.         */
  size_t          makekernel;  /* Make a kernel to create input.          */

  /* Internal */
  size_t            *channel;  /* Size of channels along each dimension.  */
  int                 domain;  /* Frequency or spatial domain conv.       */
  gal_data_t          *input;  /* Input image array.                      */
  gal_data_t         *kernel;  /* Input Kernel array.                     */
  double               *pimg;  /* Padded image array.                     */
  double               *pker;  /* Padded kernel array.                    */
  size_t                 ps0;  /* Padded size along first C axis.         */
  size_t                 ps1;  /* Padded size along second C axis.        */
  char        *freqstepsname;  /* Name of file to check frequency steps.  */
  char            *tilesname;  /* Name of file to check tiles.            */
  time_t             rawtime;  /* Starting time of the program.           */
};

#endif
