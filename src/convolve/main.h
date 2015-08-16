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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MAIN_H
#define MAIN_H



#include "mesh.h"
#include "astrthreads.h"
#include "commonparams.h"



/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astconvolve" /* Subpackage executable name. */
#define SPACK_NAME      "Convolve" /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION



#define MINGOODDIVSPEC       0.005
#define CONVFLOATINGPOINTERR 1e-10
#define COMPLEXTOREALSPEC    1  /* Spectrum of complex number.  */
#define COMPLEXTOREALPHASE   2  /* Phase of the complex number. */
#define COMPLEXTOREALREAL    3  /* Only show the real part.     */



/* User interface parameters structure */
struct uiparams
{
  char        *inputname;   /* Name of input file.           */
  char         *maskname;   /* Mask file name.               */
  char             *mhdu;   /* HDU of mask image.            */
  char       *kernelname;   /* Name of kernel file.          */
  char             *khdu;   /* HDU of input Kernel.          */
  char    *freqstepsname;   /* Name of frequency steps file. */

  int         spatialset;
  int       frequencyset;
  int        masknameset;
  int            mhduset;
  int      kernelnameset;
  int            khduset;
  int        meshsizeset;
  int            nch1set;
  int            nch2set;
  int    lastmeshfracset;
  int fullconvolutionset;
  int      makekernelset;
};




/* Processing parameters structure */
struct convolveparams
{
  struct uiparams     up;   /* user interface structure.                */
  struct commonparams cp;   /* commonparams structure.                  */
  struct meshparams   mp;   /* meshparams structure.                    */

  /* Inputs: */
  int         makekernel;   /* ==1: Make a kernel to create input.      */
  float           *input;   /* Input image array.                       */
  float          *kernel;   /* Input Kernel array.                      */
  int           anyblank;   /* If there are blank pixels in input.      */
  size_t             is0;   /* Input image size along C's first axis.   */
  size_t             is1;   /* Input image size along C's second axis.  */
  size_t             ks0;   /* Kernel size along C's first axis.        */
  size_t             ks1;   /* Kernel size along C's second axis.       */
  int         kernelflip;   /* ==1: Flip the kernel.                    */
  int         kernelnorm;   /* ==1: Normalize the kernel.               */
  int               nwcs;   /* Number of WCS headers.                   */
  struct wcsprm     *wcs;   /* WCS structure.                           */

  /* Outputs: */
  char         *meshname;   /* Name of mesh structure output.           */

  /* Operating modes: */
  int            spatial;   /* Convolve using spatial domain.           */
  int          frequency;   /* Convolve using frequency domain.         */
  int      viewfreqsteps;   /* View the frequency domain steps.         */

  /* internal: */
  time_t         rawtime;   /* Starting time of the program.            */
  double           *pimg;   /* Padded image array.                      */
  double           *pker;   /* Padded kernel array.                     */
  size_t             ps0;   /* Padded size along first C axis.          */
  size_t             ps1;   /* Padded size along second C axis.         */
};

#endif
