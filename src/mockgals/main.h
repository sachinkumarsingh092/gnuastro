/*********************************************************************
MockGals - Create mock galaxies and stars in a noisy image.
MockGals is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MAIN_H
#define MAIN_H

#include "commonparams.h"



/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astrmockgals" /* Subpackage executable name. */
#define SPACK_NAME      "MockGals"     /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION
#define LOGFILENAME     SPACK".log"



struct uiparams
{
  char          *psfname;      /* Name of PSF FITS name.             */
  char          *catname;      /* Name of catalog of parameters.     */
  int            onlypsf;      /* ==1: Only save the PSF.            */

  /* Check if all parameters are read (use .def file for
     comparison). The non optional parameters (like the catalog and
     input FITS images that come in from arguments, not options) are
     checked in the args.h files. */

  int     psffunctionset;
  int            fwhmset;
  int      moffatbetaset;
  int        psftruncset;

  int      truncationset;
  int       toleranceset;
  int      backgroundset;
  int       zeropointset;
  int            fcolset;
  int            xcolset;
  int            ycolset;
  int            rcolset;
  int            ncolset;
  int            pcolset;
  int            qcolset;
  int            mcolset;

  int          naxis1set;
  int          naxis2set;
};





struct mockgalsparams
{
  /* Other structures: */
  struct uiparams     up;      /* User interface parameters.         */
  struct commonparams cp;      /* Common parameters.                 */

  /* PSF: */
  int        psffunction;      /* PSF Moffat or Gaussian.            */
  float           psf_p1;      /* First paramr of PSF (FWHM).        */
  float           psf_p2;      /* Second param of PSF (Moffat beta). */
  float            psf_t;      /* PSF truncation radius.             */

  /* Profiles and noise */
  float       truncation;      /* Truncation radius of the profiles. */
  float        tolerance;      /* Accuracy to stop integration.      */
  float       background;      /* Sky value in the image.            */
  float        zeropoint;      /* Magnitude of zero point flux.      */
  size_t            fcol;      /* Column specifying profile function.*/
  size_t            xcol;      /* X column of profile center.        */
  size_t            ycol;      /* Y column of profile center.        */
  size_t            rcol;      /* Effective radius of profile.       */
  size_t            ncol;      /* Sersic index column of profile.    */
  size_t            pcol;      /* Position angle column of profile.  */
  size_t            qcol;      /* Axis ratio column of profile.      */
  size_t            mcol;      /* Magnitude column.                  */

  /* Output */
  size_t              s0;      /* C standard axis 0 size.            */
  size_t              s1;      /* C standard axis 1 size.            */
  char          *logname;      /* Output catalog name.               */
  int             noconv;      /* View the not convolved image.      */
  int               conv;      /* View the convolved image.          */

  /* Internal parameters: */
  time_t         rawtime;      /* Starting time of the program.      */
  size_t          psf_s0;      /* Side length of psf along axis 0.   */
  size_t          psf_s1;      /* Side length of psf along axis 0.   */
  size_t       numppcols;      /* Number of columns in the above.    */
  size_t         nummock;      /* Number of mock profiles.           */
  float             *psf;      /* Point Spread Function.             */
  double            *cat;      /* Input catalog.                     */
  size_t             cs0;      /* Number of rows in input catalog.   */
  size_t             cs1;      /* Number of columns in input catalog.*/
  double            *log;      /* Log data to be printed.            */
};

#endif
