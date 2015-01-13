/*********************************************************************
mkprof (MakeProfiles) - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
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

#include <pthread.h>

#include "commonparams.h"



/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astrmkprof" /* Subpackage executable name. */
#define SPACK_NAME      "MakeProfiles"     /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION
#define LOGFILENAME     SPACK".log"
#define PSFOUTNAME      "PSF.fits"
#define LOGNUMCOLS      3





struct built
{
  size_t              id;	/* ID of this profile in catalog. */
  float               *a;	/* Array of this profile's image. */
  double      overlapmag;	/* Magnitude in overlap.          */
  long       fpixel_m[2];	/* First overlap pixel on mock.   */
  long       lpixel_m[2];	/* Last overlap pixel on mock.    */
  long       fpixel_o[2];	/* First overlap pixel on output. */
  long       lpixel_o[2];	/* Last overlap pixel on output.  */
  struct built     *next;	/* Pointer to next element.       */
};





struct uiparams
{
  char          *psfname;      /* Name of PSF FITS name.             */
  char          *catname;      /* Name of catalog of parameters.     */

  int           tunitinp;  /* ==1: Truncation unit is in pixels.     */
                           /* ==0: It is in radial parameter.        */
  int        prepforconv;  /* Shift and expand by size of first psf. */
  size_t          xshift;  /* Shift x position of profiles.          */
  size_t          yshift;  /* Shift y position of profiles.          */

  /* Check if all parameters are read (use .def file for
     comparison). The non optional parameters (like the catalog and
     input FITS images that come in from arguments, not options) are
     checked in the args.h files. */

  int        tunitinpset;
  int       toleranceset;
  int       zeropointset;
  int          xshiftset;
  int          yshiftset;
  int     prepforconvset;

  int            fcolset;
  int            xcolset;
  int            ycolset;
  int            rcolset;
  int            ncolset;
  int            pcolset;
  int            qcolset;
  int            mcolset;
  int            tcolset;

  int          naxis1set;
  int          naxis2set;
  int      oversampleset;
};





struct mkprofparams
{
  /* Other structures: */
  struct uiparams     up;  /* User interface parameters.             */
  struct commonparams cp;  /* Common parameters.                     */

  /* Operating modes: */
  int            mginimg;  /* ==1: Build Mof and Gaus in image.      */
  int         individual;  /* ==1: Build all catalog separately.     */

  /* Profiles and noise */
  float        tolerance;  /* Accuracy to stop integration.          */
  float        zeropoint;  /* Magnitude of zero point flux.          */

  /* Catalog */
  size_t            fcol;  /* Column specifying profile function.    */
  size_t            xcol;  /* X column of profile center.            */
  size_t            ycol;  /* Y column of profile center.            */
  size_t            rcol;  /* Effective radius of profile.           */
  size_t            ncol;  /* Sersic index column of profile.        */
  size_t            pcol;  /* Position angle column of profile.      */
  size_t            qcol;  /* Axis ratio column of profile.          */
  size_t            mcol;  /* Magnitude column.                      */
  size_t            tcol;  /* Truncation of the profiles.            */

  /* Output */
  size_t              s0;  /* C standard axis 0 size.                */
  size_t              s1;  /* C standard axis 1 size.                */
  size_t      oversample;  /* Oversampling scale.                    */

  /* Internal parameters: */
  double             dos;  /* Oversampling in double type.           */
  time_t         rawtime;  /* Starting time of the program.          */
  double            *cat;  /* Input catalog.                         */
  size_t             cs0;  /* Number of rows in input catalog.       */
  size_t             cs1;  /* Number of columns in input catalog.    */
  double            *log;  /* Log data to be printed.                */
  struct built   *builtq;  /* Bottom (first) elem of build queue.    */
  pthread_cond_t  qready;  /* builtq is ready to be written.         */
  pthread_mutex_t  qlock;  /* Mutex lock to change builtq.           */
};

#endif
