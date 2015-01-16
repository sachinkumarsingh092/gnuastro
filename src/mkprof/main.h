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

#include "config.h"
#include "commonparams.h"



/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astrmkprof" /* Subpackage executable name. */
#define SPACK_NAME      "MakeProfiles"     /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION
#define LOGFILENAME     SPACK".log"
#define LOGNUMCOLS      6

#define DEGREESTORADIANS M_PI/180.0f

#define SERSICCODE      0
#define MOFFATCODE      1
#define GAUSSIANCODE    2
#define POINTCODE       3
#define MAXIMUMCODE     3

/* Log columns:

   0: ID.
   1: Overlap flux.
   2: Overlap first pixel (first axis) fpixel_o[0].
   3: Overlap first pixel (second axis) fpixel_o[1].
   4: Overlap second pixel (first axis) lpixel_o[0].
   5: Overlap second pixel (second axis) lpixel_o[1].
 */


struct builtqueue
{
  size_t               id;	/* ID of this object.             */
  float              *img;	/* Array of this profile's image. */
  long        fpixel_i[2];	/* First pixel in output image.   */
  long        lpixel_i[2];      /* Last pixel in output image.    */
  long        fpixel_o[2];	/* First pixel in this array.     */
  struct builtqueue *next;	/* Pointer to next element.       */
};





struct uiparams
{
  char          *psfname;      /* Name of PSF FITS name.               */
  char          *catname;      /* Name of catalog of parameters.       */
  int        prepforconv;  /* Shift and expand by size of first psf.   */

  /* Check if all parameters are read (use .def file for
     comparison). The non optional parameters (like the catalog and
     input FITS images that come in from arguments, not options) are
     checked in the args.h files. */

  int        tunitinpset;
  int       numrandomset;
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
  struct uiparams     up;  /* User interface parameters.               */
  struct commonparams cp;  /* Common parameters.                       */

  /* Operating modes: */
  int           psfinimg;  /* ==1: Build PSF profiles in image.        */
  int         individual;  /* ==1: Build all catalog separately.       */

  /* Profiles and noise */
  size_t       numrandom;  /* Number of radom points for integration.  */
  float        tolerance;  /* Accuracy to stop integration.            */
  float        zeropoint;  /* Magnitude of zero point flux.            */
  int           tunitinp;  /* ==1: Truncation unit is in pixels.       */
                           /* ==0: It is in radial parameter.          */
  /* Catalog */
  size_t            fcol;  /* Column specifying profile function.      */
  size_t            xcol;  /* X column of profile center.              */
  size_t            ycol;  /* Y column of profile center.              */
  size_t            rcol;  /* Effective radius of profile.             */
  size_t            ncol;  /* Sersic index column of profile.          */
  size_t            pcol;  /* Position angle column of profile.        */
  size_t            qcol;  /* Axis ratio column of profile.            */
  size_t            mcol;  /* Magnitude column.                        */
  size_t            tcol;  /* Truncation of the profiles.              */

  /* Output */
  long          naxes[2];  /* Size of the output image.                */
  long          shift[2];  /* Shift along axeses position of profiles. */
  size_t      oversample;  /* Oversampling scale.                      */

  /* Internal parameters: */
  time_t         rawtime;  /* Starting time of the program.            */
  double            *cat;  /* Input catalog.                           */
  size_t             cs0;  /* Number of rows in input catalog.         */
  size_t             cs1;  /* Number of columns in input catalog.      */
  double            *log;  /* Log data to be printed.                  */
  int          dir0file1;  /* Output is: ==0: a Dir. ==1: a file.      */
  struct builtqueue  *bq;  /* Top (last) elem of build queue.          */
  pthread_cond_t  qready;  /* bq is ready to be written.               */
  pthread_mutex_t  qlock;  /* Mutex lock to change builtq.             */
  double       halfpixel;  /* Half pixel in oversampled image.         */
};

#endif
