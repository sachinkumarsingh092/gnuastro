/*********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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

#include <pthread.h>
#include <gsl/gsl_rng.h>

#include <gnuastro/fits.h>

#include <gnuastro-internal/options.h>


/* Progarm name macros: */
#define PROGRAM_NAME   "MakeProfiles"      /* Program full name.       */
#define PROGRAM_EXEC   "astmkprof"         /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION



/* Some constants */
#define EPSREL_FOR_INTEG   2
#define DEGREESTORADIANS   M_PI/180.0


/* Modes to interpret coordinates. */
enum coord_modes
{
  MKPROF_MODE_INVALID,          /* For sanity checks.     */

  MKPROF_MODE_IMG,              /* Use image coordinates. */
  MKPROF_MODE_WCS,              /* Use WCS coordinates.   */
};



/* Types of profiles. */
enum profile_types
{
  PROFILE_INVALID,              /* Invalid (=0 by C standard). */

  PROFILE_SERSIC,               /* Sersic profile.             */
  PROFILE_MOFFAT,               /* Moffat Profile.             */
  PROFILE_GAUSSIAN,             /* Gaussian Profile.           */
  PROFILE_POINT,                /* Point profile.              */
  PROFILE_FLAT,                 /* Flat profile.               */
  PROFILE_CIRCUMFERENCE,        /* Circumference profile.      */
  PROFILE_DISTANCE,             /* Elliptical radius of pixel. */

  PROFILE_MAXIMUM_CODE,         /* Just for a sanity check.    */
};
#define MINCIRCUMWIDTH       0.5f



/* Log file:

   0: ID.
   1: Overlap magnitude.
   2: Number of accurate pixels.
   3: Fraction of accurate values.
   4: Is individual file created?   */
#define LOGNUMCOLS      5
#define LOGFILENAME     PROGRAM_EXEC".log"





struct builtqueue
{
  size_t               id;    /* ID of this object.                  */
  int               ispsf;    /* This is a PSF profile.              */
  int            overlaps;    /* ==1: Overlaps with the image.       */
  gal_data_t       *image;    /* Array of this profile's image.      */
  gal_data_t   *overlap_i;    /* Overlap tile over individual array. */
  gal_data_t   *overlap_m;    /* Overlap tile over merged array.     */
  int                func;    /* Profile's radial function.          */
  int        indivcreated;    /* ==1: an individual file is created. */
  size_t          numaccu;    /* Number of accurate pixels.          */
  double         accufrac;    /* Difference of accurate values.      */

  struct builtqueue *next;    /* Pointer to next element.            */
};





struct mkprofparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.              */
  char            *backname;  /* Name of background image file name.      */
  char             *catname;  /* Name of catalog of parameters.           */
  char             *backhdu;  /* HDU of background image.                 */
  size_t             *dsize;  /* Size of the output image.                */
  uint8_t       clearcanvas;  /* Pixels in background image set to zero.  */
  gal_data_t        *kernel;  /* Parameters to define a kernel.           */
  uint8_t        oversample;  /* Oversampling scale.                      */
  uint8_t          psfinimg;  /* ==1: Build PSF profiles in image.        */
  uint8_t        individual;  /* ==1: Build all catalog separately.       */
  uint8_t          nomerged;  /* ==1: Don't make a merged image of all.   */
  char             *typestr;  /* Type of finally merged output image.     */
  size_t          numrandom;  /* Number of radom points for integration.  */
  float           tolerance;  /* Accuracy to stop integration.            */
  uint8_t          tunitinp;  /* ==1: Truncation is in pixels, not radial.*/
  size_t             *shift;  /* Shift along axeses position of profiles. */
  uint8_t       prepforconv;  /* Shift and expand by size of first psf.   */
  float           zeropoint;  /* Magnitude of zero point flux.            */
  float         circumwidth;  /* Width of circumference (inward).         */
  uint8_t           replace;  /* Replace overlaping profile pixel values. */
  uint8_t         magatpeak;  /* Mag only for peak pixel, not all profile.*/
  uint8_t           envseed;  /* Use GSL_RNG_SEED for random seed.        */
  uint8_t              mode;  /* Coordinates in WCS or image standard.    */
  gal_list_str_t      *ccol;  /* Columns that keep coordinates.           */
  char                *fcol;  /* Column specifying profile function.      */
  char                *rcol;  /* Effective radius of profile.             */
  char                *ncol;  /* Sersic index column of profile.          */
  char                *pcol;  /* First Euler angle (X-Z-X order).         */
  char               *p2col;  /* Second Euler angle (X-Z-X order).        */
  char               *p3col;  /* Third Euler angle (X-Z-X order).         */
  char                *qcol;  /* Axis ratio1 (major/2nd dim. radius).     */
  char               *q2col;  /* Axis ratio2 (major/3rd dim. radius).     */
  char                *mcol;  /* Magnitude column.                        */
  char                *tcol;  /* Truncation of the profiles.              */
  uint8_t       mforflatpix;  /* mcol is flat pixel value (f is 4 or 5).  */
  uint8_t  mcolisbrightness;  /* mcol is total brightness, not magnitude. */
  gal_data_t         *crpix;  /* CRPIX FITS header keywords.              */
  gal_data_t         *crval;  /* CRVAL FITS header keywords.              */
  gal_data_t         *cdelt;  /* For CDELTi FITS header keywords.         */
  gal_data_t            *pc;  /* WCS PC matrix.                           */
  gal_data_t         *cunit;  /* Units of each coordinate.                */
  gal_data_t         *ctype;  /* Type of the coordinates.                 */


  /* Output */
  gal_data_t           *out;  /* Output image.                            */
  char              *outdir;  /* Output directory.                        */
  char            *basename;  /* Merged image name with no directory.     */


  /* Processing parameters: */
  size_t                num;  /* The number of profiles.                  */
  double                 *x;  /* X axis position of profile center.       */
  double                 *y;  /* Y axis position of profile center.       */
  double                 *z;  /* Z axis position of profile center.       */
  uint8_t                *f;  /* Profile function code.                   */
  float                  *r;  /* Radius of profile.                       */
  float                  *n;  /* Index of profile.                        */
  float                 *p1;  /* First Euler angle (X-Z-X order).         */
  float                 *p2;  /* Second Euler angle (X-Z-X order).        */
  float                 *p3;  /* Third Euler angle (X-Z-X order).         */
  float                 *q1;  /* Ratio of radius to second axis.          */
  float                 *q2;  /* Ratio of radius to third axis.           */
  float                  *m;  /* Magnitude of profile.                    */
  float                  *t;  /* Truncation distance.                     */
  gsl_rng              *rng;  /* Main instance of random number generator.*/
  const char      *rng_name;  /* Name of random number generator.         */
  unsigned long    rng_seed;  /* Fixed seed of random number generator.   */
  time_t            rawtime;  /* Starting time of the program.            */
  double               *cat;  /* Input catalog.                           */
  gal_data_t           *log;  /* Log data to be printed.                  */
  struct builtqueue     *bq;  /* Top (last) elem of build queue.          */
  pthread_cond_t     qready;  /* bq is ready to be written.               */
  pthread_mutex_t     qlock;  /* Mutex lock to change builtq.             */
  double          halfpixel;  /* Half pixel in oversampled image.         */
  char           *wcsheader;  /* The WCS header information for main img. */
  int            wcsnkeyrec;  /* The number of keywords in the WCS header.*/
  char       *mergedimgname;  /* Name of merged image.                    */
  int                  nwcs;  /* for WCSLIB: no. coord. representations.  */
  struct wcsprm        *wcs;  /* WCS information for this dataset.        */
  size_t               ndim;  /* Number of dimensions (for 'nomerged').   */
                              /* We can't put it in 'out' because it is   */
                              /* meaning ful there.                       */
};

#endif
