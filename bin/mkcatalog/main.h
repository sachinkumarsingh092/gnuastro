/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <gsl/gsl_rng.h>
#include <gnuastro/data.h>
#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "MakeCatalog"  /* Program full name.       */
#define PROGRAM_EXEC   "astmkcatalog" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION




/* Multiple of given number to stop searching for upper-limit magnitude. */
#define MKCATALOG_UPPERLIMIT_MINIMUM_NUM     20
#define MKCATALOG_UPPERLIMIT_MAXFAILS_MULTIP 10


/* Unit string to use if values dataset doesn't have any. */
#define MKCATALOG_NO_UNIT "input-units"



/* Intermediate/raw array elements
   ===============================

   Commonly, several high-level calculations need the same low-level
   measurements. So to avoid having to do these low-level calculations on
   each pixel multiple tiles, each thread/object will have one array of
   intermediate values which will be filled in the pass over the
   pixels. After this intermediate array is filled, and we don't need to
   pass over the pixels any more, we will use the intermediate values to
   derive the higher-level steps.

   According to the C standard, the first enum variable has a value of 0
   (int), and when none are explicitly set (with an = sign), the values of
   the subsequent enum variables are 1 larger. We want column indexs that
   also start with zero, but setting the values by hand manually as
   preprocessor macros can be buggy (repeated numbers). So defining them as
   an enum is the perfect solution. Any future column that is to be added
   (or if any are removed), we (the developers) don't have to worry.

   --------------------- POSITIONS IN FITS standard ---------------------
   Like the final output, positions in this intermediate array are also in
   the FITS standard (fastest dimension is first). */
enum objectcols
  {
    OCOL_NUMALL,         /* Number of all pixels with this label.     */
    OCOL_NUMALLXY,       /* Number of all pixels in first two dims.   */
    OCOL_NUM,            /* Number of values used in this object.     */
    OCOL_NUMXY,          /* Number of values in the first two dims.   */
    OCOL_SUM,            /* Sum of (value-sky) in object.             */
    OCOL_SUM_VAR,        /* Variance including values (not just sky). */
    OCOL_MEDIAN,         /* Median of value in object.                */
    OCOL_SIGCLIPNUM,     /* Sigma-clipped mean of this object.        */
    OCOL_SIGCLIPSTD,     /* Sigma-clipped mean of this object.        */
    OCOL_SIGCLIPMEAN,    /* Sigma-clipped mean of this object.        */
    OCOL_SIGCLIPMEDIAN,  /* Sigma-clipped mean of this object.        */
    OCOL_VX,             /* Sum of (value-sky) * x.                   */
    OCOL_VY,             /* Sum of (value-sky) * y.                   */
    OCOL_VZ,             /* Sum of (value-sky) * z.                   */
    OCOL_VXX,            /* Sum of (value-sky) * x * x.               */
    OCOL_VYY,            /* Sum of (value-sky) * y * y.               */
    OCOL_VXY,            /* Sum of (value-sky) * x * y.               */
    OCOL_SUMSKY,         /* Sum of sky value on this object.          */
    OCOL_NUMSKY,         /* Number of sky value on this object.       */
    OCOL_SUMVAR,         /* Sum of sky variance value on this object. */
    OCOL_NUMVAR,         /* Number of sky value on this object.       */
    OCOL_SUMWHT,         /* Sum of positive image pixels.             */
    OCOL_NUMWHT,         /* Number of positive pixels used for wht.   */
    OCOL_GX,             /* Geometric center of object in X.          */
    OCOL_GY,             /* Geometric center of object in Y.          */
    OCOL_GZ,             /* Geometric center of object in Z.          */
    OCOL_GXX,            /* Second order geometric variable: X*X.     */
    OCOL_GYY,            /* Second order geometric variable: Y*Y.     */
    OCOL_GXY,            /* Second order geometric variable: X*Y.     */
    OCOL_UPPERLIMIT_B,   /* Upper limit brightness.                   */
    OCOL_UPPERLIMIT_S,   /* Upper limit one-sigma value.              */
    OCOL_UPPERLIMIT_Q,   /* Quantile of object in random distribution.*/
    OCOL_UPPERLIMIT_SKEW,/* (Mean-Median)/STD of random distribution. */
    OCOL_C_NUMALL,       /* Value independent no. of pixels in clumps.*/
    OCOL_C_NUM,          /* Area of clumps in this object.            */
    OCOL_C_SUM,          /* Brightness in object clumps.              */
    OCOL_C_VX,           /* Sum of (value-sky)*x on clumps.           */
    OCOL_C_VY,           /* Sum of (value-sky)*y on obj. clumps.      */
    OCOL_C_VZ,           /* Sum of (value-sky)*z on obj. clumps.      */
    OCOL_C_GX,           /* Geometric center of clumps in object X.   */
    OCOL_C_GY,           /* Geometric center of clumps in object Y.   */
    OCOL_C_GZ,           /* Geometric center of clumps in object Z.   */
    OCOL_C_SUMWHT,       /* Sum of positive image pixels for wht.     */
    OCOL_C_NUMWHT,       /* Num of positive image pixels for wht.     */

    OCOL_NUMCOLS,        /* SHOULD BE LAST: total number of columns.  */
  };

enum clumpcols
  {
    CCOL_NUMALL,         /* Number of pixels in clump.                */
    CCOL_NUMALLXY,       /* Number of pixels in first two dims.       */
    CCOL_NUM,            /* Number of values used in clump.           */
    CCOL_NUMXY,          /* Number of values only in first two dims.  */
    CCOL_SUM,            /* River subtracted brightness.              */
    CCOL_SUM_VAR,        /* Variance including values (not just sky). */
    CCOL_MEDIAN,         /* Median of values in clump.                */
    CCOL_SIGCLIPNUM,     /* Sigma-clipped mean of this clump.         */
    CCOL_SIGCLIPSTD,     /* Sigma-clipped mean of this clump.         */
    CCOL_SIGCLIPMEAN,    /* Sigma-clipped mean of this clump.         */
    CCOL_SIGCLIPMEDIAN,  /* Sigma-clipped mean of this clump.         */
    CCOL_RIV_NUM,        /* Num river pixels around this clump.       */
    CCOL_RIV_SUM,        /* Sum of rivers around clump.               */
    CCOL_RIV_SUM_VAR,    /* Variance of sum (for error measurements). */
    CCOL_VX,             /* Sum of (value-sky) * x.                   */
    CCOL_VY,             /* Sum of (value-sky) * y.                   */
    CCOL_VZ,             /* Sum of (value-sky) * z.                   */
    CCOL_VXX,            /* Sum of flux*x*x of this clump.            */
    CCOL_VYY,            /* Sum of flux*y*y of this clump.            */
    CCOL_VXY,            /* Sum of flux*x*y of this clump.            */
    CCOL_SUMSKY,         /* Sum of sky value on this clump.           */
    CCOL_NUMSKY,         /* Number of sky value on this clump.        */
    CCOL_SUMVAR,         /* Sum of sky variance value on this clump.  */
    CCOL_NUMVAR,         /* Number of sky variance value on this clump.*/
    CCOL_SUMWHT,         /* Sum of positive image pixels for wht.     */
    CCOL_NUMWHT,         /* Num of positive image pixels for wht.     */
    CCOL_GX,             /* Geometric center of clump in X.           */
    CCOL_GY,             /* Geometric center of clump in Y.           */
    CCOL_GZ,             /* Geometric center of clump in Y.           */
    CCOL_GXX,            /* Second order geometric moment.            */
    CCOL_GYY,            /* Second order geometric moment.            */
    CCOL_GXY,            /* Second order geometric moment.            */
    CCOL_MINX,           /* Minimum X value of clump.                 */
    CCOL_MAXX,           /* Maximum X value of clump.                 */
    CCOL_MINY,           /* Minimum Y value of clump.                 */
    CCOL_MAXY,           /* Maximum Y value of clump.                 */
    CCOL_MINZ,           /* Minimum Z value of clump.                 */
    CCOL_MAXZ,           /* Maximum Z value of clump.                 */
    CCOL_UPPERLIMIT_B,   /* Upper limit brightness.                   */
    CCOL_UPPERLIMIT_S,   /* Upper limit one-sigma value.              */
    CCOL_UPPERLIMIT_Q,   /* Quantile of object in random distribution.*/
    CCOL_UPPERLIMIT_SKEW,/* (Mean-Median)/STD of random distribution. */

    CCOL_NUMCOLS,        /* SHOULD BE LAST: total number of columns.  */
  };







/* Main program parameters structure */
struct mkcatalogparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.            */
  gal_list_i32_t   *columnids;  /* The desired column codes.            */
  char           *objectsfile;  /* Input filename.                      */
  char            *valuesfile;  /* File name of objects file.           */
  char             *valueshdu;  /* HDU of objects image.                */
  char            *clumpsfile;  /* File name of objects file.           */
  char             *clumpshdu;  /* HDU of objects image.                */
  char               *skyfile;  /* File name of sky file.               */
  char                *skyhdu;  /* HDU of sky image.                    */
  char               *stdfile;  /* File name of sky STD file.           */
  char                *stdhdu;  /* HDU of sky STD image.                */

  uint8_t           clumpscat;  /* ==1: create clumps catalog.          */
  uint8_t         noclumpsort;  /* Don't sort the clumps catalog.       */
  float             zeropoint;  /* Zero-point magnitude of object.      */
  uint8_t            variance;  /* Input STD file is actually variance. */
  uint8_t        forcereadstd;  /* Read STD even if not needed.         */
  uint8_t         subtractsky;  /* ==1: subtract the Sky from values.   */
  float           sfmagnsigma;  /* Surface brightness multiple of sigma.*/
  float             sfmagarea;  /* Surface brightness area (arcsec^2).  */
  uint8_t            spectrum;  /* Object spectrum for 3D datasets.     */
  uint8_t       inbetweenints;  /* Keep rows (integer ids) with no labels. */
  double         sigmaclip[2];  /* Sigma clip column settings.          */

  char            *upmaskfile;  /* Name of upper limit mask file.       */
  char             *upmaskhdu;  /* HDU of upper limit mask file.        */
  size_t                upnum;  /* Number of upper-limit random samples.*/
  size_t             *uprange;  /* Range of random pos. around target.  */
  uint8_t             envseed;  /* Use the environment for random seed. */
  double       upsigmaclip[2];  /* Sigma clip to measure upper limit.   */
  float              upnsigma;  /* Multiple of sigma to define up-lim.  */
  int32_t       checkuplim[2];  /* Object & clump ID to check dist.     */

  /* Internal. */
  char           *relabclumps;  /* Name of new file for clump labels.   */
  time_t              rawtime;  /* Starting time of the program.        */
  gal_data_t          *values;  /* Input.                               */
  gal_data_t         *objects;  /* Object labels.                       */
  gal_data_t          *clumps;  /* Clump labels.                        */
  gal_data_t             *sky;  /* Sky.                                 */
  gal_data_t             *std;  /* Sky standard deviation.              */
  gal_data_t          *upmask;  /* Upper limit magnitude mask.          */
  float                medstd;  /* Median standard deviation value.     */
  float               cpscorr;  /* Counts-per-second correction.        */
  int32_t            *outlabs;  /* Labels in output catalog (when necessary) */
  size_t           numobjects;  /* Number of object labels in image.    */
  float               clumpsn;  /* Clump S/N threshold.                 */
  size_t            numclumps;  /* Number of clumps in image.           */
  gal_data_t      *objectcols;  /* Output columns for the objects.      */
  gal_data_t       *clumpcols;  /* Output columns for the clumps.       */
  gal_data_t           *tiles;  /* Tiles to cover each object.          */
  char            *objectsout;  /* Output objects catalog.              */
  char             *clumpsout;  /* Output clumps catalog.               */
  char            *upcheckout;  /* Name of upperlimit check table.      */
  uint8_t             *oiflag;  /* Intermediate flags for objects.      */
  uint8_t             *ciflag;  /* Intermediate flags for clumps.       */
  pthread_mutex_t       mutex;  /* Mutex to change the total numbers.   */
  size_t      clumprowsfilled;  /* No. filled clump rows at this moment.*/
  gsl_rng                *rng;  /* Main random number generator.        */
  unsigned long int  rng_seed;  /* Random number generator seed.        */
  const char        *rng_name;  /* Name of random number generator.     */
  size_t               rngmin;  /* Minimum possible value of RNG.       */
  size_t              rngdiff;  /* Difference of RNG max and min.       */
  uint8_t      uprangewarning;  /* A warning must be printed.           */
  size_t         *hostobjid_c;  /* To sort the clumps table by Obj.ID.  */
  size_t         *numclumps_c;  /* To sort the clumps table by Obj.ID.  */
  gal_data_t   *specsliceinfo;  /* Slice information for spectra.       */
  gal_data_t         *spectra;  /* Array of datasets containing spectra.*/

  char        *usedvaluesfile;  /* Ptr to final name used for values.   */
  char        *usedclumpsfile;  /* Ptr to final name used for clumps.   */
  char           *usedskyfile;  /* Ptr to final fname used for sky.     */
  char           *usedstdfile;  /* Ptr to final name used for sky std.  */

  gal_data_t          *wcs_vo;  /* Object RA-Dec flux weighted X, Y.    */
  gal_data_t          *wcs_vc;  /* Clump RA-Dec flux weighted X, Y.     */
  gal_data_t          *wcs_go;  /* Object RA-Dec geometric X,Y.         */
  gal_data_t          *wcs_gc;  /* Clump RA-Dec geometric X, Y.         */
  gal_data_t         *wcs_vcc;  /* All clumps RA-Dec flx. wht. X, Y.    */
  gal_data_t         *wcs_gcc;  /* All clumps RA-Dec geometric X, Y.    */

  char                **ctype;  /* Type of WCS axis.                    */

  uint8_t            hasblank;  /* Dataset has blank values.            */
  uint8_t              hasmag;  /* Catalog has magnitude columns.       */
  uint8_t          upperlimit;  /* Calculate upper limit magnitude.     */
};

#endif
