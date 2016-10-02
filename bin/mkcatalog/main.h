/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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

#include <gsl/gsl_rng.h>
#include <gnuastro/fits.h>

#include <commonparams.h>

/* Progarm name macros: */
#define SPACK           "astmkcatalog" /* Subpackage executable name. */
#define SPACK_NAME      "MakeCatalog"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION




/* Columns in the object and clump information tables
   ==================================================

   The information tables are not the final catalog, they keep all the
   necessary information to make the catalogs in a later step. The
   columns of that catalog can be specified by the user, while the
   information catalog colums are fixed.

   The basic idea is that through the first and second pass through
   the image, we fill in this information array. Then we build the
   catalog based on what the user has requested.

   FIRST ROW IS FLAG ROW: The object/clump IDs begin from one (1), not
   zero (0). So the first row in the oinfo, or cinfo arrays is not
   filled and we need to start from the second row (which as index
   1). So we can use the first row as a flag for the appropriate
   columns. For example, the RA and Dec need to be calculated
   together, so once they are both calculated by RA (for example),
   there is no more need to calculate them again for Dec, or vice
   versa. So in this case, we set the first row of this column to 1 so
   if it is needed again, it is not recalculated.

   NOTE: The X and Y columns should be immediately after each other,
   this is necessary for converting to RA and Dec.

   NOTE: the total area and brightness of clumps in the objects can be
   found by using the clump information table in the end.

   According to the C standard, the first enum variable has a value of
   0 (int), and when none are explicitly set (with an = sign), the
   values of the subsequent enum variables are 1 larger. We want
   column indexs that also start with zero, but setting the values by
   hand manually as preprocessor macros can be buggy (repeated
   numbers). So defining them as an enum is the perfect solution. Any
   future column that is to be added (or if any are removed), we (the
   developers) don't have to worry.
*/
enum objectcols
  {
    OAREA,               /* Area of this object.                    */
    OALLAREA,            /* Object area irrespective of threshold.  */
    ONCLUMPS,            /* Total number of clumps in this object.  */
    OBrightness,         /* Sum of (flux-sky) in object.            */
    OFlxWhtX,            /* Sum of (flux-sky)*x of this object.     */
    OFlxWhtY,            /* Sum of (flux-sky)*y of this object.     */
    OFlxWhtXX,           /* Sum of (flux-sky)*x*x of this object.   */
    OFlxWhtYY,           /* Sum of (flux-sky)*y*y of this object.   */
    OFlxWhtXY,           /* Sum of (flux-sky)*x*y of this object.   */
    OFlxWhtRA,           /* RA of (OFlxWhtX, OFlxWhtY).             */
    OFlxWhtDec,          /* Dec of (OFlxWhtX, OFlxWhtY).            */
    OAREAC,              /* Area of clumps in this object.          */
    OBrightnessC,        /* Brightness  in object clumps.           */
    OUpperBright,        /* Upper limit brightness for this object. */
    OFlxWhtCX,           /* Sum of (flux-sky)*x on object clumps.   */
    OFlxWhtCY,           /* Sum of (flux-sky)*y on obj. clumps.     */
    OFlxWhtCXX,          /* Sum of (flux-sky)*x*x on object clumps. */
    OFlxWhtCYY,          /* Sum of (flux-sky)*y*y on obj. clumps.   */
    OFlxWhtCXY,          /* Sum of (flux-sky)*x*y on obj. clumps.   */
    OPOSSHIFTX,          /* Shift in X to avoid rounding errors.    */
    OPOSSHIFTY,          /* Shift in Y to avoid rounding errors.    */
    OFlxWhtCRA,          /* RA of (OFlxWhtCX and OFlxWhtCY).        */
    OFlxWhtCDec,         /* Dec of (OFlxWhtCX and OFlxWhtCY).       */
    OSKY,                /* Sum of sky value on this object.        */
    OSTD,                /* Sum of sky STD value on this object.    */
    OPosBright,          /* Sum of positive image pixels for wht.   */
    OPosBrightC,         /* Sum of positive image pixels for wht.   */
    OGeoX,               /* Geometric center of object in X.        */
    OGeoY,               /* Geometric center of object in Y.        */
    OGeoRA,              /* RA of Geometric center of object.       */
    OGeoDec,             /* Dec of Geometric center of object.      */
    OGeoXX,              /* Second order geometric variable: X*X.   */
    OGeoYY,              /* Second order geometric variable: Y*Y.   */
    OGeoXY,              /* Second order geometric variable: X*Y.   */
    OGeoCX,              /* Geometric center of clumps in object X. */
    OGeoCY,              /* Geometric center of clumps in object Y. */
    OGeoCRA,             /* Geometric center of clumps in obj. RA.  */
    OGeoCDec,            /* Geometric center of clumps in obj. Dec. */
    OGeoCXX,             /* Second order geometric variable: X*X.   */
    OGeoCYY,             /* Second order geometric variable: Y*Y.   */
    OGeoCXY,             /* Second order geometric variable: X*Y.   */

    OCOLUMNS,            /* Keep this last: total number of columns.*/
  };

enum clumpcols
  {
    CHOSTOID,            /* ID of object hosting this clump.        */
    CINHOSTID,           /* ID of clump in host object.             */
    CAREA,               /* Area of this clump.                     */
    CALLAREA,            /* Area of clump irrespective of threshold.*/
    CFlxWhtX,            /* Sum of flux*x of this clump.            */
    CFlxWhtY,            /* Sum of flux*y of this clump.            */
    CFlxWhtXX,           /* Sum of flux*x*x of this clump.          */
    CFlxWhtYY,           /* Sum of flux*y*y of this clump.          */
    CFlxWhtXY,           /* Sum of flux*x*y of this clump.          */
    CPOSSHIFTX,          /* Shift in X to avoid rounding errors.    */
    CPOSSHIFTY,          /* Shift in Y to avoid rounding errors.    */
    CFlxWhtRA,           /* ra of (CFlxWhtX, CFlxWhtY).             */
    CFlxWhtDec,          /* Dec of (CFlxWhtX, CFlxWhtY).            */
    CBrightness,         /* River subtracted brightness.            */
    CNoRiverBrightness,  /* Sky (not river) subtracted brightness.  */
    CUpperBright,        /* Upper limit brightness for this clump.  */
    CRivAve,             /* Average value in rivers around clump.   */
    CRivArea,            /* Num river pixels around this clump.     */
    CSKY,                /* Sum of sky value on this object.        */
    CSTD,                /* Sum of sky STD value on this object.    */
    CPosBright,          /* Sum of positive image pixels for wht.   */
    CGeoX,               /* Geometric center of clump in X.         */
    CGeoY,               /* Geometric center of clump in Y.         */
    CGeoRA,              /* RA of Geometric center of clump.        */
    CGeoDec,             /* Dec of Geometric center of clump.       */
    CGeoXX,              /* Second order geometric moment.          */
    CGeoYY,              /* Second order geometric moment.          */
    CGeoXY,              /* Second order geometric moment.          */

    CCOLUMNS,            /* Keep this last: total number of columns.*/
  };





/* Codes for output columns
   ========================

   The user can ask for any columns in any order, to enable that, we
   need a fixed number ID for each column. Again (similar to the
   argument for the information array columns), it is just important
   that these numbers be separate and defining them by hand can cause
   bugs. So we are defining them as enums here. */
enum outcols
  {
    CATID,
    CATHOSTOBJID,
    CATIDINHOSTOBJ,
    CATNUMCLUMPS,
    CATAREA,
    CATCLUMPSAREA,
    CATX,
    CATY,
    CATGEOX,
    CATGEOY,
    CATCLUMPSX,
    CATCLUMPSY,
    CATCLUMPSGEOX,
    CATCLUMPSGEOY,
    CATRA,
    CATDEC,
    CATGEORA,
    CATGEODEC,
    CATCLUMPSRA,
    CATCLUMPSDEC,
    CATCLUMPSGEORA,
    CATCLUMPSGEODEC,
    CATBRIGHTNESS,
    CATCLUMPSBRIGHTNESS,
    CATNORIVERBRIGHTNESS,
    CATMAGNITUDE,
    CATMAGNITUDEERR,
    CATCLUMPSMAGNITUDE,
    CATUPPERLIMITMAG,
    CATRIVERAVE,
    CATRIVERNUM,
    CATSN,
    CATSKY,
    CATSTD,
    CATSEMIMAJOR,
    CATSEMIMINOR,
    CATPOSITIONANGLE,
    CATGEOSEMIMAJOR,
    CATGEOSEMIMINOR,
    CATGEOPOSITIONANGLE,
  };





struct uiparams
{
  char             *inputname;  /* Name of input file.               */
  char              *maskname;  /* Name of masked file.              */
  char                  *mhdu;  /* HDU of mask image.                */
  char           *objlabsname;  /* Name of object labels file.       */
  char                *objhdu;  /* HDU of object labels image.       */
  char         *clumplabsname;  /* Name of clump labels file .       */
  char              *clumphdu;  /* HDU of clump labels image.        */
  char               *skyname;  /* Sky value image file name.        */
  char                *skyhdu;  /* Sky HDU name.                     */
  char               *stdname;  /* Sky STD value image file name.    */
  char                *stdhdu;  /* Sky STD HDU name.                 */
  char            *upmaskname;  /* Name of mask for upper magnitude. */
  char             *upmaskhdu;  /* HDU Name of mask for upper mag.   */

  int             masknameset;
  int                 mhduset;
  int          objlabsnameset;
  int               objhduset;
  int        clumplabsnameset;
  int             clumphduset;
  int              skynameset;
  int               skyhduset;
  int              stdnameset;
  int               stdhduset;
  int              envseedset;

  int            zeropointset;
  int        skysubtractedset;
  int            thresholdset;

  int             intwidthset;
  int           floatwidthset;
  int            accuwidthset;
  int       floatprecisionset;
  int        accuprecisionset;
  int              nsigmagset;

  int           upmasknameset;
  int            upmaskhduset;
  int                upnumset;
  int        upsclipmultipset;
  int          upsclipaccuset;
  int             upnsigmaset;

  int                   idset;
  int            hostobjidset;
  int          idinhostobjset;
  int            numclumpsset;
  int                 areaset;
  int           clumpsareaset;
  int                    xset;
  int                    yset;
  int                 geoxset;
  int                 geoyset;
  int              clumpsxset;
  int              clumpsyset;
  int           clumpsgeoxset;
  int           clumpsgeoyset;
  int                   raset;
  int                  decset;
  int                georaset;
  int               geodecset;
  int             clumpsraset;
  int            clumpsdecset;
  int          clumpsgeoraset;
  int         clumpsgeodecset;
  int           brightnessset;
  int     clumpsbrightnessset;
  int    noriverbrightnessset;
  int            magnitudeset;
  int         magnitudeerrset;
  int      clumpsmagnitudeset;
  int        upperlimitmagset;
  int             riveraveset;
  int             rivernumset;
  int                   snset;
  int                  skyset;
  int                  stdset;
  int            semimajorset;
  int            semiminorset;
  int        positionangleset;
  int         geosemimajorset;
  int         geosemiminorset;
  int     geopositionangleset;
};





struct mkcatalogparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.         */
  struct gal_commonparams  cp;  /* Common parameters.                 */

  /* Input: */
  float                  *img;  /* Input image.                       */
  long               *objects;  /* Object labels on each pixel.       */
  long                *clumps;  /* Clump labels on each pixel.        */
  float                  *sky;  /* Sky value on each pixel.           */
  float                  *std;  /* Sky STD value on each pixel.       */
  int                    nwcs;  /* Number of WCS structures.          */
  struct wcsprm          *wcs;  /* Pointer to WCS structures.         */
  size_t                   s0;  /* Size of input (first C axis).      */
  size_t                   s1;  /* Size of input (second C axis).     */
  float             zeropoint;  /* Zeropoint magnitude of input.      */
  int           skysubtracted;  /* Input is already sky subtracted.   */
  double              nsigmag;  /* Multiple of Sky STD to report mag. */
  double            threshold;  /* Only pixels larger than this *STD. */
  int                 envseed;  /* Use GSL_RNG_SEED for random seed.  */

  /* Output: */
  char              *ocatname;  /* File name of object catalog.       */
  char              *ccatname;  /* File name of clump catalog.        */
  int                intwidth;  /* Width of integer columns.          */
  int              floatwidth;  /* Width of float columns.            */
  int               accuwidth;  /* Width of accurate float columns.   */
  int          floatprecision;  /* Floating point precision.          */
  int           accuprecision;  /* Accurate floating point precision. */

  /* Upper limit: */
  long                *upmask;  /* Mask array, only for upper liimit. */
  size_t                upnum;  /* Number of random samples.          */
  float         upsclipmultip;  /* Sigma clip multiple for std calc.  */
  float           upsclipaccu;  /* Sigma clip accuracy for std calc.  */
  float              upnsigma;  /* Multiple of sigma for upper-limit. */

  /* Operating mode: */

  /* Internal: */
  time_t              rawtime;  /* Starting time of the program.      */
  double               *oinfo;  /* Information for all the clumps.    */
  double               *cinfo;  /* Information for all the clumps.    */
  size_t           numobjects;  /* Total number of objects.           */
  size_t            numclumps;  /* Total number of clumps.            */
  struct gal_linkedlist_sll *allcolsll; /* All the input columns.     */
  size_t             *allcols;  /* Array keeping all the input cols.  */
  size_t             allncols;  /* Total number of input columns.     */
  size_t             *objcols;  /* Array of objcolsll.                */
  size_t           *clumpcols;  /* Array of clumpcolsll.              */
  size_t             objncols;  /* Num. columns in objects catalog.   */
  size_t           clumpncols;  /* Num. columns in clumps catalog.    */
  double              *objcat;  /* Output object catalog.             */
  double            *clumpcat;  /* Output clump catalog.              */
  size_t            objcurcol;  /* Current column in object catalog.  */
  size_t          clumpcurcol;  /* Current column in clump catalog.   */
  float                minstd;  /* For correction for counts/sec.     */
  float                medstd;  /* Median STD value (5 sigma mag).    */
  double              cpscorr;  /* Correction for counts/sec.         */
  double                detsn;  /* Detection S/N threshold.           */
  double              clumpsn;  /* Clumps S/N threshold.              */

  /* For going through the rows: */
  size_t               curcol;  /* Current column in the catalog.     */
  size_t           intcounter;  /* Counter of the integer columns.    */
  size_t          accucounter;  /* Counter of the integer columns.    */
  int              obj0clump1;  /* Is this an object or clump catalog.*/
  char                  *name;  /* `Object' or `Clump'.               */
  double                 *cat;  /* Either objcat or clumpcat.         */
  char              *filename;  /* Either ocatname or ccatname.       */
  size_t                  num;  /* Either numobjects or numclumps.    */
  size_t              numcols;  /* Either objncols or clumpncols.     */
  double                *info;  /* Pointer to either oinfo or cinfo.  */
  size_t                icols;  /* Either OCOLUMNS or CCOLUMNS.       */
  char                 *unitp;  /* Pointer to units array.            */
  size_t            xshiftcol;  /* Column to correct/shift positions. */
  size_t            yshiftcol;  /* Column to correct/shift positions. */
  char             line[1500];  /* Comment line.                      */
  char       description[500];  /* The description of each row.       */
  int                *intcols;  /* Indexs of integer columns.         */
  int               *accucols;  /* Indexs of accurate columns.        */
};

#endif
