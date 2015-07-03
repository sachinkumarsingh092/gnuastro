/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include "fitsarrayvv.h"
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astmkcatalog" /* Subpackage executable name. */
#define SPACK_NAME      "MakeCatalog"  /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION




/* Columns in the object and clump information tables
   ==================================================

   The information tables are not the final catalog, they keep all the
   necessary information to make the catalogs in a later step. The
   columns of that catalog can be specified by the user.

   The basic idea is that through the first or second pass through the
   image, we fill in this information array. Then we build the catalog
   based on what the user has requested.

   Note: the total area nad total flux of clumps in the objects can be
   found by using the clump information table in the end.

   The X and Y columns have to be immediately after each other. This
   is why I have used the X column indexs +1 for the Y column
   indexs. This is necessary for the WCS conversion and done so
   explicitly, so there is no chance of mistakenly putting two
   noncontiguous columns.
*/
#define OCOLUMNS     15     /* Total number of columns in object info. */
#define OAREA         0     /* Area of this object.                    */
#define ONCLUMPS      1     /* Total number of clumps in this object.  */
#define OTotFlux      2     /* Sum of (flux-sky) in object.            */
#define OFlxWhtX      3     /* Sum of (flux-sky)*x of this object.     */
#define OFlxWhtY    OFlxWhtX+1  /* Sum of (flux-sky)*y of this object. */
#define OFlxWhtRA     5     /* RA of (OFlxWhtX, OFlxWhtY).             */
#define OFlxWhtDec  OFlxWhtRA+1 /* Dec of (OFlxWhtX, OFlxWhtY).        */
#define OAREAC        7     /* Area of clumps in this object.          */
#define OTotFluxC     8     /* Sum of (flux-sky) in object clumps.     */
#define OFlxWhtCX     9     /* Sum of (flux-sky)*x on object clumps.   */
#define OFlxWhtCY   OFlxWhtCX+1 /* Sum of (flux-sky)*y on obj. clumps. */
#define OFlxWhtCRA   11     /* RA of (OFlxWhtCX and OFlxWhtCY).        */
#define OFlxWhtCDec OFlxWhtCRA+1 /* Dec of (OFlxWhtCX and OFlxWhtCY).  */
#define OSKY         13     /* Sum of sky value on this object.        */
#define OSTD         14     /* Sum of sky STD value on this object.    */

#define CCOLUMNS     12     /* Total number of columns in clump info.  */
#define CHOSTOID      0     /* ID of object hosting this clump.        */
#define CINHOSTID     1     /* ID of clump in host object.             */
#define CAREA         2     /* Area of this clump.                     */
#define CFlxWhtX      3     /* Sum of flux*x of this clump.            */
#define CFlxWhtY    CFlxWhtX+1 /* Sum of flux*y of this clump.         */
#define CFlxWhtRA     5     /* RA of (CFlxWhtX, CFlxWhtY).             */
#define CFlxWhtDec  CFlxWhtRA+1 /* Dec of (CFlxWhtX, CFlxWhtY).        */
#define CTotFlux      7     /* Sum of flux in this clump.              */
#define CAveRivFlux   8     /* Sum of flux in rivers around this clump.*/
#define CRivArea      9     /* Sum of flux in rivers around this clump.*/
#define CSKY         10     /* Sum of sky value on this object.        */
#define CSTD         11     /* Sum of sky STD value on this object.    */





/* These macros are used to identify the nature of the appropriate
   column in the output catalog.*/
#define CATID                       1
#define CATHOSTOBJID                2
#define CATIDINHOSTOBJ              3
#define CATNUMCLUMPS                4
#define CATAREA                     5
#define CATCLUMPSAREA               6
#define CATX                        7
#define CATY                        8
#define CATCLUMPSX                  9
#define CATCLUMPSY                 10
#define CATRA                      11
#define CATDEC                     12
#define CATCLUMPSRA                13
#define CATCLUMPSDEC               14
#define CATFLUX                    15
#define CATCLUMPSFLUX              16
#define CATMAGNITUDE               17
#define CATCLUMPSMAGNITUDE         18
#define CATRIVERFLUX               19
#define CATRIVERNUM                20
#define CATSN                      21
#define CATSKY                     22
#define CATSTD                     23



/* Units: */
#define CATDESCRIPTLENGTH         "%-60s"
#define CATUNITCOUNTER            "counter"
#define CATUNITFLUX               "pixel value unit"
#define CATUNITMAG                "magnitude"
#define CATUNITPIXAREA            "pixel area"
#define CATUNITPIXPOS             "pixel position"
#define CATUNITDEGREE             "degree"
#define CATUNITRATIO              "ratio"





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
  int            zeropointset;
  int        skysubtractedset;

  int             intwidthset;
  int           floatwidthset;
  int            accuwidthset;
  int       floatprecisionset;
  int        accuprecisionset;
  int              nsigmagset;

  int                   idset;
  int            hostobjidset;
  int          idinhostobjset;
  int            numclumpsset;
  int                 areaset;
  int           clumpsareaset;
  int                    xset;
  int                    yset;
  int              clumpsxset;
  int              clumpsyset;
  int                   raset;
  int                  decset;
  int             clumpsraset;
  int            clumpsdecset;
  int                 fluxset;
  int           clumpsfluxset;
  int            magnitudeset;
  int      clumpsmagnitudeset;
  int            riverfluxset;
  int             rivernumset;
  int                   snset;
  int                  skyset;
  int                  stdset;
};





struct mkcatalogparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.         */
  struct commonparams      cp;  /* Common parameters.                 */

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

  /* Output: */
  char              *ocatname;  /* File name of object catalog.       */
  char              *ccatname;  /* File name of clump catalog.        */
  int                intwidth;  /* Width of integer columns.          */
  int              floatwidth;  /* Width of float columns.            */
  int               accuwidth;  /* Width of accurate float columns.   */
  int          floatprecision;  /* Floating point precision.          */
  int           accuprecision;  /* Accurate floating point precision. */

  /* Operating mode: */

  /* Internal: */
  time_t              rawtime;  /* Starting time of the program.      */
  double               *oinfo;  /* Information for all the clumps.    */
  double               *cinfo;  /* Information for all the clumps.    */
  size_t           numobjects;  /* Total number of objects.           */
  size_t            numclumps;  /* Total number of clumps.            */
  struct sll       *allcolsll;  /* All the input columns.             */
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
  double              cpscorr;  /* Correction for counts/sec input.   */
  double               maxstd;  /* Maximum STD value (5 sigma mag).   */
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
  double                *info;  /* Either oinfo or cinfo.             */
  size_t                icols;  /* Either OCOLUMNS or CCOLUMNS.       */
  char                 *unitp;  /* Pointer to units array.            */
  char             line[1500];  /* Comment line.                      */
  char       description[500];  /* The description of each row.       */
  int                *intcols;  /* Indexs of integer columns.         */
  int               *accucols;  /* Indexs of accurate columns.        */
};

#endif
