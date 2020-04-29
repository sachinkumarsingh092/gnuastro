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
#ifndef UI_H
#define UI_H

/* For common options groups. */
#include <gnuastro-internal/options.h>





/* Option groups particular to this program. */
enum program_args_groups
{
  UI_GROUP_UPPERLIMIT = GAL_OPTIONS_GROUP_AFTER_COMMON,
  UI_GROUP_COLUMNS_IDS,
  UI_GROUP_COLUMNS_POSITION_PIXEL,
  UI_GROUP_COLUMNS_POSITION_WCS,
  UI_GROUP_COLUMNS_BRIGHTNESS,
  UI_GROUP_COLUMNS_MORPHOLOGY,
};





/* Available letters for short options:

   f g k w x y z
   E H J L O R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_CLUMPSCAT       = 'C',         /* General settings. */
  UI_KEY_VALUESFILE      = 'v',
  UI_KEY_CLUMPSFILE      = 'l',
  UI_KEY_INSKY           = 's',
  UI_KEY_INSTD           = 't',
  UI_KEY_ENVSEED         = 'e',

  UI_KEY_IDS             = 'i',         /* Catalog columns. */
  UI_KEY_HOSTOBJID       = 'j',
  UI_KEY_NUMCLUMPS       = 'c',
  UI_KEY_AREA            = 'a',
  UI_KEY_X               = 'x',
  UI_KEY_Y               = 'y',
  UI_KEY_Z               = 'z',
  UI_KEY_RA              = 'r',
  UI_KEY_DEC             = 'd',
  UI_KEY_BRIGHTNESS      = 'b',
  UI_KEY_MAGNITUDE       = 'm',
  UI_KEY_MAGNITUDEERR    = 'G',
  UI_KEY_UPPERLIMITMAG   = 'u',
  UI_KEY_SN              = 'n',
  UI_KEY_SEMIMAJOR       = 'A',
  UI_KEY_SEMIMINOR       = 'B',
  UI_KEY_AXISRATIO       = 'Q',
  UI_KEY_POSITIONANGLE   = 'p',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_VALUESHDU       = 1000,        /* General settings. */
  UI_KEY_CLUMPSHDU,
  UI_KEY_SKYHDU,
  UI_KEY_STDHDU,
  UI_KEY_WITHCLUMPS,
  UI_KEY_FORCEREADSTD,
  UI_KEY_ZEROPOINT,
  UI_KEY_SIGMACLIP,
  UI_KEY_VARIANCE,
  UI_KEY_SUBTRACTSKY,
  UI_KEY_SFMAGNSIGMA,
  UI_KEY_SFMAGAREA,
  UI_KEY_SPECTRUM,
  UI_KEY_INBETWEENINTS,
  UI_KEY_UPMASKFILE,
  UI_KEY_UPMASKHDU,
  UI_KEY_UPNUM,
  UI_KEY_UPRANGE,
  UI_KEY_UPSIGMACLIP,
  UI_KEY_UPNSIGMA,
  UI_KEY_CHECKUPLIM,
  UI_KEY_NOCLUMPSORT,

  UI_KEY_OBJID,                         /* Catalog columns. */
  UI_KEY_IDINHOSTOBJ,
  UI_KEY_AREAXY,
  UI_KEY_CLUMPSAREA,
  UI_KEY_WEIGHTAREA,
  UI_KEY_GEOAREA,
  UI_KEY_GEOAREAXY,
  UI_KEY_GEOX,
  UI_KEY_GEOY,
  UI_KEY_GEOZ,
  UI_KEY_CLUMPSX,
  UI_KEY_CLUMPSY,
  UI_KEY_CLUMPSZ,
  UI_KEY_CLUMPSGEOX,
  UI_KEY_CLUMPSGEOY,
  UI_KEY_CLUMPSGEOZ,
  UI_KEY_MINX,
  UI_KEY_MAXX,
  UI_KEY_MINY,
  UI_KEY_MAXY,
  UI_KEY_MINZ,
  UI_KEY_MAXZ,
  UI_KEY_W1,
  UI_KEY_W2,
  UI_KEY_W3,
  UI_KEY_GEOW1,
  UI_KEY_GEOW2,
  UI_KEY_GEOW3,
  UI_KEY_CLUMPSW1,
  UI_KEY_CLUMPSW2,
  UI_KEY_CLUMPSW3,
  UI_KEY_CLUMPSGEOW1,
  UI_KEY_CLUMPSGEOW2,
  UI_KEY_CLUMPSGEOW3,
  UI_KEY_BRIGHTNESSERR,
  UI_KEY_CLUMPSBRIGHTNESS,
  UI_KEY_BRIGHTNESSNORIVER,
  UI_KEY_MEAN,
  UI_KEY_MEDIAN,
  UI_KEY_CLUMPSMAGNITUDE,
  UI_KEY_UPPERLIMIT,
  UI_KEY_UPPERLIMITONESIGMA,
  UI_KEY_UPPERLIMITSIGMA,
  UI_KEY_UPPERLIMITQUANTILE,
  UI_KEY_UPPERLIMITSKEW,
  UI_KEY_RIVERAVE,
  UI_KEY_RIVERNUM,
  UI_KEY_SKY,
  UI_KEY_STD,
  UI_KEY_SIGCLIPNUMBER,
  UI_KEY_SIGCLIPMEDIAN,
  UI_KEY_SIGCLIPMEAN,
  UI_KEY_SIGCLIPSTD,
  UI_KEY_GEOSEMIMAJOR,
  UI_KEY_GEOSEMIMINOR,
  UI_KEY_GEOAXISRATIO,
  UI_KEY_GEOPOSITIONANGLE,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkcatalogparams *p);

void
ui_free_report(struct mkcatalogparams *p, struct timeval *t1);

#endif
