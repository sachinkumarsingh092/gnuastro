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
#ifndef ARGS_H
#define ARGS_H






/* Array of acceptable options. */
struct argp_option program_options[] =
  {
    /* Input options. */
    {
      "objectsfile",
      UI_KEY_OBJECTSFILE,
      "STR",
      0,
      "Image containing object/detection labels.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->objectsfile,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "objectshdu",
      UI_KEY_OBJECTSHDU,
      "STR",
      0,
      "Object image extension name or number.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->objectshdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "clumpsfile",
      UI_KEY_CLUMPSFILE,
      "STR",
      0,
      "Image containing clump labels.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->clumpsfile,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "clumpshdu",
      UI_KEY_CLUMPSHDU,
      "STR",
      0,
      "Clump image extension name or number.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->clumpshdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "skyfile",
      UI_KEY_SKYFILE,
      "STR",
      0,
      "Image containing sky values.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->skyfile,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "skyhdu",
      UI_KEY_SKYHDU,
      "STR",
      0,
      "Sky image extension name or number.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->skyhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "stdfile",
      UI_KEY_STDFILE,
      "STR",
      0,
      "Image containing sky STD values.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->stdfile,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "stdhdu",
      UI_KEY_STDHDU,
      "STR",
      0,
      "Sky image extension name or number.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->stdhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "zeropoint",
      UI_KEY_ZEROPOINT,
      "FLT",
      0,
      "Zeropoint magnitude of input dataset.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->zeropoint,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "skysubtracted",
      UI_KEY_SKYSUBTRACTED,
      0,
      0,
      "Input is already sky subtracted (for S/N).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->skysubtracted,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "threshold",
      UI_KEY_THRESHOLD,
      "FLT",
      0,
      "Use pixels more than this multiple of STD.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->threshold,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Output. */
    {
      "sfmagnsigma",
      UI_KEY_SFMAGNSIGMA,
      "FLT",
      0,
      "Surface brightness multiple of Sky STD.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->sfmagnsigma,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sfmagarea",
      UI_KEY_SFMAGAREA,
      "FLT",
      0,
      "Surface brightness area (in arcseconds^2).",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->sfmagarea,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Upper limit magnitude configurations. */
    {
      0, 0, 0, 0,
      "Upper limit magnitude:",
      ARGS_GROUP_UPPERLIMIT
    },
    {
      "upmaskfile",
      UI_KEY_UPMASKFILE,
      "STR",
      0,
      "Mask image file name only for upper limit.",
      ARGS_GROUP_UPPERLIMIT,
      &p->upmaskfile,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "upmaskhdu",
      UI_KEY_UPMASKHDU,
      "STR",
      0,
      "Mask image HDU only for upper limit.",
      ARGS_GROUP_UPPERLIMIT,
      &p->upmaskhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "upnum",
      UI_KEY_UPNUM,
      "INT",
      0,
      "Number of randomly positioned samples",
      ARGS_GROUP_UPPERLIMIT,
      &p->upnum,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "envseed",
      UI_KEY_ENVSEED,
      0,
      0,
      "Use GSL_RNG_SEED environment variable for seed.",
      ARGS_GROUP_UPPERLIMIT,
      &p->envseed,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "upsigmaclip",
      UI_KEY_UPSIGMACLIP,
      "FLT,FLT",
      0,
      "Sigma multiple and, tolerance or number.",
      ARGS_GROUP_UPPERLIMIT,
      &p->upsigmaclip,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_sigma_clip
    },
    {
      "upnsigma",
      UI_KEY_UPNSIGMA,
      "FLT",
      0,
      "Multiple of sigma to define upperlimit.",
      ARGS_GROUP_UPPERLIMIT,
      &p->upnsigma,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Output columns: IMPORTANT: Throughout MakeCatalogs, the order of
       columns is based on this. So if you make changes/additions, please
       update all the other places the columns are searched. */
    {
      0, 0, 0, 0,
      "Output catalog columns",
      ARGS_GROUP_COLUMNS
    },
    {  /* `ids' is not a unique column, it is a combination of several
          columns. */
      "ids",
      UI_KEY_IDS,
      0,
      0,
      "All IDs of objects and clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "objid",
      UI_KEY_OBJID,
      0,
      0,
      "Object label/ID.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "hostobjid",
      UI_KEY_HOSTOBJID,
      0,
      0,
      "ID of object hosting this clump.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "idinhostobj",
      UI_KEY_IDINHOSTOBJ,
      0,
      0,
      "ID of clump in host object.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "numclumps",
      UI_KEY_NUMCLUMPS,
      0,
      0,
      "Number of clumps in this object.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "area",
      UI_KEY_AREA,
      0,
      0,
      "Number of pixels in clump or object.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsarea",
      UI_KEY_CLUMPSAREA,
      0,
      0,
      "Sum of all clump areas in an object.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "weightarea",
      UI_KEY_WEIGHTAREA,
      0,
      0,
      "Area used for flux weighted positions.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "x",
      UI_KEY_X,
      0,
      0,
      "Flux weighted center in first FITS axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "y",
      UI_KEY_Y,
      0,
      0,
      "Flux weighted center in second FITS axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geox",
      UI_KEY_GEOX,
      0,
      0,
      "Geometric center in first FITS axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geoy",
      UI_KEY_GEOY,
      0,
      0,
      "Geometric center in second FITS axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsx",
      UI_KEY_CLUMPSX,
      0,
      0,
      "Flux.wht center of all clumps in obj. (X).",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsy",
      UI_KEY_CLUMPSY,
      0,
      0,
      "Flux.wht center of all clumps in obj. (Y).",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeox",
      UI_KEY_CLUMPSGEOX,
      0,
      0,
      "Geometric center of all clumps in obj. (X).",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeoy",
      UI_KEY_CLUMPSGEOY,
      0,
      0,
      "Geometric center of all clumps in obj. (Y).",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "ra",
      UI_KEY_RA,
      0,
      0,
      "Right ascension of flux weighted center.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "dec",
      UI_KEY_DEC,
      0,
      0,
      "Declination of flux weighted center.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geora",
      UI_KEY_GEORA,
      0,
      0,
      "Right ascension of geometric center.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geodec",
      UI_KEY_GEODEC,
      0,
      0,
      "Declination of geometric center.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsra",
      UI_KEY_CLUMPSRA,
      0,
      0,
      "Right ascension of f.wht center of all clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsdec",
      UI_KEY_CLUMPSDEC,
      0,
      0,
      "Declination of f.wht center of all clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeora",
      UI_KEY_CLUMPSGEORA,
      0,
      0,
      "Right ascension of geo. center of all clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeodec",
      UI_KEY_CLUMPSGEODEC,
      0,
      0,
      "Declination of geometric center of all clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "brightness",
      UI_KEY_BRIGHTNESS,
      0,
      0,
      "Brightness (sum of pixel values).",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpbrightness",
      UI_KEY_CLUMPSBRIGHTNESS,
      0,
      0,
      "Brightness of clumps in an object.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "noriverbrightness",
      UI_KEY_NORIVERBRIGHTNESS,
      0,
      0,
      "Sky (not river) subtracted clump brightness.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "magnitude",
      UI_KEY_MAGNITUDE,
      0,
      0,
      "Total magnitude of objects or clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "magnitudeerr",
      UI_KEY_MAGNITUDEERR,
      0,
      0,
      "Magnitude error of objects or clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsmagnitude",
      UI_KEY_CLUMPSMAGNITUDE,
      0,
      0,
      "Magnitude of all clumps in object.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "upperlimit",
      UI_KEY_UPPERLIMIT,
      0,
      0,
      "Upper-limit value, use other options to config.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "upperlimitmag",
      UI_KEY_UPPERLIMITMAG,
      0,
      0,
      "Upper-limit mag. use other options to config.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "riverave",
      UI_KEY_RIVERAVE,
      0,
      0,
      "Average river value surrounding a clump.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "rivernum",
      UI_KEY_RIVERNUM,
      0,
      0,
      "Number of river pixels around a clump.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "sn",
      UI_KEY_SN,
      0,
      0,
      "Signal to noise ratio of objects or clumps.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "sky",
      UI_KEY_SKY,
      0,
      0,
      "Average Sky value under this clump or object.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "std",
      UI_KEY_STD,
      0,
      0,
      "Average Sky standard deviation.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "semimajor",
      UI_KEY_SEMIMAJOR,
      0,
      0,
      "Flux weighted semi-major axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "semiminor",
      UI_KEY_SEMIMINOR,
      0,
      0,
      "Flux weighted semi-minor axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "axisratio",
      UI_KEY_AXISRATIO,
      0,
      0,
      "Flux weighted axis ratio.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "positionangle",
      UI_KEY_POSITIONANGLE,
      0,
      0,
      "Flux weighted position angle.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geosemimajor",
      UI_KEY_GEOSEMIMAJOR,
      0,
      0,
      "Geometric semi-major axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geosemiminor",
      UI_KEY_GEOSEMIMINOR,
      0,
      0,
      "Geometric semi-minor axis.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geoaxisratio",
      UI_KEY_GEOAXISRATIO,
      0,
      0,
      "Geometric axis ratio.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geopositionangle",
      UI_KEY_GEOPOSITIONANGLE,
      0,
      0,
      "Geometric position angle.",
      ARGS_GROUP_COLUMNS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },


    {0}
  };





/* Define the child argp structure
   -------------------------------

   NOTE: these parts can be left untouched.*/
struct argp
gal_options_common_child = {gal_commonopts_options,
                            gal_options_common_argp_parse,
                            NULL, NULL, NULL, NULL, NULL};

/* Use the child argp structure in list of children (only one for now). */
struct argp_child
children[]=
{
  {&gal_options_common_child, 0, NULL, 0},
  {0, 0, 0, 0}
};

/* Set all the necessary argp parameters. */
struct argp
thisargp = {program_options, parse_opt, args_doc, doc, children, NULL, NULL};
#endif
