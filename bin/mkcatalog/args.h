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
#ifndef ARGS_H
#define ARGS_H






/* Array of acceptable options. */
struct argp_option program_options[] =
  {
    /* Input options. */
    {
      "clumpsfile",
      UI_KEY_CLUMPSFILE,
      "STR",
      0,
      "Dataset containing clump labels.",
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
      "Clump labels extension name or number.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->clumpshdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "valuesfile",
      UI_KEY_VALUESFILE,
      "STR",
      0,
      "Values/brightness dataset.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->valuesfile,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "valueshdu",
      UI_KEY_VALUESHDU,
      "STR",
      0,
      "Name or number of extension containing values.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->valueshdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "insky",
      UI_KEY_INSKY,
      "STR/FLT",
      0,
      "Input Sky value or dataset.",
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
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "subtractsky",
      UI_KEY_SUBTRACTSKY,
      0,
      0,
      "Subtract the Sky dataset from the values.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->subtractsky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "instd",
      UI_KEY_INSTD,
      "STR/FLT",
      0,
      "Sky standard deviation value or dataset.",
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
      "Sky STD extension name or number.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->stdhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "variance",
      UI_KEY_VARIANCE,
      0,
      0,
      "STD input dataset is actually variance.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->variance,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "forcereadstd",
      UI_KEY_FORCEREADSTD,
      0,
      0,
      "Read STD even if no columns need it.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->forcereadstd,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
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
      "sigmaclip",
      UI_KEY_SIGMACLIP,
      "FLT,FLT",
      0,
      "Sigma-clip column multiple and tolerance.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->sigmaclip,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_sigma_clip
    },



    /* Output. */
    {
      "clumpscat",
      UI_KEY_CLUMPSCAT,
      0,
      0,
      "Make a clumps catalog also.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->clumpscat,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "noclumpsort",
      UI_KEY_NOCLUMPSORT,
      0,
      0,
      "Don't sort the clumps catalog by ID.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->noclumpsort,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
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
    {
      "spectrum",
      UI_KEY_SPECTRUM,
      0,
      0,
      "Object spectrum for cube (3D) datasets.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->spectrum,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "inbetweenints",
      UI_KEY_INBETWEENINTS,
      0,
      0,
      "Keep rows (integer ids) with no labels.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->inbetweenints,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Upper limit magnitude configurations. */
    {
      0, 0, 0, 0,
      "Upper limit magnitude settings:",
      UI_GROUP_UPPERLIMIT
    },
    {
      "upmaskfile",
      UI_KEY_UPMASKFILE,
      "STR",
      0,
      "Mask image file name only for upper limit.",
      UI_GROUP_UPPERLIMIT,
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
      UI_GROUP_UPPERLIMIT,
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
      UI_GROUP_UPPERLIMIT,
      &p->upnum,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "uprange",
      UI_KEY_UPRANGE,
      "INT,INT",
      0,
      "Range of random positions (pix) around target.",
      UI_GROUP_UPPERLIMIT,
      &p->uprange,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },
    {
      "envseed",
      UI_KEY_ENVSEED,
      0,
      0,
      "Use GSL_RNG_SEED environment variable for seed.",
      UI_GROUP_UPPERLIMIT,
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
      UI_GROUP_UPPERLIMIT,
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
      UI_GROUP_UPPERLIMIT,
      &p->upnsigma,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checkuplim",
      UI_KEY_CHECKUPLIM,
      "INT[,INT]",
      0,
      "Check random distribution for one label.",
      UI_GROUP_UPPERLIMIT,
      &p->checkuplim,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_check_upperlimit
    },





    /* ID related columns. */
    {
      0, 0, 0, 0,
      "Identifier columns",
      UI_GROUP_COLUMNS_IDS
    },
    {  /* 'ids' is not a unique column, it is a combination of several
          columns. */
      "ids",
      UI_KEY_IDS,
      0,
      0,
      "All IDs of objects and clumps.",
      UI_GROUP_COLUMNS_IDS,
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
      UI_GROUP_COLUMNS_IDS,
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
      UI_GROUP_COLUMNS_IDS,
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
      UI_GROUP_COLUMNS_IDS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },





    /* Position related columns (pixel). */
    {
      0, 0, 0, 0,
      "Positional columns (pixel)",
      UI_GROUP_COLUMNS_POSITION_PIXEL
    },
    {
      "x",
      UI_KEY_X,
      0,
      0,
      "Flux weighted center in first FITS axis.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
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
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "z",
      UI_KEY_Z,
      0,
      0,
      "Flux weighted center in third FITS axis.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
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
      UI_GROUP_COLUMNS_POSITION_PIXEL,
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
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geoz",
      UI_KEY_GEOZ,
      0,
      0,
      "Geometric center in third FITS axis.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "minx",
      UI_KEY_MINX,
      0,
      0,
      "Minimum first FITS axis position.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "maxx",
      UI_KEY_MAXX,
      0,
      0,
      "Maximum first FITS axis position.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "miny",
      UI_KEY_MINY,
      0,
      0,
      "Minimum second FITS axis position.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "maxy",
      UI_KEY_MAXY,
      0,
      0,
      "Maximum second FITS axis position.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "minz",
      UI_KEY_MINZ,
      0,
      0,
      "Minimum third FITS axis position.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "maxz",
      UI_KEY_MAXZ,
      0,
      0,
      "Maximum third FITS axis position.",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
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
      UI_GROUP_COLUMNS_POSITION_PIXEL,
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
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsz",
      UI_KEY_CLUMPSZ,
      0,
      0,
      "Flux.wht center of all clumps in obj. (Z).",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
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
      UI_GROUP_COLUMNS_POSITION_PIXEL,
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
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeoz",
      UI_KEY_CLUMPSGEOZ,
      0,
      0,
      "Geometric center of all clumps in obj. (Z).",
      UI_GROUP_COLUMNS_POSITION_PIXEL,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },





    /* Position related columns (WCS). */
    {
      0, 0, 0, 0,
      "Positional columns (WCS)",
      UI_GROUP_COLUMNS_POSITION_WCS
    },
    {
      "ra",
      UI_KEY_RA,
      0,
      0,
      "Flux weighted center right ascension.",
      UI_GROUP_COLUMNS_POSITION_WCS,
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
      "Flux weighted center declination.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "w1",
      UI_KEY_W1,
      0,
      0,
      "Flux weighted center in first WCS axis.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "w2",
      UI_KEY_W2,
      0,
      0,
      "Flux weighted center in second WCS axis.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "w3",
      UI_KEY_W3,
      0,
      0,
      "Flux weighted center in third WCS axis.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geow1",
      UI_KEY_GEOW1,
      0,
      0,
      "Geometric center in first WCS axis.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geow2",
      UI_KEY_GEOW2,
      0,
      0,
      "Geometric center in second WCS axis.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geow3",
      UI_KEY_GEOW2,
      0,
      0,
      "Geometric center in third WCS axis.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsw1",
      UI_KEY_CLUMPSW1,
      0,
      0,
      "Flux.wht center of all clumps in 1st WCS.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsw2",
      UI_KEY_CLUMPSW2,
      0,
      0,
      "Flux.wht center of all clumps in 2nd WCS.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsw3",
      UI_KEY_CLUMPSW3,
      0,
      0,
      "Flux.wht center of all clumps in 3rd WCS.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeow1",
      UI_KEY_CLUMPSGEOW1,
      0,
      0,
      "Geometric center of all clumps in 1st WCS.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeow2",
      UI_KEY_CLUMPSGEOW2,
      0,
      0,
      "Geometric center of all clumps in 2nd WCS.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "clumpsgeow3",
      UI_KEY_CLUMPSGEOW3,
      0,
      0,
      "Geometric center of all clumps in 3rd WCS.",
      UI_GROUP_COLUMNS_POSITION_WCS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },





    /* Brightness/pixel-value related columns. */
    {
      0, 0, 0, 0,
      "Brightness/magnitude related columns",
      UI_GROUP_COLUMNS_BRIGHTNESS
    },
    {
      "brightness",
      UI_KEY_BRIGHTNESS,
      0,
      0,
      "Brightness (sum of pixel values).",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "brightnesserr",
      UI_KEY_BRIGHTNESSERR,
      0,
      0,
      "Error (1-sigma) in measuring brightness.",
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "brightnessnoriver",
      UI_KEY_BRIGHTNESSNORIVER,
      0,
      0,
      "Sky (not river) subtracted clump brightness.",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "mean",
      UI_KEY_MEAN,
      0,
      0,
      "Mean of values in object/clump.",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "median",
      UI_KEY_MEDIAN,
      0,
      0,
      "Median of values in object/clump.",
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "upperlimitonesigma",
      UI_KEY_UPPERLIMITONESIGMA,
      0,
      0,
      "Upper-limit one sigma value.",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "upperlimitsigma",
      UI_KEY_UPPERLIMITSIGMA,
      0,
      0,
      "Place in random distribution (sigma multiple).",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "upperlimitquantile",
      UI_KEY_UPPERLIMITQUANTILE,
      0,
      0,
      "Quantile in random distribution (max 1).",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "upperlimitskew",
      UI_KEY_UPPERLIMITSKEW,
      0,
      0,
      "(Mean-Median)/STD of random distribution.",
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      "Sky value (per pixel).",
      UI_GROUP_COLUMNS_BRIGHTNESS,
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
      "Sky standard deviation (per pixel).",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "sigclip-number",
      UI_KEY_SIGCLIPNUMBER,
      0,
      0,
      "Number of pixels in Sigma-clipped measurement.",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "sigclip-median",
      UI_KEY_SIGCLIPMEDIAN,
      0,
      0,
      "Median after Sigma-clipping",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "sigclip-mean",
      UI_KEY_SIGCLIPMEAN,
      0,
      0,
      "Mean after Sigma-clipping",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "sigclip-std",
      UI_KEY_SIGCLIPSTD,
      0,
      0,
      "Standard deviation after Sigma-clipping",
      UI_GROUP_COLUMNS_BRIGHTNESS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },




    /* Morphology/shape related columns. */
    {
      0, 0, 0, 0,
      "Morphology/shape related columns",
      UI_GROUP_COLUMNS_MORPHOLOGY
    },
    {
      "numclumps",
      UI_KEY_NUMCLUMPS,
      0,
      0,
      "Number of clumps in this object.",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      "Number of non-blank valued pixels.",
      UI_GROUP_COLUMNS_MORPHOLOGY,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "areaxy",
      UI_KEY_AREAXY,
      0,
      0,
      "Projected area in first two dimensions.",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      "Non-blank area covered by clumps.",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      "Area used for value weighted positions.",
      UI_GROUP_COLUMNS_MORPHOLOGY,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geoarea",
      UI_KEY_GEOAREA,
      0,
      0,
      "Area labled region (irrespective of value).",
      UI_GROUP_COLUMNS_MORPHOLOGY,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_column_codes_ll
    },
    {
      "geoareaxy",
      UI_KEY_GEOAREAXY,
      0,
      0,
      "Projected geoarea in first two dimensions.",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      "RMS along major axis (in pixels).",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      "RMS along minor axis (in pixels).",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      "RMS along major axis (ignoring value).",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      "RMS along minor axis (ignoring value).",
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
      UI_GROUP_COLUMNS_MORPHOLOGY,
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
