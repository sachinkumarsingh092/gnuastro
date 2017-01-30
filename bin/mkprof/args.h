/*********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef ARGS_H
#define ARGS_H





/* Option definition array. */
struct argp_option program_options[] =
  {
    {
      "background",
      ARGS_OPTION_KEY_BACKGROUND,
      "STR",
      0,
      "A background image to make the profiles on.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->backname,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "backhdu",
      ARGS_OPTION_KEY_BACKHDU,
      "INT/STR",
      0,
      "HDU of background image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->backhdu,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      "naxis1",
      ARGS_OPTION_KEY_NAXIS1,
      "INT",
      0,
      "Number of pixels along first FITS axis.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->naxes[0],
      GAL_DATA_TYPE_LONG,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "naxis2",
      ARGS_OPTION_KEY_NAXIS2,
      "INT",
      0,
      "Number of pixels along second FITS axis.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->naxes[1],
      GAL_DATA_TYPE_LONG,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "oversample",
      ARGS_OPTION_KEY_OVERSAMPLE,
      "INT",
      0,
      "Scale of oversampling (>0 and odd).",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->oversample,
      GAL_DATA_TYPE_UCHAR,
      GAL_OPTIONS_RANGE_GT_0_ODD,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "psfinimg",
      ARGS_OPTION_KEY_PSFINIMG,
      0,
      0,
      "PSF profiles made with all in output image.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->psfinimg,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "individual",
      ARGS_OPTION_KEY_INDIVIDUAL,
      0,
      0,
      "Build all profiles separately.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->individual,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "nomerged",
      ARGS_OPTION_KEY_NOMERGED,
      0,
      0,
      "Do not create a merged image of all profiles.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->nomerged,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "type",
      ARGS_OPTION_KEY_TYPE,
      "STR",
      0,
      "uchar, short, long, longlong, float, double.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->typestr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Profiles:",
      ARGS_GROUP_PROFILES
    },
    {
      "numrandom",
      ARGS_OPTION_KEY_NUMRANDOM,
      "INT",
      0,
      "No. of random points in Monte Carlo integration.",
      ARGS_GROUP_PROFILES,
      &p->numrandom,
      GAL_DATA_TYPE_ULONG,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tolerance",
      ARGS_OPTION_KEY_TOLERANCE,
      "FLT",
      0,
      "Tolerance to switch to less accurate method.",
      ARGS_GROUP_PROFILES,
      &p->tolerance,
      GAL_DATA_TYPE_FLOAT,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tunitinp",
      ARGS_OPTION_KEY_TUNITINP,
      0,
      0,
      "Truncation is in units of pixels, not radius.",
      ARGS_GROUP_PROFILES,
      &p->tunitinp,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "xshift",
      ARGS_OPTION_KEY_XSHIFT,
      "FLT",
      0,
      "Shift profile centers and enlarge image, X axis.",
      ARGS_GROUP_PROFILES,
      &p->shift[0],
      GAL_DATA_TYPE_LONG,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "yshift",
      ARGS_OPTION_KEY_YSHIFT,
      "FLT",
      0,
      "Shift profile centers and enlarge image, Y axis.",
      ARGS_GROUP_PROFILES,
      &p->shift[1],
      GAL_DATA_TYPE_LONG,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "prepforconv",
      ARGS_OPTION_KEY_PREPFORCONV,
      0,
      0,
      "Shift and expand based on first catalog PSF.",
      ARGS_GROUP_PROFILES,
      &p->prepforconv,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "zeropoint",
      ARGS_OPTION_KEY_ZEROPOINT,
      "FLT",
      0,
      "Magnitude zero point.",
      ARGS_GROUP_PROFILES,
      &p->zeropoint,
      GAL_DATA_TYPE_FLOAT,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "circumwidth",
      ARGS_OPTION_KEY_CIRCUMWIDTH,
      "FLT",
      0,
      "Width of circumference (inward) profiles",
      ARGS_GROUP_PROFILES,
      &p->circumwidth,
      GAL_DATA_TYPE_FLOAT,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "replace",
      ARGS_OPTION_KEY_REPLACE,
      0,
      0,
      "Replace overlapping profile pixels, don't add.",
      ARGS_GROUP_PROFILES,
      &p->replace,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "magatpeak",
      ARGS_OPTION_KEY_MAGATPEAK,
      0,
      0,
      "Magnitude is for peak pixel, not full profile.",
      ARGS_GROUP_PROFILES,
      &p->magatpeak,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "envseed",
      ARGS_OPTION_KEY_ENVSEED,
      0,
      0,
      "Use GSL_RNG_SEED environment variable for seed.",
      ARGS_GROUP_PROFILES,
      &p->envseed,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Columns, by info (see `--searchin'), or number (starting from 1):",
      ARGS_GROUP_CATALOG
    },
    {
      "xcol",
      ARGS_OPTION_KEY_XCOL,
      "STR/INT",
      0,
      "Center along first FITS axis (horizontal).",
      ARGS_GROUP_CATALOG,
      &p->xcol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ycol",
      ARGS_OPTION_KEY_YCOL,
      "STR/INT",
      0,
      "Center along second FITS axis (vertical).",
      ARGS_GROUP_CATALOG,
      &p->ycol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "racol",
      ARGS_OPTION_KEY_RACOL,
      "STR/INT",
      0,
      "Center right ascension.",
      ARGS_GROUP_CATALOG,
      &p->racol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "deccol",
      ARGS_OPTION_KEY_DECCOL,
      "STR/INT",
      0,
      "Center declination.",
      ARGS_GROUP_CATALOG,
      &p->deccol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "fcol",
      ARGS_OPTION_KEY_FCOL,
      "STR/INT",
      0,
      "Sersic (1), Moffat (2), Gaussian (3), Point (4), "
      "Flat (5), Circumference (6).",
      ARGS_GROUP_CATALOG,
      &p->fcol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "rcol",
      ARGS_OPTION_KEY_RCOL,
      "STR/INT",
      0,
      "Effective radius or FWHM in pixels.",
      ARGS_GROUP_CATALOG,
      &p->rcol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ncol",
      ARGS_OPTION_KEY_NCOL,
      "STR/INT",
      0,
      "Sersic index or Moffat beta.",
      ARGS_GROUP_CATALOG,
      &p->ncol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "pcol",
      ARGS_OPTION_KEY_PCOL,
      "STR/INT",
      0,
      "Position angle.",
      ARGS_GROUP_CATALOG,
      &p->pcol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "qcol",
      ARGS_OPTION_KEY_QCOL,
      "STR/INT",
      0,
      "Axis ratio.",
      ARGS_GROUP_CATALOG,
      &p->qcol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mcol",
      ARGS_OPTION_KEY_MCOL,
      "STR/INT",
      0,
      "Magnitude.",
      ARGS_GROUP_CATALOG,
      &p->mcol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tcol",
      ARGS_OPTION_KEY_TCOL,
      "STR/INT",
      0,
      "Truncation in units of --rcol, unless --tunitinp.",
      ARGS_GROUP_CATALOG,
      &p->tcol,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mforflatpix",
      ARGS_OPTION_KEY_MFORFLATPIX,
      0,
      0,
      "mcol is flat pixel value (when fcol is 4 or 5)",
      ARGS_GROUP_CATALOG,
      &p->mforflatpix,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "WCS parameters:",
      ARGS_GROUP_WCS
    },
    {
      "crpix1",
      ARGS_OPTION_KEY_CRPIX1,
      "FLT",
      0,
      "Pixel coordinate of reference point (axis 1).",
      ARGS_GROUP_WCS,
      &p->crpix[0],
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "crpix2",
      ARGS_OPTION_KEY_CRPIX2,
      "FLT",
      0,
      "Pixel coordinate of reference point (axis 2).",
      ARGS_GROUP_WCS,
      &p->crpix[1],
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "crval1",
      ARGS_OPTION_KEY_CRVAL1,
      "FLT",
      0,
      "Right ascension at reference point (degrees).",
      ARGS_GROUP_WCS,
      &p->crval[0],
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "crval2",
      ARGS_OPTION_KEY_CRVAL2,
      "FLT",
      0,
      "Declination at reference point (degrees).",
      ARGS_GROUP_WCS,
      &p->crval[1],
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "resolution",
      ARGS_OPTION_KEY_RESOLUTION,
      "FLT",
      0,
      "Resolution of image (arcseconds/pixel).",
      ARGS_GROUP_WCS,
      &p->resolution,
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {0}
  };





/* Define the child argp structure. */
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
