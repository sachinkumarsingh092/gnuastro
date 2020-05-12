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
#ifndef ARGS_H
#define ARGS_H





/* Option definition array. */
struct argp_option program_options[] =
  {
    {
      "background",
      UI_KEY_BACKGROUND,
      "STR",
      0,
      "A background image to make the profiles on.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->backname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "backhdu",
      UI_KEY_BACKHDU,
      "INT/STR",
      0,
      "HDU of background image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->backhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "clearcanvas",
      UI_KEY_CLEARCANVAS,
      0,
      0,
      "All pixels in background image read as zero.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->clearcanvas,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "kernel",
      UI_KEY_KERNEL,
      "STR",
      0,
      "Parameters to only build one kernel.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->kernel,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_parse_kernel
    },





    {
      "mergedsize",
      UI_KEY_MERGEDSIZE,
      "INT[,INT,...]",
      0,
      "Merged image size along each dimension.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->dsize,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },
    {
      "oversample",
      UI_KEY_OVERSAMPLE,
      "INT",
      0,
      "Scale of oversampling (>0 and odd).",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->oversample,
      GAL_TYPE_UINT8,
      GAL_OPTIONS_RANGE_GT_0_ODD,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "psfinimg",
      UI_KEY_PSFINIMG,
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
      UI_KEY_INDIVIDUAL,
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
      UI_KEY_NOMERGED,
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
      0, 0, 0, 0,
      "Profiles:",
      UI_GROUP_PROFILES
    },
    {
      "mode",
      UI_KEY_MODE,
      "STR",
      0,
      "Mode of '--ccol': 'img' or 'wcs'.",
      UI_GROUP_PROFILES,
      &p->mode,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_parse_coordinate_mode
    },
    {
      "numrandom",
      UI_KEY_NUMRANDOM,
      "INT",
      0,
      "No. of random points in Monte Carlo integration.",
      UI_GROUP_PROFILES,
      &p->numrandom,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "envseed",
      UI_KEY_ENVSEED,
      0,
      0,
      "Use GSL_RNG_SEED environment variable for seed.",
      UI_GROUP_PROFILES,
      &p->envseed,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tolerance",
      UI_KEY_TOLERANCE,
      "FLT",
      0,
      "Tolerance to switch to less accurate method.",
      UI_GROUP_PROFILES,
      &p->tolerance,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tunitinp",
      UI_KEY_TUNITINP,
      0,
      0,
      "Truncation is in units of pixels, not radius.",
      UI_GROUP_PROFILES,
      &p->tunitinp,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mforflatpix",
      UI_KEY_MFORFLATPIX,
      0,
      0,
      "mcol is flat pixel value (when fcol is 5 or 6)",
      UI_GROUP_PROFILES,
      &p->mforflatpix,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mcolisbrightness",
      UI_KEY_MCOLISBRIGHTNESS,
      0,
      0,
      "mcol is total brightness, not magnitude.",
      UI_GROUP_PROFILES,
      &p->mcolisbrightness,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "shift",
      UI_KEY_SHIFT,
      "INT[, ...]",
      0,
      "Shift profile centers in output image.",
      UI_GROUP_PROFILES,
      &p->shift,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },
    {
      "prepforconv",
      UI_KEY_PREPFORCONV,
      0,
      0,
      "Shift and expand based on first catalog PSF.",
      UI_GROUP_PROFILES,
      &p->prepforconv,
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
      "Magnitude zero point.",
      UI_GROUP_PROFILES,
      &p->zeropoint,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "magatpeak",
      UI_KEY_MAGATPEAK,
      0,
      0,
      "Magnitude is for peak pixel, not full profile.",
      UI_GROUP_PROFILES,
      &p->magatpeak,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "circumwidth",
      UI_KEY_CIRCUMWIDTH,
      "FLT",
      0,
      "Width of circumference (inward) profiles",
      UI_GROUP_PROFILES,
      &p->circumwidth,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "replace",
      UI_KEY_REPLACE,
      0,
      0,
      "Replace overlapping profile pixels, don't add.",
      UI_GROUP_PROFILES,
      &p->replace,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Columns, by info (see '--searchin'), or number (starting from 1):",
      UI_GROUP_CATALOG
    },
    {
      "ccol",
      UI_KEY_CCOL,
      "STR/INT",
      0,
      "Coordinate columns (one call for each dimension).",
      UI_GROUP_CATALOG,
      &p->ccol,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "fcol",
      UI_KEY_FCOL,
      "STR/INT",
      0,
      "sersic (1), moffat (2), gaussian (3), point (4), "
      "flat (5), circumference (6), distance (7).",
      UI_GROUP_CATALOG,
      &p->fcol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "rcol",
      UI_KEY_RCOL,
      "STR/INT",
      0,
      "Effective radius or FWHM in pixels.",
      UI_GROUP_CATALOG,
      &p->rcol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ncol",
      UI_KEY_NCOL,
      "STR/INT",
      0,
      "Sersic index or Moffat beta.",
      UI_GROUP_CATALOG,
      &p->ncol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "pcol",
      UI_KEY_PCOL,
      "STR/INT",
      0,
      "Position angle (First X-Z-X Euler angle in 3D).",
      UI_GROUP_CATALOG,
      &p->pcol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "p2col",
      UI_KEY_P2COL,
      "STR/INT",
      0,
      "Second Euler angle (X-Z-X order).",
      UI_GROUP_CATALOG,
      &p->p2col,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "p3col",
      UI_KEY_P2COL,
      "STR/INT",
      0,
      "Third Euler angle (X-Z-X order).",
      UI_GROUP_CATALOG,
      &p->p3col,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "qcol",
      UI_KEY_QCOL,
      "STR/INT",
      0,
      "Axis ratio (major/dim2 in 3D).",
      UI_GROUP_CATALOG,
      &p->qcol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "q2col",
      UI_KEY_Q2COL,
      "STR/INT",
      0,
      "Axis ratio (major/dim3 in 3D).",
      UI_GROUP_CATALOG,
      &p->q2col,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mcol",
      UI_KEY_MCOL,
      "STR/INT",
      0,
      "Magnitude.",
      UI_GROUP_CATALOG,
      &p->mcol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tcol",
      UI_KEY_TCOL,
      "STR/INT",
      0,
      "Truncation in units of --rcol, unless --tunitinp.",
      UI_GROUP_CATALOG,
      &p->tcol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "WCS parameters:",
      UI_GROUP_WCS
    },
    {
      "crpix",
      UI_KEY_CRPIX,
      "FLT[, ...]",
      0,
      "Pixel coordinates of reference point.",
      UI_GROUP_WCS,
      &p->crpix,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },
    {
      "crval",
      UI_KEY_CRVAL,
      "FLT[, ...]",
      0,
      "WCS coordinates of reference point.",
      UI_GROUP_WCS,
      &p->crval,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },
    {
      "cdelt",
      UI_KEY_CDELT,
      "FLT[, ...]",
      0,
      "Resolution in each dimension.",
      UI_GROUP_WCS,
      &p->cdelt,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },
    {
      "pc",
      UI_KEY_PC,
      "FLT[, ...]",
      0,
      "WCS rotation matrix (all elements).",
      UI_GROUP_WCS,
      &p->pc,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },
    {
      "cunit",
      UI_KEY_CUNIT,
      "STR[, ... ]",
      0,
      "Units of the WCS coordinates (e.g., 'deg').",
      UI_GROUP_WCS,
      &p->cunit,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings
    },
    {
      "ctype",
      UI_KEY_CTYPE,
      "STR[, ... ]",
      0,
      "One of FITS standard WCS types.",
      UI_GROUP_WCS,
      &p->ctype,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings
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
