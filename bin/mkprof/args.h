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





/* Include necessary headers. */
#include <commonopts.h>
#include <authors-cite.h>
#include <fixedstringmacros.h>





/* Definition parameters for the argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "[BackgroundImage] Catalog";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will create a FITS image "
  "containing any number of mock astronomical profiles based on an input "
  "catalog. All the profiles will be built from the center outwards. First "
  "by Monte Carlo integration, then using the central pixel position. The "
  "tolerance level specifies when to switch to a latter.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Keys for each option.

   Available letters (-V which is used by GNU is also removed):

   a d f g j k l u v
   A E G H I J L M O Q U W Z     */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_BACKHDU_KEY         = 'B',
  ARGS_OPTION_NAXIS1_KEY          = 'x',
  ARGS_OPTION_NAXIS2_KEY          = 'y',
  ARGS_OPTION_INPUTASCANVAS_KEY   = 'C',
  ARGS_OPTION_OVERSAMPLE_KEY      = 's',
  ARGS_OPTION_INDIVIDUAL_KEY      = 'i',
  ARGS_OPTION_NOMERGED_KEY        = 'm',
  ARGS_OPTION_TYPE_KEY            = 'T',
  ARGS_OPTION_NUMRANDOM_KEY       = 'r',
  ARGS_OPTION_TOLERANCE_KEY       = 't',
  ARGS_OPTION_TUNITINP_KEY        = 'p',
  ARGS_OPTION_XSHIFT_KEY          = 'X',
  ARGS_OPTION_YSHIFT_KEY          = 'Y',
  ARGS_OPTION_PREPFORCONV_KEY     = 'c',
  ARGS_OPTION_ZEROPOINT_KEY       = 'z',
  ARGS_OPTION_CIRCUMWIDTH_KEY     = 'w',
  ARGS_OPTION_REPLACE_KEY         = 'R',
  ARGS_OPTION_ENVSEED_KEY         = 'e',
  ARGS_OPTION_MFORFLATPIX_KEY     = 'F',

  /* Only with long version. */
  ARGS_OPTION_PSFINIMG_KEY        = 1000,
  ARGS_OPTION_MAGATPEAK_KEY,
  ARGS_OPTION_XCOL_KEY,
  ARGS_OPTION_YCOL_KEY,
  ARGS_OPTION_RACOL_KEY,
  ARGS_OPTION_DECCOL_KEY,
  ARGS_OPTION_FCOL_KEY,
  ARGS_OPTION_RCOL_KEY,
  ARGS_OPTION_NCOL_KEY,
  ARGS_OPTION_PCOL_KEY,
  ARGS_OPTION_QCOL_KEY,
  ARGS_OPTION_MCOL_KEY,
  ARGS_OPTION_TCOL_KEY,
  ARGS_OPTION_CRPIX1_KEY,
  ARGS_OPTION_CRPIX2_KEY,
  ARGS_OPTION_CRVAL1_KEY,
  ARGS_OPTION_CRVAL2_KEY,
  ARGS_OPTION_RESOLUTION_KEY,
};





/* Option definition array. */
static struct argp_option options[] =
  {
    {
      0, 0, 0, 0,
      "Operating modes:",
      -1
    },





    {
      0, 0, 0, 0,
      "Input:",
      1
    },
    {
      "backhdu",
      ARGS_OPTION_BACKHDU_KEY,
      "INT/STR",
      0,
      "HDU of background image.",
      1,
      NULL, GAL_DATA_TYPE_STRING
    },



    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "naxis1",
      ARGS_OPTION_NAXIS1_KEY,
      "INT",
      0,
      "Number of pixels along first FITS axis.",
      2,
      NULL, GAL_DATA_TYPE_ULONG
    },
    {
      "naxis2",
      ARGS_OPTION_NAXIS2_KEY,
      "INT",
      0,
      "Number of pixels along second FITS axis.",
      2,
      NULL, GAL_DATA_TYPE_ULONG
    },
    {
      "inputascanvas",
      ARGS_OPTION_INPUTASCANVAS_KEY,
      0,
      0,
      "Use input image for output size and WCS.",
      2,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "oversample",
      ARGS_OPTION_OVERSAMPLE_KEY,
      "INT",
      0,
      "Scale of oversampling.",
      2,
      NULL, GAL_DATA_TYPE_ULONG
    },
    {
      "psfinimg",
      ARGS_OPTION_PSFINIMG_KEY,
      0,
      0,
      "PSF profiles made with all in output image.",
      2,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "individual",
      ARGS_OPTION_INDIVIDUAL_KEY,
      0,
      0,
      "Build all profiles separately.",
      2,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "nomerged",
      ARGS_OPTION_NOMERGED_KEY,
      0,
      0,
      "Do not create a merged image of all profiles.",
      2,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "type",
      ARGS_OPTION_TYPE_KEY,
      "STR",
      0,
      "uchar, short, long, longlong, float, double.",
      2,
      NULL, GAL_DATA_TYPE_STRING
    },





    {
      0, 0, 0, 0,
      "Profiles:",
      3
    },
    {
      "numrandom",
      ARGS_OPTION_NUMRANDOM_KEY,
      "INT",
      0,
      "No. of random points in Monte Carlo integration.",
      3,
      NULL, GAL_DATA_TYPE_ULONG
    },
    {
      "tolerance",
      ARGS_OPTION_TOLERANCE_KEY,
      "FLT",
      0,
      "Tolerance to switch to less accurate method.",
      3,
      NULL, GAL_DATA_TYPE_FLOAT
    },
    {
      "tunitinp",
      ARGS_OPTION_TUNITINP_KEY,
      0,
      0,
      "Truncation is in units of pixels, not radius.",
      3,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "xshift",
      ARGS_OPTION_XSHIFT_KEY,
      "FLT",
      0,
      "Shift profile centers and enlarge image, X axis.",
      3,
      NULL, GAL_DATA_TYPE_FLOAT
    },
    {
      "yshift",
      ARGS_OPTION_YSHIFT_KEY,
      "FLT",
      0,
      "Shift profile centers and enlarge image, Y axis.",
      3,
      NULL, GAL_DATA_TYPE_FLOAT
    },
    {
      "prepforconv",
      ARGS_OPTION_PREPFORCONV_KEY,
      0,
      0,
      "Shift and expand based on first catalog PSF.",
      3,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "zeropoint",
      ARGS_OPTION_ZEROPOINT_KEY,
      "FLT",
      0,
      "Magnitude zero point.",
      3,
      NULL, GAL_DATA_TYPE_FLOAT
    },
    {
      "circumwidth",
      ARGS_OPTION_CIRCUMWIDTH_KEY,
      "FLT",
      0,
      "Width of circumference (inward) profiles",
      3,
      NULL, GAL_DATA_TYPE_FLOAT
    },
    {
      "replace",
      ARGS_OPTION_REPLACE_KEY,
      0,
      0,
      "Replace overlapping profile pixels, don't add.",
      3,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "magatpeak",
      ARGS_OPTION_MAGATPEAK_KEY,
      0,
      0,
      "Magnitude is for peak pixel, not full profile.",
      3,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "envseed",
      ARGS_OPTION_ENVSEED_KEY,
      0,
      0,
      "Use GSL_RNG_SEED environment variable for seed.",
      3,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },





    {
      0, 0, 0, 0,
      "Catalog (column number, starting from zero):",
      4
    },
    {
      "xcol",
      ARGS_OPTION_XCOL_KEY,
      "INT/STR",
      0,
      "Center along first FITS axis (horizontal).",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "ycol",
      ARGS_OPTION_YCOL_KEY,
      "INT/STR",
      0,
      "Center along second FITS axis (vertical).",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "racol",
      ARGS_OPTION_RACOL_KEY,
      "INT/STR",
      0,
      "Center right ascension.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "deccol",
      ARGS_OPTION_DECCOL_KEY,
      "INT/STR",
      0,
      "Center declination.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "fcol",
      ARGS_OPTION_FCOL_KEY,
      "INT/STR",
      0,
      "Sersic (0), Moffat (1), Gaussian (2), Point (3),\n"
      "Flat (4), Circumference (5).",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "rcol",
      ARGS_OPTION_RCOL_KEY,
      "INT/STR",
      0,
      "Effective radius or FWHM in pixels.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "ncol",
      ARGS_OPTION_NCOL_KEY,
      "INT/STR",
      0,
      "Sersic index or Moffat beta.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "pcol",
      ARGS_OPTION_PCOL_KEY,
      "INT/STR",
      0,
      "Position angle.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "qcol",
      ARGS_OPTION_QCOL_KEY,
      "INT/STR",
      0,
      "Axis ratio.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "mcol",
      ARGS_OPTION_MCOL_KEY,
      "INT/STR",
      0,
      "Magnitude.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "tcol",
      ARGS_OPTION_TCOL_KEY,
      "INT/STR",
      0,
      "Truncation in units of --rcol, unless --tunitinp.",
      4,
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "mforflatpix",
      ARGS_OPTION_MFORFLATPIX_KEY,
      0,
      0,
      "mcol is flat pixel value (when fcol is 4 or 5)",
      4,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },





    {
      0, 0, 0, 0,
      "WCS parameters:",
      5
    },
    {
      "crpix1",
      ARGS_OPTION_CRPIX1_KEY,
      "FLT",
      0,
      "Pixel coordinate of reference point (axis 1).",
      5,
      NULL, GAL_DATA_TYPE_DOUBLE
    },
    {
      "crpix2",
      ARGS_OPTION_CRPIX2_KEY,
      "FLT",
      0,
      "Pixel coordinate of reference point (axis 2).",
      5,
      NULL, GAL_DATA_TYPE_DOUBLE
    },
    {
      "crval1",
      ARGS_OPTION_CRVAL1_KEY,
      "FLT",
      0,
      "Right ascension at reference point (degrees).",
      5,
      NULL, GAL_DATA_TYPE_DOUBLE
    },
    {
      "crval2",
      ARGS_OPTION_CRVAL2_KEY,
      "FLT",
      0,
      "Declination at reference point (degrees).",
      5,
      NULL, GAL_DATA_TYPE_DOUBLE
    },
    {
      "resolution",
      ARGS_OPTION_RESOLUTION_KEY,
      "FLT",
      0,
      "Resolution of image (arcseconds/pixel).",
      5,
      NULL, GAL_DATA_TYPE_DOUBLE
    },


    {0}
  };





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct mkprofparams *p = state->input;

  /* Pass `gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then `arg`
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign (`=`). For short "
               "options, `=` should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {

    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      gal_linkedlist_add_to_stll(&p->allargs, arg, 0);
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, options, &p->cp);
    }

  return 0;
}





/* Define the child argp structure. */
static struct argp
gal_options_common_child = {gal_commonopts_options,
                            gal_options_common_argp_parse,
                            NULL, NULL, NULL, NULL, NULL};

/* Use the child argp structure in list of children (only one for now). */
static struct argp_child
children[]=
{
  {&gal_options_common_child, 0, NULL, 0},
  {0, 0, 0, 0}
};

/* Set all the necessary argp parameters. */
static struct argp
thisargp = {options, parse_opt, args_doc, doc, children, NULL, NULL};
#endif
