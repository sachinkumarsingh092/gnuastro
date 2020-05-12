/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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
    /* Output */
    {
      "quality",
      UI_KEY_QUALITY,
      "INT",
      0,
      "Quality of output JPEG image (1 to 100).",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->quality,
      GAL_TYPE_UINT8,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "widthincm",
      UI_KEY_WIDTHINCM,
      "FLT",
      0,
      "Width in units of centimeters.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->widthincm,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "borderwidth",
      UI_KEY_BORDERWIDTH,
      "INT",
      0,
      "EPS/PDF border width in units of 1/72 inch.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->borderwidth,
      GAL_TYPE_UINT32,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "hex",
      UI_KEY_HEX,
      0,
      0,
      "Hexadecimal encoding in EPS. Default: ASCII85.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->hex,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "colormap",
      UI_KEY_COLORMAP,
      "STR[,FLT]",
      0,
      "Color map when only a single channel is given.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->colormap,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings
    },
    {
      "rgbtohsv",
      UI_KEY_RGBTOHSV,
      0,
      0,
      "Convert RGB input into HSV (in FITS output)",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->rgbtohsv,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },






    {
      0, 0, 0, 0,
      "Flux:",
      UI_GROUP_FLUX
    },
    {
      "fluxlow",
      UI_KEY_FLUXLOW,
      "FLT",
      0,
      "Lower flux truncation value.",
      UI_GROUP_FLUX,
      &p->fluxlowstr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "fluxhigh",
      UI_KEY_FLUXHIGH,
      "FLT",
      0,
      "Higher flux truncation value.",
      UI_GROUP_FLUX,
      &p->fluxhighstr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "maxbyte",
      UI_KEY_MAXBYTE,
      "INT",
      0,
      "Maximum byte value for all color channels.",
      UI_GROUP_FLUX,
      &p->maxbyte,
      GAL_TYPE_UINT8,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "forcemin",
      UI_KEY_FORCEMIN,
      0,
      0,
      "Force --fluxmin, even when smaller than minimum.",
      UI_GROUP_FLUX,
      &p->forcemin,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "forcemax",
      UI_KEY_FORCEMAX,
      0,
      0,
      "Force --fluxmax, even when larger than maximum.",
      UI_GROUP_FLUX,
      &p->forcemax,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "change",
      UI_KEY_CHANGE,
      "STR",
      0,
      "Change pixel values 'from_1:to_1,from_2:to_2'.",
      UI_GROUP_FLUX,
      &p->changestr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "changeaftertrunc",
      UI_KEY_CHANGEAFTERTRUNC,
      0,
      0,
      "First truncate then change pixel values.",
      UI_GROUP_FLUX,
      &p->changeaftertrunc,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "invert",
      UI_KEY_INVERT,
      0,
      0,
      "Invert the values in JPEG and EPS/PDF.",
      UI_GROUP_FLUX,
      &p->invert,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
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
