/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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
    /* Output */
    {
      "quality",
      ARGS_OPTION_KEY_QUALITY,
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
      ARGS_OPTION_KEY_WIDTHINCM,
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
      ARGS_OPTION_KEY_BORDERWIDTH,
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
      ARGS_OPTION_KEY_HEX,
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
      0, 0, 0, 0,
      "Flux:",
      ARGS_GROUP_FLUX
    },
    {
      "fluxlow",
      ARGS_OPTION_KEY_FLUXLOW,
      "FLT",
      0,
      "Lower flux truncation value.",
      ARGS_GROUP_FLUX,
      &p->fluxlowstr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "fluxhigh",
      ARGS_OPTION_KEY_FLUXHIGH,
      "FLT",
      0,
      "Higher flux truncation value.",
      ARGS_GROUP_FLUX,
      &p->fluxhighstr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "maxbyte",
      ARGS_OPTION_KEY_MAXBYTE,
      "INT",
      0,
      "Maximum byte value for all color channels.",
      ARGS_GROUP_FLUX,
      &p->maxbyte,
      GAL_TYPE_UINT8,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "flminbyte",
      ARGS_OPTION_KEY_FLMINBYTE,
      0,
      0,
      "Set value of fluxlow as the minimum byte value.",
      ARGS_GROUP_FLUX,
      &p->flminbyte,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "fhmaxbyte",
      ARGS_OPTION_KEY_FHMAXBYTE,
      0,
      0,
      "Set value of fluxhigh as the maximum byte value.",
      ARGS_GROUP_FLUX,
      &p->fhmaxbyte,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "change",
      ARGS_OPTION_KEY_CHANGE,
      "STR",
      0,
      "Change pixel values `from_1:to_1,from_2:to_2`.",
      ARGS_GROUP_FLUX,
      &p->changestr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "changeaftertrunc",
      ARGS_OPTION_KEY_CHANGEAFTERTRUNC,
      0,
      0,
      "First truncate then change pixel values.",
      ARGS_GROUP_FLUX,
      &p->changeaftertrunc,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "invert",
      ARGS_OPTION_KEY_INVERT,
      0,
      0,
      "Invert the values in JPEG and EPS/PDF.",
      ARGS_GROUP_FLUX,
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
