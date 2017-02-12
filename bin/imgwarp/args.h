/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (Gnuastro) package.

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

    /* Input. */
    {
      "hstartwcs",
      ARGS_OPTION_KEY_HSTARTWCS,
      "INT",
      0,
      "Header keyword number to start reading WCS.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->hstartwcs,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "hendwcs",
      ARGS_OPTION_KEY_HENDWCS,
      "INT",
      0,
      "Header keyword number to end reading WCS.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->hendwcs,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Output. */
    {
      "keepinputwcs",
      ARGS_OPTION_KEY_KEEPINPUTWCS,
      0,
      0,
      "Do not apply warp to input's WCS",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->keepinputwcs,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "maxblankfrac",
      ARGS_OPTION_KEY_MAXBLANKFRAC,
      "FLT",
      0,
      "Maximum fraction of area covered by blank.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->maxblankfrac,
      GAL_DATA_TYPE_FLOAT,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "type",
      ARGS_OPTION_KEY_TYPE,
      "STR",
      0,
      "uchar, short, long, longlong, float, double."
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->typestr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },


    {
      0, 0, 0, 0,
      "Warps:",
      ARGS_GROUP_WARPS
    },
    {
      "align",
      ARGS_OPTION_KEY_ALIGN,
      0,
      0,
      "Align the image and celestial axes."
      ARGS_GROUP_WARPS,
      &p->align,
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
      "uchar, short, long, longlong, float, double."
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->typestr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
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
