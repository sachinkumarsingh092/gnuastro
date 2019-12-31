/*********************************************************************
Warp - Warp images using projective mapping.
Warp is part of GNU Astronomy Utilities (Gnuastro) package.

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

    /* Input. */
    {
      "hstartwcs",
      UI_KEY_HSTARTWCS,
      "INT",
      0,
      "Header keyword number to start reading WCS.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->hstartwcs,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "hendwcs",
      UI_KEY_HENDWCS,
      "INT",
      0,
      "Header keyword number to end reading WCS.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->hendwcs,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Output. */
    {
      "keepwcs",
      UI_KEY_KEEPWCS,
      0,
      0,
      "Do not apply warp to input's WCS",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->keepwcs,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "coveredfrac",
      UI_KEY_COVEREDFRAC,
      "FLT",
      0,
      "Acceptable fraction of output pixel covered.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->coveredfrac,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },


    {
      0, 0, 0, 0,
      "Warps:",
      UI_GROUP_WARPS
    },
    {
      "align",
      UI_KEY_ALIGN,
      0,
      0,
      "Align the image and celestial axes.",
      UI_GROUP_WARPS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "rotate",
      UI_KEY_ROTATE,
      "FLT",
      0,
      "Rotate by the given angle in degrees.",
      UI_GROUP_WARPS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "scale",
      UI_KEY_SCALE,
      "FLT[,FLT]",
      0,
      "Scale along the given axis(es).",
      UI_GROUP_WARPS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "flip",
      UI_KEY_FLIP,
      "INT[,INT]",
      0,
      "Flip along the given axis(es).",
      UI_GROUP_WARPS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "shear",
      UI_KEY_SHEAR,
      "FLT[,FLT]",
      0,
      "Shear along the given axis(es).",
      UI_GROUP_WARPS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "translate",
      UI_KEY_TRANSLATE,
      "FLT[,FLT]",
      0,
      "Translate along the given axis(es).",
      UI_GROUP_WARPS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "project",
      UI_KEY_PROJECT,
      "FLT[,FLT]",
      0,
      "Project along the given axis(es).",
      UI_GROUP_WARPS,
      0,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "matrix",
      UI_KEY_MATRIX,
      "STR",
      0,
      "Raw transformation matrix, highest priority.",
      UI_GROUP_WARPS,
      &p->matrix,
      GAL_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_modular_warps_ll
    },
    {
      "centeroncorner",
      UI_KEY_CENTERONCORNER,
      0,
      0,
      "Center of coordinates on first pixel corner.",
      UI_GROUP_WARPS,
      &p->centeroncorner,
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
