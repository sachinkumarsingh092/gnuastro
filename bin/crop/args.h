/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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
      "mode",
      UI_KEY_MODE,
      "STR",
      0,
      "Coordinate mode `img' or `wcs'",
      GAL_OPTIONS_GROUP_INPUT,
      &p->mode,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_parse_coordinate_mode
    },
    {
      "hstartwcs",
      UI_KEY_HSTARTWCS,
      "INT",
      0,
      "Header keyword number to start reading WCS.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->hstartwcs,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "hendwcs",
      UI_KEY_HENDWCS,
      "INT",
      0,
      "Header keyword number to stop reading WCS.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->hendwcs,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "zeroisnotblank",
      UI_KEY_ZEROISNOTBLANK,
      0,
      0,
      "0.0 in float or double images are not blank.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->zeroisnotblank,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    /* Output. */
    {
      "noblank",
      UI_KEY_NOBLANK,
      0,
      0,
      "Remove parts of the crop box out of input image.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->noblank,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "suffix",
      UI_KEY_SUFFIX,
      "STR",
      0,
      "Suffix (postfix) of cropped images.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->suffix,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Crop by center",
      ARGS_GROUP_CENTER_GENERAL
    },
    {
      "checkcenter",
      UI_KEY_CHECKCENTER,
      "INT",
      0,
      "Width (in pixels) of box at center to check.",
      ARGS_GROUP_CENTER_GENERAL,
      &p->checkcenter,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_0_OR_ODD,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "width",
      UI_KEY_WIDTH,
      "FLT[,...]",
      0,
      "Width when crop is defined by its center.",
      ARGS_GROUP_CENTER_GENERAL,
      &p->width,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_parse_width_and_center
    },
    {
      "center",
      UI_KEY_CENTER,
      "FLT[,...]",
      0,
      "Central coordinates of a single crop.",
      ARGS_GROUP_CENTER_GENERAL,
      &p->center,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_parse_width_and_center
    },







    {
      0, 0, 0, 0,
      "Crop by center (when a catalog is given)",
      ARGS_GROUP_CENTER_CATALOG
    },
    {
      "catalog",
      UI_KEY_CATALOG,
      "STR",
      0,
      "Input catalog filename.",
      ARGS_GROUP_CENTER_CATALOG,
      &p->catname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "cathdu",
      UI_KEY_CATHDU,
      "STR/INT",
      0,
      "HDU of catalog, if it is a FITS table.",
      ARGS_GROUP_CENTER_CATALOG,
      &p->cathdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "namecol",
      UI_KEY_NAMECOL,
      "STR/INT",
      0,
      "Column no./info of crop filename (no suffix).",
      ARGS_GROUP_CENTER_CATALOG,
      &p->namecol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "coordcol",
      UI_KEY_COORDCOL,
      "STR/INT",
      0,
      "Column no./info containing coordinates.",
      ARGS_GROUP_CENTER_CATALOG,
      &p->coordcol,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Crop by region",
      ARGS_GROUP_REGION
    },
    {
      "section",
      UI_KEY_SECTION,
      "STR",
      0,
      "Image section string specifying crop range.",
      ARGS_GROUP_REGION,
      &p->section,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "polygon",
      UI_KEY_POLYGON,
      "STR",
      0,
      "Polygon vertices of region to crop, keep inside.",
      ARGS_GROUP_REGION,
      &p->polygon,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "outpolygon",
      UI_KEY_OUTPOLYGON,
      0,
      0,
      "Keep the polygon's outside, mask the inside.",
      ARGS_GROUP_REGION,
      &p->outpolygon,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },






    /* Operating mode */




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
