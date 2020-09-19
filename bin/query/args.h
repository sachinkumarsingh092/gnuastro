/*********************************************************************
Query - Retreive data from a remote data server.
Query is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2020, Free Software Foundation, Inc.

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
    /* Database and dataset. */
    {
      "database",
      UI_KEY_DATABASE,
      "STR",
      0,
      "Name of database (e.g., 'esa').",
      GAL_OPTIONS_GROUP_INPUT,
      &p->database,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_parse_database
    },
    {
      "query",
      UI_KEY_QUERY,
      "STR",
      0,
      "The raw query as a simple string.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->query,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
    },





    /* Query by center. */
    {
      0, 0, 0, 0,
      "Generate query by center (not compatible with '--query'):",
      UI_GROUP_BYCENTER,
    },
    {
      "dataset",
      UI_KEY_DATASET,
      "STR",
      0,
      "Name of dataset in database.",
      UI_GROUP_BYCENTER,
      &p->datasetstr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
    },
    {
      "center",
      UI_KEY_CENTER,
      "FLT[,...]",
      0,
      "Central coordinates of the query.",
      UI_GROUP_BYCENTER,
      &p->center,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },
    {
      "radius",
      UI_KEY_RADIUS,
      "FLT",
      0,
      "Radius around center to select targets.",
      UI_GROUP_BYCENTER,
      &p->radius,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },
    {
      "width",
      UI_KEY_WIDTH,
      "FLT[,FLT]",
      0,
      "Width of box to select targets.",
      UI_GROUP_BYCENTER,
      &p->width,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },
    {
      "range",
      UI_KEY_RANGE,
      "STR,FLT:FLT",
      0,
      "Range of selected targets in given column.",
      UI_GROUP_BYCENTER,
      &p->range,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_name_and_float64s
    },
    {
      "column",
      UI_KEY_COLUMN,
      "STR",
      0,
      "Column names to download from catalog.",
      UI_GROUP_BYCENTER,
      &p->columns,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
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
