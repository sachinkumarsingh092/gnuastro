/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

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
    {
      "column",
      UI_KEY_COLUMN,
      "STR",
      0,
      "Column number (counting from 1) or search string.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->columns,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "wcsfile",
      UI_KEY_WCSFILE,
      "STR",
      0,
      "File with WCS if conversion is requested.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->wcsfile,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "wcshdu",
      UI_KEY_WCSHDU,
      "STR",
      0,
      "HDU in file with WCS for conversion.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->wcshdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "catcolumnfile",
      UI_KEY_CATCOLUMNFILE,
      "STR",
      0,
      "File(s) to be concatenated by column.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->catcolumnfile,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "catcolumnhdu",
      UI_KEY_CATCOLUMNHDU,
      "STR/INT",
      0,
      "HDU/Extension(s) in catcolumnfile.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->catcolumnhdu,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "catcolumns",
      UI_KEY_CATCOLUMNS,
      "STR",
      0,
      "Columns to use in catcolumnfile.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->catcolumns,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    /* Output. */
    {
      "information",
      UI_KEY_INFORMATION,
      0,
      0,
      "Only print table and column information.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->information,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "colinfoinstdout",
      UI_KEY_COLINFOINSTDOUT,
      0,
      0,
      "Column info/metadata when printing to stdout.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->colinfoinstdout,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "catcolumnrawname",
      UI_KEY_CATCOLUMNRAWNAME,
      0,
      0,
      "Don't touch column names of --catcolumnfile.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->catcolumnrawname,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "colmetadata",
      UI_KEY_COLMETADATA,
      "STR,STR[,STR,STR]",
      0,
      "Update output metadata (name, unit, comments).",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->colmetadata,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_name_and_strings
    },





    /* Output Rows */
    {
      0, 0, 0, 0,
      "Rows in output:",
      UI_GROUP_OUTROWS
    },
    {
      "range",
      UI_KEY_RANGE,
      "STR,FLT:FLT",
      0,
      "Column, and range to limit output.",
      UI_GROUP_OUTROWS,
      &p->range,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_name_and_float64s
    },
    {
      "inpolygon",
      UI_KEY_INPOLYGON,
      "STR,STR",
      0,
      "Coord. columns that are inside '--polygon'.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->inpolygon,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings
    },
    {
      "outpolygon",
      UI_KEY_OUTPOLYGON,
      "STR,STR",
      0,
      "Coord. columns that are outside '--polygon'.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->outpolygon,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings
    },
    {
      "polygon",
      UI_KEY_POLYGON,
      "FLT:FLT[,...]",
      0,
      "Polygon for '--inpolygon' or '--outpolygon'.",
      UI_GROUP_OUTROWS,
      &p->polygon,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_colon_sep_csv
    },
    {
      "equal",
      UI_KEY_EQUAL,
      "STR,FLT[,...]",
      0,
      "Column, values to keep in output.",
      UI_GROUP_OUTROWS,
      &p->equal,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_name_and_strings
    },
    {
      "notequal",
      UI_KEY_NOTEQUAL,
      "STR,FLT[,...]",
      0,
      "Column, values to remove from output.",
      UI_GROUP_OUTROWS,
      &p->notequal,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_name_and_strings
    },
    {
      "sort",
      UI_KEY_SORT,
      "STR/INT",
      0,
      "Column name or number for sorting.",
      UI_GROUP_OUTROWS,
      &p->sort,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "descending",
      UI_KEY_DESCENDING,
      0,
      0,
      "Sort in descending order: largets first.",
      UI_GROUP_OUTROWS,
      &p->descending,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "head",
      UI_KEY_HEAD,
      "INT",
      0,
      "Only output given number of top rows.",
      UI_GROUP_OUTROWS,
      &p->head,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tail",
      UI_KEY_TAIL,
      "INT",
      0,
      "Only output given number of bottom rows.",
      UI_GROUP_OUTROWS,
      &p->tail,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "noblank",
      UI_KEY_NOBLANK,
      "STR[,STR]",
      0,
      "Remove rows with blank in given columns.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->noblank,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings
    },





    /* End. */
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
