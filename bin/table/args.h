/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

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
    {
      0, 0, 0, 0,
      "Input:",
      GAL_OPTIONS_GROUP_INPUT
    },
    {
      "column",
      ARGS_OPTION_KEY_COLUMN,
      "STR",
      0,
      "Column number (counting from 1) or search string.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->columns,
      GAL_DATA_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "searchin",
      ARGS_OPTION_KEY_SEARCHIN,
      "STR",
      0,
      "Search in column `name', `units', or `comments'.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->searchinstr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ignorecase",
      ARGS_OPTION_KEY_IGNORECASE,
      0,
      0,
      "Ignore case when matching column information.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->ignorecase,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {
      0, 0, 0, 0,
      "Output:",
      2
    },
    {
      "tabletype",
      ARGS_OPTION_KEY_TABLETYPE,
      "STR",
      0,
      "Output table type: `fits-ascii', `fits-binary'.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->tabletypestr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {
      0, 0, 0, 0,
      "Operating modes:",
      GAL_OPTIONS_GROUP_OPERATING_MODE
    },
    {
      "information",
      ARGS_OPTION_KEY_INFORMATION,
      0,
      0,
      "Only print table and column information.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->information,
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
