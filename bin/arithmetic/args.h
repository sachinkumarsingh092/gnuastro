/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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



/* Definition of program-specific options. */
struct argp_option program_options[] =
  {
    {
      "globalhdu",
      UI_KEY_GLOBALHDU,
      "STR",
      0,
      "Use this HDU for all inputs, ignore '--hdu'.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->globalhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "wcsfile",
      UI_KEY_WCSFILE,
      "STR",
      0,
      "File to use for output's WCS.",
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
      "HDU/extension to use for output's WCS.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->wcshdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      "onedasimage",
      UI_KEY_ONEDASIMAGE,
      0,
      0,
      "Write 1D outputs as an image, not a table.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->onedasimage,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "onedonstdout",
      UI_KEY_ONEDONSTDOUT,
      0,
      0,
      "Write 1D output on stdout, not in a table.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->onedonstdout,
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
