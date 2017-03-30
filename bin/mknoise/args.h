/*********************************************************************
MakeNoise - Add noise to a dataset.
MakeNoise is part of GNU Astronomy Utilities (Gnuastro) package.

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
      "stdadd",
      ARGS_OPTION_KEY_STDADD,
      "FLT",
      0,
      "Standard deviation addition constant.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->stdadd,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "background",
      ARGS_OPTION_KEY_BACKGROUND,
      "FLT",
      0,
      "Fixed background magnitude for whole input.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->background_mag,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "zeropoint",
      ARGS_OPTION_KEY_ZEROPOINT,
      "FLT",
      0,
      "Zeropoint magnitude of input.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->zeropoint,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {
      "envseed",
      ARGS_OPTION_KEY_ENVSEED,
      0,
      0,
      "Use GSL_RNG_SEED environment variable for seed",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->envseed,
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
