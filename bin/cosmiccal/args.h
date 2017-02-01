/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

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
      "redshift",
      ARGS_OPTION_KEY_REDSHIFT,
      "FLT",
      0,
      "Redshift of interest.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->redshift,
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "H0",
      ARGS_OPTION_KEY_H0,
      "FLT",
      0,
      "Current expansion rate (Hubble constant).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->H0,
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "olambda",
      ARGS_OPTION_KEY_OLAMBDA,
      "FLT",
      0,
      "Current cosmological cst. dens. per crit. dens.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->olambda,
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "omatter",
      ARGS_OPTION_KEY_OMATTER,
      "FLT",
      0,
      "Current matter density per critical density.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->omatter,
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "oradiation",
      ARGS_OPTION_KEY_ORADIATION,
      "FLT",
      0,
      "Current radiation density per critical density.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->oradiation,
      GAL_DATA_TYPE_DOUBLE,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {
      "onlyvolumne",
      ARGS_OPTION_KEY_ONLYVOLUME,
      0,
      0,
      "Only print comoving volume in Mpc^3.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->onlyvolume,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "onlyabsmagconv",
      ARGS_OPTION_KEY_ONLYABSMAGCONV,
      0,
      0,
      "Only print conversion to absolute magnitude.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->onlyabsmagconv,
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
