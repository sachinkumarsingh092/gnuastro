/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

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
      "redshift",
      UI_KEY_REDSHIFT,
      "FLT",
      0,
      "Redshift of interest.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->redshift,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "obsline",
      UI_KEY_OBSLINE,
      "STR,FLT",
      0,
      "Redshift from line and observed wavelength.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->obsline,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_parse_obsline
    },
    {
      "H0",
      UI_KEY_H0,
      "FLT",
      0,
      "Current expansion rate (Hubble constant).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->H0,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "olambda",
      UI_KEY_OLAMBDA,
      "FLT",
      0,
      "Current cosmological cst. dens. per crit. dens.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->olambda,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "omatter",
      UI_KEY_OMATTER,
      "FLT",
      0,
      "Current matter density per critical density.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->omatter,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "oradiation",
      UI_KEY_ORADIATION,
      "FLT",
      0,
      "Current radiation density per critical density.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->oradiation,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    /* Output options */
    {
      "listlines",
      UI_KEY_LISTLINES,
      0,
      0,
      "List known spectral lines.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->listlines,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Specific calculations",
      UI_GROUP_SPECIFIC
    },
    {
      "usedredshift",
      UI_KEY_USEDREDSHIFT,
      0,
      0,
      "Used redshift in this run.",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "agenow",
      UI_KEY_AGENOW,
      0,
      0,
      "Age of universe now (Ga: Giga Annum).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "criticaldensitynow",
      UI_KEY_CRITICALDENSITYNOW,
      0,
      0,
      "Critical density now (g/cm^3).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "properdistance",
      UI_KEY_PROPERDISTANCE,
      0,
      0,
      "Proper distance to z (Mpc).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "angulardimdist",
      UI_KEY_ANGULARDIMDIST,
      0,
      0,
      "Angular diameter distance (Mpc).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "arcsectandist",
      UI_KEY_ARCSECTANDIST,
      0,
      0,
      "Tangential dist. covered by 1arcsec at z (kpc).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "luminositydist",
      UI_KEY_LUMINOSITYDIST,
      0,
      0,
      "Luminosity distance to z (Mpc).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "distancemodulus",
      UI_KEY_DISTANCEMODULUS,
      0,
      0,
      "Distance modulus at z (no units).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "absmagconv",
      UI_KEY_ABSMAGCONV,
      0,
      0,
      "Conversion to absolute magnitude (no unit).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "age",
      UI_KEY_AGE,
      0,
      0,
      "Age of universe at z (Ga: Giga Annum).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "lookbacktime",
      UI_KEY_LOOKBACKTIME,
      0,
      0,
      "Look back time to z (Ga: Giga Annum).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "criticaldensity",
      UI_KEY_CRITICALDENSITY,
      0,
      0,
      "Critical density at z (g/cm^3).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "volume",
      UI_KEY_VOLUME,
      0,
      0,
      "Comoving volume (4pi str) to z (Mpc^3).",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
    },
    {
      "lineatz",
      UI_KEY_LINEATZ,
      "STR/FLT",
      0,
      "Wavelength of given line at chosen redshift",
      UI_GROUP_SPECIFIC,
      &p->specific,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value,
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
