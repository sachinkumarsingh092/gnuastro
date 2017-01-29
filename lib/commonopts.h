/*********************************************************************
Common parameters between all programs.

IMPORTANT: This header must only be included the programs, not libraries.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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
#ifndef __GAL_COMMONOPTS_H__
#define __GAL_COMMONOPTS_H__



/* Common options for all programs.

   IMPORTANT NOTE:

     This header should only be included in the programs, not the
     libraries, and in particular `options.c'. Because we want each program
     to have its own allocation of the common options structure. If it is
     included in options.c, then it will be shared between all the
     programs. */
struct argp_option gal_commonopts_options[] =
  {
    {
      0, 0, 0, 0,
      "Input:",
      GAL_OPTIONS_GROUP_INPUT
    },
    {
      "hdu",
      GAL_OPTIONS_KEY_HDU,
      "STR/INT",
      0,
      "Extension name or number of input data.",
      GAL_OPTIONS_GROUP_INPUT,
      &cp->hdu,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "searchin",
      GAL_OPTIONS_KEY_SEARCHIN,
      "STR",
      0,
      "Col. selection field: `name', `unit', `comment'.",
      GAL_OPTIONS_GROUP_INPUT,
      &cp->searchinstr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ignorecase",
      GAL_OPTIONS_KEY_IGNORECASE,
      0,
      0,
      "Ignore case when matching/searching col. info.",
      GAL_OPTIONS_GROUP_INPUT,
      &cp->ignorecase,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Output:",
      GAL_OPTIONS_GROUP_OUTPUT
    },
    {
      "output",
      GAL_OPTIONS_KEY_OUTPUT,
      "STR",
      0,
      "Output file or directory name.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->output,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "tableformat",
      GAL_OPTIONS_KEY_TABLEFORMAT,
      "STR",
      0,
      "Output table format: `fits-ascii', `fits-binary'.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->tableformatstr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "dontdelete",
      GAL_OPTIONS_KEY_DONTDELETE,
      0,
      0,
      "Don't delete output if it exists.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->dontdelete,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "keepinputdir",
      GAL_OPTIONS_KEY_KEEPINPUTDIR,
      0,
      0,
      "Keep input directory for automatic output.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->keepinputdir,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Operating modes:",
      GAL_OPTIONS_GROUP_OPERATING_MODE
    },
    {
      "quiet",
      GAL_OPTIONS_KEY_QUIET,
      0,
      0,
      "Only report errors, remain quiet about steps.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->quiet,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "numthreads",
      GAL_OPTIONS_KEY_NUMTHREADS,
      "INT",
      0,
      "Number of CPU threads to use.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->numthreads,
      GAL_DATA_TYPE_ULONG,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minmapsize",
      GAL_OPTIONS_KEY_MINMAPSIZE,
      "INT",
      0,
      "Minimum no. bytes to map arrays to hdd/ssd.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->minmapsize,
      GAL_DATA_TYPE_ULONG,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "log",
      GAL_OPTIONS_KEY_LOG,
      0,
      0,
      "No log file for programs which make one.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->log,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Internal (before control goes back to the program). */
    {
      "cite",
      GAL_OPTIONS_KEY_CITE,
      0,
      0,
      "BibTeX citation for this program.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      NULL,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "printparams",
      GAL_OPTIONS_KEY_PRINTPARAMS,
      0,
      0,
      "Print parameter values to be used and abort.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->printparams,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "config",
      GAL_OPTIONS_KEY_CONFIG,
      "STR",
      0,
      "Read file STR before continuing.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      NULL,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "setdirconf",
      GAL_OPTIONS_KEY_SETDIRCONF,
      0,
      0,
      "Set default values for this directory and abort.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->setdirconf,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "setusrconf",
      GAL_OPTIONS_KEY_SETUSRCONF,
      0,
      0,
      "Set default values for this user and abort.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->setusrconf,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "lastconfig",
      GAL_OPTIONS_KEY_LASTCONFIG,
      0,
      0,
      "Do not parse any more configuration files.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->lastconfig,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "onlyversion",
      GAL_OPTIONS_KEY_ONLYVERSION,
      "STR",
      0,
      "Only run if the program version is STR.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      NULL,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    {0}
  };

#endif
