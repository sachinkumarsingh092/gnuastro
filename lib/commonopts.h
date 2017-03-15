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
      "Select column(s) in: `name', `unit', `comment'.",
      GAL_OPTIONS_GROUP_INPUT,
      &cp->searchin,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_searchin
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




    /* Tile grid options. */
    {
      0, 0, 0, 0,
      "Tessellation (tile grid):",
      GAL_OPTIONS_GROUP_TESSELLATION
    },
    {
      "tilesize",
      GAL_OPTIONS_KEY_TILESIZE,
      "INT[,INT]",
      0,
      "Regular tile size on each dim. (FITS order).",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.tilesize,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },
    {
      "numchannels",
      GAL_OPTIONS_KEY_NUMCHANNELS,
      "INT[,..]",
      0,
      "No. of channels along each dim. (FITS order).",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.numchannels,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },
    {
      "remainderfrac",
      GAL_OPTIONS_KEY_REMAINDERFRAC,
      "FLT",
      0,
      "Fraction of remainder to split last tile.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.remainderfrac,
      GAL_DATA_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0_LT_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
    },
    {
      "workoverch",
      GAL_OPTIONS_KEY_WORKOVERCH,
      0,
      0,
      "Work (not tile) over channel edges.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.workoverch,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checktiles",
      GAL_OPTIONS_KEY_CHECKTILES,
      0,
      0,
      "Tile IDs in an image, the size of input.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.checktiles,
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
      "Output name.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->output,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "type",
      GAL_OPTIONS_KEY_TYPE,
      "STR",
      0,
      "Type of output: e.g., int16, float32, etc...",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->type,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_type
    },
    {
      "tableformat",
      GAL_OPTIONS_KEY_TABLEFORMAT,
      "STR",
      0,
      "Table format: `fits-ascii', `fits-binary'.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->tableformat,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_tableformat
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
      GAL_DATA_TYPE_SIZE_T,
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
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "log",
      GAL_OPTIONS_KEY_LOG,
      0,
      0,
      "Information about output(s) in a log file.",
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
      GAL_OPTIONS_NOT_SET,
      gal_options_print_citation
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
      "Read configuration file STR immediately.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      NULL,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_call_parse_config_file
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
      GAL_OPTIONS_NOT_SET,
      gal_options_check_version
    },



    {0}
  };

#endif
