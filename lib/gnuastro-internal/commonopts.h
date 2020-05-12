/*********************************************************************
Common parameters between all programs.

IMPORTANT: This header must only be included the programs, not libraries.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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
     libraries, and in particular 'options.c'. Because we want each program
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
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "searchin",
      GAL_OPTIONS_KEY_SEARCHIN,
      "STR",
      0,
      "Select column(s): 'name', 'unit', 'comment'.",
      GAL_OPTIONS_GROUP_INPUT,
      &cp->searchin,
      GAL_TYPE_STRING,
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
      "Ignore case in matching/searching columns.",
      GAL_OPTIONS_GROUP_INPUT,
      &cp->ignorecase,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "stdintimeout",
      GAL_OPTIONS_KEY_STDINTIMEOUT,
      "INT",
      0,
      "Micro-seconds to wait for standard input.",
      GAL_OPTIONS_GROUP_INPUT,
      &cp->stdintimeout,
      GAL_TYPE_LONG,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    /* Tile grid (tessellation) options. */
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
      "Regular tile size on dim.s (FITS order).",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.tilesize,
      GAL_TYPE_SIZE_T,
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
      "No. of channels in dim.s (FITS order).",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.numchannels,
      GAL_TYPE_STRING,
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
      GAL_TYPE_FLOAT32,
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
      "oneelempertile",
      GAL_OPTIONS_KEY_ONEELEMPERTILE,
      0,
      0,
      "Display 1 element/tile, not full input res.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->tl.oneelempertile,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "interponlyblank",
      GAL_OPTIONS_KEY_INTERPONLYBLANK,
      0,
      0,
      "Only interpolate over the blank tiles.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->interponlyblank,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "interpmetric",
      GAL_OPTIONS_KEY_INTERPMETRIC,
      "INT",
      0,
      "Interpolation metric (radial, manhattan).",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->interpmetric,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_interpmetric
    },
    {
      "interpnumngb",
      GAL_OPTIONS_KEY_INTERPNUMNGB,
      "INT",
      0,
      "No. of neighbors to use for interpolation.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &cp->interpnumngb,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
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
      "Output file name.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->output,
      GAL_TYPE_STRING,
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
      &cp->type,                /* Internally, 'cp->type' is actually an   */
      GAL_TYPE_STRING,          /* 'uint8_t', but the user gives a string. */
      GAL_OPTIONS_RANGE_GT_0,   /* So for the sanity checks to pass, we    */
      GAL_OPTIONS_NOT_MANDATORY,/* use 'GAL_TYPE_STRING' for this option.  */
      GAL_OPTIONS_NOT_SET,
      gal_options_read_type
    },
    {
      "tableformat",
      GAL_OPTIONS_KEY_TABLEFORMAT,
      "STR",
      0,
      "Table fmt: 'fits-ascii', 'fits-binary', 'txt'.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &cp->tableformat,
      GAL_TYPE_STRING,
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
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minmapsize",
      GAL_OPTIONS_KEY_MINMAPSIZE,
      "INT",
      0,
      "Minimum bytes in array to not use ram RAM.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->minmapsize,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "quietmmap",
      GAL_OPTIONS_KEY_QUIETMMAP,
      0,
      0,
      "Don't print mmap'd file's name and size.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->quietmmap,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
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
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_call_parse_config_file
    },
    {
      "checkconfig",
      GAL_OPTIONS_KEY_CHECKCONFIG,
      0,
      0,
      "List all config files and variables read.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &cp->checkconfig,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_check_config
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
      &cp->onlyversion,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_check_version
    },



    {0}
  };

#endif
