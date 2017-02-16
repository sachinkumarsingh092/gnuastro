/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

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
      ARGS_OPTION_KEY_COLUMN,
      "STR",
      0,
      "Column name or number if input is a table.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->column,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "greaterequal",
      ARGS_OPTION_KEY_GREATEREQUAL,
      "FLT",
      0,
      "Only use values greater-equal than this.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->greaterequal,
      GAL_DATA_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "lessthan",
      ARGS_OPTION_KEY_LESSTHAN,
      "FLT",
      0,
      "Only use values less than this.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->lessthan,
      GAL_DATA_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "quantrange",
      ARGS_OPTION_KEY_QUANTRANGE,
      "FLT",
      0,
      "Quantile (Q) range: from Q to 1-Q.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->quantilerange,
      GAL_DATA_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Values to print in one row",
      ARGS_GROUP_IN_ONE_ROW
    },
    {
      "number",
      ARGS_OPTION_KEY_NUMBER,
      0,
      0,
      "Print the number of data-points.",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_DATA_TYPE_INVALID,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },
    {
      "minimum",
      ARGS_OPTION_KEY_MINIMUM,
      0,
      0,
      "Print the minimum.",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },
    {
      "maximum",
      ARGS_OPTION_KEY_MAXIMUM,
      0,
      0,
      "Print the maximum.",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },
    {
      "sum",
      ARGS_OPTION_KEY_SUM,
      0,
      0,
      "Print the sum of all elements.",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },
    {
      "mean",
      ARGS_OPTION_KEY_MEAN,
      0,
      0,
      "Print the mean.",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },
    {
      "std",
      ARGS_OPTION_KEY_STD,
      0,
      0,
      "Print the standad deviation.",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },
    {
      "median",
      ARGS_OPTION_KEY_MEDIAN,
      0,
      0,
      "Print the median.",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },
    {
      "mode",
      ARGS_OPTION_KEY_MODE,
      0,
      0,
      "Print the mode (for large datasets).",
      ARGS_GROUP_IN_ONE_ROW,
      &p->toprint,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_print_in_row
    },





    {
      0, 0, 0, 0,
      "Particular calculation",
      ARGS_GROUP_PARTICULAR_STAT
    },
    {
      "asciihist",
      ARGS_OPTION_KEY_ASCIIHIST,
      0,
      0,
      "Print an ASCII histogram",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->asciihist,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "histogram",
      ARGS_OPTION_KEY_HISTOGRAM,
      0,
      0,
      "Save the histogram in output.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->histogram,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "cumulative",
      ARGS_OPTION_KEY_CUMULATIVE,
      0,
      0,
      "Save the cumulative frequency plot.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->cumulative,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sigmaclip",
      ARGS_OPTION_KEY_SIGMACLIP,
      "FLT,FLT",
      0,
      "Multiple of sigma and tolerance or number.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->sigclipstr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mirror",
      ARGS_OPTION_KEY_MIRROR,
      "FLT",
      0,
      "Quantile of mirror distribution to save.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->mirrorquant,
      GAL_DATA_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0_LT_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Histogram or Cumulative frequency plot settings",
      ARGS_GROUP_HIST_CFP
    },
    {
      "numbins",
      ARGS_OPTION_KEY_NUMBINS,
      "INT",
      0,
      "Number of bins in histogram or CFP.",
      ARGS_GROUP_HIST_CFP,
      &p->numbins,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "lowerbin",
      ARGS_OPTION_KEY_LOWERBIN,
      0,
      0,
      "Save interval lower limit, not center",
      ARGS_GROUP_HIST_CFP,
      &p->lowerbin,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "normalize",
      ARGS_OPTION_KEY_NORMALIZE,
      0,
      0,
      "Set sum of all bins to 1.",
      ARGS_GROUP_HIST_CFP,
      &p->normalize,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "onebinstart",
      ARGS_OPTION_KEY_ONEBINSTART,
      "FLT",
      0,
      "Shift bins so one bin starts on this value.",
      ARGS_GROUP_HIST_CFP,
      &p->onebinstart,
      GAL_DATA_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "maxbinone",
      ARGS_OPTION_KEY_MAXBINONE,
      0,
      0,
      "Scale such that the maximum bin has value of one.",
      ARGS_GROUP_HIST_CFP,
      &p->maxbinone,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_GT_0,
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
