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
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "refcol",
      ARGS_OPTION_KEY_REFCOL,
      "STR",
      0,
      "Reference column name or number.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->refcol,
      GAL_TYPE_STRING,
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
      GAL_TYPE_FLOAT32,
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
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "qrange",
      ARGS_OPTION_KEY_QRANGE,
      "FLT[,FLT]",
      0,
      "Quantile range: one (from Q to 1-Q) or two.",
      GAL_OPTIONS_GROUP_INPUT,
      NULL,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_read_quantile_range
    },



    /* Tessellation */
    {
      "interpolate",
      ARGS_OPTION_KEY_INTERPOLATE,
      0,
      0,
      "Interpolate over blank tiles to fill them.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &p->interpolate,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {
      0, 0, 0, 0,
      "Single value measurements",
      ARGS_GROUP_SINGLE_VALUE
    },
    {
      "number",
      ARGS_OPTION_KEY_NUMBER,
      0,
      0,
      "Number (non-blank).",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "minimum",
      ARGS_OPTION_KEY_MINIMUM,
      0,
      0,
      "Minimum.",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "maximum",
      ARGS_OPTION_KEY_MAXIMUM,
      0,
      0,
      "Maximum.",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "sum",
      ARGS_OPTION_KEY_SUM,
      0,
      0,
      "Sum.",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "mean",
      ARGS_OPTION_KEY_MEAN,
      0,
      0,
      "Mean.",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "std",
      ARGS_OPTION_KEY_STD,
      0,
      0,
      "Standad deviation.",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "median",
      ARGS_OPTION_KEY_MEDIAN,
      0,
      0,
      "Median.",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "quantile",
      ARGS_OPTION_KEY_QUANTILE,
      "FLT[,...]",
      0,
      "Quantile (multiple values acceptable).",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "quantfunc",
      ARGS_OPTION_KEY_QUANTFUNC,
      "FLT[,...]",
      0,
      "Quantile function (multiple values acceptable).",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "mode",
      ARGS_OPTION_KEY_MODE,
      0,
      0,
      "Mode (Appendix C of arXiv:1505.01664).",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "modequant",
      ARGS_OPTION_KEY_MODEQUANT,
      0,
      0,
      "Mode quantile (see --mode)",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "modesym",
      ARGS_OPTION_KEY_MODESYM,
      0,
      0,
      "Mode symmetricity (see --mode).",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "modesymvalue",
      ARGS_OPTION_KEY_MODESYMVALUE,
      0,
      0,
      "Value at mode symmetricity (see --mode).",
      ARGS_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
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
      "Print an ASCII histogram.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->asciihist,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "asciicfp",
      ARGS_OPTION_KEY_ASCIICFP,
      0,
      0,
      "Print an ASCII cumulative frequency plot.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->asciicfp,
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
      "Save the cumulative frequency plot in output.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->cumulative,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mirror",
      ARGS_OPTION_KEY_MIRROR,
      "FLT",
      0,
      "Save the histogram and CFP of the mirror dist.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->mirror,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ontile",
      ARGS_OPTION_KEY_ONTILE,
      0,
      0,
      "Single values on separate tiles, not full input.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->ontile,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sky",
      ARGS_OPTION_KEY_SKY,
      0,
      0,
      "Find the Sky and its STD over the tessellation.",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->sky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sigmaclip",
      ARGS_OPTION_KEY_SIGMACLIP,
      0,
      0,
      "Overall sigma-clipping (see `--sclipparams')",
      ARGS_GROUP_PARTICULAR_STAT,
      &p->sigmaclip,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      0, 0, 0, 0,
      "Sky and Sky STD settings",
      ARGS_GROUP_SKY
    },
    {
      "kernel",
      ARGS_OPTION_KEY_KERNEL,
      "STR",
      0,
      "File name of kernel to convolve input.",
      ARGS_GROUP_SKY,
      &p->kernelname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "khdu",
      ARGS_OPTION_KEY_KHDU,
      "STR",
      0,
      "HDU/extension name or number of kernel.",
      ARGS_GROUP_SKY,
      &p->khdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mirrordist",
      ARGS_OPTION_KEY_MIRRORDIST,
      "FLT",
      0,
      "Max. distance (error multip.) to find mode.",
      ARGS_GROUP_SKY,
      &p->mirrordist,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "modmedqdiff",
      ARGS_OPTION_KEY_MODMEDQDIFF,
      "FLT",
      0,
      "Max. mode and median quantile diff. per tile.",
      ARGS_GROUP_SKY,
      &p->modmedqdiff,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sclipparams",
      ARGS_OPTION_KEY_SCLIPPARAMS,
      "FLT,FLT",
      0,
      "Sigma clip: Multiple, and tolerance/number.",
      ARGS_GROUP_SKY,
      p->sclipparams,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_sigma_clip
    },
    {
      "smoothwidth",
      ARGS_OPTION_KEY_SMOOTHWIDTH,
      "INT",
      0,
      "Sky: flat kernel width to smooth interpolated.",
      ARGS_GROUP_SKY,
      &p->smoothwidth,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_0_OR_ODD,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checksky",
      ARGS_OPTION_KEY_CHECKSKY,
      0,
      0,
      "Store steps in `_sky_steps.fits' file.",
      ARGS_GROUP_SKY,
      &p->checksky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {
      0, 0, 0, 0,
      "Histogram and CFP settings",
      ARGS_GROUP_HIST_CFP
    },
    {
      "numbins",
      ARGS_OPTION_KEY_NUMBINS,
      "INT",
      0,
      "No. of bins in histogram or CFP tables.",
      ARGS_GROUP_HIST_CFP,
      &p->numbins,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "numasciibins",
      ARGS_OPTION_KEY_NUMASCIIBINS,
      "INT",
      0,
      "No. of bins in ASCII histogram or CFP plots.",
      ARGS_GROUP_HIST_CFP,
      &p->numasciibins,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "asciiheight",
      ARGS_OPTION_KEY_ASCIIHEIGHT,
      "INT",
      0,
      "Height of ASCII histogram or CFP plots.",
      ARGS_GROUP_HIST_CFP,
      &p->asciiheight,
      GAL_TYPE_SIZE_T,
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
    {
      "onebinstart",
      ARGS_OPTION_KEY_ONEBINSTART,
      "FLT",
      0,
      "Shift bins so one bin starts on this value.",
      ARGS_GROUP_HIST_CFP,
      &p->onebinstart,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
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
