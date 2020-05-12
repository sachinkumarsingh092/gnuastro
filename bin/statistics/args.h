/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

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






/* Array of acceptable options. */
struct argp_option program_options[] =
  {
    {
      "column",
      UI_KEY_COLUMN,
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
      UI_KEY_REFCOL,
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
      UI_KEY_GREATEREQUAL,
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
      UI_KEY_LESSTHAN,
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
      UI_KEY_QRANGE,
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
      UI_KEY_INTERPOLATE,
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
      UI_GROUP_SINGLE_VALUE
    },
    {
      "number",
      UI_KEY_NUMBER,
      0,
      0,
      "Number (non-blank).",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "minimum",
      UI_KEY_MINIMUM,
      0,
      0,
      "Minimum.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "maximum",
      UI_KEY_MAXIMUM,
      0,
      0,
      "Maximum.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "sum",
      UI_KEY_SUM,
      0,
      0,
      "Sum.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "mean",
      UI_KEY_MEAN,
      0,
      0,
      "Mean.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "std",
      UI_KEY_STD,
      0,
      0,
      "Standad deviation.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "median",
      UI_KEY_MEDIAN,
      0,
      0,
      "Median.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "quantile",
      UI_KEY_QUANTILE,
      "FLT[,...]",
      0,
      "Quantile (multiple values acceptable).",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "quantfunc",
      UI_KEY_QUANTFUNC,
      "FLT[,...]",
      0,
      "Quantile function (multiple values acceptable).",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "mode",
      UI_KEY_MODE,
      0,
      0,
      "Mode (Appendix C of arXiv:1505.01664).",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "modequant",
      UI_KEY_MODEQUANT,
      0,
      0,
      "Mode quantile (see --mode)",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "modesym",
      UI_KEY_MODESYM,
      0,
      0,
      "Mode symmetricity (see --mode).",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "modesymvalue",
      UI_KEY_MODESYMVALUE,
      0,
      0,
      "Value at mode symmetricity (see --mode).",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "sigclip-number",
      UI_KEY_SIGCLIPNUMBER,
      0,
      0,
      "Number of elements after sigma-clipping.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "sigclip-median",
      UI_KEY_SIGCLIPMEDIAN,
      0,
      0,
      "Sigma-clipped median.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "sigclip-mean",
      UI_KEY_SIGCLIPMEAN,
      0,
      0,
      "Sigma-clipped mean.",
      UI_GROUP_SINGLE_VALUE,
      &p->singlevalue,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      ui_add_to_single_value
    },
    {
      "sigclip-std",
      UI_KEY_SIGCLIPSTD,
      0,
      0,
      "Sigma-clipped standard deviation.",
      UI_GROUP_SINGLE_VALUE,
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
      UI_GROUP_PARTICULAR_STAT
    },
    {
      "asciihist",
      UI_KEY_ASCIIHIST,
      0,
      0,
      "Print an ASCII histogram.",
      UI_GROUP_PARTICULAR_STAT,
      &p->asciihist,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "asciicfp",
      UI_KEY_ASCIICFP,
      0,
      0,
      "Print an ASCII cumulative frequency plot.",
      UI_GROUP_PARTICULAR_STAT,
      &p->asciicfp,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "histogram",
      UI_KEY_HISTOGRAM,
      0,
      0,
      "Save the histogram in output.",
      UI_GROUP_PARTICULAR_STAT,
      &p->histogram,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "cumulative",
      UI_KEY_CUMULATIVE,
      0,
      0,
      "Save the cumulative frequency plot in output.",
      UI_GROUP_PARTICULAR_STAT,
      &p->cumulative,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mirror",
      UI_KEY_MIRROR,
      "FLT",
      0,
      "Save the histogram and CFP of the mirror dist.",
      UI_GROUP_PARTICULAR_STAT,
      &p->mirror,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ontile",
      UI_KEY_ONTILE,
      0,
      0,
      "Single values on separate tiles, not full input.",
      UI_GROUP_PARTICULAR_STAT,
      &p->ontile,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sky",
      UI_KEY_SKY,
      0,
      0,
      "Find the Sky and its STD over the tessellation.",
      UI_GROUP_PARTICULAR_STAT,
      &p->sky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sigmaclip",
      UI_KEY_SIGMACLIP,
      0,
      0,
      "Overall sigma-clipping (see '--sclipparams')",
      UI_GROUP_PARTICULAR_STAT,
      &p->sigmaclip,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "contour",
      UI_KEY_CONTOUR,
      "STR",
      0,
      "Contour levels, save in PGFPlots format.",
      UI_GROUP_PARTICULAR_STAT,
      &p->contour,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_float64
    },






    {
      0, 0, 0, 0,
      "Sky and Sky STD settings",
      UI_GROUP_SKY
    },
    {
      "kernel",
      UI_KEY_KERNEL,
      "STR",
      0,
      "File name of kernel to convolve input.",
      UI_GROUP_SKY,
      &p->kernelname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "khdu",
      UI_KEY_KHDU,
      "STR",
      0,
      "HDU/extension name or number of kernel.",
      UI_GROUP_SKY,
      &p->khdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "mirrordist",
      UI_KEY_MIRRORDIST,
      "FLT",
      0,
      "Max. distance (error multip.) to find mode.",
      UI_GROUP_SKY,
      &p->mirrordist,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "meanmedqdiff",
      UI_KEY_MEANMEDQDIFF,
      "FLT",
      0,
      "Max. mode and median quantile diff. per tile.",
      UI_GROUP_SKY,
      &p->meanmedqdiff,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "outliersclip",
      UI_KEY_OUTLIERSCLIP,
      "FLT,FLT",
      0,
      "Sigma-clip params for qthresh outliers.",
      UI_GROUP_SKY,
      &p->outliersclip,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_sigma_clip
    },
    {
      "outliersigma",
      UI_KEY_OUTLIERSIGMA,
      "FLT",
      0,
      "Multiple of sigma to define outliers.",
      UI_GROUP_SKY,
      &p->outliersigma,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sclipparams",
      UI_KEY_SCLIPPARAMS,
      "FLT,FLT",
      0,
      "Sigma clip: Multiple, and tolerance/number.",
      UI_GROUP_SKY,
      p->sclipparams,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_sigma_clip
    },
    {
      "smoothwidth",
      UI_KEY_SMOOTHWIDTH,
      "INT",
      0,
      "Sky: flat kernel width to smooth interpolated.",
      UI_GROUP_SKY,
      &p->smoothwidth,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_0_OR_ODD,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checksky",
      UI_KEY_CHECKSKY,
      0,
      0,
      "Store steps in '_sky_steps.fits' file.",
      UI_GROUP_SKY,
      &p->checksky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "ignoreblankintiles",
      UI_KEY_IGNOREBLANKINTILES,
      0,
      0,
      "Don't write input's blanks in the tiled output.",
      UI_GROUP_SKY,
      &p->ignoreblankintiles,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    {
      0, 0, 0, 0,
      "Histogram and CFP settings",
      UI_GROUP_HIST_CFP
    },
    {
      "numbins",
      UI_KEY_NUMBINS,
      "INT",
      0,
      "No. of bins in histogram or CFP tables.",
      UI_GROUP_HIST_CFP,
      &p->numbins,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "numasciibins",
      UI_KEY_NUMASCIIBINS,
      "INT",
      0,
      "No. of bins in ASCII histogram or CFP plots.",
      UI_GROUP_HIST_CFP,
      &p->numasciibins,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "asciiheight",
      UI_KEY_ASCIIHEIGHT,
      "INT",
      0,
      "Height of ASCII histogram or CFP plots.",
      UI_GROUP_HIST_CFP,
      &p->asciiheight,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "normalize",
      UI_KEY_NORMALIZE,
      0,
      0,
      "Set sum of all bins to 1.",
      UI_GROUP_HIST_CFP,
      &p->normalize,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "maxbinone",
      UI_KEY_MAXBINONE,
      0,
      0,
      "Scale such that the maximum bin has value of one.",
      UI_GROUP_HIST_CFP,
      &p->maxbinone,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "manualbinrange",
      UI_KEY_MANUALBINRANGE,
      0,
      0,
      "Set min/max of bins manually, not from data.",
      UI_GROUP_HIST_CFP,
      &p->manualbinrange,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "onebinstart",
      UI_KEY_ONEBINSTART,
      "FLT",
      0,
      "Shift bins so one bin starts on this value.",
      UI_GROUP_HIST_CFP,
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
