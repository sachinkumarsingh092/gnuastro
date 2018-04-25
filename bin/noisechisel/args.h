/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2018, Free Software Foundation, Inc.

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






/* Options in argp_option format. */
struct argp_option program_options[] =
  {

    /* Input options. */
    {
      "kernel",
      UI_KEY_KERNEL,
      "STR",
      0,
      "Filename of kernel to convolve with input",
      GAL_OPTIONS_GROUP_INPUT,
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
      "HDU containing kernel image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->khdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "convolved",
      UI_KEY_CONVOLVED,
      "STR",
      0,
      "Convolved image file to avoid convolution.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->convolvedname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "chdu",
      UI_KEY_CHDU,
      "STR",
      0,
      "HDU/extension of convolved image in file.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->chdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "widekernel",
      UI_KEY_WIDEKERNEL,
      "STR",
      0,
      "Filename of wider kernel for better qthresh",
      GAL_OPTIONS_GROUP_INPUT,
      &p->widekernelname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "whdu",
      UI_KEY_WHDU,
      "STR",
      0,
      "HDU containing wide kernel image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->whdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    /* Tessellation. */
    {
      "largetilesize",
      UI_KEY_LARGETILESIZE,
      "INT[,INT]",
      0,
      "Sim. to --tilesize, but for larger tiles.",
      GAL_OPTIONS_GROUP_TESSELLATION,
      &p->ltl.tilesize,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },





    /* Output options. */
    {
      "rawoutput",
      UI_KEY_RAWOUTPUT,
      0,
      0,
      "Output only detection labels & 1-elem/tile grid.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->rawoutput,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "label",
      UI_KEY_LABEL,
      0,
      0,
      "Label/count detected pixels that are connected.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->label,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    /* Detection. */
    {
      0, 0, 0, 0,
      "Detection:",
      UI_GROUP_DETECTION
    },
    {
      "mirrordist",
      UI_KEY_MIRRORDIST,
      "FLT",
      0,
      "Max. dist. (error multip.) to find mode.",
      UI_GROUP_DETECTION,
      &p->mirrordist,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "modmedqdiff",
      UI_KEY_MODMEDQDIFF,
      "FLT",
      0,
      "Max. mode and median quant diff. per tile.",
      UI_GROUP_DETECTION,
      &p->modmedqdiff,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "qthresh",
      UI_KEY_QTHRESH,
      "FLT",
      0,
      "Quantile threshold on convolved image.",
      UI_GROUP_DETECTION,
      &p->qthresh,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LT_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "qthreshtilequant",
      UI_KEY_QTHRESHTILEQUANT,
      "FLT",
      0,
      "Remove tiles at higher quantiles.",
      UI_GROUP_DETECTION,
      &p->qthreshtilequant,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "smoothwidth",
      UI_KEY_SMOOTHWIDTH,
      "INT",
      0,
      "Flat kernel width to smooth interpolated.",
      UI_GROUP_DETECTION,
      &p->smoothwidth,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_0_OR_ODD,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checkqthresh",
      UI_KEY_CHECKQTHRESH,
      0,
      0,
      "Save quantile threshold estimation in file.",
      UI_GROUP_DETECTION,
      &p->checkqthresh,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "erode",
      UI_KEY_ERODE,
      "INT",
      0,
      "Number of erosions after thresholding.",
      UI_GROUP_DETECTION,
      &p->erode,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "erodengb",
      UI_KEY_ERODENGB,
      "INT",
      0,
      "4 or 8 connectivity in erosion.",
      UI_GROUP_DETECTION,
      &p->erodengb,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "noerodequant",
      UI_KEY_NOERODEQUANT,
      "FLT",
      0,
      "Quantile for no erosion.",
      UI_GROUP_DETECTION,
      &p->noerodequant,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "opening",
      UI_KEY_OPENING,
      "INT",
      0,
      "Depth of opening after erosion.",
      UI_GROUP_DETECTION,
      &p->opening,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "openingngb",
      UI_KEY_OPENINGNGB,
      "INT",
      0,
      "4 or 8 connectivity in opening.",
      UI_GROUP_DETECTION,
      &p->openingngb,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "sigmaclip",
      UI_KEY_SIGMACLIP,
      "FLT,FLT",
      0,
      "Sigma multiple and, tolerance or number.",
      UI_GROUP_DETECTION,
      &p->sigmaclip,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_sigma_clip
    },
    {
      "minskyfrac",
      UI_KEY_MINSKYFRAC,
      "FLT",
      0,
      "Min. fraction of undetected area in tile.",
      UI_GROUP_DETECTION,
      &p->minskyfrac,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checkdetsky",
      UI_KEY_CHECKDETSKY,
      0,
      0,
      "Save Sky value estimation for pseudo-dets.",
      UI_GROUP_DETECTION,
      &p->checkdetsky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "dthresh",
      UI_KEY_DTHRESH,
      "FLT",
      0,
      "Sigma threshold for Pseudo-detections.",
      UI_GROUP_DETECTION,
      &p->dthresh,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "snminarea",
      UI_KEY_SNMINAREA,
      "INT",
      0,
      "Min. pseudo-detection area for S/N dist.",
      UI_GROUP_DETECTION,
      &p->snminarea,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checksn",
      UI_KEY_CHECKSN,
      0,
      0,
      "Save pseudo-detection S/N values to a file.",
      UI_GROUP_DETECTION,
      &p->checksn,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minnumfalse",
      UI_KEY_MINNUMFALSE,
      "INT",
      0,
      "Minimum number for S/N estimation.",
      UI_GROUP_DETECTION,
      &p->minnumfalse,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "snquant",
      UI_KEY_SNQUANT,
      "FLT",
      0,
      "Quantile in pseudo-det. to define true.",
      UI_GROUP_DETECTION,
      &p->snquant,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "detgrowquant",
      UI_KEY_DETGROWQUANT,
      "FLT",
      0,
      "Minimum quant. to expand true detections.",
      UI_GROUP_DETECTION,
      &p->detgrowquant,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "detgrowmaxholesize",
      UI_KEY_DETGROWMAXHOLESIZE,
      "INT",
      0,
      "Max. area of holes after growth to fill.",
      UI_GROUP_DETECTION,
      &p->detgrowmaxholesize,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "cleangrowndet",
      UI_KEY_CLEANGROWNDET,
      0,
      0,
      "Remove small S/N grown detections.",
      UI_GROUP_DETECTION,
      &p->cleangrowndet,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checkdetection",
      UI_KEY_CHECKDETECTION,
      0,
      0,
      "Save all the detection steps to a file.",
      UI_GROUP_DETECTION,
      &p->checkdetection,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checksky",
      UI_KEY_CHECKSKY,
      0,
      0,
      "Final sky and its STD steps in a file.",
      UI_GROUP_DETECTION,
      &p->checksky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    /* Operating mode options. */
    {
      "continueaftercheck",
      UI_KEY_CONTINUEAFTERCHECK,
      0,
      0,
      "Continue processing after checks.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->continueaftercheck,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    /* End of options. */
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
