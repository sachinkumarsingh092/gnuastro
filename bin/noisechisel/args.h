/*********************************************************************
NoiseChisel - Detect and segment signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
      "wkhdu",
      UI_KEY_WKHDU,
      "STR",
      0,
      "HDU containing wide kernel image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->wkhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "skysubtracted",
      UI_KEY_SKYSUBTRACTED,
      0,
      0,
      "Input is Sky subtracted (for error estimation).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->skysubtracted,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minskyfrac",
      UI_KEY_MINSKYFRAC,
      "FLT",
      0,
      "Min. fraction of undetected area in tile.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->minskyfrac,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minnumfalse",
      UI_KEY_MINNUMFALSE,
      "INT",
      0,
      "Minimum number for S/N estimation.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->minnumfalse,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
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
      "onlydetection",
      UI_KEY_ONLYDETECTION,
      0,
      0,
      "Stop at the end of detection.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->onlydetection,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "grownclumps",
      UI_KEY_GROWNCLUMPS,
      0,
      0,
      "Save grown clumps instead of original.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->grownclumps,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    /* Detection. */
    {
      0, 0, 0, 0,
      "Detection:",
      ARGS_GROUP_DETECTION
    },
    {
      "mirrordist",
      UI_KEY_MIRRORDIST,
      "FLT",
      0,
      "Max. dist. (error multip.) to find mode.",
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
      &p->sigmaclip,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_read_sigma_clip
    },
    {
      "checkdetsky",
      UI_KEY_CHECKDETSKY,
      0,
      0,
      "Save Sky value estimation for pseudo-dets.",
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
      &p->dthresh,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "detsnminarea",
      UI_KEY_DETSNMINAREA,
      "INT",
      0,
      "Min. pseudo-detection area for S/N dist.",
      ARGS_GROUP_DETECTION,
      &p->detsnminarea,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checkdetsn",
      UI_KEY_CHECKDETSN,
      0,
      0,
      "Save pseudo-detection S/N values to a file.",
      ARGS_GROUP_DETECTION,
      &p->checkdetsn,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "detquant",
      UI_KEY_DETQUANT,
      "FLT",
      0,
      "Quantile in pseudo-det. to define true.",
      ARGS_GROUP_DETECTION,
      &p->detquant,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0_LT_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "detgrowquant",
      UI_KEY_DETGROWQUANT,
      "FLT",
      0,
      "Minimum quant. to expand true detections.",
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
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
      ARGS_GROUP_DETECTION,
      &p->checksky,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Segmentation */
    {
      0, 0, 0, 0,
      "Segmentation:",
      ARGS_GROUP_SEGMENTATION
    },
    {
      "segsnminarea",
      UI_KEY_SEGSNMINAREA,
      "INT",
      0,
      "Minimum area of clumps for S/N estimation.",
      ARGS_GROUP_SEGMENTATION,
      &p->segsnminarea,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checkclumpsn",
      UI_KEY_CHECKCLUMPSN,
      0,
      0,
      "Save Sky clump S/N values into a file.",
      ARGS_GROUP_SEGMENTATION,
      &p->checkclumpsn,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "segquant",
      UI_KEY_SEGQUANT,
      "FLT",
      0,
      "S/N Quantile of true sky clumps.",
      ARGS_GROUP_SEGMENTATION,
      &p->segquant,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "keepmaxnearriver",
      UI_KEY_KEEPMAXNEARRIVER,
      0,
      0,
      "Keep clumps with peak touching a river.",
      ARGS_GROUP_SEGMENTATION,
      &p->keepmaxnearriver,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "gthresh",
      UI_KEY_GTHRESH,
      "FLT",
      0,
      "Multiple of STD to stop growing clumps.",
      ARGS_GROUP_SEGMENTATION,
      &p->gthresh,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minriverlength",
      UI_KEY_MINRIVERLENGTH,
      "INT",
      0,
      "Minimum len of useful grown clump rivers.",
      ARGS_GROUP_SEGMENTATION,
      &p->minriverlength,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "objbordersn",
      UI_KEY_OBJBORDERSN,
      "FLT",
      0,
      "Min. S/N for grown clumps as one object.",
      ARGS_GROUP_SEGMENTATION,
      &p->objbordersn,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checksegmentation",
      UI_KEY_CHECKSEGMENTATION,
      0,
      0,
      "Store segmentation steps in a file.",
      ARGS_GROUP_SEGMENTATION,
      &p->checksegmentation,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },


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
