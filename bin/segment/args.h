/*********************************************************************
Segment - Segment initial labels based on signal structure.
Segment is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
    /* Input options. */
    {
      "sky",
      UI_KEY_SKY,
      "STR/FLT",
      0,
      "Filename of Sky values image to subtract.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->skyname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "skyhdu",
      UI_KEY_SKYHDU,
      "STR",
      0,
      "HDU containing Sky value to subtract.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->skyhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "std",
      UI_KEY_STD,
      "STR/FLT",
      0,
      "Filename of Sky standard deviation.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->stdname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "stdhdu",
      UI_KEY_STDHDU,
      "STR",
      0,
      "HDU containing Sky standard deviation.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->stdhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "variance",
      UI_KEY_VARIANCE,
      0,
      0,
      "STD input is actually variance.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->variance,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "detection",
      UI_KEY_DETECTION,
      "STR",
      0,
      "Filename of detection image (to segment).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->detectionname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "dhdu",
      UI_KEY_DHDU,
      "STR",
      0,
      "HDU containing detection image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->dhdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "kernel",
      UI_KEY_KERNEL,
      "STR",
      0,
      "Filename of kernel to convolve with input.",
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





    /* Output. */
    {
      "rawoutput",
      UI_KEY_RAWOUTPUT,
      0,
      0,
      "Output only object and clump labels.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->rawoutput,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "onlyclumps",
      UI_KEY_ONLYCLUMPS,
      0,
      0,
      "Finish after finding true clumps.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->onlyclumps,
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
      UI_GROUP_SEGMENTATION,
      &p->grownclumps,
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




    /* Segmentation options. */
    {
      0, 0, 0, 0,
      "Segmentation:",
      UI_GROUP_SEGMENTATION
    },
    {
      "minskyfrac",
      UI_KEY_MINSKYFRAC,
      "FLT",
      0,
      "Min. fraction of undetected area in tile.",
      UI_GROUP_SEGMENTATION,
      &p->minskyfrac,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minima",
      UI_KEY_MINIMA,
      0,
      0,
      "Built internal clumps from minima.",
      UI_GROUP_SEGMENTATION,
      &p->minima,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "snminarea",
      UI_KEY_SNMINAREA,
      "INT",
      0,
      "Minimum area of clumps for S/N estimation.",
      UI_GROUP_SEGMENTATION,
      &p->snminarea,
      GAL_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checksn",
      UI_KEY_CHECKSN,
      0,
      0,
      "Save clump S/N values into a file.",
      UI_GROUP_SEGMENTATION,
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
      UI_GROUP_SEGMENTATION,
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
      "S/N Quantile of true sky clumps.",
      UI_GROUP_SEGMENTATION,
      &p->snquant,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "keepmaxnearriver",
      UI_KEY_KEEPMAXNEARRIVER,
      0,
      0,
      "Keep clumps with peak touching a river.",
      UI_GROUP_SEGMENTATION,
      &p->keepmaxnearriver,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "clumpsnthresh",
      UI_KEY_CLUMPSNTHRESH,
      "FLT",
      0,
      "S/N threshold of true clumps.",
      UI_GROUP_SEGMENTATION,
      &p->clumpsnthresh,
      GAL_TYPE_FLOAT32,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "gthresh",
      UI_KEY_GTHRESH,
      "FLT",
      0,
      "Multiple of STD to stop growing clumps.",
      UI_GROUP_SEGMENTATION,
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
      UI_GROUP_SEGMENTATION,
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
      UI_GROUP_SEGMENTATION,
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
      UI_GROUP_SEGMENTATION,
      &p->checksegmentation,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    /* Operating mode options. */





    {0}
  };





/* Define the child argp structure
   -------------------------------

   NOTE: these parts can be left untouched.*/
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
