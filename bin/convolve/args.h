/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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
    /* Inputs */
    {
      "kernel",
      UI_KEY_KERNEL,
      "STR",
      0,
      "File name of kernel for convolution.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->kernelname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
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
      "khdu",
      UI_KEY_KHDU,
      "STR",
      0,
      "HDU containing the kernel.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->khdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "kernelcolumn",
      UI_KEY_KERNELCOLUMN,
      "STR",
      0,
      "Column name or number if kernel is a table.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->kernelcolumn,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "nokernelflip",
      UI_KEY_NOKERNELFLIP,
      0,
      0,
      "Do not flip the kernel image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->nokernelflip,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "nokernelnorm",
      UI_KEY_NOKERNELNORM,
      0,
      0,
      "Do not normalize the kernel image.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->nokernelnorm,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "minsharpspec",
      UI_KEY_MINSHARPSPEC,
      "FLT",
      0,
      "Deconvolution: min spectrum of sharp img.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->minsharpspec,
      GAL_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    /* Outputs */
    {
      "checkfreqsteps",
      UI_KEY_CHECKFREQSTEPS,
      0,
      0,
      "View the steps in the frequency domain.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->checkfreqsteps,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "noedgecorrection",
      UI_KEY_NOEDGECORRECTION,
      0,
      0,
      "Do not correct the edges in the spatial domain",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->noedgecorrection,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    /* Operating mode. */
    {
      "domain",
      UI_KEY_DOMAIN,
      "STR",
      0,
      "Convolution domain: 'spatial', 'frequency'.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->domainstr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "makekernel",
      UI_KEY_MAKEKERNEL,
      "INT",
      0,
      "Make 2*INT kernel to create input image.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->makekernel,
      GAL_TYPE_SIZE_T,
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
