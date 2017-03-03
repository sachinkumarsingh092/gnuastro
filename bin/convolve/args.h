/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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
    /* Inputs */
    {
      "kernel",
      ARGS_OPTION_KEY_KERNEL,
      "STR",
      0,
      "File name of kernel for convolution.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->kernelname,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "khdu",
      ARGS_OPTION_KEY_KHDU,
      "STR",
      0,
      "HDU containing the kernel.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->khdu,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "nokernelflip",
      ARGS_OPTION_KEY_NOKERNELFLIP,
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
      ARGS_OPTION_KEY_NOKERNELNORM,
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
      ARGS_OPTION_KEY_MINSHARPSPEC,
      "FLT",
      0,
      "Deconvolution: min spectrum of sharp img.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->minsharpspec,
      GAL_DATA_TYPE_FLOAT64,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Outputs */
    {
      "checkfreqsteps",
      ARGS_OPTION_KEY_CHECKFREQSTEPS,
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





    /* Tile grid options. */
    {
      0, 0, 0, 0,
      "Tile grid (only for spatial domain):",
      ARGS_GROUP_MESH_GRID
    },
    {
      "tile",
      ARGS_OPTION_KEY_TILE,
      "INT[,INT]",
      0,
      "Size of tiles along each dim. (FITS order).",
      ARGS_GROUP_MESH_GRID,
      &p->tile,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },
    {
      "numchannels",
      ARGS_OPTION_KEY_NUMCHANNELS,
      "INT[,..]",
      0,
      "No. of channels along each dim. (FITS order).",
      ARGS_GROUP_MESH_GRID,
      &p->numchannels,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_sizes_reverse
    },
    {
      "convoverch",
      ARGS_OPTION_KEY_CONVOVERCH,
      0,
      0,
      "Convolve over channel borders.",
      ARGS_GROUP_MESH_GRID,
      &p->convoverch,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checktiles",
      ARGS_OPTION_KEY_CHECKTILES,
      0,
      0,
      "Tile IDs in an image, the size of input.",
      ARGS_GROUP_MESH_GRID,
      &p->checktiles,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Operating mode. */
    {
      "domain",
      ARGS_OPTION_KEY_DOMAIN,
      "STR",
      0,
      "Convolution domain: `spatial', `frequency'.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->domainstr,
      GAL_DATA_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "makekernel",
      ARGS_OPTION_KEY_MAKEKERNEL,
      "INT",
      0,
      "Make 2*INT kernel to create input image.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->makekernel,
      GAL_DATA_TYPE_SIZE_T,
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
