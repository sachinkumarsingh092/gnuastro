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
      GAL_DATA_TYPE_DOUBLE,
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





    /* Mesh grid options. */
    {
      0, 0, 0, 0,
      "Mesh grid (only for spatial domain):",
      ARGS_GROUP_MESH_GRID
    },
    {
      "meshsize",
      ARGS_OPTION_KEY_MESHSIZE,
      "INT",
      0,
      "Size of each mesh (tile) in the grid.",
      ARGS_GROUP_MESH_GRID,
      &p->mp.meshsize,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "nch1",
      ARGS_OPTION_KEY_NCH1,
      "INT",
      0,
      "Number of channels along first FITS axis.",
      ARGS_GROUP_MESH_GRID,
      &p->mp.nch1,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "nch2",
      ARGS_OPTION_KEY_NCH2,
      "INT",
      0,
      "Number of channels along second FITS axis.",
      ARGS_GROUP_MESH_GRID,
      &p->mp.nch2,
      GAL_DATA_TYPE_SIZE_T,
      GAL_OPTIONS_RANGE_GT_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "lastmeshfrac",
      ARGS_OPTION_KEY_LASTMESHFRAC,
      "FLT",
      0,
      "Fraction of last mesh area to add new.",
      ARGS_GROUP_MESH_GRID,
      &p->mp.lastmeshfrac,
      GAL_DATA_TYPE_FLOAT,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "lastmeshfrac",
      ARGS_OPTION_KEY_LASTMESHFRAC,
      "FLT",
      0,
      "Fraction of last mesh area to add new.",
      ARGS_GROUP_MESH_GRID,
      &p->mp.lastmeshfrac,
      GAL_DATA_TYPE_FLOAT,
      GAL_OPTIONS_RANGE_GE_0_LE_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "fullconvolution",
      ARGS_OPTION_KEY_FULLCONVOLUTION,
      0,
      0,
      "Ignore channels in spatial convolution.",
      ARGS_GROUP_MESH_GRID,
      &p->mp.fullconvolution,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "checkmesh",
      ARGS_OPTION_KEY_CHECKMESH,
      0,
      0,
      "File created to view mesh structure",
      ARGS_GROUP_MESH_GRID,
      &p->checkmesh,
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
      0,
      0,
      "Make 2*INT kernel to create input image.",
      GAL_OPTIONS_GROUP_OPERATING_MODE,
      &p->makekernel,
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
