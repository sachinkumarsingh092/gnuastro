/*********************************************************************
BuildProgram: Compile and run programs using Gnuastro's library
BuildProgram is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef ARGS_H
#define ARGS_H






/* Array of acceptable options. */
struct argp_option program_options[] =
  {
    {
      "cc",
      UI_KEY_CC,
      "STR",
      0,
      "Name of C compiler's executable.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->cc,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },

    {
      "includedir",
      UI_KEY_INCLUDE,
      "STR",
      0,
      "Directories to search for '#include's.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->include,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },

    {
      "linkdir",
      UI_KEY_LINKDIR,
      "STR",
      0,
      "Directory to search for libraries to link.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->linkdir,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },

    {
      "linklib",
      UI_KEY_LINKLIB,
      "STR",
      0,
      "Link libraries, e.g., for libgsl: '-lgsl'.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->linklib,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },

    {
      "la",
      UI_KEY_LA,
      "STR",
      0,
      "Libtool '.la' to use instead of default.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->la,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },

    {
      "tag",
      UI_KEY_TAG,
      "STR",
      0,
      "Libtool '--tag': programming language.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->tag,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },

    {
      "noenv",
      UI_KEY_NOENV,
      0,
      0,
      "No env. (e.g., LDFLAGS or CPPFLAGS) in build.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->noenv,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },





    {
      "debug",
      UI_KEY_DEBUG,
      0,
      0,
      "Debugging information in compiled binary.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->debug,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "optimize",
      UI_KEY_OPTIMIZE,
      "INT",
      0,
      "Optimization level: 0, 1, 2, 3.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->optimize,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "warning",
      UI_KEY_WARNING,
      "STR",
      0,
      "Compilation warnings on command-line.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->warning,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "onlybuild",
      UI_KEY_ONLYBUILD,
      0,
      0,
      "Don't run the built program.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->onlybuild,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "deletecompiled",
      UI_KEY_DETELECOMPILED,
      0,
      0,
      "Delete compiled program after running.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->deletecompiled,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },

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
