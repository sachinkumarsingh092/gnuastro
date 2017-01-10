/*********************************************************************
Common parameters between all programs.

IMPORTANT: This header must only be included the programs, not libraries.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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
#ifndef __GAL_COMMONOPTS_H__
#define __GAL_COMMONOPTS_H__

/* Common options for all programs.

   IMPORTANT NOTE:

     This header should only be included in the programs, not the
     libraries, and in particular `options.c'. Because we want each program
     to have its own allocation of the common options structure. If it is
     included in options.c, then it will be shared between all the
     programs.


   Free letters (-V which is used by GNU is also removed):

   a b c d e f g i j k l m n p r s t u v w x y z
   A B C E F G H I J L M O Q R T U W X Y Z

   Used numbers <= 1004

   You can use this above list to set short options for the different
   utilities.
 */
struct argp_option gal_commonopts_options[] =
  {
    /* Input/output. */
#ifndef NOT_COMMON_HDU_PARSER
    {                           /* Some utilities need to parse `hdu' them-*/
      "hdu",                    /* selves. In this case, they will define  */
      'h',                      /* `NOT_COMMON_HDU_PARSER' and set their   */
      "STR",                    /* own `hdu' option structure. */
      0,
      "Extension name or number of input data.",
      1,                        /* Input. */
      NULL, GAL_DATA_TYPE_STRING
    },
#endif
    {
      "output",
      'o',
      "STR",
      0,
      "Output file or directory name.",
      2,                        /* Output. */
      NULL, GAL_DATA_TYPE_STRING
    },
    {
      "dontdelete",
      'D',
      0,
      0,
      "Don't delete output if it exists.",
      2,                        /* Output. */
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "keepinputdir",
      'K',
      0,
      0,
      "Keep input directory for automatic output.",
      2,                        /* Output. */
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },



    /* Operating mode. */
    {
      "quiet",
      'q',
      0,
      0,
      "Only report errors, remain quiet about steps.",
      -1,                       /* Operating mode. */
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "numthreads",
      'N',
      "INT",
      0,
      "Number of CPU threads to use.",
      -1,                       /* Operating mode. */
      NULL, GAL_DATA_TYPE_ULONG
    },
    {
      "minmapsize",
      1004,
      "INT",
      0,
      "Minimum no. bytes to map arrays to hdd/ssd.",
      -1,                       /* Operating mode. */
      NULL, GAL_DATA_TYPE_ULONG
    },
    {
      "log",
      1003,
      0,
      0,
      "No log file for programs which make one.",
      -1,                       /* Operating mode. */
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },



    /* Internal (before control goes back to the program). */
    {
      "cite",
      1000,
      0,
      0,
      "BibTeX citation for this program.",
      -1,                       /* Operating mode. */
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "printparams",
      'P',
      0,
      0,
      "Print parameter values to be used and abort.",
      -1,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "setdirconf",
      'S',
      0,
      0,
      "Set default values for this directory and abort.",
      -1,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "setusrconf",
      'U',
      0,
      0,
      "Set default values for this user and abort.",
      -1,
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "onlydirconf",
      1001,
      0,
      0,
      "Only read current directory configuration file.",
      -1,                       /* Operating mode. */
      NULL, GAL_OPTIONS_NO_ARG_TYPE
    },
    {
      "onlyversion",
      1002,
      "STR",
      0,
      "Only run if the program version is this.",
      -1,                       /* Operating mode. */
      NULL, GAL_DATA_TYPE_STRING
    },



    {0}
  };

#endif
