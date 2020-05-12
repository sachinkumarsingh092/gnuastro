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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <gnuastro/data.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "BuildProgram"    /* Program full name.       */
#define PROGRAM_EXEC   "astbuildprog"    /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







/* Main program parameters structure */
struct buildprogparams
{
  /* From command-line */
  struct gal_options_common_params     cp; /* Common parameters.        */
  gal_list_str_t  *sourceargs;    /* Source files and arguments.        */
  gal_list_str_t     *include;    /* Libraries to link against.         */
  gal_list_str_t     *linkdir;    /* Libraries to link against.         */
  gal_list_str_t     *linklib;    /* Libraries to link against.         */
  char                    *la;    /* Libtool '.la' instead of default.  */
  char                    *cc;    /* C compiler to use.                 */
  uint8_t               noenv;

  char                   *tag;    /* Libtool tag (programming language).*/
  char              *optimize;    /* Optimization level.                */
  char                 *debug;    /* Keep debugging information.        */
  char               *warning;    /* Compiler warnings.                 */
  uint8_t           onlybuild;    /* Don't run the compiled program.    */
  uint8_t      deletecompiled;    /* Delete compiled program after running. */

  /* Output: */
  time_t              rawtime;  /* Starting time of the program.        */
};

#endif
