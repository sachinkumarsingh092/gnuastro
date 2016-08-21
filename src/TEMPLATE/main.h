/*********************************************************************
TEMPLATE - Source code for a blank utility for easy creation of new
           utilities, just copy and paste all the files here replacing
           TEMPLATE with the new utility name (in the source code and
           file names).

         - Add the utility name to `configure.ac' and `Makefile.am' in the
           top Gnuastro source directory.

         - Correct these top comments in all the files, don't forget the
           `astTEMPLATE.conf' and `Makefile.am' files.

TEMPLATE is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Your name <your@email>
Contributing author(s):
Copyright (C) YYYY, Free Software Foundation, Inc.

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

#include <gnuastro/fits.h>
#include <gnuastro/commonparams.h>

/* Progarm name macros: */
#define SPACK           "astTEMPLATE" /* Subpackage executable name. */
#define SPACK_NAME      "TEMPLATE"    /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_NAME") "PACKAGE_VERSION






struct uiparams
{
  char                     *inputname;  /* Name of input file.             */
};





struct TEMPLATEparams
{
  /* Other structures: */
  struct uiparams          up;  /* User interface parameters.         */
  struct gal_commonparams  cp;  /* Common parameters.                 */

  /* Input: */


  /* Internal: */
  time_t              rawtime;  /* Starting time of the program.      */
};

#endif
