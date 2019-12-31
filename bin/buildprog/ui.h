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
#ifndef UI_H
#define UI_H

/* For common options groups. */
#include <gnuastro-internal/options.h>





/* Available letters for short options:

   f i j k n p r s u v w x y z
   A B C E G H J Q R X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_CC             = 'c',
  UI_KEY_INCLUDE        = 'I',
  UI_KEY_LINKDIR        = 'L',
  UI_KEY_LINKLIB        = 'l',
  UI_KEY_ONLYBUILD      = 'b',
  UI_KEY_DEBUG          = 'g',
  UI_KEY_OPTIMIZE       = 'O',
  UI_KEY_WARNING        = 'W',
  UI_KEY_TAG            = 't',
  UI_KEY_DETELECOMPILED = 'd',
  UI_KEY_LA             = 'a',
  UI_KEY_NOENV          = 'e',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct buildprogparams *p);

void
ui_free_report(struct buildprogparams *p);

#endif
