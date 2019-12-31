/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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





/* Option groups particular to this program */
enum program_args_groups
{
  UI_GROUP_FLUX = GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Available letters for short options:

   a d e f g j k l n p r s t v y z
   E G J O Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_QUALITY             = 'u',
  UI_KEY_WIDTHINCM           = 'w',
  UI_KEY_BORDERWIDTH         = 'b',
  UI_KEY_HEX                 = 'x',
  UI_KEY_FLUXLOW             = 'L',
  UI_KEY_FLUXHIGH            = 'H',
  UI_KEY_MAXBYTE             = 'm',
  UI_KEY_FORCEMIN            = 'A',
  UI_KEY_FORCEMAX            = 'B',
  UI_KEY_CHANGE              = 'c',
  UI_KEY_CHANGEAFTERTRUNC    = 'C',
  UI_KEY_INVERT              = 'i',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_COLORMAP            = 1000,
  UI_KEY_RGBTOHSV,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct converttparams *p);

void
ui_free_report(struct converttparams *p);

#endif
