/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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





/* Option groups particular to this program. */
enum program_args_groups
{
  UI_GROUP_CENTER_GENERAL = GAL_OPTIONS_GROUP_AFTER_COMMON,
  UI_GROUP_CENTER_SINGLE,
  UI_GROUP_CENTER_CATALOG,
  UI_GROUP_REGION,
};





/* Available letters for short options:

   a d e f g i j k m r t u v y
   A B E G H J L Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_CATALOG        = 'C',
  UI_KEY_NOBLANK        = 'b',
  UI_KEY_SUFFIX         = 'p',
  UI_KEY_NAMECOL        = 'n',
  UI_KEY_SECTION        = 's',
  UI_KEY_POLYGON        = 'l',
  UI_KEY_ZEROISNOTBLANK = 'z',
  UI_KEY_MODE           = 'O',
  UI_KEY_WIDTH          = 'w',
  UI_KEY_CENTER         = 'c',
  UI_KEY_COORDCOL       = 'x',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_CATHDU         = 1000,
  UI_KEY_HSTARTWCS,
  UI_KEY_HENDWCS,
  UI_KEY_POLYGONOUT,
  UI_KEY_POLYGONSORT,
  UI_KEY_CHECKCENTER,
};








void
ui_read_check_inputs_setup(int argc, char *argv[], struct cropparams *p);

void
ui_free_report(struct cropparams *p, struct timeval *t1);

#endif
