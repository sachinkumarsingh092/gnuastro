/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

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
  UI_GROUP_OUTROWS = GAL_OPTIONS_GROUP_AFTER_COMMON,
};




/* Available letters for short options:

   a b d f g j k l m p t v x y z
   A B E G H J L O Q R X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_WCSFILE         = 'w',
  UI_KEY_WCSHDU          = 'W',
  UI_KEY_COLUMN          = 'c',
  UI_KEY_INFORMATION     = 'i',
  UI_KEY_COLINFOINSTDOUT = 'O',
  UI_KEY_RANGE           = 'r',
  UI_KEY_EQUAL           = 'e',
  UI_KEY_NOTEQUAL        = 'n',
  UI_KEY_SORT            = 's',
  UI_KEY_DESCENDING      = 'd',
  UI_KEY_HEAD            = 'H',
  UI_KEY_TAIL            = 't',
  UI_KEY_CATCOLUMN       = 'C',
  UI_KEY_CATCOLHDU       = 'u',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_POLYGON         = 1000,
  UI_KEY_INPOLYGON,
  UI_KEY_OUTPOLYGON,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct tableparams *p);

void
ui_list_select_free(struct list_select *list, int freevalue);

void
ui_free_report(struct tableparams *p);

#endif
