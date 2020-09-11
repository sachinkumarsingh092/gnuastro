/*********************************************************************
Query - Retreive data from a remote data server.
Query is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2020, Free Software Foundation, Inc.

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
  UI_GROUP_BYCENTER = GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Available letters for short options:

   a b e f g i j k m n p t u v x y z
   A B E G H J L R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_DATABASE        = 'd',
  UI_KEY_QUERY           = 'Q',
  UI_KEY_DATASET         = 's',
  UI_KEY_CENTER          = 'c',
  UI_KEY_RADIUS          = 'r',
  UI_KEY_COLUMN          = 'C',
  UI_KEY_WIDTH           = 'w',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct queryparams *p);

void
ui_free_report(struct queryparams *p, struct timeval *t1);

#endif
