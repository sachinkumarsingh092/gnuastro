/*********************************************************************
Match - A program to match catalogs and WCS warps
Match is part of GNU Astronomy Utilities (Gnuastro) package.

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



/* Group codes. */
enum program_args_groups
{
  UI_GROUP_CATALOGMATCH = GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Available letters for short options:

   b e f g i j k m n p r s t u v w x y z
   A B E G J L O Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_HDU2            = 'H',
  UI_KEY_APERTURE        = 'a',
  UI_KEY_LOGASOUTPUT     = 'l',
  UI_KEY_CCOL1           = 'c',
  UI_KEY_CCOL2           = 'C',
  UI_KEY_COORD           = 'd',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_NOTMATCHED      = 1000,
  UI_KEY_OUTCOLS,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct matchparams *p);

void
ui_free_report(struct matchparams *p, struct timeval *t1);

#endif
