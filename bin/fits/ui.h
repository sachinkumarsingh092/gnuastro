/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

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
  UI_GROUP_EXTENSION = GAL_OPTIONS_GROUP_AFTER_COMMON,
  UI_GROUP_KEYWORD,
};





/* Available letters for short options:

   b e f g i j l m x y z
   A B E G J L O W X Y
 */
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_REMOVE       = 'R',
  UI_KEY_COPY         = 'C',
  UI_KEY_CUT          = 'k',
  UI_KEY_NUMHDUS      = 'n',
  UI_KEY_PRINTALLKEYS = 'p',
  UI_KEY_ASIS         = 'a',
  UI_KEY_DELETE       = 'd',
  UI_KEY_RENAME       = 'r',
  UI_KEY_UPDATE       = 'u',
  UI_KEY_WRITE        = 'w',
  UI_KEY_COMMENT      = 'c',
  UI_KEY_HISTORY      = 'H',
  UI_KEY_DATE         = 't',
  UI_KEY_VERIFY       = 'v',
  UI_KEY_QUITONERROR  = 'Q',
  UI_KEY_DATETOSEC    = 's',


  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_TITLE        = 1000,
  UI_KEY_DATASUM,
  UI_KEY_OUTHDU,
  UI_KEY_COPYKEYS,
  UI_KEY_PRIMARYIMGHDU,
};





void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct fitsparams *p);

void
ui_free_and_report(struct fitsparams *p);

#endif
