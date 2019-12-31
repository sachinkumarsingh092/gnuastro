/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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

   a b c d e f i j k l m n p r t u v x y z
   A B C E G H J L Q R X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_GLOBALHDU       = 'g',
  UI_KEY_ONEDASIMAGE     = 'O',
  UI_KEY_ONEDONSTDOUT    = 's',
  UI_KEY_WCSFILE         = 'w',
  UI_KEY_WCSHDU          = 'W',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};


void
ui_read_check_inputs_setup(int argc, char *argv[], struct arithmeticparams *p);

size_t *
ui_read_ndim_dsize(char *filename, char *hdu, size_t *ndim);

void
freeandreport(struct arithmeticparams *p, struct timeval *t1);

#endif
