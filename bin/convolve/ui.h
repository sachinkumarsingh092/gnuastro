/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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





/* Option groups particular to this program. */
enum program_args_groups
{
  UI_GROUP_MESH_GRID = GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Available letters for short options:

   a b e f g i j l p s v w x y z
   A B E G J L O Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_KERNEL         = 'k',
  UI_KEY_KHDU           = 'u',
  UI_KEY_MINSHARPSPEC   = 'H',
  UI_KEY_CHECKFREQSTEPS = 'C',
  UI_KEY_TILESIZE       = 't',
  UI_KEY_COLUMN         = 'c',
  UI_KEY_NUMCHANNELS    = 'n',
  UI_KEY_REMAINDERFRAC  = 'r',
  UI_KEY_DOMAIN         = 'd',
  UI_KEY_MAKEKERNEL     = 'm',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_KERNELCOLUMN = 1000,
  UI_KEY_NOKERNELFLIP,
  UI_KEY_NOKERNELNORM,
  UI_KEY_NOEDGECORRECTION,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct convolveparams *p);

void
ui_free_report(struct convolveparams *p, struct timeval *t1);

#endif
