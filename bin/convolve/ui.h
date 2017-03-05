/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

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





/* Available letters for short options:

   a b e f g i j l n p s v w x y z
   A B E F G J L M O Q R W X Y Z
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_KERNEL         = 'k',
  ARGS_OPTION_KEY_KHDU           = 'u',
  ARGS_OPTION_KEY_MINSHARPSPEC   = 'H',
  ARGS_OPTION_KEY_CHECKFREQSTEPS = 'C',
  ARGS_OPTION_KEY_TILE           = 't',
  ARGS_OPTION_KEY_NUMCHANNELS    = 'c',
  ARGS_OPTION_KEY_REMAINDERFRAC  = 'r',
  ARGS_OPTION_KEY_DOMAIN         = 'd',
  ARGS_OPTION_KEY_MAKEKERNEL     = 'm',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  ARGS_OPTION_KEY_NOKERNELFLIP = 1000,
  ARGS_OPTION_KEY_NOKERNELNORM,
  ARGS_OPTION_KEY_CHECKTILES,
  ARGS_OPTION_KEY_CONVOVERCH,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct convolveparams *p);

void
ui_free_report(struct convolveparams *p, struct timeval *t1);

#endif
