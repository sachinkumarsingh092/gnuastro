/*********************************************************************
Warp - Warp images using projective mapping.
Warp is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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

   b d g i j l n u v w x y z
   A B E G H J L O Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_KEEPWCS         = 'k',
  ARGS_OPTION_KEY_COVEREDFRAC     = 'C',
  ARGS_OPTION_KEY_ALIGN           = 'a',
  ARGS_OPTION_KEY_ROTATE          = 'r',
  ARGS_OPTION_KEY_SCALE           = 's',
  ARGS_OPTION_KEY_FLIP            = 'f',
  ARGS_OPTION_KEY_SHEAR           = 'e',
  ARGS_OPTION_KEY_TRANSLATE       = 't',
  ARGS_OPTION_KEY_PROJECT         = 'p',
  ARGS_OPTION_KEY_MATRIX          = 'm',
  ARGS_OPTION_KEY_CENTERONCORNER  = 'c',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  ARGS_OPTION_KEY_HSTARTWCS       = 1000,
  ARGS_OPTION_KEY_HENDWCS,
};





/* Functions */
void
ui_read_check_inputs_setup(int argc, char *argv[], struct warpparams *p);

void
ui_free_report(struct warpparams *p, struct timeval *t1);

#endif
