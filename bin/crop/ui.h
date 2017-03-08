/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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

   e k m t u v
   A B E G H J L Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_CATALOG        = 'C',
  ARGS_OPTION_KEY_NOBLANK        = 'b',
  ARGS_OPTION_KEY_CHECKCENTER    = 'c',
  ARGS_OPTION_KEY_SUFFIX         = 'p',
  ARGS_OPTION_KEY_NAMECOL        = 'n',
  ARGS_OPTION_KEY_RACOL          = 'f',
  ARGS_OPTION_KEY_DECCOL         = 'g',
  ARGS_OPTION_KEY_RA             = 'r',
  ARGS_OPTION_KEY_DEC            = 'd',
  ARGS_OPTION_KEY_XCOL           = 'i',
  ARGS_OPTION_KEY_YCOL           = 'j',
  ARGS_OPTION_KEY_XC             = 'x',
  ARGS_OPTION_KEY_YC             = 'y',
  ARGS_OPTION_KEY_IWIDTH         = 'a',
  ARGS_OPTION_KEY_WWIDTH         = 'w',
  ARGS_OPTION_KEY_SECTION        = 's',
  ARGS_OPTION_KEY_POLYGON        = 'l',
  ARGS_OPTION_KEY_ZEROISNOTBLANK = 'z',
  ARGS_OPTION_KEY_MODE           = 'O',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  ARGS_OPTION_KEY_CATHDU         = 1000,
  ARGS_OPTION_KEY_HSTARTWCS,
  ARGS_OPTION_KEY_HENDWCS,
  ARGS_OPTION_KEY_OUTPOLYGON,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct cropparams *p);

void
ui_free_report(struct cropparams *p, struct timeval *t1);

#endif
