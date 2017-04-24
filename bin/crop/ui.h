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
  UI_KEY_CATALOG        = 'C',
  UI_KEY_NOBLANK        = 'b',
  UI_KEY_CHECKCENTER    = 'c',
  UI_KEY_SUFFIX         = 'p',
  UI_KEY_NAMECOL        = 'n',
  UI_KEY_RACOL          = 'f',
  UI_KEY_DECCOL         = 'g',
  UI_KEY_RA             = 'r',
  UI_KEY_DEC            = 'd',
  UI_KEY_XCOL           = 'i',
  UI_KEY_YCOL           = 'j',
  UI_KEY_XC             = 'x',
  UI_KEY_YC             = 'y',
  UI_KEY_IWIDTH         = 'a',
  UI_KEY_WWIDTH         = 'w',
  UI_KEY_SECTION        = 's',
  UI_KEY_POLYGON        = 'l',
  UI_KEY_ZEROISNOTBLANK = 'z',
  UI_KEY_MODE           = 'O',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_CATHDU         = 1000,
  UI_KEY_HSTARTWCS,
  UI_KEY_HENDWCS,
  UI_KEY_OUTPOLYGON,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct cropparams *p);

void
ui_free_report(struct cropparams *p, struct timeval *t1);

#endif
