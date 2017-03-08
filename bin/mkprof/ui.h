/*********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

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





/* Keys for each option.

   Available letters (-V which is used by GNU is also removed):

   a b d g j l n u v
   A E G H J L O Q W
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_BACKGROUND      = 'k',
  ARGS_OPTION_KEY_BACKHDU         = 'B',
  ARGS_OPTION_KEY_NAXIS1          = 'x',
  ARGS_OPTION_KEY_NAXIS2          = 'y',
  ARGS_OPTION_KEY_CLEARCANVAS     = 'C',
  ARGS_OPTION_KEY_OVERSAMPLE      = 's',
  ARGS_OPTION_KEY_INDIVIDUAL      = 'i',
  ARGS_OPTION_KEY_NOMERGED        = 'm',
  ARGS_OPTION_KEY_NUMRANDOM       = 'r',
  ARGS_OPTION_KEY_TOLERANCE       = 't',
  ARGS_OPTION_KEY_TUNITINP        = 'p',
  ARGS_OPTION_KEY_XSHIFT          = 'X',
  ARGS_OPTION_KEY_YSHIFT          = 'Y',
  ARGS_OPTION_KEY_PREPFORCONV     = 'c',
  ARGS_OPTION_KEY_ZEROPOINT       = 'z',
  ARGS_OPTION_KEY_CIRCUMWIDTH     = 'w',
  ARGS_OPTION_KEY_REPLACE         = 'R',
  ARGS_OPTION_KEY_ENVSEED         = 'e',
  ARGS_OPTION_KEY_MFORFLATPIX     = 'f',

  /* Only with long version. */
  ARGS_OPTION_KEY_PSFINIMG        = 1000,
  ARGS_OPTION_KEY_MAGATPEAK,
  ARGS_OPTION_KEY_XCOL,
  ARGS_OPTION_KEY_YCOL,
  ARGS_OPTION_KEY_RACOL,
  ARGS_OPTION_KEY_DECCOL,
  ARGS_OPTION_KEY_FCOL,
  ARGS_OPTION_KEY_RCOL,
  ARGS_OPTION_KEY_NCOL,
  ARGS_OPTION_KEY_PCOL,
  ARGS_OPTION_KEY_QCOL,
  ARGS_OPTION_KEY_MCOL,
  ARGS_OPTION_KEY_TCOL,
  ARGS_OPTION_KEY_CRPIX1,
  ARGS_OPTION_KEY_CRPIX2,
  ARGS_OPTION_KEY_CRVAL1,
  ARGS_OPTION_KEY_CRVAL2,
  ARGS_OPTION_KEY_RESOLUTION,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkprofparams *p);

void
ui_free_report(struct mkprofparams *p, struct timeval *t1);

#endif
