/*********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

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
  UI_GROUP_PROFILES = GAL_OPTIONS_GROUP_AFTER_COMMON,
  UI_GROUP_CATALOG,
  UI_GROUP_WCS,
};





/* Keys for each option.

   Available letters (-V which is used by GNU is also removed):

   a b d g j l n u v y
   A G H J L O Q W Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_BACKGROUND      = 'k',
  UI_KEY_BACKHDU         = 'B',
  UI_KEY_MERGEDSIZE      = 'x',
  UI_KEY_CLEARCANVAS     = 'C',
  UI_KEY_KERNEL          = 'E',
  UI_KEY_OVERSAMPLE      = 's',
  UI_KEY_INDIVIDUAL      = 'i',
  UI_KEY_NOMERGED        = 'm',
  UI_KEY_NUMRANDOM       = 'r',
  UI_KEY_TOLERANCE       = 't',
  UI_KEY_TUNITINP        = 'p',
  UI_KEY_SHIFT           = 'X',
  UI_KEY_PREPFORCONV     = 'c',
  UI_KEY_ZEROPOINT       = 'z',
  UI_KEY_CIRCUMWIDTH     = 'w',
  UI_KEY_REPLACE         = 'R',
  UI_KEY_ENVSEED         = 'e',
  UI_KEY_MFORFLATPIX     = 'f',

  /* Only with long version. */
  UI_KEY_PSFINIMG        = 1000,
  UI_KEY_MAGATPEAK,
  UI_KEY_MCOLISBRIGHTNESS,
  UI_KEY_MODE,
  UI_KEY_CCOL,
  UI_KEY_FCOL,
  UI_KEY_RCOL,
  UI_KEY_NCOL,
  UI_KEY_PCOL,
  UI_KEY_P2COL,
  UI_KEY_P3COL,
  UI_KEY_QCOL,
  UI_KEY_Q2COL,
  UI_KEY_MCOL,
  UI_KEY_TCOL,
  UI_KEY_CRPIX,
  UI_KEY_CRVAL,
  UI_KEY_CDELT,
  UI_KEY_PC,
  UI_KEY_CUNIT,
  UI_KEY_CTYPE,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkprofparams *p);

char *
ui_profile_name_write(int profile_code);

void
ui_free_report(struct mkprofparams *p, struct timeval *t1);

#endif
