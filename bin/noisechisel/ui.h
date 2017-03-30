/*********************************************************************
NoiseChisel - Detect and segment signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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

   a b f j l n u w x z
   A H J L W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_KERNEL             = 'k',
  ARGS_OPTION_KEY_SKYSUBTRACTED      = 'E',
  ARGS_OPTION_KEY_MINBFRAC           = 'B',
  ARGS_OPTION_KEY_MIRRORDIST         = 'r',
  ARGS_OPTION_KEY_MODMEDQDIFF        = 'Q',
  ARGS_OPTION_KEY_QTHRESH            = 't',
  ARGS_OPTION_KEY_ERODE              = 'e',
  ARGS_OPTION_KEY_OPENING            = 'p',
  ARGS_OPTION_KEY_SIGMACLIP          = 's',
  ARGS_OPTION_KEY_DTHRESH            = 'R',
  ARGS_OPTION_KEY_DETSNMINAREA       = 'i',
  ARGS_OPTION_KEY_DETQUANT           = 'c',
  ARGS_OPTION_KEY_DILATE             = 'd',
  ARGS_OPTION_KEY_SEGSNMINAREA       = 'm',
  ARGS_OPTION_KEY_SEGQUANT           = 'g',
  ARGS_OPTION_KEY_KEEPMAXNEARRIVER   = 'v',
  ARGS_OPTION_KEY_GTHRESH            = 'G',
  ARGS_OPTION_KEY_MINRIVERLENGTH     = 'y',
  ARGS_OPTION_KEY_OBJBORDERSN        = 'O',
  ARGS_OPTION_KEY_CONTINUEAFTERCHECK = 'C',


  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  ARGS_OPTION_KEY_KHDU              = 1000,
  ARGS_OPTION_KEY_MINNUMFALSE,
  ARGS_OPTION_KEY_ONLYDETECT,
  ARGS_OPTION_KEY_GROWNCLUMPS,
  ARGS_OPTION_KEY_SMOOTHWIDTH,
  ARGS_OPTION_KEY_CHECKQTHRESH,
  ARGS_OPTION_KEY_ERODENGB,
  ARGS_OPTION_KEY_NOERODEQUANT,
  ARGS_OPTION_KEY_OPENINGNGB,
  ARGS_OPTION_KEY_CHECKDETSKY,
  ARGS_OPTION_KEY_CHECKDETSN,
  ARGS_OPTION_KEY_CHECKDETECTION,
  ARGS_OPTION_KEY_CHECKSKY,
  ARGS_OPTION_KEY_CHECKCLUMPSN,
  ARGS_OPTION_KEY_CHECKSEGMENTATION,
};





void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct noisechiselparams *p);

void
ui_abort_after_check(struct noisechiselparams *p, char *filename,
                     char *description);

void
ui_free_report(struct noisechiselparams *p, struct timeval *t1);

#endif
