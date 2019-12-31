/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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





/* Macros. */
#define UI_NO_CONV_KERNEL_NAME "none"





/* Option groups particular to this program. */
enum program_args_groups
{
  UI_GROUP_DETECTION = GAL_OPTIONS_GROUP_AFTER_COMMON,
  UI_GROUP_SEGMENTATION,
};





/* Available letters for short options:

   a b f g i j n r u v x y z
   A E G H J O W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_LARGETILESIZE      = 'L',
  UI_KEY_KERNEL             = 'k',
  UI_KEY_WIDEKERNEL         = 'w',
  UI_KEY_MINSKYFRAC         = 'B',
  UI_KEY_MEANMEDQDIFF       = 'Q',
  UI_KEY_QTHRESH            = 't',
  UI_KEY_ERODE              = 'e',
  UI_KEY_OPENING            = 'p',
  UI_KEY_SIGMACLIP          = 's',
  UI_KEY_DTHRESH            = 'R',
  UI_KEY_SNMINAREA          = 'm',
  UI_KEY_SNQUANT            = 'c',
  UI_KEY_DETGROWQUANT       = 'd',
  UI_KEY_CONTINUEAFTERCHECK = 'C',
  UI_KEY_LABEL              = 'l',


  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_KHDU               = 1000,
  UI_KEY_CONVOLVED,
  UI_KEY_CHDU,
  UI_KEY_WHDU,
  UI_KEY_MINNUMFALSE,
  UI_KEY_SMOOTHWIDTH,
  UI_KEY_QTHRESHTILEQUANT,
  UI_KEY_OUTLIERNUM,
  UI_KEY_OUTLIERSIGMA,
  UI_KEY_OUTLIERSCLIP,
  UI_KEY_CHECKQTHRESH,
  UI_KEY_BLANKASFOREGROUND,
  UI_KEY_ERODENGB,
  UI_KEY_NOERODEQUANT,
  UI_KEY_OPENINGNGB,
  UI_KEY_SKYFRACNOBLANK,
  UI_KEY_CHECKDETSKY,
  UI_KEY_DOPENING,
  UI_KEY_DOPENINGNGB,
  UI_KEY_HOLENGB,
  UI_KEY_PSEUDOCONCOMP,
  UI_KEY_CHECKSN,
  UI_KEY_SNTHRESH,
  UI_KEY_DETGROWMAXHOLESIZE,
  UI_KEY_CLEANGROWNDET,
  UI_KEY_CHECKDETECTION,
  UI_KEY_CHECKSKY,
  UI_KEY_RAWOUTPUT,
  UI_KEY_IGNOREBLANKINTILES,
};





void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct noisechiselparams *p);

void
ui_abort_after_check(struct noisechiselparams *p, char *filename,
                     char *file2name, char *description);

void
ui_free_report(struct noisechiselparams *p, struct timeval *t1);

#endif
