/*********************************************************************
Segment - Segment initial labels based on signal structure.
Segment is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
  UI_GROUP_SEGMENTATION = GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Available letters for short options:

   a b e f g i j l n p r t u w x z
   A E H J Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_KERNEL             = 'k',
  UI_KEY_DETECTION          = 'd',
  UI_KEY_LARGETILESIZE      = 'L',
  UI_KEY_MINSKYFRAC         = 'B',
  UI_KEY_SNMINAREA          = 'm',
  UI_KEY_SNQUANT            = 'c',
  UI_KEY_KEEPMAXNEARRIVER   = 'v',
  UI_KEY_CLUMPSNTHRESH      = 's',
  UI_KEY_GTHRESH            = 'G',
  UI_KEY_MINRIVERLENGTH     = 'y',
  UI_KEY_OBJBORDERSN        = 'O',
  UI_KEY_CONTINUEAFTERCHECK = 'C',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_KHDU               = 1000,
  UI_KEY_CONVOLVED,
  UI_KEY_CHDU,
  UI_KEY_DHDU,
  UI_KEY_SKY,
  UI_KEY_SKYHDU,
  UI_KEY_STD,
  UI_KEY_STDHDU,
  UI_KEY_VARIANCE,
  UI_KEY_MINIMA,
  UI_KEY_RAWOUTPUT,
  UI_KEY_MINNUMFALSE,
  UI_KEY_ONLYCLUMPS,
  UI_KEY_GROWNCLUMPS,
  UI_KEY_CHECKSN,
  UI_KEY_CHECKSEGMENTATION,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct segmentparams *p);

void
ui_abort_after_check(struct segmentparams *p, char *filename, char *file2name,
                     char *description);

void
ui_free_report(struct segmentparams *p, struct timeval *t1);

#endif
