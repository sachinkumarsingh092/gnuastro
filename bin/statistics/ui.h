/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

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
  UI_GROUP_SINGLE_VALUE = GAL_OPTIONS_GROUP_AFTER_COMMON,
  UI_GROUP_PARTICULAR_STAT,
  UI_GROUP_SKY,
  UI_GROUP_HIST_CFP,
};





/* Available letters for short options:

   a b e f j p v w x z
   B G J L W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_COLUMN       = 'c',
  UI_KEY_REFCOL       = 'r',
  UI_KEY_GREATEREQUAL = 'g',
  UI_KEY_LESSTHAN     = 'l',
  UI_KEY_QRANGE       = 'Q',
  UI_KEY_MEAN         = 'm',
  UI_KEY_STD          = 'd',
  UI_KEY_MEDIAN       = 'E',
  UI_KEY_MODE         = 'O',
  UI_KEY_QUANTILE     = 'u',
  UI_KEY_ASCIIHIST    = 'A',
  UI_KEY_HISTOGRAM    = 'H',
  UI_KEY_CUMULATIVE   = 'C',
  UI_KEY_SIGMACLIP    = 's',
  UI_KEY_NORMALIZE    = 'n',
  UI_KEY_ONTILE       = 't',
  UI_KEY_INTERPOLATE  = 'i',
  UI_KEY_SKY          = 'y',
  UI_KEY_KERNEL       = 'k',
  UI_KEY_CONTOUR      = 'R',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_NUMBER       = 1000,
  UI_KEY_MINIMUM,
  UI_KEY_MAXIMUM,
  UI_KEY_SUM,
  UI_KEY_MODEQUANT,
  UI_KEY_MODESYM,
  UI_KEY_MODESYMVALUE,
  UI_KEY_QUANTFUNC,
  UI_KEY_ASCIICFP,
  UI_KEY_MIRROR,
  UI_KEY_NUMBINS,
  UI_KEY_NUMASCIIBINS,
  UI_KEY_ASCIIHEIGHT,
  UI_KEY_LOWERBIN,
  UI_KEY_MANUALBINRANGE,
  UI_KEY_ONEBINSTART,
  UI_KEY_MAXBINONE,
  UI_KEY_KHDU,
  UI_KEY_MIRRORDIST,
  UI_KEY_MEANMEDQDIFF,
  UI_KEY_OUTLIERNUM,
  UI_KEY_OUTLIERSIGMA,
  UI_KEY_OUTLIERSCLIP,
  UI_KEY_SMOOTHWIDTH,
  UI_KEY_CHECKSKY,
  UI_KEY_IGNOREBLANKINTILES,
  UI_KEY_SCLIPPARAMS,
  UI_KEY_SIGCLIPNUMBER,
  UI_KEY_SIGCLIPMEDIAN,
  UI_KEY_SIGCLIPMEAN,
  UI_KEY_SIGCLIPSTD,
};





/* Functions */
void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct statisticsparams *p);

void
ui_free_report(struct statisticsparams *p);

#endif
