/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

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

   a b c d e f i j k p v w x y z
   B E F G J L R W X Y Z
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_COLUMN       = 'c',
  ARGS_OPTION_KEY_REFCOL       = 'r',
  ARGS_OPTION_KEY_GREATEREQUAL = 'g',
  ARGS_OPTION_KEY_LESSTHAN     = 'l',
  ARGS_OPTION_KEY_QRANGE       = 'Q',
  ARGS_OPTION_KEY_MEAN         = 'm',
  ARGS_OPTION_KEY_STD          = 't',
  ARGS_OPTION_KEY_MEDIAN       = 'M',
  ARGS_OPTION_KEY_MODE         = 'O',
  ARGS_OPTION_KEY_QUANTILE     = 'u',
  ARGS_OPTION_KEY_ASCIIHIST    = 'A',
  ARGS_OPTION_KEY_HISTOGRAM    = 'H',
  ARGS_OPTION_KEY_CUMULATIVE   = 'C',
  ARGS_OPTION_KEY_SIGMACLIP    = 's',
  ARGS_OPTION_KEY_NORMALIZE    = 'n',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  ARGS_OPTION_KEY_NUMBER       = 1000,
  ARGS_OPTION_KEY_MINIMUM,
  ARGS_OPTION_KEY_MAXIMUM,
  ARGS_OPTION_KEY_SUM,
  ARGS_OPTION_KEY_MODEQUANT,
  ARGS_OPTION_KEY_MODESYM,
  ARGS_OPTION_KEY_MODESYMVALUE,
  ARGS_OPTION_KEY_QUANTFUNC,
  ARGS_OPTION_KEY_ASCIICFP,
  ARGS_OPTION_KEY_MIRROR,
  ARGS_OPTION_KEY_NUMBINS,
  ARGS_OPTION_KEY_NUMASCIIBINS,
  ARGS_OPTION_KEY_ASCIIHEIGHT,
  ARGS_OPTION_KEY_LOWERBIN,
  ARGS_OPTION_KEY_ONEBINSTART,
  ARGS_OPTION_KEY_MAXBINONE,
};





/* Functions */
void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct statisticsparams *p);

void
ui_free_report(struct statisticsparams *p);

#endif
