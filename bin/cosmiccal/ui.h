/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
  UI_GROUP_SPECIFIC = GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Available letters for short options:

   f j k n p t w x y
   B E J Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_REDSHIFT            = 'z',
  UI_KEY_OBSLINE             = 'O',
  UI_KEY_H0                  = 'H',
  UI_KEY_OLAMBDA             = 'l',
  UI_KEY_OMATTER             = 'm',
  UI_KEY_ORADIATION          = 'r',

  UI_KEY_USEDREDSHIFT        = 'e',
  UI_KEY_AGENOW              = 'G',
  UI_KEY_CRITICALDENSITYNOW  = 'C',
  UI_KEY_PROPERDISTANCE      = 'd',
  UI_KEY_ANGULARDIMDIST      = 'A',
  UI_KEY_ARCSECTANDIST       = 's',
  UI_KEY_LUMINOSITYDIST      = 'L',
  UI_KEY_DISTANCEMODULUS     = 'u',
  UI_KEY_ABSMAGCONV          = 'a',
  UI_KEY_AGE                 = 'g',
  UI_KEY_LOOKBACKTIME        = 'b',
  UI_KEY_CRITICALDENSITY     = 'c',
  UI_KEY_VOLUME              = 'v',
  UI_KEY_LINEATZ             = 'i',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_LISTLINES           = 1000,
};





void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct cosmiccalparams *p);

#endif
