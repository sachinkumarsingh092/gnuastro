/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

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

   b c d e f g i j k n p s t u w x y
   A B C E G J L O Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_REDSHIFT       = 'z',
  ARGS_OPTION_KEY_H0             = 'H',
  ARGS_OPTION_KEY_OLAMBDA        = 'l',
  ARGS_OPTION_KEY_OMATTER        = 'm',
  ARGS_OPTION_KEY_ORADIATION     = 'r',
  ARGS_OPTION_KEY_ONLYVOLUME     = 'v',
  ARGS_OPTION_KEY_ONLYABSMAGCONV = 'a',


  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};





void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct cosmiccalparams *p);

#endif
