/*********************************************************************
TEMPLATE - A minimal set of files and functions to define a program.
TEMPLATE is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Your Name <your@email>
Contributing author(s):
Copyright (C) YYYY, Free Software Foundation, Inc.

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

   a b c d e f g i j k l n p r s t u v w x y z
   A B C E G H J L Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_MULTIVALUE      = 'm',
  ARGS_OPTION_KEY_ONOFF           = 'O',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct TEMPLATEparams *p);

void
ui_free_report(struct TEMPLATEparams *p, struct timeval *t1);

#endif
