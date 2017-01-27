/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (Gnuastro) package.

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


/* Macros for various types of standard transformation.*/
enum standard_warps
{
  UI_WARP_INVALID,

  UI_WARP_ALIGN,
  UI_WARP_ROTATE,
  UI_WARP_SCALE,
  UI_WARP_FLIP,
  UI_WARP_SHEAR,
  UI_WARP_TRANSLATE,
  UI_WARP_PROJECT,
};


/* Functions */
void
add_to_optionwapsll(struct optionwarpsll **list, int type, char *value);

void
parse_two_values(char *str, double *v1, double *v2);

void
setparams(int argc, char *argv[], struct imgwarpparams *p);

void
freeandreport(struct imgwarpparams *p, struct timeval *t1);

#endif
