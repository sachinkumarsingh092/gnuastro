/*********************************************************************
txtarrayvv -- Convert a text file table to a C array.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef TXTARRAYVV_H
#define TXTARRAYVV_H

#include <float.h>






/* Simple macros: */
#define FMTLENGTH     50
#define DELIMITERS    " ,\t\n"
#define TXTARRAYVVLOG "txtarrayvv.log"





/* Read a text table file into a double C array. Some tables might
   have non-number elements, if so, they are put into
   replacements. Replacements is an array of strings (pointers to
   characters) similar to argv, since we need to allocate it inside
   this function you have to give its pointer, hence three
   dereferences. */
void
txttoarray(char *filename, double **array, size_t *s0, size_t *s1);

void
arraytotxt(double *array, size_t s0, size_t s1, char *comments,
	   int *int_cols, int *accu_cols, int *space, int *prec,
	   const char *filename);

#endif
