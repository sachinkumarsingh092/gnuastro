/*********************************************************************
txt -- Functions for I/O on text files.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/data.h>


/* Read the specified columns in a text file (named `filename' into an
   array of data. The desired columns should be specified by the `cols'
   array. Like a character string, the last element of this array is unique
   and will mark its end. The unique value is `-1' (which is the largest
   possible number in `size_t', since size_t is unsigned.

   When the first element of the cols array is `-1', this function will
   read all the columns from the text file and put the number of columns in
   the space pointed to by `cols' (thus replacing the `-1' value). Note
   that the number of rows is stored within each data element.*/
gal_data_t *
gal_txt_read_cols(char *filename, size_t *cols)
{
  gal_data_t *out=NULL;

  return out;
}
