/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>









/*************************************************************
 **************      Acceptable JPEG names      ***************
 *************************************************************/
int
nameisjpeg(char *name)
{
  size_t len;
  len=strlen(name);
  if (strcmp(&name[len-4], ".jpg") == 0
      || strcmp(&name[len-4], ".JPG") == 0
      || strcmp(&name[len-4], ".jpeg") == 0
      || strcmp(&name[len-5], ".JPEG") == 0
      || strcmp(&name[len-5], ".jpe") == 0
      || strcmp(&name[len-5], ".jif") == 0
      || strcmp(&name[len-5], ".jfif") == 0
      || strcmp(&name[len-5], ".jfi") == 0)
    return 1;
  else
    return 0;
}
