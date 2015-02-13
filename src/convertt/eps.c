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

#include "main.h"






/*************************************************************
 **************      Acceptable EPS names      ***************
 *************************************************************/
int
nameiseps(char *name)
{
  size_t len;
  len=strlen(name);
  if (strcmp(&name[len-3], "eps") == 0
      || strcmp(&name[len-3], "EPS") == 0
      || strcmp(&name[len-4], "epsf") == 0
      || strcmp(&name[len-4], "epsi") == 0)
    return 1;
  else
    return 0;
}





int
nameisepssuffix(char *name)
{
  if (strcmp(name, "eps") == 0 || strcmp(name, ".eps") == 0
      || strcmp(name, "EPS") == 0 || strcmp(name, ".EPS") == 0
      || strcmp(name, "epsf") == 0 || strcmp(name, ".epsf") == 0
      || strcmp(name, "epsi") == 0 || strcmp(name, ".epsi") == 0)
    return 1;
  else
    return 0;
}




















/*************************************************************
 **************       Write an EPS image        **************
 *************************************************************/
void
saveeps(struct converttparams *p)
{

}
