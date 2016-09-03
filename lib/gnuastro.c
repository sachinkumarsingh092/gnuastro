/*********************************************************************
Functions dealing with general aspects of all Gnuastro.

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

#include <string.h>

#include "gnuastro/checkset.h"
#include "gnuastro/gnuastro.h"


/* Return the version of Gnuastro. If the input argument is pointing to a
   NULL pointer, then allocate the necessary space and copy version into
   it. Otherwise, assume that the user has allocated the necessary space
   (either statically or dynamically) and just copy the PACKAGE_VERSION
   macro into the given pointer. The returned value is the pointer that the
   arument is pointing to, so both can be used.*/
char *
gal_gnuastro_version(void)
{
  char *version=NULL;
  gal_checkset_allocate_copy(PACKAGE_VERSION, &version);
  return version;
}
