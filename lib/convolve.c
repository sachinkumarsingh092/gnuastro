/*********************************************************************
convolve -- Convolve a dataset with a given kernel.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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

#include <gnuastro/tile.h>



/* Do spatial convolution on each mesh. */
static void *
gal_convolve_spatial_on_thread(void *inparm)
{
  return NULL;
}





/* Convolve a dataset with a given kernel in the spatial domain. Since
   spatial convolution can be very series of tiles arranged as an array. */
void
gal_convolve_spatial(gal_data_t *tiles, gal_data_t *kernel)
{

}
