/*********************************************************************
Interpolate - Fill blank values in a dataset
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
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/dimension.h>





/* Interpolate blank values in an array. */
void
gal_interpolate(gal_data_t *input)
{
  size_t F;
  uint8_t *flagarr;
  gal_data_t *flag;
  size_t *dinc=gal_dimension_increment(input->ndim, input->dsize);

  /* Set a value of 1 for every element that must be interpolated. */
  flag=gal_blank_flag(input);

  /* Find the proper value for each neighbor. */
  flagarr=flag->array;
  for(F=0; F<input->size; ++F)
    if(flagarr[F])
      {
        printf("to be filled: %zu\n", F);
        GAL_DIMENSION_NEIGHBOR_OP(F, input->ndim, input->dsize, 2,
                                  dinc, {printf("\tneighbor: %zu\n", nind);});
      }

  /* Clean up. */
  gal_data_free(flag);
}
