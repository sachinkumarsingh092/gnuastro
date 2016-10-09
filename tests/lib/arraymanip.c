/*********************************************************************
A test program to test array functions in Gnuastro.

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
#include <stdio.h>
#include <stdlib.h>

#include "gnuastro/array.h"
#include "gnuastro/gnuastro.h"

int
main(void)
{
  double *array;
  size_t i, size=6;

  /* Allocate the array, report an error if it wasn't allocated. */
  array=malloc(size * sizeof *array);
  if(array==NULL)
    {
      fprintf(stderr, "%zu bytes for d.\n", size);
      exit(EXIT_FAILURE);
    }

  /* Print the version of Gnuastro being used: */
  printf("Test of Gnuastro %s\n", GAL_GNUASTRO_VERSION);

  /* Fill in the test array and report its contents at the same time. */
  printf("Input array: ");
  for(i=0;i<size;++i)
    {
      array[i]=i;
      printf("%.3f, ", array[i]);
    }

  /* delete the last two characters, add a . and newline */
  printf("\b\b.\n");

  /* Now do some simple operations, and report each */
  printf("Multipyling by 2.0...\n");
  gal_array_dmultip_const(array, size, 2.0f);

  printf("Dividing by 3.0...\n");
  gal_array_ddivide_const(array, size, 3.0f);

  printf("Taking natural logarithm...\n");
  gal_array_dlog_array(array, size);

  /* Report the final array: */
  printf("Output array: ");
  for(i=0;i<size;++i)
    printf("%.3f, ", array[i]);
  printf("\b\b.\n");

  /* Cleanup and return */
  free(array);
  return EXIT_SUCCESS;
}
