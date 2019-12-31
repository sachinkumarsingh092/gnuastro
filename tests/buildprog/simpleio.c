/*********************************************************************
A test program to multithreaded building using Gnuastro's helpers.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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

#include "gnuastro/fits.h"

int
main(int argc, char *argv[])
{
  gal_data_t *image;
  char *outname="simpleio.fits";

  /* We need two arguments, note that the system also gives the executable
     name, so argc is one more than the number of arguments. */
  if(argc!=3)
    {
      fprintf(stderr, "this program only accepts two arguments");
      return EXIT_FAILURE;
    }

  /* Read the image into memory. */
  image=gal_fits_img_read(argv[1], argv[2], -1, 1);

  /* Let the user know. */
  printf("%s (hdu %s) is read into memory.\n", argv[1], argv[2]);

  /* Save the image in memory into another file. */
  gal_fits_img_write(image, outname, NULL, "BuildProgram's Simpleio");

  /* Let the user know. */
  printf("%s created.\n", outname);

  /* Clean up and return. */
  gal_data_free(image);
  return EXIT_SUCCESS;
}
