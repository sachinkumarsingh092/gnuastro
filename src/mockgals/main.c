/*********************************************************************
MockGals - Create mock galaxies and stars in a noisy image.
MockGals is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
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
#include <stdio.h>
#include <stdlib.h>

#include "timing.h"   	        /* Includes time.h and sys/time.h */

#include "main.h"
#include "mockgals.h"

#include "ui.h"		        /* needs main.h.                  */

int
main (int argc, char *argv[])
{
  struct timeval t1;
  struct mockgalsparams p={{0}, {0}, 0};

  /* Set the starting time.*/
  time(&p.rawtime);
  gettimeofday(&t1, NULL);

  /* Read the input parameters. */
  setparams(argc, argv, &p);

  /* Run Image Crop */
  mockgals(&p);

  /* Free all non-freed allocations. */
  freeandreport(&p, &t1);

  /* Return 0 for success.*/
  return 0;
}
