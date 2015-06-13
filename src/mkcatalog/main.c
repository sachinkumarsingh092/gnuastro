/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <progname.h>

#include "timing.h"   	        /* Includes time.h and sys/time.h */

#include "main.h"

#include "ui.h"		        /* needs main.h.                  */
#include "mkcatalog.h"          /* needs main.h.                  */

int
main (int argc, char *argv[])
{
  struct timeval t1;
  struct mkcatalogparams p={{0}, {0}, 0};

  /* Set the starting time. */
  time(&p.rawtime);
  gettimeofday(&t1, NULL);

  /* Set the program name (needed by non-gnu operating systems): */
  set_program_name(argv[0]);

  /* Read the input parameters. */
  setparams(argc, argv, &p);

  /* Run MakeProfiles */
  mkcatalog(&p);

  /* Free all non-freed allocations. */
  freeandreport(&p, &t1);

  /* Return successfully.*/
  return EXIT_SUCCESS;
}
