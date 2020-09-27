/*********************************************************************
Query - Retreive data from a remote data server.
Query is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2020, Free Software Foundation, Inc.

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
#include <string.h>

#include <gnuastro/wcs.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>

#include "gaia.h"
#include "query.h"




void
query_check_download(struct queryparams *p)
{
  size_t len;
  int status=0;
  char *logname;
  fitsfile *fptr;

  /* Open the FITS file and if the status value is still zero, it means
     everything worked properly. */
  fits_open_file(&fptr, p->cp.output, READONLY, &status);
  if(status==0) fits_close_file(fptr, &status);
  else
    {
      /* Add a '.log' suffix to the output filename. */
      len=strlen(p->cp.output);
      logname=gal_pointer_allocate(GAL_TYPE_STRING, len+10, 1,
                                   __func__, "logname");
      sprintf(logname, "%s.log", p->cp.output);

      /* Rename the output file to the logname file and let the user
         know. */
      rename(p->cp.output, logname);
      error(EXIT_FAILURE, 0, "the requested dataset could not be retreived! "
            "For more, please see '%s'", logname);
    }
}





void
query(struct queryparams *p)
{
  /* Download the dataset. */
  switch(p->database)
    {
    case QUERY_DATABASE_GAIA: gaia_query(p); break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to address "
            "the problem. '%d' is not a recognized database code", __func__,
            PACKAGE_BUGREPORT, p->database);
    }

  /* Make sure that the result is a readable FITS file, otherwise, abort
     with an error. */
  query_check_download(p);

  /* Let the user know that things went well. */
  if(p->cp.quiet==0)
    printf("Query output written to: %s\n", p->cp.output);
}
