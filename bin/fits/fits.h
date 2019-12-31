/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef FITS_H
#define FITS_H



enum fits_action_ids
  {
    FITS_ACTION_INVALID,        /* ==0: by C standard. */

    FITS_ACTION_DELETE,
    FITS_ACTION_RENAME,
    FITS_ACTION_UPDATE,
    FITS_ACTION_WRITE,

    FITS_ACTION_COPY,
    FITS_ACTION_REMOVE,
  };

int
fits_has_error(struct fitsparams *p, int actioncode, char *string,
               int status);

int
fits(struct fitsparams *p);

#endif
