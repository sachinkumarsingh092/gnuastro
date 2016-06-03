/*********************************************************************
Functions to that only use WCSLIB's functions, not related to FITS.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef __GAL_WCS_H__
#define __GAL_WCS_H__

#include <wcslib/wcs.h>

void
gal_wcs_xy_array_to_radec(struct wcsprm *wcs, double *xy, double *radec,
                          size_t number, size_t width);

void
gal_wcs_radec_array_to_xy(struct wcsprm *wcs, double *radec, double *xy,
                          size_t number, size_t width);

double
gal_wcs_pixel_area_arcsec2(struct wcsprm *wcs);

#endif
