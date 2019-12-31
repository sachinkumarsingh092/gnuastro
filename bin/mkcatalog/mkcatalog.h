/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#ifndef MKCATALOG_H
#define MKCATALOG_H

struct mkcatalog_passparams
{
  struct mkcatalogparams *p;    /* Main MakeCatalog paramers.           */
  double                *oi;    /* Intermediate values for objects.     */
  double                *ci;    /* Intermediate values for clumps.      */
  int32_t            object;    /* Object that is currently working on. */
  size_t        clumpsinobj;    /* The number of clumps in this object. */
  gal_data_t          *tile;    /* The tile to pass-over.               */
  int32_t             *st_o;    /* Starting pointer for object labels.  */
  int32_t             *st_c;    /* Starting pointer for clump labels.   */
  float               *st_v;    /* Starting pointer for values array.   */
  float             *st_sky;    /* Starting pointer for Sky array.      */
  float             *st_std;    /* Starting pointer for Sky STD array.  */
  size_t   start_end_inc[2];    /* Starting and ending indexs.          */
  size_t             *shift;    /* Shift coordinates.                   */
  gsl_rng              *rng;    /* Random number generator.             */
  size_t    clumpstartindex;    /* Clump starting row in final catalog. */
  gal_data_t       *up_vals;    /* Container for upper-limit values.    */
  gal_data_t      *spectrum;    /* Spectrum of each object.             */
};

void
mkcatalog_write_inputs_in_comments(struct mkcatalogparams *p,
                                   gal_list_str_t **comments, int withsky,
                                   int withstd);

void
mkcatalog(struct mkcatalogparams *p);

#endif
