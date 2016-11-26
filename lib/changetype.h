/*********************************************************************
changetype -- Changing the array of one data type to another.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef __GAL_CHANGETYPE_H__
#define __GAL_CHANGETYPE_H__

void
gal_changetype_out_is_uchar(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_char(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_ushort(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_short(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_uint(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_int(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_ulong(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_long(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_longlong(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_float(gal_data_t *in, gal_data_t *out);

void
gal_changetype_out_is_double(gal_data_t *in, gal_data_t *out);


#endif
