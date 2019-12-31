/*********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef PROFILES_H
#define PROFILES_H

double
profiles_radial_distance(struct mkonthread *mkp);

double
profiles_gaussian_total(double q);

double
profiles_gaussian(struct mkonthread *mkp);

double
profiles_moffat_alpha(double fwhm, double beta);

double
profiles_moffat_total(double alpha, double beta, double q);

double
profiles_moffat(struct mkonthread *mkp);

double
profiles_sersic_b(double n);

double
profiles_sersic_total(double n, double re, double b, double q);

double
profiles_sersic(struct mkonthread *mkp);

double
profiles_circumference(struct mkonthread *mkp);

double
profiles_flat(struct mkonthread *mkp);

#endif
