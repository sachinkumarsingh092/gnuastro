/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef CLUMPS_H
#define CLUMPS_H

/* Important sizes and values (do not change). */
#define SEGMENTNOOBJ     0
#define SEGMENTMASKED   -4
#define SEGMENTTMPCHECK -3
#define SEGMENTINIT     -2
#define SEGMENTRIVER    -1
#define INFOTABCOLS      4
#define WNGBSIZE         9	/* It is impossible to get this large! */

void
clumpsngrid(struct noisechiselparams *p);

#endif
