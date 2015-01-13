/*********************************************************************
Functions to facilitate using threads.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
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
#ifndef ASTRTHREADS_H
#define ASTRTHREADS_H

#define NONTHRDINDEX (size_t)(-1)



void
distinthreads(size_t nindexs, size_t nthrds, size_t **outthrds,
	      size_t *outthrdcols);

void
attrbarrierinit(pthread_attr_t *attr, pthread_barrier_t *b,
		size_t numthreads);

#endif
