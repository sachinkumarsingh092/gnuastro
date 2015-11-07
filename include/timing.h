/*********************************************************************
Functions to report timing in verbose mode.
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
#ifndef __GAL_TIMING_H__
#define __GAL_TIMING_H__

#include <time.h>
#include <sys/time.h>

#define GAL_TIMING_VERB_MSG_LENGTH_V     45
#define GAL_TIMING_VERB_MSG_LENGTH_T    "45"
#define GAL_TIMING_VERB_MSG_LENGTHS_2_V  65
#define GAL_TIMING_VERB_MSG_LENGTHS_2_T "65"

unsigned long int
gal_timing_time_based_rng_seed();

void
gal_timing_report(struct timeval *t1, char *jobname, size_t level);

#endif
