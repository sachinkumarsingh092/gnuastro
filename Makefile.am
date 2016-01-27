## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
## This is part of GNU Astronomy Utilities (gnuastro).
##
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
# Copyright (C) 2015, Free Software Foundation, Inc.
#
# Gnuastro is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Gnuastro is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.

## Directories to check:
## =====================
##
## Gnulib is intentionally compiled and tested prior to the gnuastro
## progams. With this configuration, the user can confirm if the
## problem is internal to gnulib or still persists after all the
## Gnulib tests have been done.
if COND_CONVERTT
  MAYBE_CONVERTT = src/convertt
endif
if COND_CONVOLVE
  MAYBE_CONVOLVE = src/convolve
endif
if COND_HEADER
  MAYBE_HEADER = src/header
endif
if COND_IMGARITH
  MAYBE_IMGARITH = src/imgarith
endif
if COND_IMGCROP
  MAYBE_IMGCROP = src/imgcrop
endif
if COND_IMGSTAT
  MAYBE_IMGSTAT = src/imgstat
endif
if COND_IMGWARP
  MAYBE_IMGWARP = src/imgwarp
endif
if COND_MKCATALOG
  MAYBE_MKCATALOG = src/mkcatalog
endif
if COND_MKNOISE
  MAYBE_MKNOISE = src/mknoise
endif
if COND_MKPROF
  MAYBE_MKPROF = src/mkprof
endif
if COND_NOISECHISEL
  MAYBE_NOISECHISEL = src/noisechisel
endif
if COND_SUBTRACTSKY
  MAYBE_SUBTRACTSKY = src/subtractsky
endif
if COND_GNULIBCHECK
  MAYBE_GNULIBCHECK = gnulib/tests
endif
SUBDIRS = gnulib/lib $(MAYBE_GNULIBCHECK) lib $(MAYBE_CONVERTT)		\
$(MAYBE_CONVOLVE) $(MAYBE_HEADER) $(MAYBE_IMGARITH) $(MAYBE_IMGCROP)	\
$(MAYBE_IMGSTAT) $(MAYBE_IMGWARP) $(MAYBE_MKCATALOG) $(MAYBE_MKNOISE)	\
$(MAYBE_MKPROF) $(MAYBE_NOISECHISEL) $(MAYBE_SUBTRACTSKY) doc tests


ACLOCAL_AMFLAGS = -I m4


## Files that are installed (and distributed)
## =========================================
dist_doc_DATA = README


## Files that are only distributed:
## ================================
EXTRA_DIST = m4/gnulib-cache.m4 include/commonargs.h	\
 include/commonparams.h include/fixedstringmacros.h	\
 include/neighbors.h
