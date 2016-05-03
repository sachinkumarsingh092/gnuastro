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


# Set the conditional directory compilation
if COND_ARITHMETIC
  MAYBE_ARITHMETIC = src/arithmetic
endif
if COND_CONVERTT
  MAYBE_CONVERTT = src/convertt
endif
if COND_CONVOLVE
  MAYBE_CONVOLVE = src/convolve
endif
if COND_COSMICCAL
  MAYBE_COSMICCAL = src/cosmiccal
endif
if COND_HEADER
  MAYBE_HEADER = src/header
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
  MAYBE_GNULIBCHECK = bootstrapped/tests
endif
SUBDIRS = bootstrapped/lib $(MAYBE_GNULIBCHECK) lib $(MAYBE_CONVERTT)   \
$(MAYBE_CONVOLVE) $(MAYBE_COSMICCAL) $(MAYBE_HEADER)                    \
$(MAYBE_ARITHMETIC) $(MAYBE_IMGCROP) $(MAYBE_IMGSTAT)                   \
$(MAYBE_IMGWARP) $(MAYBE_MKCATALOG) $(MAYBE_MKNOISE) $(MAYBE_MKPROF)    \
$(MAYBE_NOISECHISEL) $(MAYBE_SUBTRACTSKY) doc tests


## Add m4/ to the list of directories to search for m4 files
## =========================================================
ACLOCAL_AMFLAGS = -I bootstrapped/m4



## Files that are installed (and distributed)
## =========================================
dist_doc_DATA = README



## Files that are only distributed
## ===============================
EXTRA_DIST = bootstrap bootstrap.conf genauthors .dir-locals.el         \
 bootstrapped/m4/gnulib-cache.m4 bootstrapped/README .autom4te.cfg      \
 include/commonargs.h include/commonparams.h                            \
 include/fixedstringmacros.h include/neighbors.h


## Generate AUTHORS
## ================
##
## The AUTHORS file is generated automatically from the version
## control history. Note the following to cases:
##
##   - The AUTHORS file is a prerequisite for the make dist* commands
##     (all make targets including "dist"). So it is not needed when
##     we only run "make".
##
##   - AUTHORS is defined as a PHONY target so that it is called every
##     time any "make dist*" is called, even if AUTHORS already
##     exists. It is possible to build from the tar-ball (that doesn't
##     include version control), so we cannot set any dependancies for
##     it. But generally, "make dist*" is not called to often, so
##     there is no problem.
##
##   - The genauthors script will not do anything (make/update the
##     AUTHORS file) when there is no git repository.
.PHONY: AUTHORS
AUTHORS:; $(top_srcdir)/genauthors
