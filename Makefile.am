## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
## This is part of GNU Astronomy Utilities (gnuastro).
##
## Copyright (C) 2013-2015 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/
##
## gnuastro is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gnuastro is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with gnuastro. If not, see <http://www.gnu.org/licenses/>.

## Directories to check:
## =====================
if COND_MKPROF
  MAYBE_MKPROF = src/mkprof
endif
if COND_IMGCROP
  MAYBE_IMGCROP = src/imgcrop
endif
SUBDIRS = lib $(MAYBE_MKPROF) $(MAYBE_IMGCROP) doc tests





## Files that are installed (and distributed)
## =========================================
dist_doc_DATA = README





## Files that are only distributed:
## ================================
EXTRA_DIST = include/commonargs.h include/commonparams.h	\
 include/fixedstringmacros.h include/neighbors.h
