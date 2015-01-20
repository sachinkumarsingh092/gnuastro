## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
##
## This is the automake instructions to build GAU. GAU is a collection
## of utitlies for astronomical data analysis. It is free software
## released under GNU GPLv3+.
##
## Copyright (C) 2013-2014 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/

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
