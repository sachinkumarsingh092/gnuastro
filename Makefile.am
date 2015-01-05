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
if COND_MOCKGALS
  MAYBE_MOCKGALS = src/mockgals
endif
if COND_IMGCROP
  MAYBE_IMGCROP = src/imgcrop
endif
SUBDIRS = lib $(MAYBE_MOCKGALS) $(MAYBE_IMGCROP) doc tests





## Files that are installed (and distributed)
## =========================================
dist_doc_DATA = README





## Files that are only distributed:
## ================================
##
## include/: This directory keeps all the internal library headers
##   and the other header files that all the programs will need.
##
## doc/findversion.sh: This is a small script that will go through all
##   the header files of all the programs to read their version
##   numbers. It will then save them in a subpackageversions.texi file
##   that will be read by astrutils.texi to give the most uptodate
##   versions for each program.
##
## tests/basicchecks.sh: Some checks to determine if the tests should
##   be skipped (the executable was not made) or if the default file
##   should be put here for the use of the program.
##
## tests/imgcrop/makerandomcat.py: A python script which will make a
##   random catalog within the range of the input image to crop.
##
## tests/imgcrop/randcat.txt: The random catalog that was created by
##   makerandomcat.py.
EXTRA_DIST = include/arraymanip.h include/astrthreads.h			\
 include/checkset.h include/commonargs.h include/commonparams.h		\
 include/configfiles.h include/fitsarrayvv.h include/fixedstringmacros.h	\
 include/linkedlist.h include/stats.h include/timing.h			\
 include/txtarrayvv.h doc/findversions.sh				\
 tests/imgcrop/makerandomcat.py tests/imgcrop/randcat.txt		\
 tests/basicchecks.sh
