## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
##
## NOTE:
## This file will be called by the Makefile.am in the top directory,
## so this file does not need to be called directly by automake to
## make the full GNU Astronomy Utilities (Astrutils)
##
## Copyright (C) 2013-2014 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/

AM_TESTS_ENVIRONMENT = export topsrc=$(top_srcdir);

TESTS = mkprof/mosaic1.sh mkprof/mosaic2.sh mkprof/mosaic3.sh		\
mkprof/mosaic4.sh imgcrop/imgcat.sh imgcrop/wcscat.sh imgcrop/xcyc.sh	\
imgcrop/xcycnoblank.sh imgcrop/section.sh imgcrop/radec.sh

EXTRA_DIST = $(TESTS) basicchecks.sh mkprof/mkprofcat1.txt		\
mkprof/mkprofcat2.txt mkprof/mkprofcat3.txt mkprof/mkprofcat4.txt	\
imgcrop/cat.txt

CLEANFILES = *.log *.fits .astrutils/*.conf .astrutils
