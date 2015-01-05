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

TESTS = imgcrop/xcyc.sh imgcrop/xcycnoblank.sh imgcrop/section.sh	\
imgcrop/imgcat.sh imgcrop/radec.sh imgcrop/wcscat.sh

EXTRA_DIST= $(TESTS)

CLEANFILES = *.log *.fits .astrutils/*.conf .astrutils
