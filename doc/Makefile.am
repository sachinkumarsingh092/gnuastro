## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
##
## NOTE:
## This file will be called by the Makefile.am in the top directory,
## so this file does not need to be called directly by automake to
## make the full GNU Astronomy Utilities (GAU).
##
## This is the automake instructions to build GAU. GAU is a collection
## of utitlies for astronomical data analysis. It is free software
## released under GNU GPLv3+.
##
## Copyright (C) 2013-2014 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/

## Add the version numbers of the separate programs:
## ================================================= List the header
## files in this string with a space character between them, don't
## forget to add `$(top_srcdir}/src/` before them like the ones
## already there.

SHELL=/bin/bash
headers=${top_srcdir}/src/*/main.h

## Commands to make the texinfo tools:
info_TEXINFOS = astrutils.texi
astrutils_TEXINFOS = fdl.texi subpackageversions.texi formath.texi

EXTRA_DIST = findversions.sh

subpackageversions.texi: ${headers}
	./findversions.sh ${headers}
