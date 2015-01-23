## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
##
## This file will be called by the Makefile.am in the top directory,
## so this file does not need to be called directly by automake to
## make the full GNU Astronomy Utilities (Astrutils)
##
## Copyright (C) 2013-2015 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/
##
## AstrUtils is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.
##
## AstrUtils is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.

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
