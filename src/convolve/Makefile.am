## Process this file with automake to produce Makefile.inx
##
## This is part of GNU Astronomy Utilities (gnuastro). It will be
## called by the Makefile.am in the top directory.
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

# Flags:
AM_CFLAGS = -Wall -O3
AM_LDFLAGS = -lcfitsio -lwcs -lm -L$(top_builddir)/lib
AM_CPPFLAGS = -DSYSCONFIG_DIR=\"$(sysconfdir)/\" -I$(top_srcdir)/include





## ImageCrop:
bin_PROGRAMS = astconvolve

astconvolve_SOURCES = main.c main.h args.h ui.c ui.h convolve.c	\
convolve.h cite.h

astconvolve_LDADD = -lfitsarrayvv -lcheckset -ltiming -lconfigfiles





# To destribute the defaults file.
dist_sysconf_DATA = astconvolve.conf





## Generally for the man pages:
dist_man_MANS = astconvolve.1

astconvolve.1: args.h $(top_srcdir)/include/commonargs.h
	-$(HELP2MAN) --output=$@ --source="$(PACKAGE_STRING)" ./astconvolve
