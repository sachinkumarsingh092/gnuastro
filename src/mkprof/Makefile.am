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
## under the terms of the GNU General Public License as published by the
## Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.
##
## gnuastro is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with gnuastro. If not, see <http://www.gnu.org/licenses/>.


# Flags:
AM_CFLAGS = -Wall -Ofast

AM_CPPFLAGS = -DSYSCONFIG_DIR=\"$(sysconfdir)/\" -I$(top_srcdir)/include

AM_LDFLAGS = -lcfitsio -lwcs -lgsl -lgslcblas -pthread -lm	\
-L$(top_builddir)/lib





bin_PROGRAMS = astmkprof

astmkprof_SOURCES = main.c main.h ui.c ui.h mkprof.c mkprof.h args.h	\
oneprofile.c oneprofile.h profiles.c profiles.h cite.h

astmkprof_LDADD = -lconfigfiles -lfitsarrayvv -ltxtarrayvv		\
-lcheckset -llinkedlist -ltiming -lastrthreads -lbox -lstatistics	\
-larraymanip



# To destribute the defaults file.
dist_sysconf_DATA = astmkprof.conf





## Generally for the man pages:
dist_man_MANS = astmkprof.1

astmkprof.1: args.h $(top_srcdir)/include/commonargs.h
	-$(HELP2MAN) --output=$@ --source="$(PACKAGE_STRING)" ./astmkprof
