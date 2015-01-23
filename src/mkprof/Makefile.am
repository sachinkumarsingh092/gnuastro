## Process this file with automake to produce Makefile.inx
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


# Flags:
AM_CFLAGS = -Wall -Ofast

AM_CPPFLAGS = -DSYSCONFIG_DIR=\"$(sysconfdir)/\" -I$(top_srcdir)/include

AM_LDFLAGS = -lcfitsio -lwcs -lgsl -lgslcblas -pthread -lm	\
-L$(top_builddir)/lib





## ImageCrop, .astrutils is for testing
bin_PROGRAMS = astrmkprof

astrmkprof_SOURCES = main.c main.h ui.c ui.h mkprof.c mkprof.h args.h	\
oneprofile.c oneprofile.h profiles.c profiles.h cite.h

astrmkprof_LDADD = -lconfigfiles -lfitsarrayvv -ltxtarrayvv		\
-lcheckset -llinkedlist -ltiming -lastrthreads -lbox -lstatistics	\
-larraymanip



# To destribute the defaults file.
dist_sysconf_DATA = astrmkprof.conf





## Generally for the man pages:
dist_man_MANS = astrmkprof.1

astrmkprof.1: args.h $(top_srcdir)/include/commonargs.h
	-$(HELP2MAN) --output=$@ --source="$(PACKAGE_STRING)" ./astrmkprof
