## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <akhlaghi@gnu.org>
## Contributing author(s):
## Copyright (C) 2015, Free Software Foundation, Inc.
##
## Gnuastro is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Gnuastro is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.




# Flags:
AM_CFLAGS = -Wall -O3
AM_LDFLAGS = -lcfitsio -lwcs -lgsl -lgslcblas -lm	\
-L$(top_builddir)/lib
AM_CPPFLAGS = -DSYSCONFIG_DIR=\"$(sysconfdir)\" -I$(top_srcdir)/include	\
-I$(top_builddir)/gnulib/lib -I$(top_srcdir)/gnulib/lib





## ImageCrop:
bin_PROGRAMS = astmkcatalog

astmkcatalog_SOURCES = main.c main.h cite.h ui.c ui.h args.h	\
mkcatalog.c mkcatalog.h columns.c columns.h

astmkcatalog_LDADD = $(top_builddir)/gnulib/lib/libgnu.la -lconfigfiles	\
-lfitsarrayvv -lcheckset -ltiming -llinkedlist -ltxtarrayvv





# To destribute the defaults file.
dist_sysconf_DATA = astmkcatalog.conf





## Generally for the man pages:
dist_man_MANS = astmkcatalog.1

astmkcatalog.1: args.h $(top_srcdir)/include/commonargs.h
	-$(HELP2MAN) --output=$@ --source="$(PACKAGE_STRING)" ./astmkcatalog
