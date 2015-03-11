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
LDADD = gnulib/lib/libgnu.a
AM_CPPFLAGS = -DSYSCONFIG_DIR=\"$(sysconfdir)\" -I$(top_srcdir)/include	\
-I$(top_builddir)/gnulib/lib -I$(top_srcdir)/gnulib/lib


## Make libraries of the files that are shared between all the
## programs:
noinst_LTLIBRARIES = libconfigfiles.la libfitsarrayvv.la	\
libtxtarrayvv.la libcheckset.la liblinkedlist.la libtiming.la	\
libastrthreads.la libbox.la libarraymanip.la libstatistics.la	\
libspatialconvolve.la

libbox_la_SOURCES = box.c $(top_srcdir)/include/box.h
libbox_la_LIBADD = -lm

libtiming_la_SOURCES = timing.c $(top_srcdir)/include/timing.h

libcheckset_la_SOURCES = checkset.c $(top_srcdir)/include/checkset.h

liblinkedlist_la_SOURCES = linkedlist.c $(top_srcdir)/include/linkedlist.h

libstatistics_la_SOURCES = statistics.c $(top_srcdir)/include/statistics.h

libarraymanip_la_SOURCES = arraymanip.c $(top_srcdir)/include/arraymanip.h

libtxtarrayvv_la_SOURCES = txtarrayvv.c $(top_srcdir)/include/txtarrayvv.h

libastrthreads_la_SOURCES = astrthreads.c $(top_srcdir)/include/astrthreads.h
libastrthreads_la_LDFLAGS = $(LIB_PTHREAD)

libconfigfiles_la_SOURCES = configfiles.c $(top_srcdir)/include/configfiles.h

libfitsarrayvv_la_SOURCES = fitsarrayvv.c $(top_srcdir)/include/fitsarrayvv.h
libfitsarrayvv_la_LIBADD = -lwcs -lcfitsio -lpthread -lm



# Dependent libraries:
libspatialconvolve_la_SOURCES = spatialconvolve.c	\
$(top_srcdir)/include/spatialconvolve.h
libspatialconvolve_la_LDFLAGS = -lastrthreads
