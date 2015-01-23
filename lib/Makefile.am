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

# Flags for all programs:
AM_CFLAGS = -Wall -Ofast

AM_CPPFLAGS = -DSYSDATA_DIR=\"$(sysconfdir)\" -I$(top_srcdir)/include


## Make libraries of the files that are shared between all the
## programs:
noinst_LIBRARIES = libconfigfiles.a libfitsarrayvv.a libtxtarrayvv.a	\
libcheckset.a liblinkedlist.a libtiming.a libastrthreads.a libbox.a	\
libarraymanip.a libstatistics.a

libbox_a_SOURCES = box.c $(top_srcdir)/include/box.h

libtiming_a_SOURCES = timing.c $(top_srcdir)/include/timing.h

libcheckset_a_SOURCES = checkset.c $(top_srcdir)/include/checkset.h

liblinkedlist_a_SOURCES = linkedlist.c $(top_srcdir)/include/linkedlist.h

libstatistics_a_SOURCES = statistics.c $(top_srcdir)/include/statistics.h

libarraymanip_a_SOURCES = arraymanip.c $(top_srcdir)/include/arraymanip.h

libtxtarrayvv_a_SOURCES = txtarrayvv.c $(top_srcdir)/include/txtarrayvv.h

libastrthreads_a_SOURCES = astrthreads.c $(top_srcdir)/include/astrthreads.h

libconfigfiles_a_SOURCES = configfiles.c $(top_srcdir)/include/configfiles.h

libfitsarrayvv_a_SOURCES = fitsarrayvv.c $(top_srcdir)/include/fitsarrayvv.h
