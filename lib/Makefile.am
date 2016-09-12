## Process this file with automake to produce Makefile.in
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


## Necessary flags. NOTE: $(top_srcdir)/bootstrapped/lib is only necessary
## for internally compiled utilities and libraries. It must not be included
## during the tests since the bootstrapped libraries are not installed.
AM_CPPFLAGS = -I\$(top_srcdir)/bootstrapped/lib


# Define the main compiled library file
lib_LTLIBRARIES = libgnuastro.la


# Specify the library .c files
libgnuastro_la_SOURCES = array.c box.c checkset.c configfiles.c fits.c  \
  gnuastro.c linkedlist.c mesh.c mode.c polygon.c qsort.c               \
  spatialconvolve.c statistics.c threads.c timing.c txtarray.c wcs.c
libgnuastro_la_LDFLAGS = -version-info 0:0:0


# Specify the installed headers, note that we are not blindly including all
# `.h' files in the $(headersdir) directory. Some of the header files don't
# need to be installed.
headersdir=$(top_srcdir)/lib/gnuastro
pkginclude_HEADERS = $(headersdir)/array.h $(headersdir)/box.h             \
  $(headersdir)/checkset.h $(headersdir)/configfiles.h                     \
  $(headersdir)/fits.h $(headersdir)/gnuastro.h $(headersdir)/linkedlist.h \
  $(headersdir)/mesh.h $(headersdir)/mode.h $(headersdir)/polygon.h        \
  $(headersdir)/qsort.h $(headersdir)/spatialconvolve.h			   \
  $(headersdir)/statistics.h $(headersdir)/threads.h                       \
  $(headersdir)/timing.h $(headersdir)/txtarray.h $(headersdir)/wcs.h
