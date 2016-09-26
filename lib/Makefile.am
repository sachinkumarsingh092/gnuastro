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


# Define the main compiled library file and its Libtool version.
lib_LTLIBRARIES = libgnuastro.la
libgnuastro_la_LDFLAGS = -version-info $(GAL_LT_VERSION)


# Specify the library .c files
libgnuastro_la_SOURCES = array.c box.c checkset.c configfiles.c fits.c	\
  linkedlist.c mesh.c mode.c polygon.c qsort.c spatialconvolve.c	\
  statistics.c threads.c timing.c txtarray.c wcs.c


# Installed headers, note that we are not blindly including all `.h' files
# in the $(headersdir) directory. Some of the header files don't need to be
# installed.
pkginclude_HEADERS = gnuastro/gnuastro.h $(headersdir)/array.h          \
  $(headersdir)/box.h $(headersdir)/fits.h $(headersdir)/linkedlist.h   \
  $(headersdir)/mesh.h $(headersdir)/polygon.h $(headersdir)/qsort.h    \
  $(headersdir)/spatialconvolve.h $(headersdir)/statistics.h            \
  $(headersdir)/threads.h $(headersdir)/wcs.h $(headersdir)/txtarray.h


# Files to distribute in the tarball. For the headers: these are internal
# headers and don't need to be installed. Headers are only mentioned within
# the source files, and if they are not explicitly mentioned somewhere in
# the Makefile, they will not distributed, so we need to explicitly tell
# Automake to distribute them here.
headersdir=$(top_srcdir)/lib/gnuastro
EXTRA_DIST = gnuastro.pc.in gnuastro.h.in $(headersdir)/commonargs.h      \
  $(headersdir)/commonparams.h $(headersdir)/fixedstringmacros.h          \
  $(headersdir)/mode.h $(headersdir)/neighbors.h $(headersdir)/checkset.h \
  $(headersdir)/configfiles.h $(headersdir)/timing.h $(headersdir)/README


# Definitions for Gnuastro's the pkg-config file (inspired from GSL's
# Makefile.am)
pkgconfig_DATA = gnuastro.pc
pkgconfigdir = $(libdir)/pkgconfig
CLEANFILES = gnuastro.pc gnuastro/gnuastro.h


# Build Gnuastro's pkg-config file and `gnuastro.h' based on the
# information in the Makefile after the Makefile has been built.
gnuastro/gnuastro.h: Makefile gnuastro.h.in
	rm -f $@ $@.tmp
	$(MKDIR_P) gnuastro
	$(SED) -e 's|@VERSION[@]|$(VERSION)|g'                          \
	       -e 's|@HAS_PTHREAD_BARRIER[@]|$(HAS_PTHREAD_BARRIER)|g'  \
               $(srcdir)/gnuastro.h.in >> $@.tmp
	chmod a-w $@.tmp
	mv $@.tmp $@

gnuastro.pc: Makefile $(srcdir)/gnuastro.pc.in
	rm -f $@ $@.tmp
	$(SED)                                      \
	-e 's|@prefix[@]|$(prefix)|g'               \
	-e 's|@exec_prefix[@]|$(exec_prefix)|g'     \
	-e 's|@libdir[@]|$(libdir)|g'               \
	-e 's|@includedir[@]|$(includedir)|g'       \
	-e 's|@LIBS[@]|$(LIBS)|g'                   \
	-e 's|@VERSION[@]|$(VERSION)|g'             \
	'$(srcdir)/$@.in' >> $@.tmp
	chmod a-w $@.tmp
	mv $@.tmp $@
