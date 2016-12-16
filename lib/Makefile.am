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
libgnuastro_la_SOURCES = array.c box.c checkset.c configfiles.c data.c  \
  data-arithmetic-binary.c data-arithmetic-onlyint.c                    \
  data-arithmetic-other.c data-copy.c fits.c git.c linkedlist.c mesh.c  \
  mode.c polygon.c qsort.c spatialconvolve.c statistics.c table.c       \
  threads.c timing.c txt.c wcs.c





# Installed headers, note that we are not blindly including all `.h' files
# in the $(headersdir) directory. Some of the header files don't need to be
# installed.
headersdir=$(top_srcdir)/lib/gnuastro
pkginclude_HEADERS = gnuastro/config.h $(headersdir)/array.h            \
  $(headersdir)/box.h $(headersdir)/data.h $(headersdir)/fits.h         \
  $(headersdir)/git.h $(headersdir)/linkedlist.h $(headersdir)/mesh.h   \
  $(headersdir)/polygon.h $(headersdir)/qsort.h                         \
  $(headersdir)/spatialconvolve.h $(headersdir)/statistics.h            \
  $(headersdir)/table.h $(headersdir)/threads.h $(headersdir)/wcs.h     \
  $(headersdir)/txt.h





# Files to distribute in the tarball. These are internal headers and don't
# need to be installed. Headers are only mentioned within the source files,
# and if they are not explicitly mentioned somewhere in the Makefile, they
# will not distributed, so we need to explicitly tell Automake to
# distribute them here.
EXTRA_DIST = gnuastro.pc.in data-arithmetic.h data-arithmetic-binary.h	\
  data-arithmetic-onlyint.h data-arithmetic-other.h data-copy.h         \
  config.h.in checkset.h commonargs.h commonparams.h configfiles.h      \
  fixedstringmacros.h mode.h neighbors.h timing.h $(headersdir)/README





# Definitions for Gnuastro's the pkg-config file (inspired from GSL's
# Makefile.am)
pkgconfig_DATA = gnuastro.pc
pkgconfigdir = $(libdir)/pkgconfig
CLEANFILES = gnuastro.pc gnuastro/config.h





# Build Gnuastro's pkg-config file and `config.h' based on the
# information in the Makefile after the Makefile has been built.
gnuastro/config.h: Makefile config.h.in
	rm -f $@ $@.tmp
	$(MKDIR_P) gnuastro
	$(SED) -e 's|@VERSION[@]|$(VERSION)|g'                            \
	       -e 's|@HAVE_LIBGIT2[@]|$(HAVE_LIBGIT2)|g'                  \
	       -e 's|@HAVE_WCSLIB_VERSION[@]|$(HAVE_WCSLIB_VERSION)|g'    \
	       -e 's|@HAVE_PTHREAD_BARRIER[@]|$(HAVE_PTHREAD_BARRIER)|g'  \
	       -e 's|@HAVE_BIN_OP_UCHAR[@]|$(HAVE_BIN_OP_UCHAR)|g'        \
	       -e 's|@HAVE_BIN_OP_CHAR[@]|$(HAVE_BIN_OP_CHAR)|g'          \
	       -e 's|@HAVE_BIN_OP_USHORT[@]|$(HAVE_BIN_OP_USHORT)|g'      \
	       -e 's|@HAVE_BIN_OP_SHORT[@]|$(HAVE_BIN_OP_SHORT)|g'        \
	       -e 's|@HAVE_BIN_OP_UINT[@]|$(HAVE_BIN_OP_UINT)|g'          \
	       -e 's|@HAVE_BIN_OP_INT[@]|$(HAVE_BIN_OP_INT)|g'            \
	       -e 's|@HAVE_BIN_OP_ULONG[@]|$(HAVE_BIN_OP_ULONG)|g'        \
	       -e 's|@HAVE_BIN_OP_LONG[@]|$(HAVE_BIN_OP_LONG)|g'          \
	       -e 's|@HAVE_BIN_OP_LONGLONG[@]|$(HAVE_BIN_OP_LONGLONG)|g'  \
	       -e 's|@HAVE_BIN_OP_FLOAT[@]|$(HAVE_BIN_OP_FLOAT)|g'        \
	       -e 's|@HAVE_BIN_OP_DOUBLE[@]|$(HAVE_BIN_OP_DOUBLE)|g'      \
               $(srcdir)/config.h.in >> $@.tmp
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
