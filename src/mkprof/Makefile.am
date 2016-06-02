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





# Utility and its sources
bin_PROGRAMS = astmkprof

astmkprof_SOURCES = main.c main.h ui.c ui.h mkprof.c mkprof.h args.h	\
oneprofile.c oneprofile.h profiles.c profiles.h cite.h

astmkprof_LDADD = $(top_builddir)/bootstrapped/lib/libgnu.la		\
-lgalconfigfiles -lgalfitsarrayvv -lgaltxtarrayvv -lgalcheckset		\
-lgallinkedlist -lgaltiming -lgalthreads -lgalbox -lgalstatistics	\
-lgalarraymanip -lgalqsort



# To destribute the defaults file.
dist_sysconf_DATA = astmkprof.conf





## Generally for the man pages:
dist_man_MANS = astmkprof.1

## Note that in parallel builds, we want to make sure that the man page is
## only created once the executable is built, so the argument to
## bin_PROGRAMS is a prerequisite here, the $(EXEEXT) is used in
## Makefile.in.
astmkprof.1: astmkprof$(EXEEXT) args.h                                \
             $(top_srcdir)/lib/gnuastro/commonargs.h
	-$(HELP2MAN) --output=$@ ./astmkprof$(EXEEXT)                 \
	             --source="$(PACKAGE_STRING)"

# With the addition of the executable as a dependency for the man page, the
# (possibly) produced man page will not be cleaned. So we have to ask for
# its cleaning manually here.
CLEANFILES = $(builddir)/astmkprof.1
