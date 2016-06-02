## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <akhlaghi@gnu.org>
## Contributing author(s):
## Copyright (C) 2016, Free Software Foundation, Inc.
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





## Utility and its sources
bin_PROGRAMS = astcosmiccal

astcosmiccal_SOURCES = main.c main.h cite.h ui.c ui.h args.h	\
cosmiccal.c cosmiccal.h

astcosmiccal_LDADD = $(top_builddir)/bootstrapped/lib/libgnu.la	\
-lgalconfigfiles -lgalfitsarrayvv -lgalcheckset -lgaltiming	\
-lgallinkedlist





# To destribute the defaults file.
dist_sysconf_DATA = astcosmiccal.conf





## Generally for the man pages:
dist_man_MANS = astcosmiccal.1

## Note that in parallel builds, we want to make sure that the man page is
## only created once the executable is built, so the argument to
## bin_PROGRAMS is a prerequisite here, the $(EXEEXT) is used in
## Makefile.in.
astcosmiccal.1: astcosmiccal$(EXEEXT) args.h                          \
                $(top_srcdir)/lib/gnuastro/commonargs.h
	-$(HELP2MAN) --output=$@ ./astcosmiccal$(EXEEXT)              \
	             --source="$(PACKAGE_STRING)"

# With the addition of the executable as a dependency for the man page, the
# (possibly) produced man page will not be cleaned. So we have to ask for
# its cleaning manually here.
CLEANFILES = $(builddir)/astcosmiccal.1
