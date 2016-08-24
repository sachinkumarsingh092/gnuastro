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
bin_PROGRAMS = asttable

asttable_SOURCES = main.c main.h cite.h ui.c ui.h args.h	\
table.c table.h

asttable_LDADD = $(top_builddir)/bootstrapped/lib/libgnu.la             \
-lgalconfigfiles -lgalfits -lgalcheckset -lgaltiming -lgallinkedlist    \
-lgaltxtarray





# To destribute the defaults file.
# NOTE: the man page is created in doc/Makefile.am
dist_sysconf_DATA = asttable.conf
