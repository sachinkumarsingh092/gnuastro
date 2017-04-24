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


## Pre-processer flags (for Gnulib's headers). Recall that the compiled
## Gnulib library was statically linked to (copied in) Gnuastro's library.
AM_CPPFLAGS = -I\$(top_srcdir)/bootstrapped/lib



## Program definition (name, linking, sources and headers)
bin_PROGRAMS = astfits

astfits_LDADD = -lgnuastro

astfits_SOURCES = main.c ui.c extension.c fits.c keywords.c

EXTRA_DIST = main.h authors-cite.h args.h ui.h extension.c fits.h keywords.h



## The configuration file (distribute and install).
## NOTE: the man page is created in doc/Makefile.am
dist_sysconf_DATA = astfits.conf
