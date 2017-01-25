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


## Necessary flags. NOTE: $(top_srcdir)/bootstrapped/lib is only necessary
## for internally compiled utilities and libraries. It must not be included
## during the tests since the bootstrapped libraries are not installed.
AM_CPPFLAGS = -I\$(top_srcdir)/bootstrapped/lib


## Program definition (name, linking, sources and headers)
bin_PROGRAMS = astnoisechisel

astnoisechisel_LDADD = -lgnuastro

astnoisechisel_SOURCES = main.c ui.c noisechisel.c thresh.c binary.c    \
  label.c detection.c sky.c segmentation.c clumps.c

EXTRA_DIST = main.h authors-cite.h args.h ui.h noisechisel.h thresh.h   \
  binary.h label.h detection.h sky.h segmentation.h clumps.h


## The configuration file (distribute and install).
## NOTE: the man page is created in doc/Makefile.am
dist_sysconf_DATA = astnoisechisel.conf
