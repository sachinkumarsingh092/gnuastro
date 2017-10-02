## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Your Name <your@email>
## Contributing author(s):
## Copyright (C) YYYY, Free Software Foundation, Inc.
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
##
## Buildprog will also need some system-specific information that is
## gathered at compile time (for example the library installation directory
## (LIBDIR) and the executive file suffix (EXEEXT).
AM_CPPFLAGS = -I\$(top_srcdir)/bootstrapped/lib -DLIBDIR=\"$(libdir)\"  \
              -DINCLUDEDIR=\"$(includedir)\" -DEXEEXT=\"$(EXEEXT)\"



## Program definition (name, linking, sources and headers)
bin_PROGRAMS = astbuildprog

astbuildprog_LDADD = -lgnuastro

astbuildprog_SOURCES = main.c ui.c buildprog.c

EXTRA_DIST = main.h authors-cite.h args.h ui.h buildprog.h



## The configuration file (distribute and install).
## NOTE: the man page is created in doc/Makefile.am
dist_sysconf_DATA = astbuildprog.conf
