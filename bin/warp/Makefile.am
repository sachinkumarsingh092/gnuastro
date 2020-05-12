## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <mohammad@akhlaghi.org>
## Contributing author(s):
## Copyright (C) 2015-2020, Free Software Foundation, Inc.
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


## Necessary pre-processer and linker flags.
AM_LDFLAGS = -L\$(top_builddir)/lib
AM_CPPFLAGS = -I\$(top_srcdir)/bootstrapped/lib -I\$(top_srcdir)/lib

if COND_NORPATH
  MAYBE_NORPATH = $(CONFIG_LDADD)
endif


## Program definition (name, linking, sources and headers)
bin_PROGRAMS = astwarp

## Reason for linking with 'libgnu' described in 'bin/TEMPLATE/Makefile.am'.
astwarp_LDADD = $(top_builddir)/bootstrapped/lib/libgnu.la -lgnuastro \
                $(MAYBE_NORPATH)

astwarp_SOURCES = main.c ui.c warp.c

EXTRA_DIST = main.h authors-cite.h args.h ui.h warp.h



## The configuration file (distribute and install).
## NOTE: the man page is created in doc/Makefile.am
dist_sysconf_DATA = astwarp.conf
