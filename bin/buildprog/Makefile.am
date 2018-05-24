## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <mohammad@akhlaghi.org>
## Contributing author(s):
## Copyright (C) 2017-2018, Free Software Foundation, Inc.
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
##
## Buildprog will also need some system-specific information that is
## gathered at compile time (for example the library installation directory
## (LIBDIR) and the executive file suffix (EXEEXT).
AM_LDFLAGS = -L\$(top_builddir)/lib
AM_CPPFLAGS = -I\$(top_srcdir)/bootstrapped/lib -I\$(top_srcdir)/lib \
              -DLIBDIR=\"$(libdir)\" -DINCLUDEDIR=\"$(includedir)\"      \
              -DEXEEXT=\"$(EXEEXT)\"



## Program definition (name, linking, sources and headers)
bin_PROGRAMS = astbuildprog

## Reason for linking with `libgnu' described in `bin/TEMPLATE/Makefile.am'.
astbuildprog_LDADD = $(top_builddir)/bootstrapped/lib/libgnu.la -lgnuastro

astbuildprog_SOURCES = main.c ui.c buildprog.c astbuildprog.conf

EXTRA_DIST = main.h authors-cite.h args.h ui.h buildprog.h astbuildprog.conf.in



# Build the BuildProgram configuration file. The main problem is the
# possible CPPFLAGS and LDFLAGS that the user might have given at
# Gnuastro's configure time. BuildProgram also needs to know these
# directories in order to compile its programs. With the rule below, the
# directories in these variables will be written into the final
# BuildProgram configuration file.
astbuildprog.conf: $(top_srcdir)/bin/buildprog/astbuildprog.conf.in
	cp $< $@
	infoadded="no";                                \
	for i in $(CPPFLAGS); do                       \
	  if test $$infoadded = "no"; then             \
	    echo "" >> $@;                             \
	    echo "# Installation information" >> $@;   \
	    infoadded="yes";                           \
	  fi;                                          \
	  v=$$(echo $$i | sed -e 's/-I//');            \
	  echo " includedir $$v" >> $@;                \
	done;                                          \
	for i in $(LDFLAGS); do                        \
	  if test $$infoadded = "no"; then             \
	    echo "" >> $@;                             \
	    echo "# Installation information" >> $@;   \
	    infoadded="yes";                           \
	  fi;                                          \
	  v=$$(echo $$i | sed -e 's/-L//');            \
	  echo " linkdir $$v" >> $@;                   \
	done



## The configuration file (clean, distribute and install).
## NOTE: the man page is created in doc/Makefile.am
CLEANFILES = astbuildprog.conf
dist_sysconf_DATA = astbuildprog.conf
