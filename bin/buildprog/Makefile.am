## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <mohammad@akhlaghi.org>
## Contributing author(s):
## Copyright (C) 2017-2020, Free Software Foundation, Inc.
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

if COND_NORPATH
  MAYBE_NORPATH = $(CONFIG_LDADD)
endif


## Program definition (name, linking, sources and headers)
bin_PROGRAMS = astbuildprog

## Reason for linking with 'libgnu' described in 'bin/TEMPLATE/Makefile.am'.
astbuildprog_LDADD = $(top_builddir)/bootstrapped/lib/libgnu.la -lgnuastro \
                     $(MAYBE_NORPATH)

# Basic program sources.
astbuildprog_SOURCES = main.c ui.c buildprog.c

# Extra files that must be distributed in the tarball.
EXTRA_DIST = main.h authors-cite.h args.h ui.h buildprog.h astbuildprog.conf.in



# Build the BuildProgram configuration file. The main problem is the
# possible CPPFLAGS and LDFLAGS that the user might have given at
# Gnuastro's configure time. BuildProgram also needs to know these
# directories in order to compile its programs. With the rule below, the
# directories in these variables will be written into the final
# BuildProgram configuration file.
#
# The user might give non-'-I' arguments to CPPFLAGS and non-'-L' arguments
# to LDFLAGS (see bug #54346). But here, we only want these two options,
# and nothing else. To make the check, we'll see if '$i' (the input string)
# is the same as '$v' or not. If the string starts with '-I' or '-L', it
# will change between the two (the '-I' or '-L' have been removed), so it
# shouldn't be the same.
astbuildprog.conf: $(top_srcdir)/bin/buildprog/astbuildprog.conf.in
	cat $< > $@
	infoadded="no";                                \
	for i in $(CPPFLAGS); do                       \
	  v=$$(echo $$i | sed -e 's/^-I//');           \
	  if test $$i != $$v; then                     \
	    if test $$infoadded = "no"; then           \
	      echo "" >> $@;                           \
	      echo "# Installation information" >> $@; \
	      infoadded="yes";                         \
	    fi;                                        \
	    echo " includedir $$v" >> $@;              \
	  fi;                                          \
	done;                                          \
	for i in $(LDFLAGS); do                        \
	  v=$$(echo $$i | sed -e 's/^-L//');           \
	  if test $$i != $$v; then                     \
	    if test $$infoadded = "no"; then           \
	      echo "" >> $@;                           \
	      echo "# Installation information" >> $@; \
	      infoadded="yes";                         \
	    fi;                                        \
	    echo " linkdir $$v" >> $@;                 \
	  fi;                                          \
	done
	echo "# Build-time LDFLAGS: $(LDFLAGS)" >> $@
	echo "# Build-time CPPFLAGS: $(CPPFLAGS)" >> $@



## The configuration file in BuildProgram must be built on the users's
## system (see above). The built configuration file is needed during 'make
## check'. But the (built) configuration file should NOT be distributed in a
## tarball. Since all built files are distributed, the only way to avoid
## the configuration file being distributed it so use a 'nodist_' prefix.
##
## NOTE: the man page is created in doc/Makefile.am
nodist_sysconf_DATA = astbuildprog.conf
check_DATA = astbuildprog.conf
CLEANFILES = astbuildprog.conf
