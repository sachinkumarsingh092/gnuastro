## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
## This is part of GNU Astronomy Utilities (gnuastro).
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





## Sources to be created before everything else
## ============================================
##
## These target(s) will be created before anything else when running
## 'make', 'make check', or 'make install', see the "Built sources" section
## of the AUTOMAKE manual and the git-version-gen script comments for more
## information. Note that AUTHORS depends on '$(top_srcdir)/.version', and
## '$(top_srcdir)/.version' depends on $(top_srcdir)/configure'. So any
## time the version is updated, both of these will be re-built. But during
## usual utility and library development, '$(top_srcdir)/configure' is not
## commonly updated, so this will not slow down the process.
BUILT_SOURCES = $(top_srcdir)/AUTHORS



## Directories to check:
## =====================
##
## Gnulib is intentionally compiled and tested prior to the gnuastro
## progams. With this configuration, the user can confirm if the
## problem is internal to gnulib or still persists after all the
## Gnulib tests have been done.
if COND_ARITHMETIC
  MAYBE_ARITHMETIC = bin/arithmetic
endif
if COND_BUILDPROG
  MAYBE_BUILDPROG = bin/buildprog
endif
if COND_CONVERTT
  MAYBE_CONVERTT = bin/convertt
endif
if COND_CONVOLVE
  MAYBE_CONVOLVE = bin/convolve
endif
if COND_COSMICCAL
  MAYBE_COSMICCAL = bin/cosmiccal
endif
if COND_CROP
  MAYBE_CROP = bin/crop
endif
if COND_FITS
  MAYBE_FITS = bin/fits
endif
if COND_MATCH
  MAYBE_MATCH = bin/match
endif
if COND_MKCATALOG
  MAYBE_MKCATALOG = bin/mkcatalog
endif
if COND_MKNOISE
  MAYBE_MKNOISE = bin/mknoise
endif
if COND_MKPROF
  MAYBE_MKPROF = bin/mkprof
endif
if COND_NOISECHISEL
  MAYBE_NOISECHISEL = bin/noisechisel
endif
if COND_SEGMENT
  MAYBE_SEGMENT = bin/segment
endif
if COND_STATISTICS
  MAYBE_STATISTICS = bin/statistics
endif
if COND_TABLE
  MAYBE_TABLE = bin/table
endif
#if COND_TEMPLATE
#  MAYBE_TEMPLATE = bin/TEMPLATE
#endif
if COND_WARP
  MAYBE_WARP = bin/warp
endif
if COND_GNULIBCHECK
  MAYBE_GNULIBCHECK = bootstrapped/tests
endif



## Subdirectories to build
## =======================
##
## Note that by default 'COND_TEMPLATE' is not set in configure, it is
## commented, and exists only as a template for you to copy and paste to
## name your new utility. The same rule is applied here (in the 'if'
## conditions above). When 'MAYBE_TEMPLATE' is not defined, then Make will
## see it as a blank string and igonore it, so there is no problem with
## having an uncommented 'MAYBE_TEMPLATE' as a value in 'SUBDIRS'.
SUBDIRS = bootstrapped/lib $(MAYBE_GNULIBCHECK) lib $(MAYBE_ARITHMETIC)    \
  $(MAYBE_BUILDPROG) $(MAYBE_CONVERTT) $(MAYBE_CONVOLVE)                   \
  $(MAYBE_COSMICCAL) $(MAYBE_CROP) $(MAYBE_FITS) $(MAYBE_MATCH)            \
  $(MAYBE_MKCATALOG) $(MAYBE_MKNOISE) $(MAYBE_MKPROF) $(MAYBE_NOISECHISEL) \
  $(MAYBE_SEGMENT) $(MAYBE_STATISTICS) $(MAYBE_TABLE) $(MAYBE_TEMPLATE)    \
  $(MAYBE_WARP) bin/script doc tests



## Add m4/ to the list of directories to search for m4 files
## =========================================================
ACLOCAL_AMFLAGS = -I bootstrapped/m4



## Files that are installed (and distributed)
## =========================================
dist_doc_DATA = README



# Installed system configuration files
# ====================================
dist_sysconf_DATA = bin/gnuastro.conf



## Files that are only distributed
## ===============================
##
## Note that 'COPYING' (containing the GNU GPL) is included in the
## distribution tarball by default in Automake, but not other license
## files, so we have to manually add 'COPYING.FDL'.
EXTRA_DIST = COPYING.FDL genauthors .dir-locals.el .version \
             developer-build bootstrapped/README .autom4te.cfg



## Print the completion messages
## =============================
##
## After the jobs in all the subdirectories are done, print the following
## messages to confirm that everything is complete and guide the users on
## what they should do next.
##
## When running ./configure, the user can opt-out of these messages using
## the GUIDEMESSAGE variable that is set when they run ./configure with the
## '--distable-guide-message' option.
##
## Note that all-local is a prerequisite of 'make check' too, so we will
## only print its message when 'make' was called with no options. Make will
## set MAKECMDGOALS to blank if there are no arguments, however, the way
## Automake works, its value is set to "all-am".
all-local:
	@if [ x$(MAKECMDGOALS) = x"all-am" ] && [ x$(GUIDEMESSAGE) = xyes ]; then    \
	 echo;                                                                       \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo "Gnuastro $(VERSION), was successfully built.";                        \
	 echo "Please check the build on your system by running:";                   \
	 echo;                                                                       \
	 echo "    make check -j8";                                                  \
	 echo;                                                                       \
	 echo "(You can change the 8 to the number of CPU threads available.)";      \
	 echo "(The following \"Leaving directory\" notices can be ignored.)";       \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo;                                                                       \
	fi

## We cannot do the writing persmission test at configure time. Here is a
## quote from the generated ./configure script: The prefix variables...
##
##   "... are left unexpanded so users can "make install exec_prefix=/foo"
##   and all the variables that are supposed to be based on exec_prefix by
##   default will actually change."
##
## So actually, they can also do that after 'make check'. Usually
## professional users would want to do such a thing, so they can ignore the
## NOTE. This note is mostly for beginners and it is not written to convey
## that this is the ONLY solution.
check-local:
	@if [ x$(GUIDEMESSAGE) = xyes ]; then                                        \
	 echo;                                                                       \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo "Your build of Gnuastro $(VERSION) didn't fail on any tests.";         \
	 echo "To install Gnuastro, please run the command(s) below:";               \
	 echo;                                                                       \
	 if [ ! -w $(prefix) ] ; then                                                \
	  echo "(NOTE: As the user $$USER, you don't have writing permissions";      \
	  echo "in \"$(prefix)\". To install Gnuastro there, you must get ";         \
	  echo "privileges first, as described below. If you can't, configure ";     \
	  echo "Gnuastro again, for example: \"./configure --prefix ~/.local\".)";   \
	  echo;                                                                      \
	  echo "    sudo make install";                                              \
	  echo "or";                                                                 \
	  echo "    su";                                                             \
	 fi;                                                                         \
	 echo "    make install";                                                    \
	 if [ ! -w $(prefix) ] ; then                                                \
	  echo "    exit";                                                           \
	 fi;                                                                         \
	 echo;                                                                       \
	 echo "(The following \"Leaving directory\" notices can be ignored.)";       \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo;                                                                       \
	fi

## Note that the '\' characters in the GNU head here are not printed on the
## command line. So we have to consider them. The ASCII GNU head is taken
## from: https://www.gnu.org/graphics/gnu-ascii.html
install-exec-local:
	@if [ x$(GUIDEMESSAGE) = xyes ]; then                                        \
	 echo;                                                                       \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo "      ,           , ";                                                \
	 echo "     /             \ ";                                               \
	 echo "    ((__-^^-,-^^-__))       Congratulations!";                        \
	 echo "     \`-_---' \`---_-'        GNU Astronomy Utilities (Gnuastro),";   \
	 echo "      \`--|o\` 'o|--'         Version $(VERSION)";                    \
	 echo "         \  \`  / ";                                                  \
	 echo "          ): :(             Successfully installed on this system.";  \
	 echo "          :o_o: ";                                                    \
	 echo "           \"-\" ";                                                   \
	 echo;                                                                       \
	 echo "More information                     Command to run ";                \
	 echo "----------------                     -------------- ";                \
	 echo "Complete official Gnuastro book:    ' info gnuastro           '";     \
	 echo "Entertaining tutorials:             ' info gnuastro Tutorials '";     \
	 echo "Dedicated help mailing list:        ' info help-gnuastro      '";     \
	 echo "Instructions for reporting bugs:    ' info bug-gnuastro       '";     \
	 echo "Effectively use Info:               ' info info               '";     \
	 echo;                                                                       \
	 echo;                                                                       \
	 echo "Environment variables to check: ";                                    \
	 echo "  - '$(prefix)/bin' in PATH.";                                        \
	 echo "  - '$(prefix)/lib' in LD_LIBRARY_PATH.";                             \
	 echo "(for an intro, run 'info gnuastro \"Installation directory\"')";      \
	 echo;                                                                       \
	 echo;                                                                       \
	 echo "To stay up to date with future releases, please subscribe to: ";      \
	 echo "    https://lists.gnu.org/mailman/listinfo/info-gnuastro";            \
	 echo;                                                                       \
	 echo;                                                                       \
	 echo "(The following lines can be ignored.)";                               \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo;                                                                       \
	fi





## $(top_srcdir)/.version
## ======================
##
## This file is created from the $(VERSION) variable which was defined by
## the 'git-version-gen' script (located at address below), which is run
## when the '$(top_srcdir)/configure' script is being built by Autoconf.
##
##    $(top_srcdir)/bootstrapped/build-aux/git-version-gen
##
## Note that contrary to what is proposed by 'git-version-gen', here the
## creation of '$(top_srcdir)/.version' depends on the
## '$(top_srcdir)/configure' script. Therefore, anytime a the VERSION
## variable is updated there, '$(top_srcdir)/.version' is also
## updated. During development, of the main functionality of Gnuastro
## (utilities and libraries), the '$(top_srcdir/configure' script is rarely
## updated, so '$(top_srcdir)/.version' will not be rebuilt and thus it
## won't harm the speed of tests during development.
$(top_srcdir)/.version: $(top_srcdir)/configure
	echo $(VERSION) > $@-t && mv $@-t $@





## $(top_srcdir)/AUTHORS
## =====================
##
## This file is generated automatically from the version controlled
## history. Note the following:
##
##   - '$(top_srcdir)/AUTHORS' is defined as a BUILT_SOURCES variable, so
##     it is the first thing that is built (even in multi-threaded runs).
##
##   - '$(top_srcdir)/AUTHORS' is updated only when
##     '$(top_srcdir)/.version' is updated. '$(top_srcdir)/.version' its
##     self is only updated when '$(top_srcdir)/configure' is updated. But
##     generally, the '$(top_srcdir)/configure' script is not updated
##     regularly during development and outside of version control.
##
##   - The '$(top_srcdir)/genauthors' script will not do anything
##     (make/update the '$(top_srcdir)AUTHORS' file) when there is no git
##     repository. '$(top_srcdir)/AUTHORS' is only necessary when building
##     a tarball distribution.
$(top_srcdir)/AUTHORS: $(top_srcdir)/.version
	$(top_srcdir)/genauthors $(top_srcdir)





## Run when building a distribution
## ================================
##
## These targets will be created when building a (tarball)
## distribution. Note that AUTHORS depends on '.version'.
dist-hook: $(top_srcdir)/AUTHORS
	echo $(VERSION) > $(distdir)/.tarball-version




## Nice joke
## =========
##
## This joke is taken from the "The art of unix programming", and a quote
## by Stuart Feldman (creator of Make): "One of the older Unix jokes is
## "Make love" which results in "Don’t know how to make love"."
love:
	@echo "Don't know how to make love!"
