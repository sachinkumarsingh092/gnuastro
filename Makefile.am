## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
## This is part of GNU Astronomy Utilities (gnuastro).
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



## Sources to be run before everything else
## ========================================
##
## The .version file is managed by the by the following file to
## automatically generate the version:
##
##    bootstrapped/build-aux/git-version-gen
##
## we want it to be created before anything else after running `make',
## `make check', or `make install', see the "Built sources" section of the
## AUTOMAKE manual and the git-version-gen script comments for more
## information.
BUILT_SOURCES = $(top_srcdir)/.version
$(top_srcdir)/.version:; echo $(VERSION) > $@-t && mv $@-t $@
dist-hook:; echo $(VERSION) > $(distdir)/.tarball-version



## Directories to check:
## =====================
##
## Gnulib is intentionally compiled and tested prior to the gnuastro
## progams. With this configuration, the user can confirm if the
## problem is internal to gnulib or still persists after all the
## Gnulib tests have been done.
if COND_ARITHMETIC
  MAYBE_ARITHMETIC = src/arithmetic
endif
if COND_CONVERTT
  MAYBE_CONVERTT = src/convertt
endif
if COND_CONVOLVE
  MAYBE_CONVOLVE = src/convolve
endif
if COND_COSMICCAL
  MAYBE_COSMICCAL = src/cosmiccal
endif
if COND_HEADER
  MAYBE_HEADER = src/header
endif
if COND_IMGCROP
  MAYBE_IMGCROP = src/imgcrop
endif
if COND_IMGSTAT
  MAYBE_IMGSTAT = src/imgstat
endif
if COND_IMGWARP
  MAYBE_IMGWARP = src/imgwarp
endif
if COND_MKCATALOG
  MAYBE_MKCATALOG = src/mkcatalog
endif
if COND_MKNOISE
  MAYBE_MKNOISE = src/mknoise
endif
if COND_MKPROF
  MAYBE_MKPROF = src/mkprof
endif
if COND_NOISECHISEL
  MAYBE_NOISECHISEL = src/noisechisel
endif
if COND_SUBTRACTSKY
  MAYBE_SUBTRACTSKY = src/subtractsky
endif
if COND_GNULIBCHECK
  MAYBE_GNULIBCHECK = bootstrapped/tests
endif
SUBDIRS = bootstrapped/lib $(MAYBE_GNULIBCHECK) lib $(MAYBE_CONVERTT)   \
$(MAYBE_CONVOLVE) $(MAYBE_COSMICCAL) $(MAYBE_HEADER)                    \
$(MAYBE_ARITHMETIC) $(MAYBE_IMGCROP) $(MAYBE_IMGSTAT)                   \
$(MAYBE_IMGWARP) $(MAYBE_MKCATALOG) $(MAYBE_MKNOISE) $(MAYBE_MKPROF)    \
$(MAYBE_NOISECHISEL) $(MAYBE_SUBTRACTSKY) doc tests



## Add m4/ to the list of directories to search for m4 files
## =========================================================
ACLOCAL_AMFLAGS = -I bootstrapped/m4



## Files that are installed (and distributed)
## =========================================
dist_doc_DATA = README



## Files that are only distributed
## ===============================
EXTRA_DIST = bootstrap bootstrap.conf genauthors .dir-locals.el         \
 $(top_srcdir)/.version bootstrapped/m4/gnulib-cache.m4                 \
 bootstrapped/README .autom4te.cfg include/commonargs.h                 \
 include/commonparams.h include/fixedstringmacros.h include/neighbors.h



## Print the completion messages
## =============================
##
## After the jobs in all the subdirectories are done, print the following
## messages to confirm that everything is complete and guide the users on
## what they should do next.
##
## Note that all-local is a prerequisite of `make check' too, so we will
## only print its message when `make' was called with no options. Make will
## set MAKECMDGOALS to blank if there are no arguments, however, the way
## Automake works, its value is set to "all-am".
all-local:
	@if [ x$(MAKECMDGOALS) = x"all-am" ]; then                                   \
	 echo;                                                                       \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo "Gnuastro $(VERSION), was successfully built.";                        \
	 echo "Please check the build on your system by running:";                   \
	 echo;                                                                       \
	 echo "    make check";                                                      \
	 echo;                                                                       \
	 echo "(The following \"Leaving directory\" notices can be ignored.)";       \
	 echo "==================================================================="; \
	 echo "==================================================================="; \
	 echo;                                                                       \
	fi

check-local:
	@echo
	@echo "==================================================================="
	@echo "==================================================================="
	@echo "All checks succeeded for the built Gnuastro $(VERSION)."
	@echo "To install Gnuastro, please run:"
	@echo
	@if [ ! -w $(prefix) ] ; then                                                \
	  echo "(NOTE: As the user $$USER, you don't have writing permissions";      \
	  echo "in \"$(prefix)\", therefore the \"su\" command below should be run"; \
	  echo "first. Configure Gnuastro again if you can't run \"su\". This";      \
	  echo "is an example that will work: \"./configure --prefix ~/.local\".)";  \
	  echo;                                                                      \
	  echo "    su";                                                             \
	fi
	@echo "    make install"
	@echo
	@echo "(The following \"Leaving directory\" notices can be ignored.)"
	@echo "==================================================================="
	@echo "==================================================================="
	@echo

install-exec-local:
	@echo
	@echo "==================================================================="
	@echo "==================================================================="
	@echo "Gnuastro $(VERSION) is now successfully installed."
	@echo
	@echo "                  !!!! Congratulations !!! "
	@echo
	@echo "The complete official Gnuastro book is available by running:"
	@echo "    info gnuastro"
	@echo
	@echo "We have prepared some entertaining tutorials, please run:"
	@echo "    info gnuastro Tutorials"
	@echo
	@echo "There is a dedicated mailing list for helping you use Gnuastro:"
	@echo "    help-gnuastro@gnu.org"
	@echo
	@echo "If you have found a bug, please run the following command. It "
	@echo "contains instructions for the best ways of sharing them so we"
	@echo "can fix it as soon as possible."
	@echo "    info bug-gnuastro"
	@echo
	@echo "To stay up to date with future releases, please subscribe to "
	@echo "the Gnuastro announcement mailing list at:"
	@echo "    https://lists.gnu.org/mailman/listinfo/info-gnuastro"
	@echo
	@echo
	@echo "(The following lines can be ignored.)"
	@echo "==================================================================="
	@echo "==================================================================="
	@echo



## Generate AUTHORS for distributions
## ==================================
##
## The AUTHORS file is generated automatically from the version
## control history. Note the following two cases:
##
##   - The AUTHORS file is a prerequisite for the make dist* commands (all
##     make targets including "dist") in the Makefile.in that Automake
##     builds from this file. So it is not needed when we only run "make".
##
##   - AUTHORS is defined as a double colon rule so that it is called every
##     time any "make dist*" is called, even if AUTHORS already exists. It
##     is possible for it to be called in a non version controlled source,
##     in that case, it will not be updated (or created if it doesn't
##     already exist, see `genauthors'). But generally, "make dist*" is not
##     called too often and outside of version control.
##
##   - The genauthors script will not do anything (make/update the
##     AUTHORS file) when there is no git repository.
AUTHORS::; $(top_srcdir)/genauthors
