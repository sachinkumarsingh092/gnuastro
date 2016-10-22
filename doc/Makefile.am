## Process this file with automake to produce Makefile.inx
##
## This is part of GNU Astronomy Utilities (gnuastro). It will be
## called by the Makefile.am in the top directory.
##
## Original author:
##     Mohammad Akhlaghi <akhlaghi@gnu.org>
## Contributing author(s):
## Copyright (C) 2015, Free Software Foundation, Inc.
##
## GNU Astronomy Utilities is free software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## GNU Astronomy Utilities is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty
## of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GNU Astronomy Utilities. If not, see
## <http://www.gnu.org/licenses/>.


## Set the shell
SHELL=/bin/bash

BUILT_SOURCES = $(srcdir)/authors.texi


## fdl.texi (The GNU Free documentation license in Texinfo format) is
## bootstrapped from Gnulib. Therefore it is not in this
## directory. However, it is needed by gnuastro.texi to build the
## documentation. It is the job of AM_MAKEINFOFLAGS, to pass options to the
## documentation build programs. BOOTSTRPDOC is defined in `configure.ac'.
bootstrpdoc = $(top_srcdir)/bootstrapped/doc
AM_MAKEINFOFLAGS = -I $(bootstrpdoc) -I $(srcdir)


# This is just a temporary work-around since Automake does not pass the
# AM_MAKEINFOFLAGS to texi2dvi (which builds DVI and PDF outputs). This was
# reported in Automake's bug 23599. Note that a space is necessary between
# the `-I' and directory for texi2dvi.
TEXI2DVI = texi2dvi -I $(bootstrpdoc) -I $(srcdir)


## Commands to make the texinfo tools.
info_TEXINFOS = gnuastro.texi
gnuastro_TEXINFOS = $(bootstrpdoc)/fdl.texi formath.texi	\
                    $(srcdir)/authors.texi


## Files not predefined by Automake, and not in dependencies that must
## be included in the final tar-ball distribution.
EXTRA_DIST = genauthors README


## We want to build the authors.texi file only when we are building in the
## version controlled source. For the non version controlled souce,
## deleting this file is like somehow deleting formath.texi, it will not be
## rebuilt and result in an error.
$(srcdir)/authors.texi: $(top_srcdir)/configure
	$(top_srcdir)/doc/genauthors $(top_srcdir)

## Images:
infognuastrodir=$(infodir)/gnuastro-figures/
dist_infognuastro_DATA = $(top_srcdir)/doc/gnuastro-figures/*










## Man pages: to keep things managable, sort the utilities alphabetically.
dist_man_MANS = man/astarithmetic.1 man/astconvertt.1 man/astconvolve.1 \
man/astcosmiccal.1 man/astheader.1 man/astimgcrop.1 man/astimgstat.1    \
man/astimgwarp.1 man/astmkcatalog.1 man/astmknoise.1 man/astmkprof.1    \
man/astnoisechisel.1 man/astsubtractsky.1 man/asttable.1


## See if help2man is present or not. When help2man doesn't exist, we don't
## want to overwhelm the user with any commands, so we just let them know
## that the distributed man pages will be used.
if COND_HASHELP2MAN
  MAYBE_HELP2MAN = help2man --output=$@ --source="$(PACKAGE_STRING)"
else
  MAYBE_HELP2MAN = @echo "Using distributed man page for"
endif


# Build the `man' directory and then put all the man pages in
# it. Unfortunately as far as I know, pattern rules are not portable in all
# implementations of Make, so we have to list all the utilities manually.
toputildir=$(top_builddir)/bin
ALLMANSDEP = $(top_srcdir)/lib/commonargs.h
man/astarithmetic.1: $(top_srcdir)/bin/arithmetic/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/arithmetic/astarithmetic
man/astconvertt.1: $(top_srcdir)/bin/convertt/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/convertt/astconvertt
man/astconvolve.1: $(top_srcdir)/bin/convolve/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/convolve/astconvolve
man/astcosmiccal.1: $(top_srcdir)/bin/cosmiccal/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/cosmiccal/astcosmiccal
man/astheader.1: $(top_srcdir)/bin/header/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/header/astheader
man/astimgcrop.1: $(top_srcdir)/bin/imgcrop/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/imgcrop/astimgcrop
man/astimgstat.1: $(top_srcdir)/bin/imgstat/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/imgstat/astimgstat
man/astimgwarp.1: $(top_srcdir)/bin/imgwarp/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/imgwarp/astimgwarp
man/astmkcatalog.1: $(top_srcdir)/bin/mkcatalog/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/mkcatalog/astmkcatalog
man/astmknoise.1: $(top_srcdir)/bin/mknoise/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/mknoise/astmknoise
man/astmkprof.1: $(top_srcdir)/bin/mkprof/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/mkprof/astmkprof
man/astnoisechisel.1: $(top_srcdir)/bin/noisechisel/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/noisechisel/astnoisechisel
man/astsubtractsky.1: $(top_srcdir)/bin/subtractsky/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/subtractsky/astsubtractsky
man/asttable.1: $(top_srcdir)/bin/table/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) $(toputildir)/table/asttable
#man/astTEMPLATE.1: $(top_srcdir)/bin/TABLE/args.h  $(ALLMANSDEP)
#	$(MAYBE_HELP2MAN) $(toputildir)/TEMPLATE/astTEMPLATE
