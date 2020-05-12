## Process this file with automake to produce Makefile.inx
##
## This is part of GNU Astronomy Utilities (gnuastro). It will be
## called by the Makefile.am in the top directory.
##
## Original author:
##     Mohammad Akhlaghi <mohammad@akhlaghi.org>
## Contributing author(s):
## Copyright (C) 2015-2020, Free Software Foundation, Inc.
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
SHELL=/bin/sh

BUILT_SOURCES = $(srcdir)/authors.texi


## fdl.texi (The GNU Free documentation license in Texinfo format) is
## bootstrapped from Gnulib. Therefore it is not in this
## directory. However, it is needed by gnuastro.texi to build the
## documentation. It is the job of AM_MAKEINFOFLAGS, to pass options to the
## documentation build programs. BOOTSTRPDOC is defined in 'configure.ac'.
bootstrpdoc = $(top_srcdir)/bootstrapped/doc
AM_MAKEINFOFLAGS = -I $(bootstrpdoc) -I $(srcdir)


# This is just a temporary work-around since Automake does not pass the
# AM_MAKEINFOFLAGS to texi2dvi (which builds DVI and PDF outputs). This was
# reported in Automake's bug 23599. Note that a space is necessary between
# the '-I' and directory for texi2dvi.
TEXI2DVI = texi2dvi -I $(bootstrpdoc) -I $(srcdir)


## Commands to make the texinfo tools.
info_TEXINFOS = gnuastro.texi
gnuastro_TEXINFOS = $(bootstrpdoc)/fdl.texi $(bootstrpdoc)/gpl-3.0.texi \
                    $(srcdir)/authors.texi formath.texi


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










## Man pages to build
## ==================
##
## Installing all the programs is not mandatory in Gnuastro. If a program
## is not built, we shouldn't build its man-page either.
if COND_ARITHMETIC
  MAYBE_ARITHMETIC_MAN = man/astarithmetic.1
endif
if COND_BUILDPROG
  MAYBE_BUILDPROG_MAN = man/astbuildprog.1
endif
if COND_CONVERTT
  MAYBE_CONVERTT_MAN = man/astconvertt.1
endif
if COND_CONVOLVE
  MAYBE_CONVOLVE_MAN = man/astconvolve.1
endif
if COND_COSMICCAL
  MAYBE_COSMICCAL_MAN = man/astcosmiccal.1
endif
if COND_CROP
  MAYBE_CROP_MAN = man/astcrop.1
endif
if COND_FITS
  MAYBE_FITS_MAN = man/astfits.1
endif
if COND_MATCH
  MAYBE_MATCH_MAN = man/astmatch.1
endif
if COND_MKCATALOG
  MAYBE_MKCATALOG_MAN = man/astmkcatalog.1
endif
if COND_MKNOISE
  MAYBE_MKNOISE_MAN = man/astmknoise.1
endif
if COND_MKPROF
  MAYBE_MKPROF_MAN = man/astmkprof.1
endif
if COND_NOISECHISEL
  MAYBE_NOISECHISEL_MAN = man/astnoisechisel.1
endif
if COND_SEGMENT
  MAYBE_SEGMENT_MAN = man/astsegment.1
endif
if COND_STATISTICS
  MAYBE_STATISTICS_MAN = man/aststatistics.1
endif
if COND_TABLE
  MAYBE_TABLE_MAN = man/asttable.1
endif
if COND_WARP
  MAYBE_WARP_MAN = man/astwarp.1
endif
#if COND_TEMPLATE
#  MAYBE_TEMPLATE_MAN = man/astTEMPLATE.1
#endif
dist_man_MANS = $(MAYBE_ARITHMETIC_MAN) $(MAYBE_BUILDPROG_MAN)          \
  $(MAYBE_CONVERTT_MAN) $(MAYBE_CONVOLVE_MAN) $(MAYBE_COSMICCAL_MAN)    \
  $(MAYBE_CROP_MAN) $(MAYBE_FITS_MAN) $(MAYBE_MATCH_MAN)                \
  $(MAYBE_MKCATALOG_MAN) $(MAYBE_MKNOISE_MAN) $(MAYBE_MKPROF_MAN)       \
  $(MAYBE_NOISECHISEL_MAN) $(MAYBE_SEGMENT_MAN) $(MAYBE_STATISTICS_MAN) \
  $(MAYBE_TABLE_MAN) $(MAYBE_WARP_MAN) man/astscript-sort-by-night.1


## See if help2man is present or not. When help2man doesn't exist, we don't
## want to overwhelm the user with any commands, so we just let them know
## that the distributed man pages will be used.
if COND_HASHELP2MAN
  MAYBE_HELP2MAN = help2man --no-discard-stderr --output=$@ \
                   --source="$(PACKAGE_STRING)"
else
  MAYBE_HELP2MAN = @echo "Using distributed man page for"
endif


# Build the 'man' directory and then put all the man pages in
# it. Unfortunately as far as I know, pattern rules are not portable in all
# implementations of Make, so we have to list all the utilities manually.
toputildir=$(top_builddir)/bin
ALLMANSDEP = $(top_srcdir)/lib/gnuastro-internal/options.h
man/astarithmetic.1: $(top_srcdir)/bin/arithmetic/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "arithmetic operations on images and numbers" \
	                  --libtool $(toputildir)/arithmetic/astarithmetic

man/astbuildprog.1: $(top_srcdir)/bin/buildprog/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "compile, link with Gnuastro library and its dependencies, and run a C program"      \
	                  --libtool $(toputildir)/buildprog/astbuildprog

man/astconvertt.1: $(top_srcdir)/bin/convertt/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "convert known data types to each other"      \
	                  --libtool $(toputildir)/convertt/astconvertt

man/astconvolve.1: $(top_srcdir)/bin/convolve/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "convolve an image with a given kernel"       \
	                  --libtool $(toputildir)/convolve/astconvolve

man/astcosmiccal.1: $(top_srcdir)/bin/cosmiccal/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "estimate cosmological values"                \
	                  --libtool $(toputildir)/cosmiccal/astcosmiccal

man/astcrop.1: $(top_srcdir)/bin/crop/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "crop regions of a dataset"                   \
	                  --libtool $(toputildir)/crop/astcrop

man/astfits.1: $(top_srcdir)/bin/fits/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "view and manipulate FITS headers"            \
	                  --libtool $(toputildir)/fits/astfits

man/astmatch.1: $(top_srcdir)/bin/match/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "match catalogs by 1D or 2D positions"        \
	                  --libtool $(toputildir)/match/astmatch

man/astmkcatalog.1: $(top_srcdir)/bin/mkcatalog/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "Make a catalog from labeled input images"    \
	                  --libtool $(toputildir)/mkcatalog/astmkcatalog

man/astmknoise.1: $(top_srcdir)/bin/mknoise/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "add noise to an image"                       \
	                  --libtool $(toputildir)/mknoise/astmknoise

man/astmkprof.1: $(top_srcdir)/bin/mkprof/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "make standard 2D profiles on an image"       \
	                  --libtool $(toputildir)/mkprof/astmkprof

man/astnoisechisel.1: $(top_srcdir)/bin/noisechisel/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "detect signal in a noisy image"              \
	                  --libtool $(toputildir)/noisechisel/astnoisechisel

man/astscript-sort-by-night.1: $(top_srcdir)/bin/script/sort-by-night.in   \
                               $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "Sort input FITS files by night"              \
	                  --libtool $(toputildir)/script/astscript-sort-by-night

man/astsegment.1: $(top_srcdir)/bin/segment/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "segmentation based on signal structure"      \
	                  --libtool $(toputildir)/segment/astsegment

man/aststatistics.1: $(top_srcdir)/bin/statistics/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "calculate statistics of a dataset"           \
	                  --libtool $(toputildir)/statistics/aststatistics

man/asttable.1: $(top_srcdir)/bin/table/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "read and write ASCII and FITS tables"        \
	                  --libtool $(toputildir)/table/asttable

#man/astTEMPLATE.1: $(top_srcdir)/bin/TEMPLATE/args.h  $(ALLMANSDEP)
#	$(MAYBE_HELP2MAN) -n "put a short description here"                \
	                  $(toputildir)/TEMPLATE/astTEMPLATE

man/astwarp.1: $(top_srcdir)/bin/warp/args.h  $(ALLMANSDEP)
	$(MAYBE_HELP2MAN) -n "warp (transform) input dataset"              \
	                  --libtool $(toputildir)/warp/astwarp
