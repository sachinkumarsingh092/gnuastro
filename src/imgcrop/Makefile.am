## Process this file with automake to produce Makefile.inx
##
## Hand written file: used as input into Autotools.
##
## NOTE:
## This file will be called by the Makefile.am in the top directory,
## so this file does not need to be called directly by automake to
## make the full GNU Astronomy Utilities (Astrutils)
##
## Copyright (C) 2013-2014 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/


# Flags:
AM_CFLAGS = -Wall -Ofast

AM_CPPFLAGS = -DSYSDATA_DIR=\"$(pkgdatadir)/\" -I$(top_srcdir)/include

AM_LDFLAGS = -lcfitsio -lwcs -pthread -lm -L$(top_builddir)/lib





## ImageCrop, .astrutils is for testing
bin_PROGRAMS = astrimgcrop 

astrimgcrop_SOURCES = main.c main.h ui.c ui.h args.h imgcrop.c	\
imgcrop.h wcsmode.c wcsmode.h crop.c crop.h

astrimgcrop_LDADD = -ldefaults -lfitsarrayvv -ltxtarrayvv -lcheckset	\
-llinkedlist -ltiming -lastrthreads



# To destribute the defaults file.
dist_pkgdata_DATA = astrimgcrop.def





## Generally for the man pages:
dist_man_MANS = astrimgcrop.1

astrimgcrop.1: args.h $(top_srcdir)/include/commonargs.h 
	-$(HELP2MAN) --output=$@ --source="$(PACKAGE_STRING)" ./astrimgcrop
