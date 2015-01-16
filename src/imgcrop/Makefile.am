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

AM_CPPFLAGS = -DSYSCONFIG_DIR=\"$(sysconfdir)/\" -I$(top_srcdir)/include

AM_LDFLAGS = -lcfitsio -lwcs -pthread -lm -L$(top_builddir)/lib





## ImageCrop:
bin_PROGRAMS = astrimgcrop

astrimgcrop_SOURCES = main.c main.h ui.c ui.h args.h imgcrop.c	\
imgcrop.h wcsmode.c wcsmode.h crop.c crop.h cite.h

astrimgcrop_LDADD = -lconfigfiles -lfitsarrayvv -ltxtarrayvv	\
-lcheckset -llinkedlist -ltiming -lastrthreads -lbox





# To destribute the defaults file.
dist_sysconf_DATA = astrimgcrop.conf





## Generally for the man pages:
dist_man_MANS = astrimgcrop.1

astrimgcrop.1: args.h $(top_srcdir)/include/commonargs.h
	-$(HELP2MAN) --output=$@ --source="$(PACKAGE_STRING)" ./astrimgcrop
