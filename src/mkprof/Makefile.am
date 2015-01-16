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

AM_LDFLAGS = -lcfitsio -lwcs -lgsl -lgslcblas -pthread -lm	\
-L$(top_builddir)/lib





## ImageCrop, .astrutils is for testing
bin_PROGRAMS = astrmkprof

astrmkprof_SOURCES = main.c main.h ui.c ui.h mkprof.c mkprof.h args.h	\
oneprofile.c oneprofile.h profiles.c profiles.h cite.h

astrmkprof_LDADD = -lconfigfiles -lfitsarrayvv -ltxtarrayvv		\
-lcheckset -llinkedlist -ltiming -lastrthreads -lbox -lstatistics	\
-larraymanip



# To destribute the defaults file.
dist_sysconf_DATA = astrmkprof.conf





## Generally for the man pages:
dist_man_MANS = astrmkprof.1

astrmkprof.1: args.h $(top_srcdir)/include/commonargs.h
	-$(HELP2MAN) --output=$@ --source="$(PACKAGE_STRING)" ./astrmkprof
