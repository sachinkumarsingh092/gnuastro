## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
##
## NOTE:
## This file will be called by the Makefile.am in the top directory,
## so this file does not need to be called directly by automake to
## make the full GNU Astronomy Utilities (AstrUtils).
##
## Copyright (C) 2013-2014 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/


# Flags for all programs:
AM_CFLAGS = -Wall -Ofast
AM_CPPFLAGS = -DSYSDATA_DIR=\"$(pkgdatadir)\" -I$(top_srcdir)/include


## Make libraries of the files that are shared between all the
## programs:
noinst_LIBRARIES = libdefaults.a libfitsarrayvv.a libtxtarrayvv.a	\
libcheckset.a liblinkedlist.a libtiming.a libastrthreads.a

libastrthreads_a_SOURCES = astrthreads.c

libdefaults_a_SOURCES = defaults.c

libfitsarrayvv_a_SOURCES = fitsarrayvv.c

libtxtarrayvv_a_SOURCES = txtarrayvv.c

libcheckset_a_SOURCES = checkset.c

liblinkedlist_a_SOURCES = linkedlist.c

libtiming_a_SOURCES = timing.c
