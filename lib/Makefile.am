## Process this file with automake to produce Makefile.in
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





## Make libraries of the files that are shared between all the
## programs:
headersdir=$(top_srcdir)/lib/gnuastro
noinst_LTLIBRARIES = libgalconfigfiles.la libgalfits.la			  \
libgaltxtarrayvv.la libgalcheckset.la libgallinkedlist.la libgaltiming.la \
libgalthreads.la libgalbox.la libgalarraymanip.la libgalstatistics.la	  \
libgalspatialconvolve.la libgalmesh.la libgalmode.la libgalqsort.la	  \
libgalpolygon.la libgalwcs.la

libgalbox_la_SOURCES = box.c $(headersdir)/box.h

libgalwcs_la_SOURCES = wcs.c $(headersdir)/wcs.h

libgalfits_la_SOURCES = fits.c $(headersdir)/fits.h

libgalmesh_la_SOURCES = mesh.c $(headersdir)/mesh.h

libgalmode_la_SOURCES = mode.c $(headersdir)/mode.h

libgalqsort_la_SOURCES = qsort.c $(headersdir)/qsort.h

libgaltiming_la_SOURCES = timing.c $(headersdir)/timing.h

libgalpolygon_la_SOURCES = polygon.c $(headersdir)/polygon.h

libgalthreads_la_SOURCES = threads.c $(headersdir)/threads.h

libgalcheckset_la_SOURCES = checkset.c $(headersdir)/checkset.h

libgallinkedlist_la_SOURCES = linkedlist.c $(headersdir)/linkedlist.h

libgalstatistics_la_SOURCES = statistics.c $(headersdir)/statistics.h

libgalarraymanip_la_SOURCES = arraymanip.c $(headersdir)/arraymanip.h

libgaltxtarrayvv_la_SOURCES = txtarrayvv.c $(headersdir)/txtarrayvv.h

libgalconfigfiles_la_SOURCES = configfiles.c $(headersdir)/configfiles.h



# Dependent libraries:
libgalspatialconvolve_la_SOURCES = spatialconvolve.c	\
	$(headersdir)/spatialconvolve.h
libgalspatialconvolve_la_LDFLAGS = -lgalthreads -lgalbox
