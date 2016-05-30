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

BUILT_SOURCES = authors.texi


## fdl.texi (The GNU Free documentation license in Texinfo format) is
## bootstrapped from Gnulib. Therefore it is not in this
## directory. However, it is needed by gnuastro.texi to build the
## documentation. It is the job of AM_MAKEINFOFLAGS, to pass options to the
## documentation build programs. BOOTSTRPDOC is defined in `configure.ac'.
bootstrpdoc = $(top_srcdir)/bootstrapped/doc
AM_MAKEINFOFLAGS = -I $(bootstrpdoc)


# This is just a temporary work-around since Automake does not pass the
# AM_MAKEINFOFLAGS to texi2dvi (which builds DVI and PDF outputs). This was
# reported in Automake's bug 23599.
TEXI2DVI = texi2dvi -I $(bootstrpdoc)


## Commands to make the texinfo tools.
info_TEXINFOS = gnuastro.texi
gnuastro_TEXINFOS = $(bootstrpdoc)/fdl.texi formath.texi authors.texi


## Files not predefined by Automake, and not in dependencies that must
## be included in the final tar-ball distribution.
EXTRA_DIST = genauthors README


## We want to build the authors.texi file only when we are building in the
## version controlled source. For the non version controlled souce,
## deleting this file is like somehow deleting formath.texi, it will not be
## rebuilt and result in an error.
authors.texi: $(top_srcdir)/configure
	$(top_srcdir)/doc/genauthors $(top_srcdir)


## Images:
infognuastrodir=$(infodir)/gnuastro-figures/
dist_infognuastro_DATA = $(top_srcdir)/doc/gnuastro-figures/*
