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

## Add the version numbers of the separate programs:
## ================================================= List the header
## files in this string with a space character between them, don't
## forget to add `$(top_srcdir}/src/` before them like the ones
## already there.

SHELL=/bin/bash


## Commands to make the texinfo tools:
info_TEXINFOS = gnuastro.texi
gnuastro_TEXINFOS = fdl.texi authors.texi utilityversions.texi formath.texi


## Files not predefined by Automake, and not in dependencies that must
## be included in the final tar-ball distribution.
EXTRA_DIST = findversions.sh README


## Generate a list of Gnuastro authors from repo history:
authors.texi:
        # Only update this file if a .git directory is present in the
        # parent directory. Otherwise, this is being run from a
        # tar-ball and asuthors.texi cannot be updated.
        #
        # In case there is no git repository, then just touch
        # authors.texi. If one already exists, only its date will be
        # unchanged. However, if it doesn't already exist a blank file
        # will be created (which will allow us to generate the
        # documentation).
	if [ -d ../.git ]; then                                          \
	  git shortlog --numbered --summary --email                      \
	    | sed -e 's/</ /' -e 's/>/ /' -e 's/@/@@/' -e "s/è/@'e/"     \
            | awk '{printf "%s %s (%s, %s)@*\n", $$2, $$3, $$4, $$1}'    \
            > $@;                                                        \
        else touch $@;                                                   \
	fi


# Find the version of each utility:
headers=${top_srcdir}/src/*/main.h
utilityversions.texi: ${headers}
	./findversions.sh ${headers}


## Images:
infognuastrodir=${infodir}/gnuastro-figures/
dist_infognuastro_DATA = ${top_srcdir}/doc/gnuastro-figures/*
