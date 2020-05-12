## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <mohammad@akhlaghi.org>
## Contributing author(s):
## Copyright (C) 2019-2020, Free Software Foundation, Inc.
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





## List of programs (scripts in this directory) to install under the
## 'prefix/bin' directory ('bin_SCRIPTS'), files necessary to distribute
## with the tarball ('EXTRA_DIST') and output files (to be cleaned with
## 'make clean').
bin_SCRIPTS = astscript-sort-by-night

EXTRA_DIST = sort-by-night.in

CLEANFILES = $(bin_SCRIPTS)





## Command to do basic substitutions (anything surrounded by an '@').
do_subst = sed -e 's,[@]VERSION[@],$(VERSION),g' \
               -e 's,[@]SCRIPT_NAME[@],$@,g'





## Rules to build the scripts
astscript-sort-by-night: sort-by-night.in Makefile
	$(do_subst) < $(srcdir)/sort-by-night.in > $@
	chmod +x $@
