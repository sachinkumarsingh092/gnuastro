# This is not actually a test of any of the programs, it just brings
# in all the configuration files into a locally created .gnuastro
# directory for all the tests to use. It is part of the GNU Astronomy
# Utilities (Gnuastro).
#
# Original author:
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
# Copyright (C) 2016-2020, Free Software Foundation, Inc.
#
# Gnuastro is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# Gnuastro is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.




# Make the .gnuastro directory
# ----------------------------
#
# This directory will keep the default configuration files for all the
# programs. If it already exists, delete it. 'mkdir_p' is the equivalent to
# GNU's 'mkdir -p' which might not exist on some systems. It comes from
# Autoconf's tests and is exported to the test shell scripts from the
# 'tests/Makefile.am' file.
$mkdir_p .gnuastro





# Common options for all programs
# -------------------------------
#
# Copy the common options while adding the following optios only for make
# check.
#
#   - 'lastconfig' will make sure that the program stop searching for
#     configuration files after this one.
#
#   - Log files are not necessary during tests, they are mainly used for
#     reporting extra information about successful tests. Failed messages
#     will be printed on the command-line not in a log-file. So to keep
#     this directory clean, we'll ask the programs to not generate any.
cat > addedoptions.txt <<EOF


# Added only for "make check":
 lastconfig            1
 log                   0
EOF
cat $topsrc/bin/gnuastro.conf addedoptions.txt > .gnuastro/gnuastro.conf
rm addedoptions.txt





# Bring utility configuration files
# ---------------------------------
#
# Each utility's configuration file is copied in the 'tests' directory for
# easy readability. Note that some programs may need to build their
# configuration files during compilation. Hence, their configuration files
# are in the build directory, not the source directory.
for prog in arithmetic buildprog convertt convolve cosmiccal crop fits    \
                       match mkcatalog mknoise mkprof noisechisel segment \
                       statistics table warp
do
    if test -f $topsrc/bin/$prog/ast$prog.conf; then
        ctopdir=$topsrc
    else
        ctopdir=$topbuild
    fi
    cp $ctopdir/bin/$prog/*.conf .gnuastro/
done
