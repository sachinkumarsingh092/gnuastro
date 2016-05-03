# This is not actually a test of any of the programs, it just brings
# in all the configuration files into a locally created .gnuastro
# directory for all the tests to use. It is part of the GNU Astronomy
# Utilities (Gnuastro).
#
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
# Copyright (C) 2016, Free Software Foundation, Inc.
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
# This directory will keep all the configuration files for this run of
# `make check'.
if [ ! -d .gnuastro ]; then
    mkdir .gnuastro
fi





# Add options for `make check'
# ----------------------------
#
# It might be necessary to add specific options to all the programs during
# `make check'. Therefore we have defined the file `addedoptions.txt' to
# keep these extra options and append them to the configuration file in the
# source directory of the utility.
#
#   - The onlydirconf option is added so the utilities don't go looking
#     into the user's home and system wide directories (which might contain
#     configuration files from older versions). If the option names have
#     changed or an option has been removed, such sitations will cause a
#     failed test.
cat > addedoptions.txt <<EOF


# Added only for "make check":
onlydirconf           1
EOF





# Bring utility configuration files
# ---------------------------------
#
# Each utility's configuration file is read and appended with the
# addedoptions.txt file to create the configuration file which will be used
# by `make check'.
for prog in arithmetic convertt convolve cosmiccal header imgcrop \
            imgstat imgwarp mkcatalog mknoise mkprof noisechisel  \
            subtractsky
do

    # Copy the configuration file from the utility source and add the
    # options added here.
    cat $topsrc/src/$prog/ast$prog.conf addedoptions.txt          \
        > .gnuastro/ast$prog.conf

done
