#! /bin/sh
# Script to rebuild and test a utility during development.
#
# During the development of Gnuastro, you often make changes in a utility
# then need to rebuild and run it with special arguments and options (to
# test your work). This script is made to facilitate this process. It will
# take these steps:
#
#   1. Delete an existing utility executable.
#
#   2. Run Make on all Gnuastro.
#
#   3. If Make was successful, then go into the output directory, and run
#      utility with the given arguments and options. The executable is run
#      within an output directory (possibly different from the source or
#      build directories) so if you need to make lots of temporary test
#      files, there they won't get mixed up with non-output files.
#
# Combined with the `tmpfs-config-make', this script can be used to greatly
# simplify the development process. After running that script once, for
# subsequent builds during your development, you can run this script from
# the top source directory (by running `./tests/during-dev.sh', or giving
# this to the `compile' command in Emacs). Note that you have to set the
# first few variables (directories, utility name, arguments and options)
# manually before each major development activity.
#
# This file will be changed alot during development. So please make sure
# you have left the variable values empty before making commit so each
# developer is guided to set them for their own case. Alternatively you can
# reset this file to its version controlled status before making your
# commit:
#
#     git checkout -- tests/during-dev.sh
#
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
# Copyright (C) 2016, Free Software Foundation, Inc.
#
# Gnuastro is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# Gnuastro is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with Gnuastro. If not, see <http://www.gnu.org/licenses/>.




# SET INPUT PARAMETERS
# ====================

# Set the basic test directories. If you are building over the source
# directory, then set `builddir' to `./'. If you want the outputs to be in
# the top source directory, set it to `./'. Since 'build' is the assumed
# symbolic link in `tmpfs-config-make', it is also assumed in the version
# controlled version of this script.
builddir=build
outdir=


# Set the utility name, along with its arguments and options.
utilname=
arguments=
options=





# RUN THE PROCEDURES
# ==================

# First, make sure the variables are set. Note that arguments and options
# are not absolutly vital! so they are not checked here. The utility will
# warn and halt if it needs them.
if [ x$outdir = x ];   then echo "outdir is not set.";   exit 1; fi
if [ x$builddir = x ]; then echo "builddir is not set."; exit 1; fi
if [ x$utilname = x ]; then echo "utilname is not set."; exit 1; fi


# If builddir is relative, then append the current directory to make it
# absolute. This is done because we will be going into the output directory
# for executing the utility and we need to know the absolute address of the
# top build directory.
if [ ! "${builddir:0:1}" = "/" ]; then
   builddir=$(pwd)/$builddir
fi


# Set the utility's executable file name
utility=$builddir/src/$utilname/ast$utilname


# If the utility is already built, then remove the executable.
if [ -f $utility ]; then rm $utility; fi


# Make Gnuastro (note that during development, it is sometimes necessary to
# edit/rebuild the libraries too). If Make is successful, then change to
# the output directory and run the utility with the given arguments and
# options.
curdir=$(pwd)
if make -C $builddir; then
    cd $outdir
    $utility $arguments $options
fi
