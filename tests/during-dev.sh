#! /bin/bash
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
# Combined with the 'developer-build', this script can be used to greatly
# simplify the development process. After running that script once, for
# subsequent builds during your development, you can run this script from
# the top source directory (by running './tests/during-dev.sh', or giving
# this to the 'compile' command in Emacs). Note that you have to set the
# first few variables (directories, utility name, arguments and options)
# manually before each major development activity.
#
# This file will be changed alot during development. So please make sure
# you have left the variable values empty before committing. This will
# remind each developer running this script to set the values based on
# their particular work. Alternatively you can reset this file to its
# version controlled status before committing your main work with the
# command below:
#
#     git checkout -- tests/during-dev.sh
#
# This file can also be used as a model to write a test for the work you
# have done (to be checked with 'make check'). Just copy and paste an
# existing test from the utility and replace the last few lines based on
# this file.
#
# Original author:
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
# Copyright (C) 2016-2020, Free Software Foundation, Inc.
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
# directory, then set 'builddir' to './'. If you want the outputs to be in
# the top source directory, set it to './'. Since 'build' is the assumed
# symbolic link in 'developer-build', it is also assumed in the version
# controlled version of this script. Note, if your directory names have
# space characters in them, quote the full value
numjobs=8
builddir=build
outdir=



# Set the utility name, along with its arguments and options. NOTE, for
# multiple arguments and options, please put them all between quotation
# marks so the space characters are included. If there are spaces in the
# values of the options or arguments, quote them two times (once for this
# script, and once for the utility. In such cases it might be easier to
# just add the argument/option to the final script that runs the utility
# rather than these variables.
utilname=
arguments=
options=



# RUN THE PROCEDURES
# ==================

# Stop the script if there are any errors.
set -e


# First, make sure the variables are set. Note that arguments and options
# are not absolutly vital! so they are not checked here. The utility will
# warn and halt if it needs them.
if [ x"$outdir"   = x ]; then echo "outdir is not set.";   exit 1; fi
if [ x"$numjobs"  = x ]; then echo "numjobs is not set.";  exit 1; fi
if [ x"$utilname" = x ]; then echo "utilname is not set."; exit 1; fi
if [ x"$builddir" = x ]; then echo "builddir is not set."; exit 1; fi


# Make sure 'utilname' doesn't start with 'ast' (a common mistake).
astprefix="${utilname:0:3}"
if [ x"$astprefix" = x"ast" ]; then
    echo "'utilname' must not start with 'ast'."; exit 1;
fi


# If builddir is relative, then append the current directory to make it
# absolute. This is done because we will be going into the output directory
# for executing the utility and we need to know the absolute address of the
# top build directory.
srcdir=$(pwd)
if [ ! "${builddir:0:1}" = "/" ]; then
   builddir="$srcdir/$builddir"
fi


# Set the utility's executable file name
longprefix="${utilname:0:6}"
if [ x"$longprefix" = x"script" ]; then
    execdir="script"
else
    execdir="$utilname"
fi
utility="$builddir/bin/$execdir/ast$utilname"


# If the utility is already built, then remove the executable so it is
# definitely remade.
if [ -f "$utility" ]; then rm "$utility"; fi


# Make Gnuastro (note that during development, it is sometimes necessary to
# edit/rebuild the libraries too). If Make is successful, then change to
# the output directory and run the utility with the given arguments and
# options.
#
# Before actually running put a copy of the configuration file in the
# output directory and also add the onlydirconf option so user or system
# wide configuration files don't interfere.
if make -j$numjobs -C "$builddir"; then

    # Change to the output directory.
    cd "$outdir"

    # Make the .gnuastro directory if it doesn't exist.
    if [ ! -d .gnuastro ]; then
       mkdir .gnuastro
    fi

    # Put a copy of this utility's configuration file there and add the
    # onlydirconf option. We are first printing an empty line just in case
    # the last line in the configuration file doesn't actualy end with a
    # new line (in which case the appended string will be added to the end
    # of the last line).
    if [ $utilname = buildprog ]; then
        extraopts="--la=$builddir/lib/libgnuastro.la"
        extraopts="$extraopts -I$srcdir/lib -L$builddir/lib"
        topconfdir="$builddir"
    else
        topconfdir="$srcdir"
    fi

    # Copy the configuration file(s).
    cfiles="$srcdir/bin/gnuastro.conf"
    if [ x"$longprefix" != x"script" ]; then
        cfiles="$cfiles $topconfdir/bin/$utilname/ast$utilname.conf"
    fi
    cp $cfiles .gnuastro/

    # Append 'lastconfig' option to 'gnuastro.conf', so the program doesn't
    # go into the system headers.
    echo ""               >> .gnuastro/gnuastro.conf
    echo " lastconfig 1"  >> .gnuastro/gnuastro.conf

    # Run the built utility with the given arguments and options.
    "$utility" $arguments $options $extraopts

    # Clean up.
    rm -rf .gnuastro
fi
