#! /bin/bash

# This shell script will convert all the images into all the necessary
# formats that Texinfo accepts so the manual can be created in any
# format. This script is part of the manual of GNU Astronomy Utilities
# (Gnuastro).
#
# There is one input argument for this script:
#
#   1) The address of the folder keeping the images to convert.
#
#
# Original author:
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
# Copyright (C) 2015-2019, Free Software Foundation, Inc.
#
# Gnuastro is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Gnuastro is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.


# Initial settings:
set -o nounset                  # Stop if a variable is not set.
set -o errexit                  # Stop if a program returns false.




# There should only be two types of files here:
#
#     png: Files that were originally rasterized.
#          Converts to eps and txt.
#
#     eps: vector graphics files.
#          Converts to pdf, png and txt.
FIGDIR=$1
for FILENAME in $FIGDIR/*
do

    # Get the necessary information from the file name.
    NODIR=$(basename "$FILENAME")
    BASENAME="${NODIR%.*}"
    EXTENSION="${NODIR##*.}"


    # When the inputs are PNG, they should only be conerted to
    # EPS. Later, all the necessary files will be converted to txt
    # too.
    if [ $EXTENSION = "png" ]; then
        OUTNAME=$FIGDIR/$BASENAME".eps"
        if [ ! -f $OUTNAME ]; then
            convert $FILENAME $OUTNAME
        fi


    # Initially there will only be EPS and PNG files. But this script
    # might be run after that. After the initial run, all initial PNG
    # files will also have an EPS file. So running this script will
    # add a PDF file for them! We don't want that. The best way to
    # avoid this is to check if a PNG file has already been created
    # for this EPS file. Only if it hasn't then go along and make the
    # PDF and EPS files.
    elif [ $EXTENSION = "eps" ]; then
        OUTNAME=$FIGDIR$BASENAME".png"
        if [ ! -f $OUTNAME ]; then
            epspdf $FILENAME
            convert -density 150 $FILENAME $OUTNAME
        fi
    else
        continue
    fi


    # Write the text file:
    TXTNAME=$FIGDIR$BASENAME".txt"
    echo $FILENAME > $TXTNAME

done
