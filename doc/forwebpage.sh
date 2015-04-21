#! /usr/bin/bash

# This shell script will first run gendocs.sh on the Texinfo source
# which creates the standrard GNU webpage style documentation
# page. Then it does some modifications to the HTMLs produced:
#
# 1. Add a CSS sytlesheet.
# 2. Add MathJax.
# 3. Add a "GNU Astronomy Utilities manual" title on the top.
# 4. Add a link to read in different formats.
# 5. If MathJax was used in a page, add a link to Javascript.html
# 6. Add a link to the main Gnuastro webpage on the bottom.
# 7. Remove the border from cartouche in the HTML.
#
# After making all these changes, it copies all the necessry files in
# in the folder that will be linked to the GNU servers and it will run
# CVS to upload them there.
#
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
# Copyright (C) 2015, Free Software Foundation, Inc.
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
# along with gnuastro. If not, see <http://www.gnu.org/licenses/>.

umask 002

if [ -d ./manual ]; then rm -rf ./manual; fi

# Run gendocs.sh to generate all the files:
./gendocs.sh --email bug-gnuastro@gnu.org gnuastro \
             "GNU Astronomy Utilities manual"
rm gnuastro.aux gnuastro.cp gnuastro.cps gnuastro.fn gnuastro.ky \
   gnuastro.log gnuastro.pg gnuastro.toc gnuastro.tp gnuastro.vr


# Copy the two necessary files in the manual directory:
cp javascript.html style.css ./manual/


# Make the proper corrections to the HTML files:
echo
echo %%%%% Correcting the HTMLs %%%%%
thismonth=$(date +"%B %Y")
if [ -f tmp.html ]; then rm tmp.html; fi

for file in manual/gnuastro.html manual/html_node/*.html
do
    if grep -q '\\(\|$$' "$file";
    then hasjavascript="yes"
    else hasjavascript="no"
    fi

    if [ "$file" != manual/gnuastro.html ] && [ "$file" != manual/html_node/index.html ];
    then addtitle="yes"
    else addtitle="no"
    fi

    while read -r line; do

        # Actions that must be done before a given line:
        if [ "$line" = "</head>" ]; then
            if [ "$file" = manual/gnuastro.html ]; then
                echo "<link rel=\"stylesheet\" type=\"text/css\" href=\"./style.css\">" >> tmp.html
            else
                echo "<link rel=\"stylesheet\" type=\"text/css\" href=\"../style.css\">" >> tmp.html
            fi
            if [ $hasjavascript = "yes" ]; then
                echo "<script type=\"text/javascript\"" >> tmp.html
                echo "src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\">" >> tmp.html
                echo "</script>" >> tmp.html
                echo "" >> tmp.html
            fi
        fi

        if [ "$line" = "</body>" ]; then
            echo "<hr>" >> tmp.html
            echo "<p><a href=\"http://www.gnu.org/software/gnuastro/manual\">Read in other formats</a>." >> tmp.html
            if [ $hasjavascript = "yes" ]; then
                echo "<br><a href=\"http://www.gnu.org/software/gnuastro/manual/javascript.html\" rel=\"jslicense\">JavaScript license information</a>" >> tmp.html
            fi
            echo "<br><a href=\"http://www.gnu.org/software/gnuastro\">GNU Astronomy Utilities</a> manual, $thismonth.</p>" >> tmp.html
        fi

        #Do any modifications to the line:
        temp1=${line:0:5}
        if [ "$temp1" = "<tabl" ]; then
            temp2=$(echo "$line" | sed -e "s/class=\"cartouche\" border=\"1\"/class=\"cartouche\"/")
            line=$temp2
        fi

        # Add the line:
        echo $line >> tmp.html

        # Actions after a given line:
        if [ $addtitle = "yes" ] &&  [ "$temp1" = "<body" ]; then
            echo "<h2>GNU Astronomy Utilities manual</h2>" >> tmp.html
        fi
    done < "$file"
    mv tmp.html "$file"
done


# Copy the generated files to the proper directory. We do not want to
# copy the full directory there because the CVS information will be
# removed.
cp -R manual/* www/gnuastro/manual/
rm -rf ./manual


# Run CVS to upload the page to the GNU server
cd www/gnuastro/
cvs commit
