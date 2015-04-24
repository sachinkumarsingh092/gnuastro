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


# Get the current month to put under all the pages.
thismonth=$(date +"%B %Y")


# Make the proper corrections to the HTML files:
echo
echo %%%%% Correcting the HTMLs %%%%%


# Correct the address of the `(dir)' links on the top pages of both
# HTML outputs. In the ./manual/gnuastro.html, it is `dir.html#top'
# which should be change to index.html. In
# ./manual/html_node/index.html, it is `../dir/index.html' which
# should become ../index.html
cat ./manual/gnuastro.html | sed s/dir\.html\#Top/index.html/g > tmp.txt
mv tmp.txt ./manual/gnuastro.html
cat ./manual/html_node/index.html | sed -e 's/\/dir\//\//g' > tmp.txt
mv tmp.txt ./manual/html_node/index.html


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
                cssbase="./"
                jsbase="../"
            else
                cssbase="../"
                jsbase="../../"
            fi
            echo "<link rel=\"stylesheet\" type=\"text/css\" href=\""$cssbase"style.css\">" >> tmp.html
            if [ $hasjavascript = "yes" ]; then
                echo "<script type=\"text/javascript\"" >> tmp.html
                echo "src=\""$jsbase"MathJax/MathJax.js?config=TeX-AMS-MML_HTMLorMML\">" >> tmp.html
                echo "</script>" >> tmp.html
                echo "" >> tmp.html
            fi
        fi

        if [ "$line" = "</body>" ]; then

            # Put a blank line to finish the main body of the page.
            echo "<hr>" >> tmp.html

            # Add a notice in case Javascript is disabled.
            if [ $hasjavascript = "yes" ]; then
                echo "" >> tmp.html
                echo "<noscript>" >> tmp.html
                echo "<table class=\"cartouche\" style=\"background-color:#FF9912;\"><tr><td><p>" >> tmp.html
                echo "<b>Warning:</b> This page uses <a href=\"http://www.mathjax.org/\">MathJax</a>" >> tmp.html
                echo "to render TeX equations. MathJax requires JavaScript for the rendering." >> tmp.html
                echo "However, scripts are disabled." >> tmp.html
                echo "<br /><br />" >> tmp.html
                echo "To see the equations, you can either use <a href=\"https://www.gnu.org/software/librejs/\">LibreJS</a>" >> tmp.html
                echo "to allow trusted scripts, or get the full manual in <a href=\""$cssbase"gnuastro.pdf\">PDF</a>." >> tmp.html
                echo "</p></td></tr></table>" >> tmp.html
                echo "</noscript>" >> tmp.html
                echo "" >> tmp.html
            fi
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
        if [ "$temp1" = "<body" ]; then

            # Add a title.
            if [ $addtitle = "yes" ]; then
            echo "<h2>GNU Astronomy Utilities manual</h2>" >> tmp.html
            fi
        fi
    done < "$file"
    mv tmp.html "$file"
done
echo %%%%% DONE %%%%%

# Copy the generated files to the proper directory. We do not want to
# copy the full directory there because the CVS information will be
# removed.
cp -R manual/* www/gnuastro/manual/
rm -rf ./manual


# Run CVS to upload the page to the GNU server
cd www/gnuastro/
cvs commit -m "Update."
