#! /usr/bin/bash

umask 002

if [ -d ./manual ]; then rm -rf ./manual; fi

# Run gendocs.sh to generate all the files:
./gendocs.sh --email bug-gnuastro@gnu.org gnuastro \
             "GNU Astronomy Utilities manual"
rm gnuastro.aux gnuastro.cp gnuastro.cps gnuastro.fn gnuastro.ky \
   gnuastro.log gnuastro.pg gnuastro.toc gnuastro.tp gnuastro.vr



# Make the proper corrections to the HTML files:
echo
echo
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
        if [ $hasjavascript = "yes" ] && [ "$line" = "</head>" ]; then
            echo "<script type=\"text/javascript\"" >> tmp.html
            echo "src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\">" >> tmp.html
            echo "</script>" >> tmp.html
            echo "" >> tmp.html
        fi

        if [ "$line" = "</body>" ]; then
            echo "<hr>" >> tmp.html
            echo "<p><a href=\"http://www.gnu.org/software/gnuastro/manual\">Read in other formats</a>." >> tmp.html
            if [ $hasjavascript = "yes" ]; then
                echo "<br><a href=\"http://www.gnu.org/software/gnuastro/manual/javascript.html\" rel=\"jslicense\">JavaScript license information</a>" >> tmp.html
            fi
            echo "<br><a href=\"http://www.gnu.org/software/gnuastro\">GNU Astronomy Utilities</a> manual, $thismonth.</p>" >> tmp.html
        fi

        # Add the line:
        echo $line >> tmp.html

        # Actions after a given line:
        temporary=${line:0:5}
        if [ $addtitle = "yes" ] &&  [ "$temporary" = "<body" ]; then
            echo "<h2>GNU Astronomy Utilities manual</h2>" >> tmp.html
        fi
    done < "$file"
    mv tmp.html "$file"
done


# Copy the generated files to the proper directory
cp -R javascript.html manual/* ../www/gnuastro/manual/
rm -rf ./manual


# Run CVS to upload the page to the GNU server
cd ../www/gnuastro/
cvs commit
