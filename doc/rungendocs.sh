./gendocs.sh --email bug-gnuastro@gnu.org gnuastro "GNU Astronomy Utilities manual"

cp -R manual/* ../www/gnuastro/manual/
rm -rf ./manual

rm gnuastro.aux gnuastro.cp gnuastro.cps gnuastro.fn gnuastro.ky \
   gnuastro.log gnuastro.pg gnuastro.toc gnuastro.tp gnuastro.vr

cd ../www/gnuastro/
cvs commit
