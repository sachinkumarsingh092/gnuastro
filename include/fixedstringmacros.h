/*********************************************************************
Fixed strings for use in all utilities.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef COPYRIGHT_H
#define COPYRIGHT_H

#define SHORTCOPYRIGHT "Copyright (C) 2015, Free Software Foundation, Inc."


#define SHORTLICENSE   "License GPLv3+: GNU General public license version 3 or later."


#define COPYRIGHT SHORTCOPYRIGHT"\n"SHORTLICENSE"\n"			\
  "This is free software: you are free to change and redistribute it.\n" \
  "There is NO WARRANTY, to the extent permitted by law."		\


#define TOPHELPINFO   "\n"SPACK_NAME" is part of "PACKAGE_STRING".\n"



/* This is fixed for all the packages. */
#define MOREHELPINFO   "\nFor more information, please run any of the "	\
  "following commands. They will respectively show you the `Invoking "	\
  SPACK_NAME"' subsection, the complete `"SPACK_NAME"' section, or the "\
  "full "PACKAGE_NAME" manual. In particular the first contains "       \
  "a very complete explanation of each option.\n\n"                       \
  "     info "SPACK"\n\n"                                                 \
  "     info "SPACK_NAME"\n\n"						\
  "     info "PACKAGE_TARNAME"\n\n"                                     \
  "If you couldn't find your answer in the manual, you can get "        \
  "direct help from experienced Gnuastro users and developers. "        \
  "For more information, please run:\n"                                 \
  "     info help-gnuastro\n\n"                                         \
  SPACK_NAME" options:"                                                 \



#define GNUASTROBIBTEX                                                  \
  "@ARTICLE{noisechisel,\n"                                             \
  "   author = {{Akhlaghi}, M. and {Ichikawa}, T.},\n"                  \
  "    title = \"{Noise-based Detection and Segmentation of Nebulous Objects}\",\n" \
  "  journal = {\apjs},\n"                                              \
  "archivePrefix = \"arXiv\",\n"                                        \
  "   eprint = {1505.01664},\n"                                         \
  " primaryClass = \"astro-ph.IM\",\n"                                  \
  " keywords = {galaxies: irregular, galaxies: photometry, "            \
  "galaxies: structure, methods: data analysis, "                       \
  "techniques: image processing, techniques: photometric},\n"           \
  "     year = 2015,\n"                                                 \
  "    month = sep,\n"                                                  \
  "   volume = 220,\n"                                                  \
  "      eid = {1},\n"                                                  \
  "    pages = {1},\n"                                                  \
  "      doi = {10.1088/0067-0049/220/1/1},\n"                          \
  "   adsurl = {http://adsabs.harvard.edu/abs/2015ApJS..220....1A},\n"  \
  "  adsnote = {Provided by the SAO/NASA Astrophysics Data System}\n"   \
  "}"


/* This can be used in the end of error messages related to option
   values. */
#define HOWTOCHECKVALUES                                                \
  " You can check all the input values with the `--printparams' "       \
  "(-P) option."

#endif
