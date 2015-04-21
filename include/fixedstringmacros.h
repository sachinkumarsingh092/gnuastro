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
  "full "PACKAGE_NAME" manual.\n"                                       \
  "      info "SPACK"\n"                                                \
  "      info "SPACK_NAME"\n"						\
  "      info "PACKAGE_TARNAME"\n"                                      \
  "If you couldn't find your answer in the manual, you can contact "    \
  "`help-gnuastro@gnu.org' for direct help from experienced Gnuastro "  \
  "users and developers.\n"                                             \
  "\nOptions:"								\



#define ASTRUTILSBIBTEX 						\
    "@ARTICLE{gnuastro,\n"						\
    "   author = {{Akhlaghi}, M. and {Ichikawa}, T.},\n"		\
    "    title = \"{Noise based detection and segmentation of nebulous" \
    " objects}\"\n"							\
    "  journal = {\\apjs},\n"						\
    " keywords = {galaxies: irregular, galaxies: photometry, galaxies: " \
    "structure, methods: data analysis, techniques: image processing, "	\
    "techniques: photometric},\n"					\
    "     year = 2015,\n"						\
    "    month = XXX,\n"						\
    "   volume = XXX,\n"						\
    "      eid = {XX},\n"						\
    "    pages = {XX},\n"						\
    "      doi = {XXXXXXXXX},\n"					\
    "   adsurl = {XXXXXXXXXXXXXXXXXXXXXX},\n"				\
    "}"									\

#endif
