/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef AUTHORS_CITE_H
#define AUTHORS_CITE_H

/* When any specific citation is necessary, please add its BibTeX (from ADS
   hopefully) to this variable along with a title decribing what this
   paper/book does for the progarm in a short line. In the following line
   put a row of '-' with the same length and then put the BibTeX.

   This macro will be used in 'gal_options_print_citation' function of
   'lib/options.c' (from the top Gnuastro source code directory). */

#define PROGRAM_BIBTEX                                                  \
  "Description of MakeCatalog\n"                                        \
  "--------------------------\n"                                        \
  "@ARTICLE{makecatalog,\n"                                             \
  "       author = {{Akhlaghi}, Mohammad},\n"                           \
  "        title = \"{Separating Detection and Catalog Production}\",\n" \
  "      journal = {ASPC},\n"                                           \
  "         year = \"2019\",\n"                                         \
  "        month = \"Oct\",\n"                                          \
  "       volume = {521},\n"                                            \
  "        pages = {299},\n"                                            \
  "archivePrefix = {arXiv},\n"                                          \
  "       eprint = {1611.06387},\n"                                     \
  " primaryClass = {astro-ph.IM},\n"                                    \
  "       adsurl = {https://ui.adsabs.harvard.edu/abs/2019ASPC..521..299A},\n" \
  "      adsnote = {Provided by the SAO/NASA Astrophysics Data System}\n" \
  "}\n"

#define PROGRAM_AUTHORS "Mohammad Akhlaghi"

#endif
