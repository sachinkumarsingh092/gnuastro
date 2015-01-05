/*********************************************************************
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef COPYRIGHT_H
#define COPYRIGHT_H

#define SHORTCOPYRIGHT "Copyright (C) 2013-2015, Free Software Foundation."
#define SHORTLICENSE   "License GPLv3+: GNU General public license version 3 or later."


#define COPYRIGHT SHORTCOPYRIGHT"\n"SHORTLICENSE"\n"			\
  "This is free software: you are free to change and redistribute it.\n" \
  "There is NO WARRANTY, to the extent permitted by law."		\

#define TOPHELPINFO   "\n"SPACK_NAME" is part of "PACKAGE_STRING".\n"

/* This is fixed for all the packages. */
#define MOREHELPINFO   "\nFor more information, please run any of the "	\
  "following commands. They will respectively show you the `invoking "	\
  SPACK_NAME"` section of the manual, the full "SPACK_NAME" section "	\
  "of the manual, and the full "PACKAGE_NAME" manual.\n\n"		\
  "  info "SPACK"\n"							\
  "  info "SPACK_NAME"\n"						\
  "  info "PACKAGE_TARNAME"\n\n"					\
  SPACK_NAME" was configured on this system at "CONFIGDATE", "		\
  CONFIGTIME".\n\n"							\
  "Options:"								\

#endif
