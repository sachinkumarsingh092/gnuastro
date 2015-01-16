/*********************************************************************
mkprof (MakeProfiles) - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
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
#ifndef CITE_H
#define CITE_H

#define MKPROFBIBTEX ""

#define PRINTCITEABORT {						\
    printf("\nWe hope %s has been useful for your research.\n"		\
	   "Citations are vital for the continued work on %s.\n"	\
	   "Thank you for citing it in your research paper.\n"		\
	   "\nPlease cite as \"%s\":\n\n%s\n\n%s",			\
	   SPACK_NAME, SPACK_NAME, SPACK_STRING, ASTRUTILSBIBTEX,	\
	   MKPROFBIBTEX);						\
    exit(EXIT_SUCCESS);							\
}

#endif
