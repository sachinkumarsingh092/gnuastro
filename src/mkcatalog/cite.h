/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef CITE_H
#define CITE_H

#define HEADERBIBTEX ""

#define PRINTCITEABORT {						\
    printf("\nWe hope %s has been useful for your research.\n"		\
	   "Citations are vital for the continued work on %s.\n"	\
	   "Thank you for citing it in your research paper.\n"		\
	   "\nPlease cite as \"%s\":\n\n%s\n\n%s",			\
	   SPACK_NAME, SPACK_NAME, SPACK_STRING, GNUASTROBIBTEX,	\
	   HEADERBIBTEX);						\
    exit(EXIT_SUCCESS);							\
}

#endif
