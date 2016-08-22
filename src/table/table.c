/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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
#include <config.h>

#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>

#include "main.h"




/* Print all the column information: */
void
printinfo(struct tableparams *p)
{
  size_t i;
  char *typestring=NULL;

  printf("%s (hdu: %s)\n", p->up.fitsname, p->cp.hdu);
  printf("Number of rows: %lu\n", p->nrows);
  printf("Column information\n");
  printf("------------------\n");
  for(i=1;i<=p->ncols;++i)
    {
      switch(p->typecode[i])
        {
        case TBIT:
          typestring="bit";
          break;
        case TBYTE:
          typestring="byte";
          break;
        case TLOGICAL:
          typestring="logicals";
          break;
        case TSTRING:
          typestring="string";
          break;
        case TSHORT:
          typestring="short";
          break;
        case TLONG:
          typestring="long";
          break;
        case TLONGLONG:
          typestring="longlong";
          break;
        case TFLOAT:
          typestring="float";
          break;
        case TDOUBLE:
          typestring="double";
          break;
        case TCOMPLEX:
          typestring="complex";
          break;
        case TDBLCOMPLEX:
          typestring="dblcomplex";
          break;
        case TSBYTE:
          typestring="signed byte";
          break;
        case TUINT:
          typestring="unsigned int";
          break;
        case TUSHORT:
          typestring="unsigned short";
          break;
        default:
          error(EXIT_FAILURE, 0, "%d (from TFORM%lu='%c') is not a "
                "recognized CFITSIO datatype.",
                p->typecode[i], i, p->tform[i][0]);
        }
      printf("%-3lu %-25s %-20s %s\n", i, p->ttype[i], typestring,
             p->tunit[i]);
    }

}




/* Top level function */
void
table(struct tableparams *p)
{
  if(p->information) printinfo(p);
}
