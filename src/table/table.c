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



/**************************************************************/
/***************        Input table         *******************/
/**************************************************************/
/* Set the format string for each column: */
void
setformatstring(struct tableparams *p, size_t outcolid)
{
  struct uiparams *up=&p->up;
  char width[10], accu[10], *type=NULL;
  struct outcolumn *ocol=&p->ocols[outcolid];


  switch(ocol->datatype)
    {
    case TBIT:
      error(EXIT_FAILURE, 0, "Table doesn't print TBIT data type "
            "currently, please contact us at %s so we can implement "
            "it.", PACKAGE_BUGREPORT);

    case TBYTE:
      type="u";
      sprintf(width, "%lu", up->sintwidth);
      accu[0]='\0';
      break;

    case TLOGICAL: case TSBYTE:
      type="d";
      sprintf(width, "%lu", up->sintwidth);
      accu[0]='\0';
      break;

    case TSTRING:
      type="s";
      if(up->sintwidth) sprintf(width, "%lu", up->sintwidth);
      else width[0]='\0';
      accu[0]='\0';
      break;

    case TSHORT:
      type="d";
      sprintf(width, "%lu", up->sintwidth);
      accu[0]='\0';
      break;

    case TLONG:
      type="ld";
      sprintf(width, "%lu", up->lintwidth);
      accu[0]='\0';
      break;

    case TLONGLONG:
      type="ld";
      sprintf(width, "%lu", up->lintwidth);
      accu[0]='\0';
      break;

    case TFLOAT:
      type = up->feg=='f' ? "f" : ( up->feg=='e' ? "e" : "g");
      sprintf(width, "%lu", up->floatwidth);
      sprintf(accu, ".%lu", up->floatprecision);
      break;

    case TDOUBLE:
      type = up->feg=='f' ? "f" : ( up->feg=='e' ? "e" : "g");
      sprintf(width, "%lu", up->doublewidth);
      sprintf(accu, ".%lu", up->doubleprecision);
      break;

    case TCOMPLEX:
      error(EXIT_FAILURE, 0, "Table doesn't print TCOMPLEX data type "
            "currently, please contact us at %s so we can implement "
            "it.", PACKAGE_BUGREPORT);
      break;

    case TDBLCOMPLEX:
      error(EXIT_FAILURE, 0, "Table doesn't print TDBLCOMPLEX data type "
            "currently, please contact us at %s so we can implement "
            "it.", PACKAGE_BUGREPORT);
      break;

    case TINT:
      type="d";
      sprintf(width, "%lu", up->sintwidth);
      accu[0]='\0';
      break;

    case TUINT:
      type="u";
      sprintf(width, "%lu", up->sintwidth);
      accu[0]='\0';
      break;

    case TUSHORT:
      type="u";
      sprintf(width, "%lu", up->sintwidth);
      accu[0]='\0';
      break;

    case TULONG:
      type="lu";
      sprintf(width, "%lu", up->lintwidth);
      accu[0]='\0';
      break;

    default:
      error(EXIT_FAILURE, 0, "datatype value of %d not recognized in "
            "gal_fits_datatype_alloc", ocol->datatype);
    }

  /* Put the type, width and accu into the format string for this
     column: */
  sprintf(ocol->fmt, "%%-%s%s%s", width, accu, type);
}





/* Read all the input columns */
void
readinputcols(struct tableparams *p)
{
  size_t i;
  void *nulval;
  struct outcolumn *col;
  int datatype, status=0;

  /* Get the contents of each table column: */
  for(i=0;i<p->nocols;++i)
    {
      /* Variables for simple reading */
      col=&p->ocols[i];
      datatype=col->datatype;

      /* Allocate the blank value and array to keep the actual of this
         column. */
      nulval=gal_fits_datatype_blank(datatype);
      col->data=gal_fits_datatype_alloc(p->nrows, datatype);

      /* Call CFITSIO to read the column information. */
      fits_read_col(p->fitsptr, datatype, col->inindex+1, 1, 1,
                    p->nrows, nulval, col->data, &col->anynul,
                    &status);

      /* Free the space allocated for the blank value, we don't need it any
         more: it is an internal macro to Gnuastro (see `fits.h'), */
      free(nulval);

      /* Set the format string to print the column values. */
      setformatstring(p, i);
    }
}




















/**************************************************************/
/***************       Output table         *******************/
/**************************************************************/
void
printoutput(struct tableparams *p)
{
  size_t i, row;
  struct outcolumn *ocols=p->ocols;

  for(row=0;row<p->nrows;++row)
    {
      for(i=0;i<p->nocols;++i)
        switch(ocols[i].datatype)
          {
          case TBIT:
            error(EXIT_FAILURE, 0, "Table doesn't print TBIT data type "
                  "currently, please contact us at %s so we can implement "
                  "it.", PACKAGE_BUGREPORT);

          case TBYTE:
            printf(ocols[i].fmt, ((unsigned char *)ocols[i].data)[row]);
            break;

          case TLOGICAL: case TSBYTE:
            printf(ocols[i].fmt, ((char *)ocols[i].data)[row]);
            break;

          case TSTRING:
            printf(ocols[i].fmt, ((char **)ocols[i].data)[row]);
            break;

          case TSHORT:
            printf(ocols[i].fmt, ((short *)ocols[i].data)[row]);
            break;

          case TLONG:
            printf(ocols[i].fmt, ((long *)ocols[i].data)[row]);
            break;

          case TLONGLONG:
            printf(ocols[i].fmt, ((LONGLONG *)ocols[i].data)[row]);
            break;

          case TFLOAT:
            printf(ocols[i].fmt, ((float *)ocols[i].data)[row]);
            break;

          case TDOUBLE:
            printf(ocols[i].fmt, ((double *)ocols[i].data)[row]);
            break;

          case TCOMPLEX:
            error(EXIT_FAILURE, 0, "Table doesn't print TCOMPLEX data type "
                  "currently, please contact us at %s so we can implement "
                  "it.", PACKAGE_BUGREPORT);
            break;

          case TDBLCOMPLEX:
            error(EXIT_FAILURE, 0, "Table doesn't print TDBLCOMPLEX data "
                  "type currently, please contact us at %s so we can "
                  "implement it.", PACKAGE_BUGREPORT);
            break;

          case TINT:
            printf(ocols[i].fmt, ((char *)ocols[i].data)[row]);
            break;

          case TUINT:
            printf(ocols[i].fmt, ((unsigned int *)ocols[i].data)[row]);
            break;

          case TUSHORT:
            printf(ocols[i].fmt, ((unsigned short *)ocols[i].data)[row]);
            break;

          case TULONG:
            printf(ocols[i].fmt, ((unsigned long *)ocols[i].data)[row]);
            break;

          default:
            error(EXIT_FAILURE, 0, "datatype value of %d not recognized in "
                  "printoutput", ocols[i].datatype);
          }
      printf("\n");
    }
}



















/**************************************************************/
/***************       Top function         *******************/
/**************************************************************/
void
table(struct tableparams *p)
{
  readinputcols(p);

  printoutput(p);
}
