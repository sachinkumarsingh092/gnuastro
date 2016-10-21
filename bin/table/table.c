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
#include <unistd.h>

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
      sprintf(width, "%zu", up->sintwidth);
      accu[0]='\0';
      break;

    case TLOGICAL: case TSBYTE:
      type="d";
      sprintf(width, "%zu", up->sintwidth);
      accu[0]='\0';
      break;

    case TSTRING:
      type="s";
      if(up->sintwidth) sprintf(width, "%zu", up->sintwidth);
      else width[0]='\0';
      accu[0]='\0';
      break;

    case TSHORT:
      type="d";
      sprintf(width, "%zu", up->sintwidth);
      accu[0]='\0';
      break;

    case TLONG:
      type="ld";
      sprintf(width, "%zu", up->lintwidth);
      accu[0]='\0';
      break;

    case TLONGLONG:
      type="ld";
      sprintf(width, "%zu", up->lintwidth);
      accu[0]='\0';
      break;

    case TFLOAT:
      type = up->feg=='f' ? "f" : ( up->feg=='e' ? "e" : "g");
      sprintf(width, "%zu", up->floatwidth);
      sprintf(accu, ".%zu", up->floatprecision);
      break;

    case TDOUBLE:
      type = up->feg=='f' ? "f" : ( up->feg=='e' ? "e" : "g");
      sprintf(width, "%zu", up->doublewidth);
      sprintf(accu, ".%zu", up->doubleprecision);
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
      sprintf(width, "%zu", up->sintwidth);
      accu[0]='\0';
      break;

    case TUINT:
      type="u";
      sprintf(width, "%zu", up->sintwidth);
      accu[0]='\0';
      break;

    case TUSHORT:
      type="u";
      sprintf(width, "%zu", up->sintwidth);
      accu[0]='\0';
      break;

    case TULONG:
      type="lu";
      sprintf(width, "%zu", up->lintwidth);
      accu[0]='\0';
      break;

    default:
      error(EXIT_FAILURE, 0, "datatype value of %d not recognized in "
            "setformatstring (table.c)", ocol->datatype);
    }

  /* Put the type, width and accu into the format string for this
     column: */
  sprintf(ocol->fmt, "%%-%s%s%s", width, accu, type);
}





/* Read all the input columns */
void
readinputcols(struct tableparams *p)
{
  double *colfromtxt;
  struct outcolumn *col;
  int datatype, status=0;
  size_t i, j, nrows=p->nrows, incols=p->up.ncols;

  /* Get the contents of each table column: */
  for(i=0;i<p->nocols;++i)
    {
      /* Variables for simple reading */
      col=&p->ocols[i];

      datatype=col->datatype;

      /* Allocate the blank value for this column. Note that we will also
         need the blankvalue for a text file when outputing to a FITS. */
      col->nulval=gal_fits_datatype_blank(datatype);

      /* Read the input column. */
      if(p->fitsptr)
        {
          /* Allocate space for the data in this column */
          col->data=gal_fits_datatype_alloc(nrows, datatype);

          /* Call CFITSIO to read the column information. */
          fits_read_col(p->fitsptr, datatype, col->inindex+1, 1, 1,
                        nrows, col->nulval, col->data, &col->anynul,
                        &status);
        }
      else
        {
          /* This is a text file, read by Gnuastro's current txtarray
             library. This library currently only reads a 2D table into a
             2D array of type double. The important thing here is that the
             array is row-contiguous. But here we want column contiguous
             data. So we allocate an array to only put this column's values
             in.*/
          errno=0;
          colfromtxt=col->data=malloc(nrows * col->esize);
          if(col->data==NULL)
            error(EXIT_FAILURE, errno, "%zu bytes for col->data",
                  nrows * col->esize);

          for(j=0;j<nrows;++j)
            colfromtxt[j]=p->up.txtarray[ j * incols + col->inindex ];
        }

      /* Set the format string to print the column values if the output is
         to be printed as text (either in a text file or to the standard
         output. */
      if(p->outputtotxt || p->outputtostdout)
        setformatstring(p, i);
    }
}




















/**************************************************************/
/***************       Output table         *******************/
/**************************************************************/
void
saveouttofits(struct tableparams *p)
{
  size_t i;
  int status=0;
  fitsfile *fptr;
  struct uiparams *up=&p->up;
  char **ttype, **tform, **tunit;
  struct outcolumn *ocols=p->ocols;

  /* Allocate the information arrays for CFITSIO. */
  errno=0;
  ttype=malloc(p->nocols*sizeof *ttype);
  if(ttype==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for ttype",
          p->nocols*sizeof *ttype);
  errno=0;
  tform=malloc(p->nocols*sizeof *tform);
  if(tform==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for tform",
          p->nocols*sizeof *tform);
  errno=0;
  tunit=malloc(p->nocols*sizeof *tunit);
  if(tunit==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for tunit",
          p->nocols*sizeof *tunit);

  /* Fill in the information arrays: */
  for(i=0;i<p->nocols;++i)
    {
      tform[i]=up->ttstr[ ocols[i].inindex ];
      ttype[i]=up->tname[ ocols[i].inindex ];
      tunit[i]=up->tunit[ ocols[i].inindex ];
    }

  /* Open the output FITS file. */
  if(access(p->cp.output,F_OK) != -1 )
    fits_open_file(&fptr, p->cp.output, READWRITE, &status);
  else
    fits_create_file(&fptr, p->cp.output, &status);

  /* Create a new table extension */
  fits_create_tbl(fptr, p->fitstabletype, p->nrows, p->nocols, ttype,
                  tform, tunit, "Table", &status);

  /* Write this column's data into the FITS file. */
  for(i=0;i<p->nocols;++i)
    fits_write_colnull(fptr, ocols[i].datatype, i+1, 1, 1, p->nrows,
                       ocols[i].data, ocols[i].nulval, &status);

  /* Include the ending comments and close the file. */
  gal_fits_write_keys_version(fptr, NULL, SPACK_STRING);
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);

  /* Clean up */
  free(ttype);
  free(tform);
  free(tunit);
}





void
printoutput(struct tableparams *p)
{
  FILE *fp;
  size_t i, row;
  struct outcolumn *ocols=p->ocols;

  /* Determine the output stream and open the file for writing if its a
     file. */
  if(p->outputtotxt)
    {
      errno=0;
      fp=fopen(p->cp.output, "w");
      if(fp==NULL)
        error(EXIT_FAILURE, errno, "%s", p->cp.output);
    }
  else
    fp=stdout;

  /* Print each column and each row: */
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
            fprintf(fp, ocols[i].fmt, ((unsigned char *)ocols[i].data)[row]);
            break;

          case TLOGICAL: case TSBYTE:
            fprintf(fp, ocols[i].fmt, ((char *)ocols[i].data)[row]);
            break;

          case TSTRING:
            fprintf(fp, ocols[i].fmt, ((char **)ocols[i].data)[row]);
            break;

          case TSHORT:
            fprintf(fp, ocols[i].fmt, ((short *)ocols[i].data)[row]);
            break;

          case TLONG:
            fprintf(fp, ocols[i].fmt, ((long *)ocols[i].data)[row]);
            break;

          case TLONGLONG:
            fprintf(fp, ocols[i].fmt, ((LONGLONG *)ocols[i].data)[row]);
            break;

          case TFLOAT:
            fprintf(fp, ocols[i].fmt, ((float *)ocols[i].data)[row]);
            break;

          case TDOUBLE:
            fprintf(fp, ocols[i].fmt, ((double *)ocols[i].data)[row]);
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
            fprintf(fp, ocols[i].fmt, ((char *)ocols[i].data)[row]);
            break;

          case TUINT:
            fprintf(fp, ocols[i].fmt, ((unsigned int *)ocols[i].data)[row]);
            break;

          case TUSHORT:
            fprintf(fp, ocols[i].fmt, ((unsigned short *)ocols[i].data)[row]);
            break;

          case TULONG:
            fprintf(fp, ocols[i].fmt, ((unsigned long *)ocols[i].data)[row]);
            break;

          default:
            error(EXIT_FAILURE, 0, "datatype value of %d not recognized in "
                  "printoutput", ocols[i].datatype);
          }
      fprintf(fp, "\n");
    }

  /* If we printed to a file, then close it. */
  if(p->outputtotxt)
    {
      errno=0;
      if( fclose(fp) == EOF )
        error(EXIT_FAILURE, errno, "%s", p->cp.output);
    }
}



















/**************************************************************/
/***************       Top function         *******************/
/**************************************************************/
void
table(struct tableparams *p)
{
  readinputcols(p);

  if(p->outputtofits)
    saveouttofits(p);
  else
    printoutput(p);
}
