/*********************************************************************
Functions to convert a FITS array to a C array and vice versa.
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
#ifndef FITSMATRIX_H
#define FITSMATRIX_H

#include <math.h>
#include <float.h>
#include <stdint.h>

#include <fitsio.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>
#include <wcslib/wcs.h>

#define FITSBYTEBLANK     UCHAR_MAX	/* 0 is often meaningful here! */
#define FITSSHORTBLANK    INT16_MIN
#define FITSLONGBLANK     INT32_MIN
#define FITSLLONGBLANK    INT64_MIN
#define FITSFLOATBLANK    NAN




/*

For some reason, CFITSIO does not use the standard stdint fixed size
types! It uses the subjective 'short', 'int' and 'long' variables
which can differ in size from system to system!!!!!!!!!!!!!!!

In the 32bit systems that 'long' was 32 bits or 4 bytes, has passed
but the names have stuck! The FITS standard defines LONG_IMG as a
32bit signed type, but CFITSIO converts it to a local 'long' which is
64 bits on a modern (64 bit) system!!!! This is simply absurd and very
confusing!!!! It should have stuck to the standard, not the name of
the variable!

Because of this we have to stick to this wrong convention too.

 */


/*************************************************************
 ******************         Basic          *******************
 *************************************************************/
void
fitsioerror(int status, char *message);

int
nameisfits(char *name);

int
nameisfitssuffix(char *name);

void
numhduensions(char *filename, int *numhdu);





/*************************************************************
 ******************         Header          ******************
 *************************************************************/
/* To create a linked list of headers. */
struct fitsheaderll
{
  int                 kfree;   /* ==1, keyname will be freed.          */
  int                 vfree;   /* ==1, value will be freed.            */
  int                 cfree;   /* ==1, comment will be freed.          */
  int              datatype;   /* Data type of the keyword             */
  char             *keyname;   /* Name of keyword.                     */
  void               *value;   /* Pointer to the value of the keyword. */
  char             *comment;   /* Comment for the keyword.             */
  char                *unit;   /* Units of the keyword.                */
  struct fitsheaderll *next;   /* Pointer to the next element.         */
};

void
add_to_fitsheaderll(struct fitsheaderll **list, int datatype,
                    char *keyname, int kfree, void *value, int vfree,
                    char *comment, int cfree, char *unit);

void
add_to_fitsheaderllend(struct fitsheaderll **list, int datatype,
		       char *keyname, int kfree, void *value, int vfree,
		       char *comment, int cfree, char *unit);
void
filenameinkeywords(char *keynamebase, char *filename,
		   struct fitsheaderll **list);

void
addwcstoheader(fitsfile *fptr, char *wcsheader, int nkeyrec);

void
updatekeys(fitsfile *fptr, struct fitsheaderll **keylist);

void
copyrightandend(fitsfile *fptr, struct fitsheaderll *headers,
                char *spack_string);





/*************************************************************
 ******************        Read/Write        *****************
 *************************************************************/
void *
bitpixblank(int bitpix);

void
convertblank(void *array, int bitpix, size_t size, void *value);

int
bitpixtodtype(int bitpix);

void
imgbitpixsize(fitsfile *fptr, int *bitpix, long *naxis);

void
readfitshdu(char *filename, char *hdu, int desiredtype, fitsfile **outfptr);

void *
bitpixalloc(size_t size, int bitpix);

void
changetype(void *in, int inbitpix, size_t size, size_t numblank,
           void **out, int outbitpix);

void
readwcs(fitsfile *fptr, int *nwcs, struct wcsprm **wcs);

void
readfitswcs(char *filename, char *hdu, int *nwcs, struct wcsprm **wcs);

size_t
fitsimgtoarray(char *filename, char *hdu, int *bitpix, void **array,
               size_t *s0, size_t *s1);

void
arraytofitsimg(char *filename, char *hdu, int bitpix, void *array,
	       size_t s0, size_t s1, size_t numblank, struct wcsprm *wcs,
	       struct fitsheaderll *headers, char *spack_string);

void
atofcorrectwcs(char *filename, char *hdu, int bitpix, void *array,
	       size_t s0, size_t s1, char *wcsheader, int wcsnkeyrec,
	       double *crpix, char *spack_string);


#endif
