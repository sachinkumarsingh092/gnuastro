/*********************************************************************
Functions to convert a FITS array to a C array and vice versa.
This is part of GNU Astronomy Utilities (AstrUtils) package.

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
#ifndef FITSMATRIX_H
#define FITSMATRIX_H

#include <math.h>
#include <float.h>
#include <stdint.h>

#include <fitsio.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>
#include <wcslib/wcs.h>

#define FITSBYTENUL     UINT8_MAX	/* 0 is often meaningful here! */
#define FITSSHORTNUL    INT16_MIN
#define FITSLONGNUL     INT32_MIN
#define FITSLLONGNUL    INT64_MIN
#define FITSFLOATNUL    NAN




/*************************************************************
 ******************         Basic          *******************
 *************************************************************/
void
fitsioerror(int status, char *message);

int
nameisfits(char *name);

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
add_to_fitsheaderllend(struct fitsheaderll **list, int datatype,
		       char *keyname, int kfree, void *value, int vfree,
		       char *comment, int cfree, char *unit);
void
filenameinkeywords(char *keynamebase, char *filename,
		   struct fitsheaderll **list);
void
addwcstoheader(fitsfile *fptr, struct wcsprm *wcs);

void
updatekeys(fitsfile *fptr, struct fitsheaderll **keylist);

void
copyrightandend(fitsfile *fptr, char *spack_string);





/*************************************************************
 ******************        Read/Write        *****************
 *************************************************************/
void *
bitpixnull(int bitpix);

void
convertnul(void *array, int bitpix, size_t size, void *value);

int
bitpixtodtype(int bitpix);

void
imgbitpixsize(fitsfile *fptr, int *bitpix, long *naxis);

void
readfitshdu(char *filename, char *hdu, int desiredtype, fitsfile **outfptr);

void *
bitpixalloc(size_t size, int bitpix);

void
changetype(void *in, int inbitpix, size_t size, void **out, int outbitpix);

void
readwcs(fitsfile *fptr, int *nwcs, struct wcsprm **wcs);

int
fitsimgtoarray(char *filename, char *hdu, void *bitnul, int *bitpix,
	       void **array, size_t *s0, size_t *s1);

void
arraytofitsimg(char *filename, char *hdu, int bitpix, void *array,
	       size_t s0, size_t s1, struct wcsprm *wcs,
	       char *spack_string);


#endif
