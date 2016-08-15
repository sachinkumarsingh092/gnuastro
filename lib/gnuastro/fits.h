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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef __GAL_FITS_H__
#define __GAL_FITS_H__

#include <math.h>
#include <float.h>
#include <stdint.h>

#include <fitsio.h>
#include <wcslib/wcs.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>

#define GAL_FITS_STRING_BLANK   NULL
#define GAL_FITS_BYTE_BLANK     UCHAR_MAX /* 0 is often meaningful here! */
#define GAL_FITS_SHORT_BLANK    INT16_MIN
#define GAL_FITS_LONG_BLANK     INT32_MIN
#define GAL_FITS_LLONG_BLANK    INT64_MIN
#define GAL_FITS_FLOAT_BLANK    NAN





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
gal_fits_io_error(int status, char *message);

int
gal_fits_name_is_fits(char *name);

int
gal_fits_name_is_fits_suffix(char *name);

void
gal_fits_num_hdus(char *filename, int *numhdu);





/*************************************************************
 ******************         Header          ******************
 *************************************************************/
/* To create a linked list of headers. */
struct gal_fits_header_ll
{
  int                 kfree;   /* ==1, keyname will be freed.          */
  int                 vfree;   /* ==1, value will be freed.            */
  int                 cfree;   /* ==1, comment will be freed.          */
  int              datatype;   /* Data type of the keyword             */
  char             *keyname;   /* Name of keyword.                     */
  void               *value;   /* Pointer to the value of the keyword. */
  char             *comment;   /* Comment for the keyword.             */
  char                *unit;   /* Units of the keyword.                */
  struct gal_fits_header_ll *next;   /* Pointer to the next element.         */
};





struct gal_fits_read_header_keys
{
  int            status;
  char         *keyname;
  int          datatype;
  char  str[FLEN_VALUE];
  unsigned char       u;
  short               s;
  long                l;
  LONGLONG            L;
  float               f;
  double              d;
};





void
gal_fits_read_keywords(char *filename, char *hdu,
                       struct gal_fits_read_header_keys *out,
                       size_t num);

void
gal_fits_add_to_fits_header_ll(struct gal_fits_header_ll **list,
                               int datatype, char *keyname, int kfree,
                               void *value, int vfree, char *comment,
                               int cfree, char *unit);

void
gal_fits_add_to_fits_header_ll_end(struct gal_fits_header_ll **list,
                                   int datatype, char *keyname, int kfree,
                                   void *value, int vfree, char *comment,
                                   int cfree, char *unit);

void
gal_fits_file_name_in_keywords(char *keynamebase, char *filename,
                               struct gal_fits_header_ll **list);

void
gal_fits_add_wcs_to_header(fitsfile *fptr, char *wcsheader, int nkeyrec);

void
gal_fits_update_keys(fitsfile *fptr,
                     struct gal_fits_header_ll **keylist);

void
gal_fits_copyright_end(fitsfile *fptr,
                       struct gal_fits_header_ll *headers,
                       char *spack_string);





/*************************************************************
 ******************        Read/Write        *****************
 *************************************************************/
void *
gal_fits_bitpix_blank(int bitpix);

void
gal_fits_convert_blank(void *array, int bitpix, size_t size, void *value);

int
gal_fits_bitpix_to_dtype(int bitpix);

void
gal_fits_img_bitpix_size(fitsfile *fptr, int *bitpix, long *naxis);

void
gal_fits_read_hdu(char *filename, char *hdu, int desiredtype,
                  fitsfile **outfptr);

void *
gal_fits_bitpix_alloc(size_t size, int bitpix);

void
gal_fits_change_type(void *in, int inbitpix, size_t size, int anyblank,
                     void **out, int outbitpix);

void
gal_fits_read_wcs_from_pointer(fitsfile *fptr, int *nwcs,
                               struct wcsprm **wcs,
                               size_t hstart, size_t hend);

void
gal_fits_read_wcs(char *filename, char *hdu, size_t hstartwcs,
                  size_t hendwcs, int *nwcs, struct wcsprm **wcs);

int
gal_fits_hdu_to_array(char *filename, char *hdu, int *bitpix,
                      void **array, size_t *s0, size_t *s1);

void
gal_fits_array_to_file(char *filename, char *hdu, int bitpix,
                       void *array, size_t s0, size_t s1, int anyblank,
                       struct wcsprm *wcs,
                       struct gal_fits_header_ll *headers,
                       char *spack_string);

void
gal_fits_atof_correct_wcs(char *filename, char *hdu, int bitpix,
                          void *array, size_t s0, size_t s1,
                          char *wcsheader, int wcsnkeyrec,
                          double *crpix, char *spack_string);





/**************************************************************/
/**********          Check prepare file            ************/
/**************************************************************/
void
gal_fits_file_or_ext_name(char *inputname, char *inhdu, int othernameset,
                          char **othername, char *ohdu, int ohduset,
                          char *type);

void
gal_fits_set_mask_name(char *inputname, char **maskname, char *inhdu,
                       char *mhdu);

void
gal_fits_file_to_double(char *inputname, char *maskname, char *inhdu,
                        char *mhdu, double **img, int *inbitpix,
                        int *anyblank, size_t *ins0, size_t *ins1);

void
gal_fits_file_to_float(char *inputname, char *maskname, char *inhdu,
                       char *mhdu, float **img, int *inbitpix,
                       int *anyblank, size_t *ins0, size_t *ins1);

void
gal_fits_file_to_long(char *inputname, char *inhdu, long **img,
                      int *inbitpix, int *anyblank, size_t *ins0,
                      size_t *ins1);

void
gal_fits_prep_float_kernel(char *inputname, char *inhdu, float **kernel,
                           size_t *ins0, size_t *ins1);

#endif
