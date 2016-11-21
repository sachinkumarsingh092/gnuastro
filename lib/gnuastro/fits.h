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

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <limits.h>
#include <fitsio.h>
#include <wcslib/wcs.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>
#include <gsl/gsl_complex.h>

#include <gnuastro/data.h>

/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */



/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */





/*************************************************************
 ******************         Basic          *******************
 *************************************************************/
void
gal_fits_io_error(int status, char *message);

int
gal_fits_name_is_fits(char *name);

int
gal_fits_suffix_is_fits(char *suffix);





/*************************************************************
 ******************         Header          ******************
 *************************************************************/
/* To create a linked list of headers. */
struct gal_fits_key_ll
{
  int                    kfree;   /* ==1, free keyword name.   */
  int                    vfree;   /* ==1, free keyword value.  */
  int                    cfree;   /* ==1, free comment.        */
  int                 datatype;   /* Keyword value datatype.   */
  char                *keyname;   /* Keyword Name.             */
  void                  *value;   /* Keyword value.            */
  char                *comment;   /* Keyword comment.          */
  char                   *unit;   /* Keyword unit.             */
  struct gal_fits_key_ll *next;   /* Pointer next keyword.     */
};





struct gal_fits_key
{
  int            status;        /* CFITSIO status.        */
  char         *keyname;        /* Name of keyword.       */
  int          datatype;        /* Type of keyword value. */
  char  str[FLEN_VALUE];        /* String value.          */
  unsigned char       u;        /* Byte value.            */
  short               s;        /* Short integer value.   */
  long                l;        /* Long integer value.    */
  LONGLONG            L;        /* Long Long value.       */
  float               f;        /* Float value.           */
  double              d;        /* Double value.          */
};





void
gal_fits_read_keywords(char *filename, char *hdu, struct gal_fits_key *out,
                       size_t num);

void
gal_fits_add_to_key_ll(struct gal_fits_key_ll **list, int datatype,
                       char *keyname, int kfree, void *value, int vfree,
                       char *comment, int cfree, char *unit);

void
gal_fits_add_to_key_ll_end(struct gal_fits_key_ll **list, int datatype,
                           char *keyname, int kfree, void *value, int vfree,
                           char *comment, int cfree, char *unit);

void
gal_fits_file_name_in_keywords(char *keynamebase, char *filename,
                               struct gal_fits_key_ll **list);

void
gal_fits_add_wcs_to_header(fitsfile *fptr, char *wcsheader, int nkeyrec);

void
gal_fits_update_keys(fitsfile *fptr, struct gal_fits_key_ll **keylist);

void
gal_fits_write_keys_version(fitsfile *fptr, struct gal_fits_key_ll *headers,
                            char *spack_string);





/*************************************************************
 ******************        Read/Write        *****************
 *************************************************************/
void *
gal_fits_datatype_blank(int datatype);

void
gal_fits_convert_blank(void *array, int bitpix, size_t size, void *value);

void
gal_fits_blank_to_value(void *array, int datatype, size_t size, void *value);

int
gal_fits_bitpix_to_datatype(int bitpix);

void
gal_fits_img_bitpix_size(fitsfile *fptr, int *bitpix, long *naxis);

void
gal_fits_read_hdu(char *filename, char *hdu, unsigned char img0_tab1,
                  fitsfile **outfptr);

void *
gal_fits_datatype_alloc(size_t size, int datatype);

size_t
gal_fits_datatype_size(int datatype);

void
gal_fits_change_type(void *in, int inbitpix, size_t size, int anyblank,
                     void **out, int outbitpix);

void
gal_fits_num_hdus(char *filename, int *numhdu);

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
                       struct wcsprm *wcs, struct gal_fits_key_ll *headers,
                       char *spack_string);

void
gal_fits_atof_correct_wcs(char *filename, char *hdu, int bitpix,
                          void *array, size_t s0, size_t s1,
                          char *wcsheader, int wcsnkeyrec,
                          double *crpix, char *spack_string);





/**************************************************************/
/**********                  Table                 ************/
/**************************************************************/
int
gal_fits_tform_to_datatype(char tform);

void
gal_fits_table_size(fitsfile *fitsptr, size_t *nrows, size_t *ncols);

int
gal_fits_table_type(fitsfile *fptr);





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



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_FITS_H__ */
