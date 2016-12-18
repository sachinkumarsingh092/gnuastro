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

/* When we are within Gnuastro's building process, `IN_GNUASTRO_BUILD' is
   defined. In the build process, installation information (in particular
   `GAL_CONFIG_HAVE_WCSLIB_VERION' that we need in `fits.c') is kept in
   `config.h'. When building a user's programs, this information is kept in
   `gnuastro/config.h'. Note that all `.c' files must start with the
   inclusion of `config.h' and that `gnuastro/config.h' is only created at
   installation time (not present during the building of Gnuastro).*/
#ifndef IN_GNUASTRO_BUILD
#include <gnuastro/config.h>
#endif

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <math.h>
#include <float.h>

#include <fitsio.h>
#include <wcslib/wcs.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>

#include <gnuastro/data.h>
#include <gnuastro/table.h>
#include <gnuastro/linkedlist.h>

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





/* To create a linked list of headers. */
struct gal_fits_key_ll
{
  int                    kfree;   /* ==1, free keyword name.   */
  int                    vfree;   /* ==1, free keyword value.  */
  int                    cfree;   /* ==1, free comment.        */
  int                     type;   /* Keyword value type.       */
  char                *keyname;   /* Keyword Name.             */
  void                  *value;   /* Keyword value.            */
  char                *comment;   /* Keyword comment.          */
  char                   *unit;   /* Keyword unit.             */
  struct gal_fits_key_ll *next;   /* Pointer next keyword.     */
};





/*************************************************************
 **************        Reporting errors:       ***************
 *************************************************************/
void
gal_fits_io_error(int status, char *message);





/*************************************************************
 **************           FITS names           ***************
 *************************************************************/
int
gal_fits_name_is_fits(char *name);

int
gal_fits_suffix_is_fits(char *suffix);

void
gal_fits_file_or_ext_name(char *inputname, char *inhdu, int othernameset,
                          char **othername, char *ohdu, int ohduset,
                          char *type);



/*************************************************************
 **************           Type codes           ***************
 *************************************************************/
int
gal_fits_bitpix_to_type(int bitpix);

int
gal_fits_type_to_bitpix(int type);

int
gal_fits_tform_to_type(char tform);

int
gal_fits_type_to_datatype(int type);

int
gal_fits_datatype_to_type(int type);





/*************************************************************
 **************        Get information         ***************
 *************************************************************/
void
gal_fits_num_hdus(char *filename, int *numhdu);

void
gal_fits_img_info(fitsfile *fptr, int *type, size_t *ndim, long **dsize);





/**************************************************************/
/**********                  HDU                   ************/
/**************************************************************/

fitsfile *
gal_fits_read_hdu(char *filename, char *hdu, unsigned char img0_tab1);





/**************************************************************/
/**********            Header keywords             ************/
/**************************************************************/
void
gal_fits_read_keywords_fptr(fitsfile *fptr, gal_data_t *keysll,
                            int readcomment, int readunit);

void
gal_fits_read_keywords(char *filename, char *hdu, gal_data_t *keysll,
                       int readcomment, int readunit);

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
 ***********       Read WCS from FITS pointer      ***********
 *************************************************************/
void
gal_fits_read_wcs_from_pointer(fitsfile *fptr, int *nwcs, struct wcsprm **wcs,
                               size_t hstartwcs, size_t hendwcs);

void
gal_fits_read_wcs(char *filename, char *hdu, size_t hstartwcs,
                  size_t hendwcs, int *nwcs, struct wcsprm **wcs);





/*************************************************************
 ******************     Array functions      *****************
 *************************************************************/
gal_data_t *
gal_fits_read_img_hdu(char *filename, char *hdu, char *maskname,
                      char *mhdu, size_t minmapsize);

gal_data_t *
gal_fits_read_to_type(char *inputname, char *inhdu, char *maskname,
                      char *mhdu, int type, size_t minmapsize);

gal_data_t *
gal_fits_read_float_kernel(char *inputname, char *inhdu, float **outkernel,
                           size_t *ins0, size_t *ins1);

fitsfile *
gal_fits_write_img_fitsptr(gal_data_t *data, char *filename, char *extname);

void
gal_fits_write_img(gal_data_t *data, char *filename, char *extname,
                   struct gal_fits_key_ll *headers, char *spack_string);

void
gal_fits_write_img_update_crpix(gal_data_t *data, char *filename,
                                char *extname,
                                struct gal_fits_key_ll *headers,
                                double *crpix, char *spack_string);





/**************************************************************/
/**********                  Table                 ************/
/**************************************************************/
void
gal_fits_table_size(fitsfile *fitsptr, size_t *nrows, size_t *ncols);

int
gal_fits_table_type(fitsfile *fptr);

gal_data_t *
gal_fits_table_info(char *filename, char *hdu, size_t *numcols,
                    int *tabletype);

gal_data_t *
gal_fits_table_read(char *filename, char *hdu, gal_data_t *colinfo,
                    struct gal_linkedlist_sll *indexll, int minmapsize);

void
gal_fits_table_write(gal_data_t *cols, char *comments, int tabletype,
                     char *filename, int dontdelete);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_FITS_H__ */
