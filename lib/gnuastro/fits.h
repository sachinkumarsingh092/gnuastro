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



/* Macros. */
#define GAL_FITS_MAX_NDIM 999




/* To create a linked list of headers. */
struct gal_fits_key_ll
{
  int                    kfree;   /* ==1, free keyword name.   */
  int                    vfree;   /* ==1, free keyword value.  */
  int                    cfree;   /* ==1, free comment.        */
  uint8_t                 type;   /* Keyword value type.       */
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

char *
gal_fits_name_save_as_string(char *filename, char *hdu);



/*************************************************************
 **************           Type codes           ***************
 *************************************************************/
int
gal_fits_bitpix_to_type(int bitpix);

int
gal_fits_type_to_bitpix(uint8_t type);

char
gal_fits_type_to_bin_tform(uint8_t type);

int
gal_fits_type_to_datatype(uint8_t type);

int
gal_fits_datatype_to_type(int datatype, int is_table_column);




/**************************************************************/
/**********                  HDU                   ************/
/**************************************************************/
fitsfile *
gal_fits_open_to_write(char *filename);

void
gal_fits_hdu_num(char *filename, int *numhdu);

int
gal_fits_hdu_type(char *filename, char *hdu);

fitsfile *
gal_fits_hdu_open(char *filename, char *hdu, int iomode);

fitsfile *
gal_fits_hdu_open_type(char *filename, char *hdu, unsigned char img0_tab1);




/**************************************************************/
/**********            Header keywords             ************/
/**************************************************************/
void
gal_fits_key_clean_str_value(char *string);

void
gal_fits_key_read_from_ptr(fitsfile *fptr, gal_data_t *keysll,
                           int readcomment, int readunit);

void
gal_fits_key_read(char *filename, char *hdu, gal_data_t *keysll,
                  int readcomment, int readunit);

void
gal_fits_key_add_to_ll(struct gal_fits_key_ll **list, uint8_t type,
                       char *keyname, int kfree, void *value, int vfree,
                       char *comment, int cfree, char *unit);

void
gal_fits_key_add_to_ll_end(struct gal_fits_key_ll **list, uint8_t type,
                           char *keyname, int kfree, void *value, int vfree,
                           char *comment, int cfree, char *unit);

void
gal_fits_key_write_filename(char *keynamebase, char *filename,
                            struct gal_fits_key_ll **list);

void
gal_fits_key_write_wcsstr(fitsfile *fptr, char *wcsstr, int nkeyrec);

void
gal_fits_key_write(fitsfile *fptr, struct gal_fits_key_ll **keylist);

void
gal_fits_key_write_version(fitsfile *fptr, struct gal_fits_key_ll *headers,
                           char *program_string);





/*************************************************************
 ******************     Array functions      *****************
 *************************************************************/
void
gal_fits_img_info(fitsfile *fptr, int *type, size_t *ndim, size_t **dsize);

gal_data_t *
gal_fits_img_read(char *filename, char *hdu, size_t minmapsize);

gal_data_t *
gal_fits_img_read_to_type(char *inputname, char *inhdu, uint8_t type,
                          size_t minmapsize);

gal_data_t *
gal_fits_img_read_kernel(char *filename, char *hdu, size_t minmapsize);

fitsfile *
gal_fits_img_write_to_ptr(gal_data_t *data, char *filename);

void
gal_fits_img_write(gal_data_t *data, char *filename,
                   struct gal_fits_key_ll *headers, char *program_string);

void
gal_fits_img_write_to_type(gal_data_t *data, char *filename,
                           struct gal_fits_key_ll *headers,
                           char *program_string, int type);

void
gal_fits_img_write_corr_wcs_str(gal_data_t *data, char *filename,
                                char *wcsheader, int nkeyrec, double *crpix,
                                struct gal_fits_key_ll *headers,
                                char *program_string);





/**************************************************************/
/**********                  Table                 ************/
/**************************************************************/
void
gal_fits_tab_size(fitsfile *fitsptr, size_t *nrows, size_t *ncols);

int
gal_fits_tab_type(fitsfile *fptr);

gal_data_t *
gal_fits_tab_info(char *filename, char *hdu, size_t *numcols,
                    size_t *numrows, int *tabletype);

gal_data_t *
gal_fits_tab_read(char *filename, char *hdu, size_t numrows,
                    gal_data_t *colinfo, struct gal_linkedlist_sll *indexll,
                    int minmapsize);

void
gal_fits_tab_write(gal_data_t *cols, struct gal_linkedlist_stll *comments,
                   int tabletype, char *filename, int dontdelete);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_FITS_H__ */
