/*********************************************************************
Functions to that only use WCSLIB's functions, not related to FITS.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#ifndef __GAL_WCS_H__
#define __GAL_WCS_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <fitsio.h>
#include <wcslib/wcs.h>

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
 ***********               Read WCS                ***********
 *************************************************************/
struct wcsprm *
gal_wcs_read_fitsptr(fitsfile *fptr, size_t hstartwcs, size_t hendwcs,
                     int *nwcs);

struct wcsprm *
gal_wcs_read(char *filename, char *hdu, size_t hstartwcs,
             size_t hendwcs, int *nwcs);





/**************************************************************/
/**********              Utilities                 ************/
/**************************************************************/
struct wcsprm *
gal_wcs_copy(struct wcsprm *wcs);

void
gal_wcs_remove_dimension(struct wcsprm *wcs, size_t fitsdim);

void
gal_wcs_on_tile(gal_data_t *tile);

double *
gal_wcs_warp_matrix(struct wcsprm *wcs);

void
gal_wcs_decompose_pc_cdelt(struct wcsprm *wcs);

double
gal_wcs_angular_distance_deg(double r1, double d1, double r2, double d2);

double *
gal_wcs_pixel_scale(struct wcsprm *wcs);

double
gal_wcs_pixel_area_arcsec2(struct wcsprm *wcs);





/**************************************************************/
/**********              Conversion                ************/
/**************************************************************/
gal_data_t *
gal_wcs_world_to_img(gal_data_t *coords, struct wcsprm *wcs, int inplace);

gal_data_t *
gal_wcs_img_to_world(gal_data_t *coords, struct wcsprm *wcs, int inplace);





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_WCS_H__ */
