/*********************************************************************
tile -- work with tesselations over a host dataset.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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
#ifndef __GAL_TILE_H__
#define __GAL_TILE_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
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








/***********************************************************************/
/**************             About block           **********************/
/***********************************************************************/
gal_data_t *
gal_tile_block(gal_data_t *input);

void
gal_tile_block_start_coord(gal_data_t *tile, size_t *start_coord);

void *
gal_tile_block_start_end(gal_data_t *tile, gal_data_t *work,
                         size_t *start_end);

size_t
gal_tile_block_increment(gal_data_t *block, size_t *tsize,
                         size_t num_increment);

void
gal_tile_block_check_tiles(gal_data_t *tiles, char *filename,
                           char *program_name);




/***********************************************************************/
/**************           Tile full dataset         ********************/
/***********************************************************************/
size_t *
gal_tile_all_sanity_check(char *filename, char *hdu, gal_data_t *input,
                          size_t *tile, size_t *numchannels);

size_t
gal_tile_all_position(gal_data_t *input, size_t *regular,
                      float remainderfrac, gal_data_t **out, size_t multiple);

void
gal_tile_all_position_two_layers(gal_data_t *input, size_t *channel_size,
                                 size_t *tile_size, float remainderfrac,
                                 gal_data_t **channels, gal_data_t **tiles);




__END_C_DECLS    /* From C++ preparations */

#endif
