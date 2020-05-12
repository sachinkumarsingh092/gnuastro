/*********************************************************************
tile -- work with tesselations over a host dataset.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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
#include <gnuastro/fits.h>
#include <gnuastro/dimension.h>

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
/**************              Single tile              ******************/
/***********************************************************************/
void
gal_tile_start_coord(gal_data_t *tile, size_t *start_coord);

void
gal_tile_start_end_coord(gal_data_t *tile, size_t *start_end, int rel_block);

void *
gal_tile_start_end_ind_inclusive(gal_data_t *tile, gal_data_t *work,
                                 size_t *start_end_inc);




/***********************************************************************/
/**************           Series of tiles             ******************/
/***********************************************************************/
gal_data_t *
gal_tile_series_from_minmax(gal_data_t *block, size_t *minmax, size_t number);





/***********************************************************************/
/**************           Allocated block         **********************/
/***********************************************************************/
gal_data_t *
gal_tile_block(gal_data_t *tile);

size_t
gal_tile_block_increment(gal_data_t *block, size_t *tsize,
                         size_t num_increment, size_t *coord);

gal_data_t *
gal_tile_block_write_const_value(gal_data_t *tilevalues, gal_data_t *tilesll,
                                 int withblank, int initialize);

gal_data_t *
gal_tile_block_check_tiles(gal_data_t *tiles);

void *
gal_tile_block_relative_to_other(gal_data_t *tile, gal_data_t *other);

void
gal_tile_block_blank_flag(gal_data_t *tile_ll, size_t numthreads);





/***********************************************************************/
/**************           Tile full dataset         ********************/
/***********************************************************************/

struct gal_tile_two_layer_params
{
  /* Inputs */
  size_t             *tilesize; /* Tile size along each dim. (C order).   */
  size_t          *numchannels; /* Channel no. along each dim. (C order). */
  float          remainderfrac; /* Frac. of remainers in each dim to cut. */
  uint8_t           workoverch; /* Convolve over channel borders.         */
  uint8_t           checktiles; /* Tile IDs in an img, the size of input. */
  uint8_t       oneelempertile; /* Only use one element for each tile.    */

  /* Internal parameters. */
  size_t                  ndim; /* The number of dimensions.              */
  size_t              tottiles; /* Total number of tiles in all dim.      */
  size_t          tottilesinch; /* Number of tiles in one channel.        */
  size_t           totchannels; /* Total number of channels in all dim.   */
  size_t          *channelsize; /* Size of channels along each dimension. */
  size_t             *numtiles; /* Tile no. in each dim. over-all.        */
  size_t         *numtilesinch; /* Tile no. in each dim. on one channel.  */
  char          *tilecheckname; /* Name of file to check tiles.           */
  size_t          *permutation; /* Tile pos. in memory --> pos. overall.  */
  size_t           *firsttsize; /* See 'gal_tile_full_regular_first'.     */

  /* Actual tile and channel data structures. */
  gal_data_t            *tiles; /* Tiles array (also linked with 'next'). */
  gal_data_t         *channels; /* Channels array (linked with 'next').   */
};


size_t *
gal_tile_full(gal_data_t *input, size_t *regular,
              float remainderfrac, gal_data_t **out, size_t multiple,
              size_t **firsttsize);

void
gal_tile_full_sanity_check(char *filename, char *hdu, gal_data_t *input,
                           struct gal_tile_two_layer_params *tl);

void
gal_tile_full_two_layers(gal_data_t *input,
                         struct gal_tile_two_layer_params *tl);

void
gal_tile_full_permutation(struct gal_tile_two_layer_params *tl);

void
gal_tile_full_values_write(gal_data_t *tilevalues,
                           struct gal_tile_two_layer_params *tl,
                           int withblank, char *filename,
                           gal_fits_list_key_t *keys, char *program_string);

gal_data_t *
gal_tile_full_values_smooth(gal_data_t *tilevalues,
                            struct gal_tile_two_layer_params *tl,
                            size_t width, size_t numthreads);

size_t
gal_tile_full_id_from_coord(struct gal_tile_two_layer_params *tl,
                            size_t *coord);

void
gal_tile_full_blank_flag(gal_data_t *tile_ll, size_t numthreads);

void
gal_tile_full_free_contents(struct gal_tile_two_layer_params *tl);





/***********************************************************************/
/**************           Function-like macros        ******************/
/***********************************************************************/
/* Useful when the input and other types are already known. We want this to
   be self-sufficient (and be possible to call it independent of
   'GAL_TILE_PARSE_OPERATE'), so some variables (basic definitions) that
   are already defined in 'GAL_TILE_PARSE_OPERATE' re-defined here. */
#define GAL_TILE_PO_OISET(IT, OT, IN, OTHER, PARSE_OTHER, CHECK_BLANK, OP) { \
    IT *i=IN->array;                                                    \
    gal_data_t *tpo_other=OTHER; /* 'OTHER' may be NULL. */             \
    gal_data_t *tpo_oblock = OTHER ? gal_tile_block(OTHER) : NULL;      \
                                                                        \
    size_t tpo_s_e_i_junk[2]={0,0};                                     \
    IT b, *tpo_st=NULL, *tpo_f=i+IN->size;                              \
    size_t tpo_i_increment=0, tpo_num_i_inc=1;                          \
    size_t tpo_o_increment=0, tpo_num_o_inc=1;                          \
    int tpo_parse_other=(OTHER && PARSE_OTHER);                         \
    gal_data_t *tpo_iblock = gal_tile_block(IN);                        \
    OT *tpo_ost=NULL, *o = tpo_other ? tpo_other->array : NULL;         \
    int tpo_hasblank = CHECK_BLANK ? gal_blank_present(IN, 0) : 0;      \
    size_t tpo_s_e_i[2]={0,tpo_iblock->size-1}; /* -1: this is INCLUSIVE */ \
                                                                        \
                                                                        \
    /* A small sanity check: if 'OTHER' is given, and it is a block, */ \
    /* then it must have the same size as 'IN's block. On the other  */ \
    /* hand, when 'OTHER' is a tile, its must have 'IN's size.       */ \
    if( tpo_parse_other )                                               \
      {                                                                 \
        if( OTHER==tpo_oblock )    /* 'OTHER' is a block. */            \
          {                                                             \
            if( gal_dimension_is_different(tpo_iblock, tpo_oblock) )    \
              {                                                         \
                /* 'error' function, is a GNU extension, see above. */  \
                fprintf(stderr, "GAL_TILE_PO_OISET: when "              \
                        "'PARSE_OTHER' is non-zero, the allocated "     \
                        "block size of 'IN' and 'OTHER' must be "       \
                        "equal, but they are not: %zu and %zu "         \
                        "elements respectively\n", tpo_iblock->size,    \
                        tpo_oblock->size);                              \
                exit(EXIT_FAILURE);                                     \
              }                                                         \
          }                                                             \
        else                                                            \
          if( gal_dimension_is_different(IN, OTHER) )                   \
            {                                                           \
              /* The 'error' function, is a GNU extension and this */   \
              /* is a header, not a library which the user has to  */   \
              /* compile every time (on their own system).         */   \
              fprintf(stderr, "GAL_TILE_PO_OISET: when "                \
                      "'PARSE_OTHER' is non-zero, the sizes of 'IN' "   \
                      "and 'OTHER' must be equal (in all "              \
                      "dimensions), but they are not: %zu and %zu "     \
                      "elements respectively\n", IN->size,              \
                      tpo_other->size);                                 \
              exit(EXIT_FAILURE);                                       \
            }                                                           \
      }                                                                 \
                                                                        \
                                                                        \
    /* Write the blank value for the input type into 'b'. */            \
    gal_blank_write(&b, tpo_iblock->type);                              \
                                                                        \
                                                                        \
    /* If this is a tile, not a full block, then we need to set the  */ \
    /* starting pointers ('tpo_st' and 'tpo_ost'). The latter needs  */ \
    /* special attention: if it is a block, then we will use the     */ \
    /* the same starting element as the input tile. If 'OTHER' is a  */ \
    /* tile, then use its own starting position (recall that we have */ \
    /* already made sure that 'IN' and 'OTHER' have the same size.   */ \
    if(IN!=tpo_iblock)                                                  \
      {                                                                 \
        tpo_st = gal_tile_start_end_ind_inclusive(IN, tpo_iblock,       \
                                                  tpo_s_e_i);           \
        if( tpo_parse_other )                                           \
          tpo_ost = ( OTHER==tpo_oblock                                 \
                      ? ( (OT *)(tpo_oblock->array)                     \
                          + ( tpo_st - (IT *)(tpo_iblock->array) ) )    \
                      : gal_tile_start_end_ind_inclusive(tpo_other,     \
                                                         tpo_oblock,    \
                                                         tpo_s_e_i_junk) ); \
      }                                                                 \
                                                                        \
                                                                        \
    /* Go over contiguous patches of memory. */                         \
    while( tpo_s_e_i[0] + tpo_i_increment <= tpo_s_e_i[1] )             \
      {                                                                 \
        /* If we are on a tile, reset 'i' and 'o'. */                   \
        if(IN!=tpo_iblock)                                              \
          {                                                             \
            tpo_f = ( ( i = tpo_st + tpo_i_increment )                  \
                      + IN->dsize[IN->ndim-1] );                        \
            if(tpo_parse_other) o = tpo_ost + tpo_o_increment;          \
          }                                                             \
                                                                        \
        /* Do the operation depending the nature of the blank value. */ \
        /* Recall that for integer types, the blank value must be    */ \
        /* checked with '=='. But for floats, the blank value can be */ \
        /* a NaN. Recall that a NAN will fail any comparison         */ \
        /* including '=='. So when the blank value is not equal to   */ \
        /* itself, then it is floating point and is a NAN. In that   */ \
        /* case, the only way to check if a data element is blank or */ \
        /* not is to see if the value of each element is equal to    */ \
        /* itself or not. */                                            \
        if(tpo_hasblank)                                                \
          {                                                             \
            if(b==b)                                                    \
              do{if(*i!=b)  {OP;} if(tpo_parse_other) ++o;} while(++i<tpo_f);\
            else                                                        \
              do{if(*i==*i) {OP;} if(tpo_parse_other) ++o;} while(++i<tpo_f);\
          }                                                             \
        else                                                            \
          do    {           {OP;} if(tpo_parse_other) ++o;} while(++i<tpo_f);\
                                                                        \
        /* Set the incrementation. On a fully allocated iblock (when */ \
        /* 'IN==tpo_iblock'), we have already gone through the whole */ \
        /* array, so we'll set the incrementation to the size of the */ \
        /* whole block. This will stop the 'while' loop above. On a  */ \
        /* tile, we need to increment to the next contiguous patch   */ \
        /* of memory to continue parsing this tile. */                  \
        tpo_i_increment += ( IN==tpo_iblock                             \
                             ? tpo_iblock->size                         \
                             : gal_tile_block_increment(tpo_iblock,     \
                                                        IN->dsize,      \
                                                        tpo_num_i_inc++, \
                                                        NULL) );        \
                                                                        \
        /* Similarly, increment the other array if necessary. Like   */ \
        /* the above, when 'OTHER' is a full block, we'll just use   */ \
        /* the same increment as 'IN'. Otherwise, when 'OTHER' is a  */ \
        /* tile, calculate its increment based on its own block.     */ \
        if(tpo_parse_other)                                             \
          {                                                             \
            if(OTHER==tpo_oblock) tpo_o_increment=tpo_i_increment;      \
            else                                                        \
              tpo_o_increment += gal_tile_block_increment(tpo_oblock,   \
                                                          tpo_other->dsize, \
                                                          tpo_num_o_inc++, \
                                                          NULL);        \
          }                                                             \
      }                                                                 \
                                                                        \
    /* This is done in case the caller doesn't need 'o' to avoid */     \
    /* compiler warnings. */                                            \
    o = o ? o+0 : NULL;                                                 \
  }


#define GAL_TILE_PO_OSET(OT, IN, OTHER, PARSE_OTHER, CHECK_BLANK, OP) { \
  switch(tpo_iblock->type)                                              \
    {                                                                   \
    case GAL_TYPE_UINT8:                                                \
      GAL_TILE_PO_OISET(uint8_t,  OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_INT8:                                                 \
      GAL_TILE_PO_OISET(int8_t,   OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_UINT16:                                               \
      GAL_TILE_PO_OISET(uint16_t, OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_INT16:                                                \
      GAL_TILE_PO_OISET(int16_t,  OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_UINT32:                                               \
      GAL_TILE_PO_OISET(uint32_t, OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_INT32:                                                \
      GAL_TILE_PO_OISET(int32_t,  OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_UINT64:                                               \
      GAL_TILE_PO_OISET(uint64_t, OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_INT64:                                                \
      GAL_TILE_PO_OISET(int64_t,  OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_FLOAT32:                                              \
      GAL_TILE_PO_OISET(float,    OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      GAL_TILE_PO_OISET(double,   OT,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
      break;                                                            \
    default:                                                            \
      { /* 'error' function might not be available for the user. */     \
        fprintf(stderr, "GAL_TILE_PO_OSET: type code %d not recognized",\
                tpo_iblock->type);                                      \
        exit(EXIT_FAILURE);                                             \
      }                                                                 \
    }                                                                   \
  }


/* Parse over a region of memory (can be an n-dimensional tile or a fully
   allocated block of memory) and do a certain operation. If 'OTHER' is not
   NULL, this macro will also parse it at the same time . Note that OTHER
   must either have only one element (for the whole input) or have exactly
   the same number of elements as the input (one value for one
   pixel/element of the input). See the documentation for more on this
   macro and some examples. */
#define GAL_TILE_PARSE_OPERATE(IN, OTHER, PARSE_OTHER, CHECK_BLANK, OP) { \
    gal_data_t *tpo_iblock = gal_tile_block(IN);                        \
    gal_data_t *tpo_oblock = OTHER ? gal_tile_block(OTHER) : NULL;      \
                                                                        \
    /* First set the OTHER type. */                                     \
    if(OTHER)                                                           \
      switch(tpo_oblock->type)                                          \
        {                                                               \
        case GAL_TYPE_UINT8:                                            \
          GAL_TILE_PO_OSET(uint8_t, IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_INT8:                                             \
          GAL_TILE_PO_OSET(int8_t,  IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_UINT16:                                           \
          GAL_TILE_PO_OSET(uint16_t,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_INT16:                                            \
          GAL_TILE_PO_OSET(int16_t, IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_UINT32:                                           \
          GAL_TILE_PO_OSET(uint32_t,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_INT32:                                            \
          GAL_TILE_PO_OSET(int32_t, IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_UINT64:                                           \
          GAL_TILE_PO_OSET(uint64_t,IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_INT64:                                            \
          GAL_TILE_PO_OSET(int64_t, IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_FLOAT32:                                          \
          GAL_TILE_PO_OSET(float,   IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        case GAL_TYPE_FLOAT64:                                          \
          GAL_TILE_PO_OSET(double,  IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
          break;                                                        \
        default:                                                        \
          {                                                             \
            fprintf(stderr, "type code %d not recognized in "           \
                    "'GAL_TILE_PARSE_OPERATE'", tpo_oblock->type);      \
            exit(EXIT_FAILURE);                                         \
          }                                                             \
        }                                                               \
    else                                                                \
      /* When 'OTHER==NULL', its type is irrelevant, we'll just use */  \
      /*'int' as a place holder. */                                     \
      GAL_TILE_PO_OSET(int,         IN,OTHER,PARSE_OTHER,CHECK_BLANK,OP); \
  }



__END_C_DECLS    /* From C++ preparations */

#endif
