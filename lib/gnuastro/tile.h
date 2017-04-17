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
/**************           Allocated block         **********************/
/***********************************************************************/
gal_data_t *
gal_tile_block(gal_data_t *input);

size_t
gal_tile_block_increment(gal_data_t *block, size_t *tsize,
                         size_t num_increment, size_t *coord);

gal_data_t *
gal_tile_block_write_const_value(gal_data_t *tilevalues, gal_data_t *tilesll,
                                 int initialize);

gal_data_t *
gal_tile_block_check_tiles(gal_data_t *tiles);

void *
gal_tile_block_relative_to_other(gal_data_t *tile, gal_data_t *other);





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
  size_t           *firsttsize; /* See `gal_tile_full_regular_first'.     */

  /* Actual tile and channel data structures. */
  gal_data_t            *tiles; /* Tiles array (also linked with `next'). */
  gal_data_t         *channels; /* Channels array (linked with `next').   */
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
                           char *filename, char *program_string);

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
#define GAL_TILE_PO_OISET(IT, OT, OP, IN, OUT, PARSE_OUT, CHECK_BLANK) { \
    gal_data_t *out_w=OUT; /* Since `OUT' may be NULL. */               \
    OT *ost, *o = out_w ? out_w->array : NULL;                          \
    IT b=0, *st=NULL, *i=IN->array, *f=i+IN->size;                      \
                                                                        \
    /* Write the blank value for the input type into `b'). Note that */ \
    /* a tile doesn't necessarily have to have a type. */               \
    if(hasblank) gal_blank_write(&b, iblock->type);                     \
                                                                        \
    /* If this is a tile, not a full block. */                          \
    if(IN!=iblock)                                                      \
      {                                                                 \
        st=gal_tile_start_end_ind_inclusive(IN, iblock, s_e_ind);       \
        if(parse_out)                                                   \
          ost = (OT *)(oblock->array) + ( st - (IT *)(iblock->array) ); \
      }                                                                 \
                                                                        \
    /* Go over contiguous patches of memory. */                         \
    while( s_e_ind[0] + increment <= s_e_ind[1] )                       \
      {                                                                 \
                                                                        \
        /* If we are on a tile, reset `i' and `f'. Also, if there   */  \
        /* is more than one element in `OUT', then set that. Note   */  \
        /* that it is upto the caller to increment `o' in `OP'.     */  \
        if(IN!=iblock)                                                  \
          {                                                             \
            f = ( i = st + increment ) + IN->dsize[IN->ndim-1];         \
            if(parse_out) o = ost + increment;                          \
          }                                                             \
                                                                        \
        /* Do the operation depending the nature of the blank value. */ \
        /* Recall that for integer types, the blank value must be    */ \
        /* checked with `=='. But for floats, the blank value can be */ \
        /* a NaN. Recall that a NAN will fail any comparison         */ \
        /* including `=='. So when the blank value is not equal to   */ \
        /* itself, then it is floating point and is a NAN. In that   */ \
        /* case, the only way to check if a data element is blank or */ \
        /* not is to see if the value of each element is equal to    */ \
        /* itself or not. */                                            \
        if(hasblank)                                                    \
          {                                                             \
            if(b==b) do{if(*i!=b) {OP;}if(parse_out) ++o;} while(++i<f); \
            else     do{if(*i==*i){OP;}if(parse_out) ++o;} while(++i<f); \
          }                                                             \
        else         do          {{OP;}if(parse_out) ++o;} while(++i<f); \
                                                                        \
                                                                        \
        /* Set the incrementation. On a fully allocated iblock (when */ \
        /* `IN==iblock'), we have already gone through the whole   */   \
        /* array, so we'll set the incrementation to the size of the */ \
        /* while block which will stop the `while' loop above. On a  */ \
        /* tile, we need to increment to the next contiguous patch   */ \
        /* of memory to continue parsing this tile. */                  \
        increment += ( IN==iblock ? iblock->size                        \
                       : gal_tile_block_increment(iblock, IN->dsize,    \
                                                  num_increment++,      \
                                                  NULL) );              \
      }                                                                 \
                                                                        \
    /* This is done in case the caller doesn't need `o' to avoid */     \
    /* compiler warnings. */                                            \
    o = o ? o+0 : NULL;                                                 \
  }


#define GAL_TILE_PO_OSET(OT, OP, IN, OUT, PARSE_OUT, CHECK_BLANK) {     \
  switch(iblock->type)                                                  \
    {                                                                   \
    case GAL_TYPE_UINT8:                                                \
      GAL_TILE_PO_OISET(uint8_t,  OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_INT8:                                                 \
      GAL_TILE_PO_OISET(int8_t,   OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_UINT16:                                               \
      GAL_TILE_PO_OISET(uint16_t, OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_INT16:                                                \
      GAL_TILE_PO_OISET(int16_t,  OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_UINT32:                                               \
      GAL_TILE_PO_OISET(uint32_t, OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_INT32:                                                \
      GAL_TILE_PO_OISET(int32_t,  OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_UINT64:                                               \
      GAL_TILE_PO_OISET(uint64_t, OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_INT64:                                                \
      GAL_TILE_PO_OISET(int64_t,  OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_FLOAT32:                                              \
      GAL_TILE_PO_OISET(float,    OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    case GAL_TYPE_FLOAT64:                                              \
      GAL_TILE_PO_OISET(double,   OT,OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
      break;                                                            \
    default:                                                            \
      {                                                                 \
        fprintf(stderr, "type code %d not recognized in "               \
              "`GAL_TILE_PO_OSET'", iblock->type);                      \
        exit(EXIT_FAILURE);                                             \
      }                                                                 \
    }                                                                   \
  }


/* Parse over a region of memory (can be an n-dimensional tile or a fully
   allocated block of memory) and do a certain operation. If `OUT' is not
   NULL, it can also save the output of the operation into the `gal_data_t'
   that it points to. Note that OUT must either have only one element (for
   the whole input) or have exactly the same number of elements as the
   input (one value for one pixel/element of the input). The input
   arguments are:

       `OP': Operator: this can be any number of C experssions. This macro
             is going to define an `IT *i' variable which will increment
             over each element of the input array/tile. It will also define
             an `OT *o' which you can use to access the output. If
             `PARSE_OUT' is non-zero, then `o' will also be incremented to
             the same index element but in the output array. You can use
             these along with any other variable you define before this
             macro to process the input and store it in the output. Note
             that `i' and `o' will be incremented once in here, so don't
             increment them in `OP'.

       `IN': Input `gal_data_t', can be a tile or an allocated block of
             memory.

       `OUT': Output `gal_data_t'. It can be NULL. In that case, `o' will
             be NULL and should not be used. If given,

       `PARSE_OUT': Parse the output along with the input. When this is
             non-zero, then the `o' pointer (described in `OP') will be
             incremented to cover the same region of the input in the
             output.

       `CHECK_BLANK': If it is non-zero, then the input will be checked for
             blank values.

   For `OP' you have access to some other variables besides `i' and `o':

      `i': Pointer to this element of input to parse.

      `o': Pointer to corresponding element in output.

      `b': blank value in input't type.

   See `lib/statistics.c' for some example applications of this function.

  -----------
  TIPS/TRICKS
  -----------

  You can use a given tile on a dataset that it was not initialized with
  (but has the same size). You can do that like this:

     void *tarray;
     gal_data_t *tblock;

     tarray=tile->array;
     tblock=tile->block;
     tile->array=gal_tile_block_relative_to_other(tile, THE_OTHER_DATASET);
     tile->block=THE_OTHER_DATASET;  <<-- `tile->block' must corrected
                                          AFTER `tile->array'.
     GAL_TILE_PARSE_OPERATE(...)

     tile->array=tarray;
     tile->block=tblock;

  You can work on one other array while working on the one that tile
  actually points to. To do that, you can make a fake tile and pass that as
  the `OUT' argument to this macro. While the name is `OUT' it can be used
  in any manner you like. To do this, you can take the following steps. In
  this example, `tile' is the actual tile that you have.

     gal_data_t *faketile;

     // These can be done outside a loop.
     faketile=gal_data_alloc(NULL, OTHER_DATASET->type, 1, &dsize,
                            NULL, 0, -1, NULL, NULL, NULL);
     free(faketile->array);
     free(faketile->dsize);
     faketile->block=OTHER_DATASET;
     faketile->ndim=OTHER_DATASET->ndim;

     // These can be done in the loop over tiles.
     faketile->size=tile->size;
     faketile->dsize=tile->dsize;
     faketile->array=gal_tile_block_relative_to_other(tile, OTHER_DATASET);

     // Do your processing....
     GAL_TILE_PARSE_OPERATE(PROCESSING, tile, faketile, 1, 1);

     // Clean up.
     bintile->array=NULL;
     bintile->dsize=NULL;
     gal_data_free(bintile);
*/
#define GAL_TILE_PARSE_OPERATE(OP, IN, OUT, PARSE_OUT, CHECK_BLANK) {   \
    int parse_out=(OUT && PARSE_OUT);                                   \
    size_t increment=0, num_increment=1;                                \
    gal_data_t *iblock = gal_tile_block(IN);                            \
    gal_data_t *oblock = OUT ? gal_tile_block(OUT) : NULL;              \
    int hasblank = CHECK_BLANK ? gal_blank_present(IN) : 0;             \
    size_t s_e_ind[2]={0,iblock->size-1}; /* -1: this is INCLUSIVE */   \
                                                                        \
    /* A small sanity check. */                                         \
    if( parse_out && gal_data_dsize_is_different(iblock, oblock) )      \
      {                                                                 \
        /* The `error' function, is a GNU extension. */                 \
        fprintf(stderr, "when `PARSE_OUT' is non-zero, the "            \
                "allocated block size of the input and output of "      \
                "`GAL_TILE_PARSE_OPERATE' must be equal, but they are " \
                "not: %zu and %zu elements respectively)",              \
                iblock->size, oblock->size);                            \
        exit(EXIT_FAILURE);                                             \
      }                                                                 \
                                                                        \
    /* First set the OUTPUT type. */                                    \
    if(OUT)                                                             \
      switch(oblock->type)                                              \
        {                                                               \
        case GAL_TYPE_UINT8:                                            \
          GAL_TILE_PO_OSET(uint8_t,  OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_INT8:                                             \
          GAL_TILE_PO_OSET(int8_t,   OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_UINT16:                                           \
          GAL_TILE_PO_OSET(uint16_t, OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_INT16:                                            \
          GAL_TILE_PO_OSET(int16_t,  OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_UINT32:                                           \
          GAL_TILE_PO_OSET(uint32_t, OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_INT32:                                            \
          GAL_TILE_PO_OSET(int32_t,  OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_UINT64:                                           \
          GAL_TILE_PO_OSET(uint64_t, OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_INT64:                                            \
          GAL_TILE_PO_OSET(int64_t,  OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_FLOAT32:                                          \
          GAL_TILE_PO_OSET(float,    OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        case GAL_TYPE_FLOAT64:                                          \
          GAL_TILE_PO_OSET(double,   OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
          break;                                                        \
        default:                                                        \
          {                                                             \
            fprintf(stderr, "type code %d not recognized in "           \
                    "`GAL_TILE_PARSE_OPERATE'", oblock->type);          \
            exit(EXIT_FAILURE);                                         \
          }                                                             \
        }                                                               \
    else                                                                \
      /* When `OUT==NULL', its type is irrelevant, we'll just use */    \
      /*`int' as a place holder. */                                     \
      GAL_TILE_PO_OSET(int,          OP,IN,OUT,PARSE_OUT,CHECK_BLANK);  \
  }



__END_C_DECLS    /* From C++ preparations */

#endif
