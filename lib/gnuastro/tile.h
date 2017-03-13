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
gal_tile_block_check_tiles(gal_data_t *tiles);




/***********************************************************************/
/**************           Tile full dataset         ********************/
/***********************************************************************/
size_t *
gal_tile_all_sanity_check(char *filename, char *hdu, gal_data_t *input,
                          size_t *tile, size_t *numchannels);

size_t *
gal_tile_all_position(gal_data_t *input, size_t *regular,
                      float remainderfrac, gal_data_t **out, size_t multiple);

size_t *
gal_tile_all_position_two_layers(gal_data_t *input, size_t *channel_size,
                                 size_t *tile_size, float remainderfrac,
                                 gal_data_t **channels, gal_data_t **tiles);




/*********************************************************************/
/********************         On threads          ********************/
/*********************************************************************/
struct gal_tile_thread_param
{
  size_t            id; /* Id of this thread.                            */
  void         *params; /* Input structure for higher-level settings.    */
  size_t       *indexs; /* Indexes of actions to be done in this thread. */
  pthread_barrier_t *b; /* Pointer the barrier for all threads.          */
};

void
gal_tile_function_on_threads(gal_data_t *tiles, void *(*meshfunc)(void *),
                             size_t numthreads, void *caller_params);





/***********************************************************************/
/**************           Function-like macros        ******************/
/***********************************************************************/
/* Parse over a region of memory (can be an n-dimensional tile or a fully
   allocated block of memory) and do a certain operation. If `OUT' is not
   NULL, it can also save the output into the `gal_data_t' that it points
   to. Note that it must either have only one element (for the whole input)
   or have exactly the same number of elements as the input (one value for
   one pixel/element of the input). The input arguments are:

       `IT': input type (C type) for example `int32_t' or `float'.

       `OT': output type (C type) for example `int32_t' or `float'.

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
             be NULL and should not be used.

       `PARSE_OUT': Parse the output along with the input. When this is
             non-zero, then the `o' pointer (described in `OP') will be
             incremented to cover the same region of the input in the
             output.

   See `lib/statistics.c' for some example applications of this function.
*/
#define GAL_TILE_PARSE_OPERATE(IT, OT, OP, IN, OUT, PARSE_OUT) {        \
    size_t increment=0, num_increment=1;                                \
    OT *ost, *o = OUT ? OUT->array : NULL;                              \
    gal_data_t *iblock=gal_tile_block(IN);                              \
    IT b, *st, *i=IN->array, *f=i+IN->size;                             \
    gal_data_t *oblock = OUT ? gal_tile_block(OUT) : NULL;              \
    int hasblank=gal_blank_present(IN), parse_out=(OUT && PARSE_OUT);   \
    size_t s_e_ind[2]={0,iblock->size-1}; /* -1: this is INCLUSIVE */   \
                                                                        \
    /* Write the blank value for the input type into `b'). Note that */ \
    /* a tile doesn't necessarily have to have a type. */               \
    gal_blank_write(&b, iblock->type);                                  \
                                                                        \
    /* A small sanity check. */                                         \
    if( parse_out && gal_data_dsize_is_different(iblock, oblock) )      \
      error(EXIT_FAILURE, 0, "when `PARSE_OUT' is non-zero, the "       \
            "allocated block size of the input and output of "          \
            "`GAL_TILE_PARSE_OPERATE' must be equal, but they are "     \
            "not: %zu and %zu elements respectively)", iblock->size,    \
            oblock->size);                                              \
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
            if(b==b) do if(*i!=b)  {OP; if(parse_out) ++o;} while(++i<f); \
            else     do if(*i==*i) {OP; if(parse_out) ++o;} while(++i<f); \
          }                                                             \
        else         do            {OP; if(parse_out) ++o;} while(++i<f); \
                                                                        \
        /* Set the incrementation. On a fully allocated iblock (when */ \
        /* `IN==iblock'), we have already gone through the whole     */ \
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





__END_C_DECLS    /* From C++ preparations */

#endif
