/*********************************************************************
Interpolate - Fill blank values in a dataset
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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/data.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>
#include <gnuastro/interpolate.h>
#include <gnuastro/permutation.h>






/*********************************************************************/
/********************      Nearest neighbor       ********************/
/*********************************************************************/
/* We want the flags to be powers of two so we can use bit-wise checks. */
#define INTERPOLATE_FLAGS_NO       0
#define INTERPOLATE_FLAGS_CHECKED  1<<0
#define INTERPOLATE_FLAGS_BLANK    1<<1





/* Paramters for interpolation on threads. */
struct interpolate_params
{
  gal_data_t                    *input;
  gal_data_t                      *out;
  gal_data_t                   *blanks;
  size_t                  numneighbors;
  uint8_t                *thread_flags;
  void                       *ngb_vals;
  int                        onlyblank;
  struct gal_tile_two_layer_params *tl;
};





/* Run the interpolation on many threads. */
static void *
interpolate_close_neighbors_on_thread(void *in_prm)
{
  /* Variables that others depend on. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct interpolate_params *prm=(struct interpolate_params *)(tprm->params);
  struct gal_tile_two_layer_params *tl=prm->tl;
  int correct_index=(tl && tl->totchannels>1 && tl->workoverch==0);
  gal_data_t *input=prm->input;

  /* Rest of variables. */
  float pdist;
  gal_data_t *nearest, *single;
  size_t ngb_counter, dist, *dinc;
  struct gal_linkedlist_tosll *lQ, *sQ;
  uint8_t *b, *bf, *bb, type=input->type;
  void *values=input->array, *out=prm->out->array;
  size_t pind, twidth=gal_data_sizeof(input->type);
  size_t i, index, fullind, chstart, ndim=input->ndim;
  size_t size = (correct_index ? tl->tottilesinch : input->size);
  size_t *icoord=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);
  size_t *ncoord=gal_data_malloc_array(GAL_DATA_TYPE_SIZE_T, ndim);
  size_t *dsize = (correct_index ? tl->numtilesinch : input->dsize);
  uint8_t *fullflag=&prm->thread_flags[tprm->id*input->size], *flag=fullflag;
  void *nv=gal_data_ptr_increment(prm->ngb_vals, tprm->id*prm->numneighbors,
                                  input->type);

  /* Initializations. */
  bb=prm->blanks->array;
  bf=(b=fullflag)+input->size;
  dinc=gal_dimension_increment(ndim, dsize);
  do *b = *bb++ ? INTERPOLATE_FLAGS_BLANK : 0; while(++b<bf);


  /* Put the nearest neighbor values into a structure for easy processing
     later. */
  nearest=gal_data_alloc(nv, input->type, 1, &prm->numneighbors, NULL, 0,
                         0, NULL, NULL, NULL);


  /* Go over all the points given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {
      /* For easy reading. */
      fullind=tprm->indexs[i];


      /* If the caller only wanted to interpolate over blank values and
         this value is not blank (we know from the flags), then just set
         the output value at this element to the input value and go to the
         next element. */
      if(prm->onlyblank && !(fullflag[fullind] & INTERPOLATE_FLAGS_BLANK) )
        {
          /* It should be `input->array', not `values'! Because when
             channels are treated separately, `values' is going to
             change. */
          memcpy(gal_data_ptr_increment(out,          fullind, type),
                 gal_data_ptr_increment(input->array, fullind, type),
                 twidth);
          continue;
        }


      /* Correct the index (if necessary). When the values come from a
         tiled dataset, the caller might want to interpolate the values of
         each channel separately (not mix values from different
         channels). In such a case, the tiles of each channel (and their
         values in `input' are contiguous. So we need to correct
         `tprm->indexs[i]' (which is the index over the whole tessellation,
         including all channels). */
      if(correct_index)
        {
          /* Index of this tile in its channel. */
          index = fullind % tl->tottilesinch;

          /* Index of the first tile in this channel. */
          chstart = (fullind / tl->tottilesinch) * tl->tottilesinch;

          /* Correct the values and flag pointers so we can only work in
             this channel. */
          values = gal_data_ptr_increment(input->array, chstart, type);
          flag = gal_data_ptr_increment(fullflag, chstart,
                                        GAL_DATA_TYPE_UINT8);
        }
      else
        index=fullind;


      /* Reset all checked bits in the flags array to 0. */
      ngb_counter=0;
      bf=(b=flag)+size;
      do *b &= ~(INTERPOLATE_FLAGS_CHECKED); while(++b<bf);


      /* Get the coordinates of this pixel (to be interpolated). */
      gal_dimension_index_to_coord(index, ndim, dsize, icoord);


      /* Start parsing the neighbors. We will use a two-way ordered linked
         list structure. To start from the nearest and go out to the
         farthest. */
      lQ=sQ=NULL;
      gal_linkedlist_add_to_tosll_end(&lQ, &sQ, index, 0.0f);
      while(sQ)
        {
          /* Pop-out (p) an index from the queue: */
          gal_linkedlist_pop_from_tosll_start(&lQ, &sQ, &pind, &pdist);

          /* If this isn't a blank value (recall that `index' was the first
             element added to the list which might be blank). */
          if( !(flag[pind] & INTERPOLATE_FLAGS_BLANK) )
            {
              /* Copy the value into the `nv' array. */
              memcpy(gal_data_ptr_increment(nv, ngb_counter, type),
                     gal_data_ptr_increment(values, pind, type),
                     twidth);

              /* If we have filled all the elements, break out. */
              if(++ngb_counter>=prm->numneighbors) break;
            }

          /* Go over all the neighbors of this popped pixel and add them to
             the list of neighbors to be checked. */
          GAL_DIMENSION_NEIGHBOR_OP(pind, ndim, dsize, 1, dinc,
           {
             /* Only look at neighbors that have not been checked. VERY
                IMPORTANT: we must not check for blank values here,
                otherwise we won't be able to parse over extended blank
                regions. */
             if( !(flag[nind] & INTERPOLATE_FLAGS_CHECKED) )
               {
                 /* Get the coordinates of this neighbor. */
                 gal_dimension_index_to_coord(nind, ndim, dsize, ncoord);

                 /* Distance of this neighbor to the one to be filled. */
                 dist=gal_dimension_dist_manhattan(icoord, ncoord, ndim);

                 /* Add this neighbor to the list. */
                 gal_linkedlist_add_to_tosll_end(&lQ, &sQ, nind, dist);

                 /* Flag this neighbor as checked. */
                 flag[nind] |= INTERPOLATE_FLAGS_CHECKED;
               }
           } );

          /* If there are no more meshes to add to the queue, then this
             shows, there were not enough points for
             interpolation. Normally, this loop should only be exited
             through the `currentnum>=numnearest' check above. */
          if(sQ==NULL)
            error(EXIT_FAILURE, 0, "only %zu neighbors found while you had "
                  "asked to use %zu neighbors for close neighbor "
                  "interpolation", ngb_counter, prm->numneighbors);
        }

      /* Calculate the median of the values and write it in. */
      single=gal_statistics_median(nearest, 1);
      memcpy(gal_data_ptr_increment(out, fullind, type), single->array,
             twidth);

      /* Clean up. */
      gal_data_free(single);
    }


  /* Clean up, wait for all the other threads to finish and return. */
  free(icoord);
  free(ncoord);
  nearest->array=NULL;
  gal_data_free(nearest);
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Interpolate blank values in an array. If the `tl!=NULL', then it is
   assumed that the tile values correspond to given tessellation. Such that
   `input[i]' corresponds to `tiles[i]' in the tessellation. */
gal_data_t *
gal_interpolate_close_neighbors(gal_data_t *input,
                                struct gal_tile_two_layer_params *tl,
                                size_t numneighbors, size_t numthreads,
                                int onlyblank)
{
  struct interpolate_params prm;
  int permute=(tl && tl->totchannels>1 && tl->workoverch);


  /* Initialize the constant parameters. */
  prm.tl=tl;
  prm.input=input;
  prm.onlyblank=onlyblank;
  prm.numneighbors=numneighbors;


  /* Flag the blank values. */
  prm.blanks=gal_blank_flag(input);


  /* Allocate space for the output. */
  prm.out=gal_data_alloc(NULL, input->type, input->ndim, input->dsize,
                         input->wcs, 0, input->minmapsize, "INTERPOLATED",
                         input->unit, NULL);


  /* If the input is from a tile structure and the user has asked to ignore
     channels, then re-order the values. */
  if(permute)
    {
      /* Prepare the permutation (if necessary/not already defined). */
      gal_tile_full_permutation(tl);

      /* Re-order values to ignore channels (if necessary). */
      gal_permutation_apply(input, tl->permutation);
      gal_permutation_apply(prm.blanks, tl->permutation);
    }


  /* Allocate space for all the flag values of all the threads here (memory
     in each thread is limited) and this is cleaner. */
  prm.thread_flags=gal_data_malloc_array(GAL_DATA_TYPE_UINT8,
                                         numthreads*input->size);


  /* Allocate space for the list of neighbor values in each thread. */
  prm.ngb_vals=gal_data_malloc_array(input->type, numthreads*numneighbors);


  /* Spin off the threads. */
  gal_threads_spin_off(interpolate_close_neighbors_on_thread, &prm,
                       input->size, numthreads);


  /* If the values were permuted for the interpolation, then re-order the
     values back to their original location (so they correspond to their
     tile indexs. */
  if(permute)
    {
      gal_permutation_apply_inverse(input, tl->permutation);
      gal_permutation_apply_inverse(prm.out, tl->permutation);
    }


  /* Clean up and return. */
  gal_data_free(prm.blanks);
  free(prm.thread_flags);
  free(prm.ngb_vals);
  return prm.out;
}
