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

#include <gnuastro-internal/checkset.h>




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
  size_t                           num;
  gal_data_t                      *out;
  gal_data_t                   *blanks;
  size_t                  numneighbors;
  uint8_t                *thread_flags;
  int                        onlyblank;
  struct gal_linkedlist_vll  *ngb_vals;
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
  int correct_index=(tl && tl->totchannels>1 && !tl->workoverch);
  gal_data_t *input=prm->input;

  /* Rest of variables. */
  void *nv;
  float pdist;
  uint8_t *b, *bf, *bb;
  struct gal_linkedlist_vll *tvll;
  struct gal_linkedlist_tosll *lQ, *sQ;
  size_t ngb_counter, dist, pind, *dinc;
  size_t i, index, fullind, chstart=0, ndim=input->ndim;
  gal_data_t *median, *tin, *tout, *tnear, *nearest=NULL;
  size_t *icoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);
  size_t *ncoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);
  size_t size = (correct_index ? tl->tottilesinch : input->size);
  size_t *dsize = (correct_index ? tl->numtilesinch : input->dsize);
  uint8_t *fullflag=&prm->thread_flags[tprm->id*input->size], *flag=fullflag;


  /* Initializations. */
  bb=prm->blanks->array;
  bf=(b=fullflag)+input->size;
  dinc=gal_dimension_increment(ndim, dsize);
  do *b = *bb++ ? INTERPOLATE_FLAGS_BLANK : 0; while(++b<bf);


  /* Put the allocated space to keep the neighbor values into a structure
     for easy processing. */
  tin=input;
  for(tvll=prm->ngb_vals; tvll!=NULL; tvll=tvll->next)
    {
      nv=gal_data_ptr_increment(tvll->v, tprm->id*prm->numneighbors,
                                input->type);
      gal_data_add_to_ll(&nearest, nv, tin->type, 1, &prm->numneighbors,
                         NULL, 0, -1, NULL, NULL, NULL);
      tin=tin->next;
    }
  gal_data_reverse_ll(&nearest);


  /* Go over all the points given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. */
      fullind=tprm->indexs[i];


      /* If the caller only wanted to interpolate over blank values and
         this value is not blank (we know from the flags), then just set
         the output value at this element to the input value and go to the
         next element. */
      if(prm->onlyblank && !(fullflag[fullind] & INTERPOLATE_FLAGS_BLANK) )
        {
          tin=input;
          for(tout=prm->out; tout!=NULL; tout=tout->next)
            {
              memcpy(gal_data_ptr_increment(tout->array, fullind, tin->type),
                     gal_data_ptr_increment(tin->array,  fullind, tin->type),
                     gal_type_sizeof(tin->type));
              tin=tin->next;
            }
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

          /* Set the channel's starting pointer for the flags. */
          flag = gal_data_ptr_increment(fullflag, chstart, GAL_TYPE_UINT8);
        }
      else
        {
          chstart=0;
          index=fullind;
        }


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

          /* If this isn't a blank value then add its values to the list of
             neighbor values. Note that we didn't check whether the values
             were blank or not when adding this pixel to the queue. */
          if( !(flag[pind] & INTERPOLATE_FLAGS_BLANK) )
            {
              tin=input;
              for(tnear=nearest; tnear!=NULL; tnear=tnear->next)
                {
                  memcpy(gal_data_ptr_increment(tnear->array, ngb_counter,
                                                tin->type),
                         gal_data_ptr_increment(tin->array, chstart+pind,
                                                tin->type),
                         gal_type_sizeof(tin->type));
                  tin=tin->next;
                }

              /* If we have filled all the elements clean up the linked
                 list and break out. */
              if(++ngb_counter>=prm->numneighbors)
                {
                  if(lQ) gal_linkedlist_tosll_free(lQ);
                  break;
                }
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

      /* Calculate the median of the values and write it in the output. */
      tout=prm->out;
      for(tnear=nearest; tnear!=NULL; tnear=tnear->next)
        {
          /* Find the median and copy it. */
          median=gal_statistics_median(tnear, 1);
          memcpy(gal_data_ptr_increment(tout->array, fullind, tout->type),
                 median->array, gal_type_sizeof(tout->type));

          /* Clean up and go to next array. */
          gal_data_free(median);
          tout=tout->next;
        }
    }


  /* Clean up. */
  for(tnear=nearest; tnear!=NULL; tnear=tnear->next) tnear->array=NULL;
  gal_data_free_ll(nearest);
  free(icoord);
  free(ncoord);
  free(dinc);


  /* Wait for all the other threads to finish and return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* When no interpolation is needed, then we can just copy the input into
   the output. */
static gal_data_t *
interpolate_copy_input(gal_data_t *input, int aslinkedlist)
{
  gal_data_t *tin, *tout;

  /* Make a copy of the first input. */
  tout=gal_data_copy(input);
  tout->next=NULL;

  /* If we have a linked list, copy each element. */
  if(aslinkedlist)
    {
      /* Note that we have already copied the first input. */
      for(tin=input->next; tin!=NULL; tin=tin->next)
        {
          /* Copy this dataset (will also copy flags). */
          tout->next=gal_data_copy(tin);
          tout=tout->next;
        }

      /* Output is the reverse of the input, so reverse it. */
      gal_data_reverse_ll(&tout);
    }

  /* Return the copied list. */
  return tout;
}





/* Interpolate blank values in an array. If the `tl!=NULL', then it is
   assumed that the tile values correspond to given tessellation. Such that
   `input[i]' corresponds to `tiles[i]' in the tessellation. */
gal_data_t *
gal_interpolate_close_neighbors(gal_data_t *input,
                                struct gal_tile_two_layer_params *tl,
                                size_t numneighbors, size_t numthreads,
                                int onlyblank, int aslinkedlist)
{
  gal_data_t *tin, *tout;
  struct interpolate_params prm;
  size_t ngbvnum=numthreads*numneighbors;
  int permute=(tl && tl->totchannels>1 && tl->workoverch);


  /* If there are no blank values in the array we should only fill blank
     values, then simply copy the input and abort. */
  if( (input->flag | GAL_DATA_FLAG_BLANK_CH)     /* Zero bit is meaningful.*/
      && !(input->flag | GAL_DATA_FLAG_HASBLANK) /* There are no blanks.   */
      && onlyblank )                             /* Only interpolate blank.*/
    return interpolate_copy_input(input, aslinkedlist);


  /* Initialize the constant parameters. */
  prm.tl           = tl;
  prm.ngb_vals     = NULL;
  prm.input        = input;
  prm.onlyblank    = onlyblank;
  prm.numneighbors = numneighbors;
  prm.num          = aslinkedlist ? gal_data_num_in_ll(input) : 1;


  /* Flag the blank values. */
  prm.blanks=gal_blank_flag(input);


  /* If the input is from a tile structure and the user has asked to ignore
     channels, then re-order the values. */
  if(permute)
    {
      /* Prepare the permutation (if necessary/not already defined). */
      gal_tile_full_permutation(tl);

      /* Re-order values to ignore channels (if necessary). */
      gal_permutation_apply(input, tl->permutation);
      gal_permutation_apply(prm.blanks, tl->permutation);

      /* If this is a linked list, then permute remaining nodes. */
      if(aslinkedlist)
        for(tin=input->next; tin!=NULL; tin=tin->next)
          gal_permutation_apply(tin, tl->permutation);
    }


  /* Allocate space for the (first) output. */
  prm.out=gal_data_alloc(NULL, input->type, input->ndim, input->dsize,
                         input->wcs, 0, input->minmapsize, NULL,
                         input->unit, NULL);
  gal_linkedlist_add_to_vll(&prm.ngb_vals,
                            gal_data_malloc_array(input->type, ngbvnum));


  /* If we are given a list of datasets, make the necessary
     allocations. The reason we are doing this after a check of
     `aslinkedlist' is that the `input' might have a `next' element, but
     the caller might not have called `aslinkedlist'. */
  prm.out->next=NULL;
  if(aslinkedlist)
    for(tin=input->next; tin!=NULL; tin=tin->next)
      {
        /* A small sanity check. */
        if( gal_data_dsize_is_different(input, tin) )
          error(EXIT_FAILURE, 0, "The other datasets in the list must "
                "have the same dimension and size in "
                "`gal_interpolate_close_neighbors'");

        /* Allocate the output array for this node. */
        gal_data_add_to_ll(&prm.out, NULL, tin->type, tin->ndim, tin->dsize,
                           tin->wcs, 0, tin->minmapsize, NULL, tin->unit,
                           NULL);

        /* Allocate the space for the neighbor values of this input. */
        gal_linkedlist_add_to_vll(&prm.ngb_vals,
                                  gal_data_malloc_array(tin->type, ngbvnum));
      }
  gal_data_reverse_ll(&prm.out);
  gal_linkedlist_reverse_vll(&prm.ngb_vals);


  /* Allocate space for all the flag values of all the threads here (memory
     in each thread is limited) and this is cleaner. */
  prm.thread_flags=gal_data_malloc_array(GAL_TYPE_UINT8,
                                         numthreads*input->size);


  /* Spin off the threads. */
  gal_threads_spin_off(interpolate_close_neighbors_on_thread, &prm,
                       input->size, numthreads);


  /* If the values were permuted for the interpolation, then re-order the
     values back to their original location (so they correspond to their
     tile indexs. */
  if(permute)
    {
      gal_permutation_apply_inverse(input, tl->permutation);
      for(tout=prm.out; tout!=NULL; tout=tout->next)
        gal_permutation_apply_inverse(tout, tl->permutation);
    }


  /* The interpolated array doesn't have blank values. So set the blank
     flag to 0 and set the use-zero to 1. */
  for(tout=prm.out; tout!=NULL; tout=tout->next)
    {
      tout->flag |= GAL_DATA_FLAG_BLANK_CH;
      tout->flag &= ~GAL_DATA_FLAG_HASBLANK;
    }


  /* Clean up and return. */
  free(prm.thread_flags);
  gal_data_free(prm.blanks);
  gal_linkedlist_free_vll(prm.ngb_vals, 1);
  return prm.out;
}
