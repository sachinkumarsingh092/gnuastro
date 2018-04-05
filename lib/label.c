/*********************************************************************
label -- Work on labeled (integer valued) datasets.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018, Free Software Foundation, Inc.

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

#include <gnuastro/list.h>
#include <gnuastro/qsort.h>
#include <gnuastro/label.h>
#include <gnuastro/dimension.h>






/****************************************************************
 *****************   Over segmentation       ********************
 ****************************************************************/
/* Over-segment the region specified by its indexs into peaks and their
   respective regions (clumps). This is very similar to the immersion
   method of Vincent & Soille(1991), but here, we will not separate the
   image into layers, instead, we will work based on the ordered flux
   values. If a certain pixel (at a certain level) has no neighbors, it is
   a local maximum and will be assigned a new label. If it has a labeled
   neighbor, it will take that label and if there is more than one
   neighboring labeled region that pixel will be a `river` pixel.

   DON'T FORGET: SET THE FLAGS FOR CONV EQUAL TO INPUT IN SEGMENT.

*/
size_t
gal_label_oversegment(gal_data_t *input, gal_data_t *indexs,
                      gal_data_t *label, size_t *topinds)
{
  size_t ndim=input->ndim;

  float *arr=input->array;
  gal_list_sizet_t *Q=NULL, *cleanup=NULL;
  size_t *a, *af, ind, *dsize=input->dsize;
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  int32_t n1, nlab, rlab, curlab=1, *clabel=label->array;

  /*********************************************
   For checks and debugging:*
  gal_data_t *crop;
  size_t extcount=1;
  int32_t *cr, *crf;
  size_t checkdsize[2]={10,10};
  size_t checkstart[2]={50,145};
  char *filename="clumpbuild.fits";
  size_t checkstartind=gal_dimension_coord_to_index(2, dsize, checkstart);
  gal_data_t *tile=gal_data_alloc(gal_data_ptr_increment(arr, checkstartind,
                                                         input->type),
                                  GAL_TYPE_INVALID, 2, checkdsize,
                                  NULL, 0, 0, NULL, NULL, NULL);
  tile->block=input;
  gal_checkset_writable_remove(filename, 0, 0);
  if(p->cp.numthreads!=1)
    error(EXIT_FAILURE, 0, "in the debugging mode of `clumps_oversegment' "
          "only one thread must be used");
  crop=gal_data_copy(tile);
  gal_fits_img_write(crop, filename, NULL, PROGRAM_NAME);
  gal_data_free(crop);
  printf("blank: %u\nriver: %u\ntmpcheck: %u\ninit: %u\n",
         (int32_t)GAL_BLANK_INT32, (int32_t)GAL_LABEL_RIVER,
         (int32_t)GAL_LABEL_TMPCHECK, (int32_t)GAL_LABEL_INIT);
  tile->array=gal_tile_block_relative_to_other(tile, clabel);
  tile->block=clabel;
  **********************************************/


  /* If the size of the indexs is zero, then this function is pointless. */
  if(indexs->size==0) return 0;


  /* Sort the given indexs based on their flux (`gal_qsort_index_arr' is
     defined as static in `gnuastro/qsort.h') */
  gal_qsort_index_arr=input->array;
  qsort(indexs->array, indexs->size, sizeof(size_t),
        gal_qsort_index_float_decreasing);


  /* Initialize the region we want to over-segment. */
  af=(a=indexs->array)+indexs->size;
  do clabel[*a]=GAL_LABEL_INIT; while(++a<af);


  /* Go over all the given indexs and pull out the clumps. */
  af=(a=indexs->array)+indexs->size;
  do
    /* When regions of a constant flux or masked regions exist, some later
       indexs (although they have same flux) will be filled before hand. If
       they are done, there is no need to do them again. */
    if(clabel[*a]==GAL_LABEL_INIT)
      {
        /* It might happen where one or multiple regions of the pixels
           under study have the same flux. So two equal valued pixels of
           two separate (but equal flux) regions will fall immediately
           after each other in the sorted list of indexs and we have to
           account for this.

           Therefore, if we see that the next pixel in the index list has
           the same flux as this one, it does not guarantee that it should
           be given the same label. Similar to the breadth first search
           algorithm for finding connected components, we will search all
           the neighbours and the neighbours of those neighbours that have
           the same flux of this pixel to see if they touch any label or
           not and to finally give them all the same label. */
        if( (a+1)<af && arr[*a]==arr[*(a+1)] )
          {
            /* Label of first neighbor found. */
            n1=0;

            /* A small sanity check. */
            if(Q!=NULL || cleanup!=NULL)
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so "
                    "we can fix this problem. `Q' and `cleanup' should be "
                    "NULL but while checking the equal flux regions they "
                    "aren't", __func__, PACKAGE_BUGREPORT);

            /* Add this pixel to a queue. */
            gal_list_sizet_add(&Q, *a);
            gal_list_sizet_add(&cleanup, *a);
            clabel[*a] = GAL_LABEL_TMPCHECK;

            /* Find all the pixels that have the same flux and are
               connected. */
            while(Q!=NULL)
              {
                /* Pop an element from the queue. */
                ind=gal_list_sizet_pop(&Q);

                /* Look at the neighbors and see if we already have a
                   label. */
                GAL_DIMENSION_NEIGHBOR_OP(ind, ndim, dsize, ndim, dinc,
                   {
                     /* If it is already decided to be a river, then stop
                        looking at the neighbors. */
                     if(n1!=GAL_LABEL_RIVER)
                       {
                         /* For easy reading. */
                         nlab=clabel[ nind ];

                         /* This neighbor's label isn't zero. */
                         if(nlab)
                           {
                             /* If this neighbor has not been labeled yet
                                and has an equal flux, add it to the queue
                                to expand the studied region.*/
                             if( nlab==GAL_LABEL_INIT && arr[nind]==arr[*a] )
                               {
                                 clabel[nind]=GAL_LABEL_TMPCHECK;
                                 gal_list_sizet_add(&Q, nind);
                                 gal_list_sizet_add(&cleanup, nind);
                               }
                             else
                               n1=( nlab>0

                                    /* If this neighbor has a positive
                                       nlab, it belongs to another object,
                                       so if `n1' has not been set for the
                                       whole region (n1==0), put `nlab'
                                       into `n1'. If `n1' has been set and
                                       is different from `nlab' then this
                                       whole equal flux region should be a
                                       wide river because it is connecting
                                       two connected regions.*/
                                    ? ( n1
                                        ? (n1==nlab ? n1 : GAL_LABEL_RIVER)
                                        : nlab )

                                    /* If the data has blank pixels (recall
                                       that blank in int32 is negative),
                                       see if the neighbor is blank and if
                                       so, set the label to a river. Since
                                       the flag checking can be done
                                       outside this loop, for datasets with
                                       no blank element this last step will
                                       be completley ignored. */
                                    : ( ( (input->flag
                                           & GAL_DATA_FLAG_HASBLANK)
                                          && nlab==GAL_BLANK_INT32 )
                                        ? GAL_LABEL_RIVER : n1 ) );
                           }

                         /* If this neigbour has a label of zero, then we
                            are on the edge of the indexed region (the
                            neighbor is not in the initial list of pixels
                            to segment). When over-segmenting the noise and
                            the detections, `clabel' is zero for the parts
                            of the image that we are not interested in
                            here. */
                         else clabel[*a]=GAL_LABEL_RIVER;
                       }
                   } );
              }

            /* Set the label that is to be given to this equal flux
               region. If `n1' was set to any value, then that label should
               be used for the whole region. Otherwise, this is a new
               label, see the case for a non-flat region. */
            if(n1) rlab = n1;
            else
              {
                rlab = curlab++;
                if( topinds )              /* This is a local maximum of   */
                  topinds[rlab]=*a;        /* this region, save its index. */
              }

            /* Give the same label to the whole connected equal flux
               region, except those that might have been on the side of
               the image and were a river pixel. */
            while(cleanup!=NULL)
              {
                ind=gal_list_sizet_pop(&cleanup);
                /* If it was on the sides of the image, it has been
                   changed to a river pixel. */
                if( clabel[ ind ]==GAL_LABEL_TMPCHECK ) clabel[ ind ]=rlab;
              }
          }

        /* The flux of this pixel is not the same as the next sorted
           flux, so simply find the label for this object. */
        else
          {
            /* `n1' is the label of the first labeled neighbor found, so
               we'll initialize it to zero. */
            n1=0;

            /* Go over all the fully connected neighbors of this pixel and
               see if all the neighbors (with maximum connectivity: the
               number of dimensions) that have a non-macro value belong to
               one label or not. If the pixel is neighboured by more than
               one label, set it as a river pixel. Also if it is touching a
               zero valued pixel (which does not belong to this object),
               set it as a river pixel.*/
            GAL_DIMENSION_NEIGHBOR_OP(*a, ndim, dsize, ndim, dinc,
               {
                 /* When `n1' has already been set as a river, there is no
                    point in looking at the other neighbors. */
                 if(n1!=GAL_LABEL_RIVER)
                   {
                     /* For easy reading. */
                     nlab=clabel[ nind ];

                     /* If this neighbor is on a non-processing label, then
                        set the first neighbor accordingly. Note that we
                        also want the zero valued neighbors (detections if
                        working on sky, and sky if working on detection):
                        we want rivers between the two domains. */
                     n1 = ( nlab

                            /* nlab is non-zero. */
                            ? ( nlab>0

                                /* Neighbor has a meaningful label, so
                                   check with any previously found labeled
                                   neighbors. */
                                ? ( n1
                                    ? ( nlab==n1 ? n1 : GAL_LABEL_RIVER )
                                    : nlab )

                                /* If the dataset has blank values and this
                                   neighbor is blank, then the pixel should
                                   be a river. Note that the blank checking
                                   can be optimized out, so if the input
                                   doesn't have blank values,
                                   `nlab==GAL_BLANK_INT32' will never be
                                   checked. */
                                : ( (input->flag & GAL_DATA_FLAG_HASBLANK)
                                    && nlab==GAL_BLANK_INT32
                                    ? GAL_LABEL_RIVER : n1 ) )

                            /* `nlab==0' (the neighbor lies in the other
                               domain (sky or detections). To avoid the
                               different domains touching, this pixel
                               should be a river. */
                            : GAL_LABEL_RIVER );
                   }
               });

            /* Either assign a new label to this pixel, or give it the one
               of its neighbors. If n1 equals zero, then this is a new
               peak, and a new label should be created.  But if n1!=0, it
               is either a river pixel (has more than one labeled neighbor
               and has been set to `GAL_LABEL_RIVER' before) or all its
               neighbors have the same label. In both such cases, rlab
               should be set to n1.*/
            if(n1) rlab = n1;
            else
              {
                rlab = curlab++;
                if( topinds )
                  topinds[ rlab ]=*a;
              }

            /* Put the found label in the pixel. */
            clabel[ *a ] = rlab;
          }

        /*********************************************
         For checks and debugging:
        if(    *a / dsize[1] >= checkstart[0]
            && *a / dsize[1] <  checkstart[0] + checkdsize[0]
            && *a % dsize[1] >= checkstart[1]
            && *a % dsize[1] <  checkstart[1] + checkdsize[1] )
          {
            printf("%zu (%zu: %zu, %zu): %u\n", ++extcount, *a,
                   (*a%dsize[1])-checkstart[1], (*a/dsize[1])-checkstart[0],
                   clabel[*a]);
            crop=gal_data_copy(tile);
            crf=(cr=crop->array)+crop->size;
            do if(*cr==GAL_LABEL_RIVER) *cr=0; while(++cr<crf);
            gal_fits_img_write(crop, filename, NULL, PROGRAM_NAME);
            gal_data_free(crop);
          }
        **********************************************/
      }
  while(++a<af);

  /*********************************************
   For checks and debugging:
  tile->array=NULL;
  gal_data_free(tile);
  printf("Total number of clumps: %u\n", curlab-1);
  **********************************************/

  /* Clean up. */
  free(dinc);

  /* Return the total number of clumps. */
  return curlab-1;
}





/* Grow the given labels over the list of indexs. In over-segmentation
   (watershed algorithm), the indexs have to be sorted by pixel value and
   local maximums (that aren't touching any already labeled pixel) get a
   separate label. However, here the final number of labels will not
   change. All pixels that aren't directly touching a labeled pixel just
   get pushed back to the start of the loop, and the loop iterates until
   its size doesn't change any more. This is because in a generic scenario
   some of the indexed pixels might not be reachable. As a result, it is
   not necessary for the given indexs to be sorted here.

   This function looks for positive-valued neighbors and will label a pixel
   if it touches one. Therefore, it is is very important that only
   pixels/labels that are intended for growth have positive values before
   calling this function. Any non-positive (zero or negative) value will be
   ignored by this function. Thus, it is recommended that while filling in
   the `indexs' array values, you initialize all those pixels with
   `GAL_LABEL_INIT', and set non-labeled pixels that you don't want to grow
   to 0.

   This function will write into both the input datasets. After this
   function, some of the non-positive `labels' pixels will have a new
   positive label and the number of useful elements in `indexs' will have
   decreased. Those indexs that couldn't be labeled will remain inside
   `indexs'. If `withrivers' is non-zero, then pixels that are immediately
   touching more than one positive value will be given a `GAL_LABEL_RIVER'
   label.

   Note that the `indexs->array' is not re-allocated to its new size at the
   end. But since `indexs->dsize[0]' and `indexs->size' have new values,
   the extra elements just won't be used until they are ultimately freed by
   `gal_data_free'.

   Input:

     labels: The labels array that must be operated on. The pixels that
             must be "grown" must have the value `GAL_LABEL_INIT' (negative).

     sorted_indexs: The indexs of the pixels that must be grown.

     withrivers: as described above.

     connectivity: connectivity to define neighbors for growth.  */
void
gal_label_grow_indexs(gal_data_t *labels, gal_data_t *indexs, int withrivers,
                      int connectivity)
{
  int searchngb;
  size_t *iarray=indexs->array;
  int32_t n1, nlab, *olabel=labels->array;
  size_t *s, *sf, thisround, ninds=indexs->size;
  size_t *dinc=gal_dimension_increment(labels->ndim, labels->dsize);

  /* Some basic sanity checks: */
  if(labels->type!=GAL_TYPE_INT32)
    error(EXIT_FAILURE, 0, "%s: `labels' has to have type of int32_t",
          __func__);
  if(indexs->type!=GAL_TYPE_SIZE_T)
    error(EXIT_FAILURE, 0, "%s: `indexs' must be `size_t' type but it is "
          "`%s'", __func__, gal_type_name(indexs->type,1));
  if(indexs->ndim!=1)
    error(EXIT_FAILURE, 0, "%s: `indexs' has to be a 1D array, but it is "
          "%zuD", __func__, indexs->ndim);

  /* The basic idea is this: after growing, not all the blank pixels are
     necessarily filled, for example the pixels might belong to two regions
     above the growth threshold. So the pixels in between them (which are
     below the threshold will not ever be able to get a label). Therefore,
     the safest way we can terminate the loop of growing the objects is to
     stop it when the number of pixels left to fill in this round
     (thisround) equals the number of blanks.

     To start the loop, we set `thisround' to one more than the number of
     indexed pixels. Note that it will be corrected immediately after the
     loop has started, it is just important to pass the `while'. */
  thisround=ninds+1;
  while( thisround > ninds )
    {
      /* `thisround' will keep the number of pixels to be inspected in this
         round. `ninds' will count the number of pixels left without a
         label by the end of this round. Since `ninds' comes from the
         previous loop (or outside, for the first round) it has to be saved
         in `thisround' to begin counting a fresh. */
      thisround=ninds;
      ninds=0;

      /* Go over all the available indexs. NOTE: while the `indexs->array'
         pointer remains unchanged, `indexs->size' can/will change (get
         smaller) in every loop. */
      sf = (s=indexs->array) + indexs->size;
      do
        {
          /* We'll begin by assuming the nearest neighbor of this pixel
             has no label (has a value of 0). */
          n1=0;

          /* Check the neighbors of this pixel. Note that since this
             macro has multiple loops within it, we can't use
             break. We'll use the `searchngb' variable instead. */
          searchngb=1;
          GAL_DIMENSION_NEIGHBOR_OP(*s, labels->ndim, labels->dsize,
            connectivity, dinc,
            {
              if(searchngb)
                {
                  /* For easy reading. */
                  nlab = olabel[nind];

                  /* This neighbor's label is meaningful. */
                  if(nlab>0)                   /* This is a real label. */
                    {
                      if(n1)       /* A prev. ngb label has been found. */
                        {
                          if( n1 != nlab )    /* Different label from   */
                            {    /* prevously found ngb for this pixel. */
                              n1=GAL_LABEL_RIVER;
                              searchngb=0;
                            }
                        }
                      else
                        {   /* This is the first labeld neighbor found. */
                          n1=nlab;

                          /* If we want to completely fill in the region
                             (`withrivers==0'), then there is no point in
                             looking in other neighbors, the first
                             neighbor we find, is the one we'll use. */
                          if(!withrivers) searchngb=0;
                        }
                    }
                }
            } );

          /* The loop over neighbors (above) finishes with three
             possibilities:

             n1==0                    --> No labeled neighbor was found.
             n1==GAL_LABEL_RIVER      --> Connecting two labeled regions.
             n1>0                     --> Only has one neighbouring label.

             The first one means that no neighbors were found and this
             pixel should be kept for the next loop (we'll be growing the
             objects pixel-layer by pixel-layer). In the other two cases,
             we just need to write in the value of `n1'. */
          if(n1)
            {
              /* Set the label. */
              olabel[*s]=n1;

              /* If this pixel is a river (can only happen when
                 `withrivers' is zero), keep it in the loop, because we
                 want the `indexs' dataset to contain all non-positive
                 (non-labeled) pixels, including rivers. */
              if(n1==GAL_LABEL_RIVER)
                iarray[ ninds++ ] = *s;
            }
          else
            iarray[ ninds++ ] = *s;

          /* Correct the size of the `indexs' dataset. */
          indexs->size = indexs->dsize[0] = ninds;
        }
      while(++s<sf);
    }

  /* Clean up. */
  free(dinc);
}
