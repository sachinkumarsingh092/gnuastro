/*********************************************************************
NoiseChisel - Detect and segment signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/qsort.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "ui.h"
#include "clumps.h"
#include "threshold.h"










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
   neighboring labeled region that pixel will be a `river` pixel. */
void
clumps_oversegment(struct clumps_thread_params *cltprm)
{
  struct noisechiselparams *p=cltprm->clprm->p;
  size_t ndim=p->input->ndim;

  float *arr=p->conv->array;
  gal_data_t *indexs=cltprm->indexs;
  gal_list_sizet_t *Q=NULL, *cleanup=NULL;
  size_t *a, *af, ind, *dsize=p->input->dsize;
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  int32_t n1, nlab, rlab, curlab=1, *clabel=p->clabel->array;

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
                                                         p->conv->type),
                                  GAL_TYPE_INVALID, 2, checkdsize,
                                  NULL, 0, 0, NULL, NULL, NULL);
  tile->block=p->conv;
  gal_checkset_writable_remove(filename, 0, 0);
  if(p->cp.numthreads!=1)
    error(EXIT_FAILURE, 0, "in the debugging mode of `clumps_oversegment' "
          "only one thread must be used");
  crop=gal_data_copy(tile);
  gal_fits_img_write(crop, filename, NULL, PROGRAM_NAME);
  gal_data_free(crop);
  printf("blank: %u\nriver: %u\ntmpcheck: %u\ninit: %u\nmaxlab: %u\n",
         (int32_t)GAL_BLANK_INT32, (int32_t)CLUMPS_RIVER,
         (int32_t)CLUMPS_TMPCHECK, (int32_t)CLUMPS_INIT,
         (int32_t)CLUMPS_MAXLAB);
  tile->array=gal_tile_block_relative_to_other(tile, p->clabel);
  tile->block=p->clabel;
  **********************************************/


  /* If the size of the indexs is zero, then this function is pointless. */
  if(indexs->size==0) { cltprm->numinitclumps=0; return; }


  /* Sort the given indexs based on their flux (`gal_qsort_index_arr' is
     defined as static in `gnuastro/qsort.h') */
  gal_qsort_index_arr=p->conv->array;
  qsort(indexs->array, indexs->size, sizeof(size_t),
        gal_qsort_index_float_decreasing);


  /* Initialize the region we want to over-segment. */
  af=(a=indexs->array)+indexs->size; do clabel[*a]=CLUMPS_INIT; while(++a<af);


  /* Go over all the given indexs and pull out the clumps. */
  af=(a=indexs->array)+indexs->size;
  do
    /* When regions of a constant flux or masked regions exist, some later
       indexs (although they have same flux) will be filled before hand. If
       they are done, there is no need to do them again. */
    if(clabel[*a]==CLUMPS_INIT)
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
            clabel[*a] = CLUMPS_TMPCHECK;

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
                     if(n1!=CLUMPS_RIVER)
                       {
                         /* For easy reading. */
                         nlab=clabel[ nind ];

                         /* This neighbor's label isn't zero. */
                         if(nlab)
                           {
                             /* If this neighbor has not been labeled yet
                                and has an equal flux, add it to the queue
                                to expand the studied region.*/
                             if( nlab==CLUMPS_INIT && arr[nind]==arr[*a] )
                               {
                                 clabel[nind]=CLUMPS_TMPCHECK;
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
                                        ? (n1==nlab ? n1 : CLUMPS_RIVER)
                                        : nlab )

                                    /* If the data has blank pixels (recall
                                       that blank in int32 is negative),
                                       see if the neighbor is blank and if
                                       so, set the label to a river. Since
                                       the flag checking can be done
                                       outside this loop, for datasets with
                                       no blank element this last step will
                                       be completley ignored. */
                                    : ( ( (p->input->flag
                                           & GAL_DATA_FLAG_HASBLANK)
                                          && nlab==GAL_BLANK_INT32 )
                                        ? CLUMPS_RIVER : n1 ) );
                           }

                         /* If this neigbour has a label of zero, then we
                            are on the edge of the indexed region (the
                            neighbor is not in the initial list of pixels
                            to segment). When over-segmenting the noise and
                            the detections, `clabel' is zero for the parts
                            of the image that we are not interested in
                            here. */
                         else clabel[*a]=CLUMPS_RIVER;
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
                if( cltprm->topinds )       /* This is a local maximum of  */
                  cltprm->topinds[rlab]=*a; /* this region, save its index.*/
              }

            /* Give the same label to the whole connected equal flux
               region, except those that might have been on the side of
               the image and were a river pixel. */
            while(cleanup!=NULL)
              {
                ind=gal_list_sizet_pop(&cleanup);
                /* If it was on the sides of the image, it has been
                   changed to a river pixel. */
                if( clabel[ ind ]==CLUMPS_TMPCHECK ) clabel[ ind ]=rlab;
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
               number of dimensions) that have a non-macro value (less than
               CLUMPS_MAXLAB) belong to one label or not. If the pixel is
               neighboured by more than one label, set it as a river
               pixel. Also if it is touching a zero valued pixel (which
               does not belong to this object), set it as a river pixel.*/
            GAL_DIMENSION_NEIGHBOR_OP(*a, ndim, dsize, ndim, dinc,
               {
                 /* When `n1' has already been set as a river, there is no
                    point in looking at the other neighbors. */
                 if(n1!=CLUMPS_RIVER)
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
                                    ? ( nlab==n1 ? n1 : CLUMPS_RIVER )
                                    : nlab )

                                /* If the dataset has blank values and this
                                   neighbor is blank, then the pixel should
                                   be a river. Note that the blank checking
                                   can be optimized out, so if the input
                                   doesn't have blank values,
                                   `nlab==GAL_BLANK_INT32' will never be
                                   checked. */
                                : ( (p->input->flag & GAL_DATA_FLAG_HASBLANK)
                                    && nlab==GAL_BLANK_INT32
                                    ? CLUMPS_RIVER : n1 ) )

                            /* `nlab==0' (the neighbor lies in the other
                               domain (sky or detections). To avoid the
                               different domains touching, this pixel
                               should be a river. */
                            : CLUMPS_RIVER );
                   }
               });

            /* Either assign a new label to this pixel, or give it the one
               of its neighbors. If n1 equals zero, then this is a new
               peak, and a new label should be created.  But if n1!=0, it
               is either a river pixel (has more than one labeled neighbor
               and has been set to `CLUMPS_RIVER' before) or all its
               neighbors have the same label. In both such cases, rlab
               should be set to n1.*/
            if(n1) rlab = n1;
            else
              {
                rlab = curlab++;
                if( cltprm->topinds )
                  cltprm->topinds[ rlab ]=*a;
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
            do if(*cr==CLUMPS_RIVER) *cr=0; while(++cr<crf);
            gal_fits_img_write(crop, filename, NULL, PROGRAM_NAME);
            gal_data_free(crop);
          }
        **********************************************/
      }
  while(++a<af);

  /* Save the total number of clumps. */
  cltprm->numinitclumps=curlab-1;

  /* Set all the river pixels to zero. Note that this is only necessary for
     the detected clumps. When finding clumps over the Sky, we will be
     going over the full tile and removing rivers after this function. This
     is because, we set the borders of the tile to a river value and didn't
     include them in the list of indexs. */
  if(cltprm->clprm->sky0_det1)
    {
      af=(a=indexs->array)+indexs->size;
      do if( clabel[*a]==CLUMPS_RIVER ) clabel[*a]=CLUMPS_INIT; while(++a<af);
    }

  /*********************************************
   For checks and debugging:
  tile->array=NULL;
  gal_data_free(tile);
  printf("Total number of clumps: %u\n", curlab-1);
  **********************************************/

  /* Clean up. */
  free(dinc);
}




















/**********************************************************************/
/*****************              Grow clumps           *****************/
/**********************************************************************/
/* Make the preparations for the intiial growing the clumps to identify
   objects: a single standard deviation for the whole object and preparing
   the labels (because the growth is going to happen on the `olabel'
   image. */
void
clumps_grow_prepare_initial(struct clumps_thread_params *cltprm)
{
  gal_data_t *indexs=cltprm->indexs;
  gal_data_t *input=cltprm->clprm->p->input;
  struct noisechiselparams *p=cltprm->clprm->p;

  size_t *s, *sf, *dsize=input->dsize;
  size_t ndiffuse=0, coord[2], *dindexs;
  double wcoord[2]={0.0f,0.0f}, brightness=0.0f;
  float glimit, *imgss=input->array, *std=p->std->array;
  int32_t *olabel=p->olabel->array, *clabel=p->clabel->array;


  /* Find the flux weighted center (meaningful only for positive valued
     pixels). */
  sf=(s=indexs->array)+indexs->size;
  do
    if( imgss[ *s ] > 0.0f )
      {
        brightness += imgss[ *s ];
        wcoord[0]  += imgss[ *s ] * (*s/dsize[1]);
        wcoord[1]  += imgss[ *s ] * (*s%dsize[1]);
      }
  while(++s<sf);


  /* Calculate the center, if no pixels were positive, use the
     geometric center (irrespective of flux). */
  if(brightness==0.0f)
    {
      sf=(s=indexs->array)+indexs->size;
      do
        {
          wcoord[0] += *s / dsize[1];
          wcoord[1] += *s % dsize[1];
        }
      while(++s<sf);
      brightness = indexs->size;
    }


  /* Convert floatint point coordinates to FITS integers. */
  coord[0] = GAL_DIMENSION_FLT_TO_INT(wcoord[0]);
  coord[1] = GAL_DIMENSION_FLT_TO_INT(wcoord[1]);


  /* Find the growth limit  */
  cltprm->std = std[ gal_tile_full_id_from_coord(&p->cp.tl, coord) ];
  glimit = p->gthresh * cltprm->std;


  /* Allocate space to keep the diffuse indexs over this detection. We need
     to keep the actual indexs since it is our only connection to the
     object at this stage: we are also going to re-label the pixels to
     grow. For most astronomical objects, the major part of the detection
     area is going to be diffuse flux, so we will just allocate the same
     size as `indexs' array (the `dsize' will be corrected after getting
     the exact number. */
  cltprm->diffuseindexs=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1,
                                       cltprm->indexs->dsize, NULL, 0,
                                       p->cp.minmapsize, NULL, NULL, NULL);
  dindexs=cltprm->diffuseindexs->array;
  sf=(s=indexs->array)+indexs->size;
  do
    {
      olabel[*s] = clabel[*s];
      if( clabel[*s]==CLUMPS_INIT )
        if( imgss[*s]>glimit ) dindexs[ ndiffuse++ ] = *s;
    }
  while(++s<sf);


  /* Correct the sizes of the `diffuseindexs' data structure. */
  cltprm->diffuseindexs->size = cltprm->diffuseindexs->dsize[0] = ndiffuse;
}





/* Add all the remaining pixels in the detection (below the growth
   threshold, or those that were not touching). Note that initially
   `diffuseindexs' was filled with the pixels that are above the growth
   threshold. That was necessary for identifying the objects. Now that we
   have identified the objects and labeled them, we want to add the
   remaining diffuse pixels to it too before doing the final growth.

   Note that the most efficient way is just to re-fill the `diffuseindexs'
   array instead of adding the pixels below the threshold and sorting them
   afterwards.*/
void
clumps_grow_prepare_final(struct clumps_thread_params *cltprm)
{
  size_t ndiffuse=0;
  size_t *dindexs=cltprm->diffuseindexs->array;
  int32_t *olabel=cltprm->clprm->p->olabel->array;
  size_t *s=cltprm->indexs->array, *sf=s+cltprm->indexs->size;

  /* Recall that we initially allocated `diffuseindexs' to have the same
     size as the indexs. So there is no problem if there are more pixels in
     this final round compared to the initial round. */
  do
    if( olabel[*s] < 0 )
      dindexs[ ndiffuse++ ] = *s;
  while(++s<sf);

  /* Correct the sizes of the `diffuseindexs' data structure. */
  cltprm->diffuseindexs->size = cltprm->diffuseindexs->dsize[0] = ndiffuse;
}





/* Grow the true clumps over the diffuse regions of a detection. Note that
   unlike before, were river pixels would get a separate label for them
   selves, here, they don't, they just get set back to SEGMENTINIT. This is
   because some of the pixels that lie immediately between two labeled
   regions might not be in the blankinds array (they were below the
   threshold). So we have to find river pixels later on, after the growth
   is done independently.

   This function is going to be used before identifying objects and also
   after it (to completely fill in the diffuse area). The distinguishing
   point between these two steps is the presence of rivers, so you can use
   the `withrivers' argument. */
void
clumps_grow(struct clumps_thread_params *cltprm, int withrivers)
{
  struct noisechiselparams *p=cltprm->clprm->p;
  gal_data_t *diffuseindexs=cltprm->diffuseindexs;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  int searchngb;
  size_t *diarray=cltprm->diffuseindexs->array;
  int32_t n1, nlab, *olabel=p->olabel->array;
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  size_t *s, *sf, thisround, ndiffuse=diffuseindexs->size;

  /* The basic idea is this: after growing, not all the blank pixels are
     necessarily filled, for example the pixels might belong to two regions
     above the growth threshold. So the pixels in between them (which are
     below the threshold will not ever be able to get a label). Therefore,
     the safest way we can terminate the loop of growing the objects is to
     stop it when the number of pixels left to fill in this round
     (thisround) equals the number of blanks.

     To start the loop, we set `thisround' to one more than the number of
     diffuse pixels. Note that it will be corrected immediately after the
     loop has started, it is just important to pass the `while'. */
  thisround=ndiffuse+1;
  while( thisround > ndiffuse )
    {
      /* `thisround' will keep the number of pixels to be inspected in this
         round. `ndiffuse' will count the number of pixels left without a
         label by the end of this round. Since `ndiffuse' comes from the
         previous loop (or outside, for the first round) it has to be saved
         in `thisround' to begin counting a fresh. */
      thisround=ndiffuse;
      ndiffuse=0;

      /* Go over all the available indexs. */
      sf=(s=diffuseindexs->array)+diffuseindexs->size;
      do
        {
          /* We'll begin by assuming the nearest neighbor of this pixel has
             no label (has a value of 0). */
          n1=0;

          /* Check the very closest neighbors of each pixel (4-connectivity
             in a 2D image). Note that since this macro has multiple loops
             within it, we can't use break. We'll use a variable instead. */
          searchngb=1;
          GAL_DIMENSION_NEIGHBOR_OP(*s, ndim, dsize, 1, dinc,
            {
              if(searchngb)
                {
                  /* For easy reading. */
                  nlab = olabel[nind];

                  /* This neighbor's label is meaningful. */
                  if(nlab>0)                      /* This is a real label. */
                    {
                      if(n1)  /* A prev. neighboring label has been found. */
                        {
                          if( n1 != nlab )       /* Different label from   */
                            {  /* prevously found neighbor for this pixel. */
                              n1=CLUMPS_RIVER;
                              searchngb=0;
                            }
                        }
                      else
                        {      /* This is the first labeld neighbor found. */
                          n1=nlab;

                          /* If we want to completely fill in the region
                             (`withrivers==0'), then there is no point in
                             looking in other neighbors, the first neighbor
                             we find is the one we'll use. */
                          if(!withrivers) searchngb=0;
                        }
                    }
                }
            } );

          /* The loop above over neighbors finishes with three
             possibilities:

               n1==0                    --> No labeled neighbor was found.
               n1==CLUMPS_RIVER         --> Connecting two labeled regions.
               n1>0                     --> Only has one neighbouring label.

             The first one means that no neighbors were found and this
             pixel should be kept for the next loop (we'll be growing the
             objects pixel-layer by pixel-layer). In the other two cases,
             we just need to write in the value of `n1'. */
          if(n1)
            {
              olabel[*s]=n1;
              if(withrivers && n1==CLUMPS_RIVER)   /* To keep rivers in */
                diarray[ ndiffuse++ ] = *s;        /* the diffuse list. */
            }
          else
            diarray[ ndiffuse++ ] = *s;

          /* Correct the size of the `diffuseindexs' dataset. */
          diffuseindexs->size = diffuseindexs->dsize[0] = ndiffuse;
        }
      while(++s<sf);
    }

  /* Clean up. */
  free(dinc);
}


















/**********************************************************************/
/*****************             S/N threshold          *****************/
/**********************************************************************/
/* In this function we want to find the general information for each clump
   in an over-segmented labeled array. The signal in each clump is the
   average signal inside it subtracted by the average signal in the river
   pixels around it. So this function will go over all the pixels in the
   object (already found in deblendclumps()) and add them appropriately.

   The output is an array of size cltprm->numinitial*INFO_NCOLS. as listed
   below.*/
enum infocols
  {
    INFO_X,              /* Flux weighted X center col, 0 by C std. */
    INFO_Y,              /* Flux weighted Y center col.             */
    INFO_NFF,            /* Number of non-negative pixels (for X,Y).*/
    INFO_INFLUX,         /* Tatal flux within clump.                */
    INFO_INAREA,         /* Tatal area within clump.                */
    INFO_RIVFLUX,        /* Tatal flux within rivers around clump.  */
    INFO_RIVAREA,        /* Tatal area within rivers around clump.  */
    INFO_INSTD,          /* Standard deviation at clump center.     */

    INFO_NCOLS,          /* Total number of columns.                */
  };
static void
clumps_get_raw_info(struct clumps_thread_params *cltprm)
{
  struct noisechiselparams *p=cltprm->clprm->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  size_t i, *a, *af, ii, coord[2];
  double *row, *info=cltprm->info->array;
  size_t nngb=gal_dimension_num_neighbors(ndim);
  float *arr=p->input->array, *std=p->std->array;
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  int32_t lab, nlab, *ngblabs, *clabel=p->clabel->array;

  /* Allocate the array to keep the neighbor labels of river pixels. */
  ngblabs=gal_data_malloc_array(GAL_TYPE_INT32, nngb, __func__, "ngblabs");

  /* Go over all the pixels in this region. */
  af=(a=cltprm->indexs->array)+cltprm->indexs->size;
  do
    if( !isnan(arr[ *a ]) )
      {
        /* This pixel belongs to a clump. */
        if( clabel[ *a ]>0 )
          {
            lab=clabel[*a];
            ++info[ lab * INFO_NCOLS + INFO_INAREA ];
            info[   lab * INFO_NCOLS + INFO_INFLUX ] += arr[*a];
            if( arr[*a]>0.0f )
              {
                info[ lab * INFO_NCOLS + INFO_NFF ] += arr[*a];
                info[ lab * INFO_NCOLS + INFO_X ] += arr[*a] * (*a/dsize[1]);
                info[ lab * INFO_NCOLS + INFO_Y ] += arr[*a] * (*a%dsize[1]);
              }
          }

        /* This pixel belongs to a river (has a value of zero and isn't
           blank). */
        else
          {
            /* We are on a river pixel. So the value of this pixel has to
               be added to any of the clumps in touches. But since it might
               touch a labeled region more than once, we use `ngblabs' to
               keep track of which label we have already added its value
               to. `ii` is the number of different labels this river pixel
               has already been considered for. `ngblabs' will keep the list
               labels. */
            ii=0;
            memset(ngblabs, 0, nngb*sizeof *ngblabs);

            /* Look into the 8-connected neighbors (recall that a
               connectivity of `ndim' means all pixels touching it (even on
               one vertice). */
            GAL_DIMENSION_NEIGHBOR_OP(*a, ndim, dsize, ndim, dinc, {
                /* This neighbor's label. */
                nlab=clabel[ nind ];

                /* We only want those neighbors that are not rivers (>0) or
                   any of the flag values. */
                if(nlab>0)
                  {
                    /* Go over all already checked labels and make sure
                       this clump hasn't already been considered. */
                    for(i=0;i<ii;++i) if(ngblabs[i]==nlab) break;

                    /* This neighbor clump hasn't been considered yet: */
                    if(i==ii)
                      {
                        ngblabs[ii++] = nlab;
                        ++info[ nlab * INFO_NCOLS + INFO_RIVAREA ];
                        info[   nlab * INFO_NCOLS + INFO_RIVFLUX ] += arr[*a];
                      }
                  }
              } );
          }
      }
  while(++a<af);


  /* Do the final preparations. All the calculations are only necessary for
     the clumps that satisfy the minimum area. So there is no need to waste
     time on the smaller ones. */
  for(lab=1; lab<=cltprm->numinitclumps; ++lab)
    {
      row = &info [ lab * INFO_NCOLS ];
      if ( row[INFO_INAREA] > p->segsnminarea )
        {
          /* Especially over the undetected regions, it might happen that
             none of the pixels were positive. In that case, set the total
             area of the clump to zero so it is no longer considered.*/
          if( row[INFO_NFF]==0.0f ) row[INFO_INAREA]=0;
          else
            {
              coord[0]=GAL_DIMENSION_FLT_TO_INT(row[INFO_X]/row[INFO_NFF]);
              coord[1]=GAL_DIMENSION_FLT_TO_INT(row[INFO_Y]/row[INFO_NFF]);
              row[INFO_INSTD] = std[ gal_tile_full_id_from_coord(&p->cp.tl,
                                                                 coord) ];
              /* For a check
              printf("---------\n");
              printf("\t%f --> %zu\n", row[INFO_Y]/row[INFO_NFF], coord[1]);
              printf("\t%f --> %zu\n", row[INFO_X]/row[INFO_NFF], coord[0]);
              printf("%u: (%zu, %zu): %.3f\n", lab, coord[1]+1,
                     coord[0]+1, row[INFO_INSTD]);
              */
            }
        }
    }

  /* Clean up. */
  free(dinc);
  free(ngblabs);
}




/* Make an S/N table for the clumps in a given region. */
void
clumps_make_sn_table(struct clumps_thread_params *cltprm)
{
  struct noisechiselparams *p=cltprm->clprm->p;
  size_t tablen=cltprm->numinitclumps+1;

  float *snarr;
  int32_t *indarr=NULL;
  double I, O, Ni, var, *row;
  int sky0_det1=cltprm->clprm->sky0_det1;
  size_t i, ind, counter=0, infodsize[2]={tablen, INFO_NCOLS};


  /* If there were no initial clumps, then ignore this function. */
  if(cltprm->numinitclumps==0) { cltprm->sn=NULL; return; }


  /* Allocate the arrays to keep the final S/N table (and possibly S/N
     index) for this object or tile. */
  cltprm->sn        = &cltprm->clprm->sn[ cltprm->id ];
  cltprm->sn->ndim  = 1;                        /* Depends on `cltprm->sn' */
  cltprm->sn->type  = GAL_TYPE_FLOAT32;
  cltprm->sn->dsize = gal_data_malloc_array(GAL_TYPE_SIZE_T, 1, __func__,
                                            "cltprm->sn->dsize");
  cltprm->sn->array = gal_data_malloc_array(cltprm->sn->type, tablen,
                                            __func__, "cltprm->sn->array");
  cltprm->sn->size  = cltprm->sn->dsize[0] = tablen;       /* After dsize. */
  if( cltprm->clprm->snind )
    {
      cltprm->snind        = &cltprm->clprm->snind [ cltprm->id ];
      cltprm->snind->ndim  = 1;              /* Depends on `cltprm->snind' */
      cltprm->snind->type  = GAL_TYPE_INT32;
      cltprm->snind->dsize = gal_data_malloc_array(GAL_TYPE_SIZE_T, 1,
                                                   __func__,
                                                   "cltprm->snind->dsize");
      cltprm->snind->size  = cltprm->snind->dsize[0]=tablen;/* After dsize */
      cltprm->snind->array = gal_data_malloc_array(cltprm->snind->type,
                                                   tablen, __func__,
                                                   "cltprm->snind->array");
    }
  else cltprm->snind=NULL;


  /* Allocate the array to keep the raw information of each clump. Note the
     `+1' in `infodsize', this is because the labels begin with 1 and we
     want each label to have one row on the same label.*/
  cltprm->info=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, infodsize,
                              NULL, 1, p->cp.minmapsize, NULL, NULL, NULL);


  /* First get the raw information necessary for making the S/N table. */
  clumps_get_raw_info(cltprm);


  /* Calculate the signal to noise ratio for successful clumps */
  snarr=cltprm->sn->array;
  if(cltprm->snind) indarr=cltprm->snind->array;
  for(i=1;i<tablen;++i)
    {
      /* For readability. */
      row = &( ((double *)(cltprm->info->array))[ i * INFO_NCOLS ] );
      Ni  = row[ INFO_INAREA ];
      I   = row[ INFO_INFLUX ]  / row[ INFO_INAREA ];
      O   = row[ INFO_RIVFLUX ] / row[ INFO_RIVAREA ];


      /* If the inner flux is smaller than the outer flux (happens only in
         noise cases) or the area is smaller than the minimum area to
         calculate signal-to-noise, then set the S/N of this segment to
         zero. */
      if( I>O && Ni>p->segsnminarea )
        {
          /* Here we have done sky subtraction once. However, if the sky
             was already subtracted (informed by the user), then the
             varience should be multiplied by 2.  */
          var = ( (p->skysubtracted ? 2.0f : 1.0f)
                  * row[INFO_INSTD] * row[INFO_INSTD] );

          /* Calculate the Signal to noise ratio, if we are on the noise
             regions, we don't care about the IDs of the clumps anymore, so
             store the Signal to noise ratios contiguously (for easy
             sorting and etc). Note that counter will always be smaller and
             equal to i. */
          ind = sky0_det1 ? i : counter++;
          if(cltprm->snind) indarr[ind]=i;
          snarr[ind]=( sqrt(Ni/p->cpscorr) * (I-O)
                       / sqrt( (I>0?I:-1*I) + (O>0?O:-1*O) + var ) );
        }
      else
        {
          /* Only over detections, we should put a NaN when the S/N  */
          if(sky0_det1)
            {
              snarr[i]=NAN;
              if(cltprm->snind) indarr[i]=i;
            }
        }
    }


  /* If we are in Sky mode, the sizes have to be corrected */
  if(sky0_det1==0)
    {
      cltprm->sn->dsize[0] = cltprm->sn->size = counter;
      if(cltprm->snind) cltprm->snind->dsize[0] = cltprm->snind->size=counter;
    }


  /* Clean up. */
  gal_data_free(cltprm->info);
}





/* Correct the labels of the clumps that will be used in determining the
   S/N threshold for true clumps.   */
static void
clumps_correct_sky_labels_for_check(struct clumps_thread_params *cltprm,
                                    gal_data_t *tile)
{
  gal_data_t *newinds;
  int32_t *ninds, curlab, *l, *lf;
  size_t len=cltprm->numinitclumps+1;
  struct noisechiselparams *p=cltprm->clprm->p;

  /* A small sanity check. */
  if(gal_tile_block(tile)!=p->clabel)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to address "
          "the problem. `tile->block' must point to the `clabel' dataset",
          __func__, PACKAGE_BUGREPORT);


  /* Allocate a dataset with the new indexs, note that it will need to have
     one element for each initial label (the excluded clumps need to be set
     to zero). So we also need to clear the allocated space. */
  newinds=gal_data_alloc(NULL, p->clabel->type, 1, &len, NULL, 0,
                         p->cp.minmapsize, NULL, NULL, NULL);


  /* Get the next available label for these clumps. If more than one thread
     was used, we are first going to lock the mutex (so no other thread
     changes these values), we will then read the shared number for this
     thread to use, then update the shared number and finally, unlock the
     mutex so other threads can do the same when they get to this point. */
  if(p->cp.numthreads>1) pthread_mutex_lock(&cltprm->clprm->labmutex);
  curlab        = p->numclumps+1;   /* Note that counting begins from 1. */
  p->numclumps += cltprm->snind->size;
  if(p->cp.numthreads>1) pthread_mutex_unlock(&cltprm->clprm->labmutex);


  /* Initialize the newinds array to CLUMPS_INIT (which be used as a new
     label for all the clumps that must be removed. */
  lf = (l=newinds->array) + newinds->size; do *l++=CLUMPS_INIT; while(l<lf);


  /* The new indexs array has been initialized to zero. So we just need to
     go over the labels in `cltprm->sninds' and give them a value of
     `curlab++'. */
  ninds=newinds->array;
  lf = (l=cltprm->snind->array) + cltprm->snind->size;
  do { ninds[*l]=curlab++; *l=ninds[*l]; } while(++l<lf);



  /* Go over this tile and correct the values. */
  GAL_TILE_PARSE_OPERATE( tile, NULL, 0, 1,
                          {if(*i>0) *i=ninds[ *(int32_t *)i ];} );

  /* Clean up. */
  gal_data_free(newinds);
}





static void *
clumps_find_make_sn_table(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct clumps_params *clprm=(struct clumps_params *)(tprm->params);
  struct noisechiselparams *p=clprm->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  void *tarray;
  double numdet;
  gal_data_t *tile, *tblock, *tmp;
  uint8_t *binary=p->binary->array;
  struct clumps_thread_params cltprm;
  size_t i, c, ind, tind, num, numsky, *indarr;
  size_t *scoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                       "scoord");
  size_t *icoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                       "icoord");


  /* Initialize the parameters for this thread. */
  cltprm.clprm   = clprm;
  cltprm.topinds = NULL;


  /* Go over all the tiles/detections given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* IDs. */
      cltprm.id = tind  = tprm->indexs[i];
      tile = &p->ltl.tiles[tind];


      /* Change the tile's pointers to the binary image (which has 1 for
         detected pixels and 0 for un-detected regions). */
      tarray=tile->array;
      tblock=tile->block;
      tile->array = gal_tile_block_relative_to_other(tile, p->binary);
      tile->block = p->binary;


      /* Get the number of usable elements in this tile (note that tiles
         can have blank pixels), so we can't simply use `tile->size'. */
      if(tile->flag & GAL_DATA_FLAG_HASBLANK)
        {
          tmp=gal_statistics_number(tile);
          num=*((size_t *)(tmp->array));
          gal_data_free(tmp);
        }
      else num=tile->size;


      /* Find the number of detected pixels over this tile. Since this is
         the binary image, this is just the sum of all the pixels. */
      tmp=gal_statistics_sum(tile);
      numdet=*((double *)(tmp->array));
      gal_data_free(tmp);


      /* See if this tile should be used or not (has enough undetected
         pixels). Note that it might happen that some tiles are fully
         blank. In such cases, it is important to first check the number of
         detected pixels. */
      numsky=num-numdet;
      if( num && (float)numsky/(float)num > p->minskyfrac )
        {
          /* Add the indexs of all undetected pixels in this tile into an
             array. */
          cltprm.indexs=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &numsky,
                                       NULL, 0, p->cp.minmapsize, NULL, NULL,
                                       NULL);

          /* Change the tile's block to the clump labels dataset (because
             we'll need to set the labels of the rivers on the edge of the
             tile here). */
          tile->array = gal_tile_block_relative_to_other(tile, p->clabel);
          tile->block = p->clabel;


          /* We need to set all the pixels on the edge of the tile to
             rivers and not include them in the list of indexs to set
             clumps. To do that, we need this tile's starting
             coordinates. */
          gal_dimension_index_to_coord(gal_data_ptr_dist(p->clabel->array,
                                                         tile->array,
                                                         p->clabel->type),
                                       ndim, dsize, scoord);


          /* Add the index of every sky element to the array of
             indexs. Note that since we know the array is always of type
             `int32_t', we can call the `GAL_TILE_PO_OISET' macro to avoid
             having to deal with multiple possible types in
             `GAL_TILE_PARSE_OPERATE'. Since the OUT macro-variable is
             NULL, the `int' is just a place-holder, it will not be
             used. */
          c=0;
          indarr=cltprm.indexs->array;
          GAL_TILE_PO_OISET(int32_t, int, tile, NULL, 0, 1, {
              /* This pixel's index over all the image. */
              ind = (int32_t *)i - (int32_t *)(p->clabel->array);
              gal_dimension_index_to_coord(ind, ndim, dsize, icoord);

              /* If the pixel is on the tile edge, set it as river and
                 don't include it in the indexs. */
              if( icoord[0]==scoord[0]
                  || icoord[0]==scoord[0]+tile->dsize[0]-1
                  || icoord[1]==scoord[1]
                  || icoord[1]==scoord[1]+tile->dsize[1]-1 )
                *(int32_t *)i=CLUMPS_RIVER;

              /* This pixel is not on the edge, check if it had a value of
                 `0' in the binary image (is not detected) then add it to
                 the list of indexs (note that the binary image also
                 contains the blank pixels, so only sky regions have a
                 value of 0 in the binary image). */
              else if( binary[ind]==0 )
                {
                  /*
                  if(c!=cltprm.indexs->size)
                    {
                      if(cltprm.id==282) *i+=2;
                  */
                      indarr[c++]=gal_data_ptr_dist(p->clabel->array, i,
                                                    p->clabel->type);
                  /*
                    }
                  else
                    if(cltprm.id==282)
                      {
                        int32_t *clabel=p->clabel->array;
                        size_t kjd=gal_data_ptr_dist(p->clabel->array, i,
                                                     p->clabel->type);
                        printf("%zu, %zu: %u\n", kjd%dsize[1]+1,
                               kjd/dsize[1]+1, clabel[kjd]);
                      }
                  */
                }
            });


          /* Correct the number of indexs. */
          cltprm.indexs->size=cltprm.indexs->dsize[0]=c;


          /* Generate the clumps over this region. */
          clumps_oversegment(&cltprm);


          /* Set all river pixels to CLUMPS_INIT (to be distinguishable
             from the detected regions). */
          GAL_TILE_PO_OISET( int32_t, int, tile, NULL, 0, 1,
                             {if(*i==CLUMPS_RIVER) *i=CLUMPS_INIT;} );


          /* For a check, the step variable will be set. */
          if(clprm->step==1)
            { gal_data_free(cltprm.indexs); continue; }


          /* Make the clump S/N table. */
          clumps_make_sn_table(&cltprm);


          /* If the user wanted to check the steps, remove the clumps that
             weren't used from the `clabel' image (they have been already
             excluded from the table). */
          if(cltprm.snind)
            clumps_correct_sky_labels_for_check(&cltprm, tile);


          /* Clean up. */
          gal_data_free(cltprm.indexs);
        }

      /* Reset the tile's pointers back to what they were. */
      tile->array=tarray;
      tile->block=tblock;
    }

  /* Clean up. */
  free(scoord);
  free(icoord);

  /* Wait for the all the threads to finish and return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Find the true clump signal to noise value from the clumps in the sky
   region.

   Each thread will find the useful signal to noise values for the tiles
   that have been assigned to it. It will then store the pointer to the S/N
   table into the sntablearr array (with the size of the number of
   meshs). If no clumps could be found in a mesh, then
   sntablearr[i]=NULL. Otherwise, it points to an array of the useful S/N
   values in that clump. Note that we don't care about the order of S/N
   values any more! There is also an accompanying array to keep the number
   of elements in the final S/N array of each mesh: numclumpsarr.

   Using these two arrays, after all the threads are finished, we can
   concatenate all the S/N values into one array and send it to the main
   findsnthresh function in thresh.c. */
void
clumps_true_find_sn_thresh(struct noisechiselparams *p)
{
  char *msg;
  struct timeval t1;
  size_t i, j, c, numsn=0;
  struct clumps_params clprm;
  gal_list_str_t *comments=NULL;
  gal_data_t *sn, *snind, *quant, *claborig;

  /* Get starting time for later reporting if necessary. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);


  /* Initialize/allocate the clump parameters structure, Note that the S/N
     indexs are also needed when we want to check the segmentation steps
     (they are used to correct the indexs in the final output). */
  clprm.p=p;
  clprm.sky0_det1=0;
  clprm.sn=gal_data_array_calloc(p->ltl.tottiles);
  clprm.snind = ( p->checksegmentation || p->checkclumpsn
                  ? gal_data_array_calloc(p->ltl.tottiles) : NULL );


  /* If the user wants to check the steps of get an S/N table, then we need
     a unique label for each clump. But in each region, the labels start
     from 1. So we need a central place to keep the next available
     label. Since `p->numclumps' is not used yet, we will use it here. When
     multiple threads are used, we will need a mutex to make sure that only
     one thread can change this central variable at every one moment. */
  if(p->checksegmentation || p->checkclumpsn)
    {
      p->numclumps=0;
      if( p->cp.numthreads > 1 ) pthread_mutex_init(&clprm.labmutex, NULL);
    }


  /* Spin off the threads to start the work. Note that several steps are
     done on each tile within a thread. So if the user wants to check
     steps, we need to break out of the processing get an over-all output,
     then reset the input and call it again. So it will be slower, but its
     is natural, since the user is testing to find the correct combination
     of parameters for later use. */
  if(p->segmentationname)
    {
      /* Necessary initializations. */
      clprm.step=1;
      claborig=p->clabel;
      p->clabel=gal_data_copy(claborig);

      /* Do each step. */
      while(clprm.step<3)
        {
          /* Reset the temporary copy of clabel back to its original. */
          if(clprm.step>1)
            memcpy(p->clabel->array, claborig->array,
                   claborig->size*gal_type_sizeof(claborig->type));

          /* Do this step. */
          gal_threads_spin_off(clumps_find_make_sn_table, &clprm,
                               p->ltl.tottiles, p->cp.numthreads);

          /* Set the extension name. */
          switch(clprm.step)
            {
            case 1: p->clabel->name = "SKY_CLUMPS_ALL";  break;
            case 2: p->clabel->name = "SKY_CLUMPS_FOR_SN";   break;
            default:
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so "
                    "we can address the issue. The value %d is not valid for "
                    "clprm.step", __func__, PACKAGE_BUGREPORT, clprm.step);
            }

          /* Write the demonstration array into the check image. The
             default values are hard to view, so we'll make a copy of the
             demo, set all Sky regions to blank and all clump macro values
             to zero. */
          gal_fits_img_write(p->clabel, p->segmentationname, NULL,
                             PROGRAM_NAME);

          /* Increment the step counter. */
          ++clprm.step;
        }

      /* Clean up (we don't need the original any more). */
      gal_data_free(claborig);
      p->clabel->name=NULL;
    }
  else
    {
      clprm.step=0;
      gal_threads_spin_off(clumps_find_make_sn_table, &clprm,
                           p->ltl.tottiles, p->cp.numthreads);
    }


  /* Destroy the mutex if it was initialized. */
  if( p->cp.numthreads>1 && (p->checksegmentation || p->checkclumpsn) )
    pthread_mutex_destroy(&clprm.labmutex);


  /* Find the total number of S/N values we have for all the clumps. */
  for(i=0;i<p->ltl.tottiles;++i)
    if(clprm.sn[i].ndim)  /* Only on tiles were an S/N was calculated. */
      numsn+=clprm.sn[i].size;
  if( numsn < p->minnumfalse )
    error(EXIT_FAILURE, 0, "only %zu clumps could be identified in the "
          "undetected regions. This is less than %zu (value to "
          "`--minnumfalse' option). Please either decrease this value or "
          "other options to change prior processing steps", numsn,
          p->minnumfalse);


  /* Allocate the space to keep all the S/N values. */
  sn=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &numsn, NULL, 0,
                    p->cp.minmapsize, "CLUMP_S/N", "ratio",
                    "Signal-to-noise ratio");
  snind = ( p->checkclumpsn
            ? gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &numsn, NULL, 0,
                             p->cp.minmapsize, "CLUMP_ID", "counter",
                             "Unique ID for this clump.")
            : NULL );


  /* Copy the S/N values of all the clumps into the unified array. */
  c=0;
  for(i=0;i<p->ltl.tottiles;++i)
    if(clprm.sn[i].ndim)
      for(j=0;j<clprm.sn[i].size;++j)
        {
          ((float *)(sn->array))[c] = ((float *)(clprm.sn[i].array))[j];
          if(snind)
            ((int32_t *)(snind->array))[c] =
              ((int32_t *)(clprm.snind[i].array))[j];
          ++c;
        }


  /* The S/N array of sky clumps is desiged to have no blank values, so set
     the flags accordingly to avoid a redundant blank search. */
  sn->flag = GAL_DATA_FLAG_BLANK_CH | GAL_DATA_FLAG_HASBLANK;


  /* If the user wanted to see the S/N table, then save it. */
  if(p->checkclumpsn)
    {
      /* Make the comments, then write the table and free the comments. */
      if(p->cp.numthreads>1)
        gal_list_str_add(&comments, "NOTE: In multi-threaded mode, clump "
                         "IDs differ in each run and are not sorted.", 1);
      gal_list_str_add(&comments, "See also: `SKY_CLUMPS_FOR_SN' HDU of "
                       "output with `--checksegmentation'.", 1);
      gal_list_str_add(&comments, "S/N of clumps over undetected regions.",
                       1);
      threshold_write_sn_table(p, sn, snind, p->clumpsn_s_name, comments);
      gal_list_str_free(comments, 1);
    }


  /* Find the desired quantile from the full S/N distribution. */
  quant = gal_statistics_quantile(sn, p->segquant, 1);
  p->clumpsnthresh = *((float *)(quant->array));
  if(!p->cp.quiet)
    {
      asprintf(&msg, "Clump S/N: %.2f (%.3f quant of %zu).",
               p->clumpsnthresh, p->segquant, sn->size);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* Clean up. */
  gal_data_free(sn);
  gal_data_free(snind);
  gal_data_free(quant);
  gal_data_array_free(clprm.sn, p->ltl.tottiles, 1);
  gal_data_array_free(clprm.snind, p->ltl.tottiles, 1);
}


















/***********************************************************************/
/*****************           Over detections           *****************/
/***********************************************************************/
/* Put the indexs of each labeled region into an array of `gal_data_t's
   (where each element is a dataset containing the respective label's
   indexs). */
gal_data_t *
clumps_det_label_indexs(struct noisechiselparams *p)
{
  size_t i, *areas;
  int32_t *a, *l, *lf;
  gal_data_t *labindexs=gal_data_array_calloc(p->numdetections+1);

  /* Find the area in each detected objects (to see how much space we need
     to allocate). */
  areas=gal_data_calloc_array(GAL_TYPE_SIZE_T, p->numdetections+1, __func__,
                              "areas");
  if(p->input->flag & GAL_DATA_FLAG_HASBLANK)
    {
      lf=(l=p->olabel->array)+p->olabel->size; /* Blank pixels have a      */
      do if(*l>0) ++areas[*l]; while(++l<lf);  /* negative value in int32. */
    }
  else
    {
      lf=(l=p->olabel->array)+p->olabel->size; do ++areas[*l]; while(++l<lf);
      areas[0]=0;
    }

  /* For a check.
  for(i=0;i<p->numdetections+1;++i)
    printf("detection %zu: %zu\n", i, areas[i]);
  exit(0);
  */

  /* Allocate/Initialize the dataset containing the indexs of each
     object. We don't want the labels of the non-detected regions
     (areas[0]). So we'll set that to zero.*/
  for(i=1;i<p->numdetections+1;++i)
    gal_data_initialize(&labindexs[i], NULL, GAL_TYPE_SIZE_T, 1,
                        &areas[i], NULL, 0, p->cp.minmapsize, NULL, NULL,
                        NULL);

  /* Put the indexs into each dataset. We will use the areas array again,
     but this time, use it as a counter. */
  memset(areas, 0, (p->numdetections+1)*sizeof *areas);
  lf=(a=l=p->olabel->array)+p->olabel->size;
  do
    if(*l>0)  /* No undetected regions (*l==0), or blank (<0) */
      ((size_t *)(labindexs[*l].array))[ areas[*l]++ ] = l-a;
  while(++l<lf);

  /* Clean up and return. */
  free(areas);
  return labindexs;
}





/* Only keep true clumps over detections. */
void
clumps_det_keep_true_relabel(struct clumps_thread_params *cltprm)
{
  struct noisechiselparams *p=cltprm->clprm->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  int istouching;
  size_t i, *s, *sf, *dinc;
  float *sn = cltprm->sn ? cltprm->sn->array : NULL;
  int32_t *l, *lf, *newlabs, curlab=1, *clabel=p->clabel->array;

  /* If there were no clumps over the detection, then just set the number
     of true clumps to zero, otherwise, see which ones should be
     removed. */
  if(cltprm->sn)
    {
      /* Allocate the necessary arrays. */
      newlabs=gal_data_malloc_array(GAL_TYPE_INT32,
                                    cltprm->numinitclumps+1, __func__,
                                    "newlabs");
      dinc=gal_dimension_increment(ndim, dsize);

      /* Initialize the new labels with CLUMPS_INIT (so the diffuse area
         can be distinguished from the clumps). */
      lf=(l=newlabs)+cltprm->numinitclumps+1;
      do *l++=CLUMPS_INIT; while(l<lf);

      /* Set the new labels. Here we will also be removing clumps with a peak
         that touches a river pixel. */
      if(p->keepmaxnearriver)
        {
          for(i=1;i<cltprm->numinitclumps+1;++i)
            if( sn[i] > p->clumpsnthresh ) newlabs[i]=curlab++;
        }
      else
        {
          for(i=1;i<cltprm->numinitclumps+1;++i)
            {
              /* Check if all the neighbors of this top element are
                 touching a river or not. */
              istouching=0;
              GAL_DIMENSION_NEIGHBOR_OP(cltprm->topinds[i], ndim, dsize,
                                        ndim, dinc,
                                        {
                                          if(clabel[nind]==0)
                                            istouching=1;
                                        });

              /* If the peak isn't touching a river, then check its S/N and
                 if that is also good, give it a new label. */
              if( !istouching && sn[i] > p->clumpsnthresh )
                newlabs[i]=curlab++;
            }
        }

      /* Correct the clump labels. Note that the non-clumpy regions over the
         detections (rivers) have already been initialized to CLUMPS_INIT
         (which is negative). So we'll just need to correct the ones with a
         value larger than 0. */
      sf=(s=cltprm->indexs->array)+cltprm->indexs->size;
      do if(clabel[*s]>0) clabel[*s] = newlabs[ clabel[*s] ]; while(++s<sf);

      /* Save the total number of true clumps in this detection. */
      cltprm->numtrueclumps=curlab-1;

      /* Clean up. */
      free(dinc);
      free(newlabs);
    }
  else cltprm->numtrueclumps=0;
}
