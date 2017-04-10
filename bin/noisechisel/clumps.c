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
#include <gnuastro/linkedlist.h>

#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "ui.h"










/**********************************************************************/
/*****************     Basic structures and macros    *****************/
/**********************************************************************/
/* Parameters for all threads. */
struct clumps_params
{
  int               sky0_det1;  /* If working on the Sky or Detections.    */
  gal_data_t              *sn;  /* Array of clump S/N tables.              */
  gal_data_t           *snind;  /* Array of clump S/N index (for check).   */
  struct noisechiselparams *p;  /* Pointer to main NoiseChisel parameters. */
};

/* Parameters for one thread. */
struct clumps_thread_params
{
  size_t                   id;  /* ID of this detection/tile over tile.    */
  size_t             *topinds;  /* Indexs of all local maxima.             */
  size_t            numclumps;  /* The number of clumps found in this run. */
  gal_data_t          *indexs;  /* Array containing indexs of this object. */
  gal_data_t            *info;  /* Information for all clumps.             */
  gal_data_t              *sn;  /* Signal-to-noise ratio for these clumps. */
  gal_data_t           *snind;  /* Index of S/N for these clumps.          */
  struct clumps_params *clprm;  /* Pointer to main structure.              */
};

/* Constants for the clump over-segmentation. */
#define CLUMPS_RIVER     UINT32_MAX-2
#define CLUMPS_TMPCHECK  UINT32_MAX-3
#define CLUMPS_INIT      UINT32_MAX-4
#define CLUMPS_MAXLAB    UINT32_MAX-5   /* Largest clump label (unsigned). */

















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
  size_t *a, *af, ind, *dsize=p->input->dsize;
  struct gal_linkedlist_sll *Q=NULL, *cleanup=NULL;
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  uint32_t n1, nlab, rlab, curlab=1, *clabel=p->clabel->array;

  /*********************************************
   For checks and debugging:*
  gal_data_t *crop;
  size_t extcount=1;
  uint32_t *cr, *crf;
  size_t checkdsize[2]={10,10};
  size_t checkstart[2]={50,145};
  char *filename="clumpbuild.fits";
  size_t checkstartind=gal_dimension_coord_to_index(2, dsize, checkstart);
  gal_data_t *tile=gal_data_alloc(gal_data_ptr_increment(arr, checkstartind,
                                                         p->conv->type),
                                  GAL_TYPE_INVALID, 2, checkdsize,
                                  NULL, 0, 0, NULL, NULL, NULL);
  tile->block=p->conv;
  gal_checkset_check_remove_file(filename, 0, 0);
  if(p->cp.numthreads!=1)
    error(EXIT_FAILURE, 0, "in the debugging mode of `clumps_oversegment' "
          "only one thread must be used");
  crop=gal_data_copy(tile);
  gal_fits_img_write(crop, filename, NULL, PROGRAM_STRING);
  gal_data_free(crop);
  printf("blank: %u\nriver: %u\ntmpcheck: %u\ninit: %u\nmaxlab: %u\n",
         (uint32_t)GAL_BLANK_UINT32, (uint32_t)CLUMPS_RIVER,
         (uint32_t)CLUMPS_TMPCHECK, (uint32_t)CLUMPS_INIT,
         (uint32_t)CLUMPS_MAXLAB);
  tile->array=gal_tile_block_relative_to_other(tile, p->clabel);
  tile->block=p->clabel;
  **********************************************/


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
              error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we "
                    "can fix this problem. In `clumps_oversegment', `Q' and "
                    "`cleanup' should be NULL but while checking the equal"
                    "flux regions they aren't", PACKAGE_BUGREPORT);

            /* Add this pixel to a queue. */
            gal_linkedlist_add_to_sll(&Q, *a);
            gal_linkedlist_add_to_sll(&cleanup, *a);
            clabel[*a] = CLUMPS_TMPCHECK;

            /* Find all the pixels that have the same flux and are
               connected. */
            while(Q!=NULL)
              {
                /* Pop an element from the queue. */
                gal_linkedlist_pop_from_sll(&Q, &ind);

                /* Look at the neighbors and see if we already have a
                   label. */
                GAL_DIMENSION_NEIGHBOR_OP(ind, ndim, dsize, ndim, dinc,
                   {
                     /* For easy reading. */
                     nlab=clabel[ nind ];

                     /* This neighbor is within the search region. */
                     if(nlab)
                       {
                         /* If this neighbour has not been labeled yet
                            and has an equal flux, add it to the queue
                            to expand the studied region.*/
                         if( nlab==CLUMPS_INIT && arr[nind]==arr[*a] )
                           {
                             clabel[nind]=CLUMPS_TMPCHECK;
                             gal_linkedlist_add_to_sll(&cleanup, nind);
                             gal_linkedlist_add_to_sll(&Q, nind);
                           }

                         /* If this neighbour has a positive nlab, it
                            belongs to another object, so if `n1' has not
                            been set for the whole region, put `nlab' into
                            `n1'. If `n1' has been set and is different
                            from `nlab' then this whole equal flux region
                            should be a wide river because it is connecting
                            two connected regions.*/
                         else if(nlab<CLUMPS_MAXLAB)
                           {
                             if(n1==0)            n1 = nlab;
                             else if(nlab!=n1)    n1 = CLUMPS_RIVER;
                           }
                       }

                     /* If this neigbour has a label of zero, then we are
                        on the edge of the indexed region (the neighbor is
                        not in the initial list of pixels to segment). When
                        over-segmenting the noise and the detections,
                        `clabel' is zero for the parts of the image that we
                        are not interested in here. */
                     else
                       clabel[*a]=CLUMPS_RIVER;
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
                if( cltprm->topinds)        /* This is a local maximum of  */
                  cltprm->topinds[rlab]=*a; /* this region, save its index.*/
              }

            /* Give the same label to the whole connected equal flux
               region, except those that might have been on the side of
               the image and were a river pixel. */
            while(cleanup!=NULL)
              {
                gal_linkedlist_pop_from_sll(&cleanup, &ind);
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
               see if all the neighbors (connectivity equals the number of
               dimensions) that have a non-macro value (less than
               CLUMPS_MAXLAB) belong to one label or not. If the pixel is
               neighboured by more than one label, set it as a river
               pixel. Also if it is touching a zero valued pixel (which
               does not belong to this object), set it as a river pixel.*/
            GAL_DIMENSION_NEIGHBOR_OP(*a, ndim, dsize, ndim, dinc,
               {
                 nlab=clabel[ nind ];
                 if(nlab<CLUMPS_MAXLAB)
                   {
                     if(nlab==0)           n1=CLUMPS_RIVER;
                     else
                       {
                         if(n1)
                           {
                             if(nlab!=n1)  n1=CLUMPS_RIVER;
                           }
                         else              n1=nlab;
                       }
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
            gal_fits_img_write(crop, filename, NULL, PROGRAM_STRING);
            gal_data_free(crop);
          }
        **********************************************/
      }
  while(++a<af);

  /* Save the total number of clumps. */
  cltprm->numclumps=curlab-1;


  /* Set all the river pixels to zero, this is only necessary for the
     clumps over detections, not the sky clumps (they will be converted
     over all the image. */
  if(cltprm->clprm->sky0_det1)
    {
      af=(a=indexs->array)+indexs->size;
      do if( clabel[*a]==CLUMPS_RIVER ) clabel[*a]=0; while(++a<af);
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
/*****************             S/N threshold          *****************/
/**********************************************************************/
/* In this function we want to find the general information for each clump
   in an over-segmented labeled array. The signal in each clump is the
   average signal inside it subtracted by the average signal in the river
   pixels around it. So this function will go over all the pixels in the
   object (already found in deblendclumps()) and add them appropriately.

   The output is an array of size numclumps*INFO_NCOLS. as listed
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
  uint32_t lab, nlab, *ngblabs, *clabel=p->clabel->array;

  /* Allocate the array to keep the neighbor labels of river pixels. */
  ngblabs=gal_data_malloc_array(GAL_TYPE_UINT32, nngb);

  /* Go over all the pixels in this region. */
  af=(a=cltprm->indexs->array)+cltprm->indexs->size;
  do
    if( !isnan(arr[ *a ]) )
      {
        /* This pixel belongs to a clump. */
        if( clabel[ *a ] )
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
                if(nlab && nlab<CLUMPS_MAXLAB)
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
  for(lab=1; lab<=cltprm->numclumps; ++lab)
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
static void
clumps_make_sn_table(struct clumps_thread_params *cltprm)
{
  struct noisechiselparams *p=cltprm->clprm->p;
  size_t tablen=cltprm->numclumps+1;

  float *snarr;
  uint32_t *indarr=NULL;
  double I, O, Ni, var, *row;
  int sky0_det1=cltprm->clprm->sky0_det1;
  size_t i, ind, counter=0, infodsize[2]={tablen, INFO_NCOLS};


  /* Allocate the arrays to keep the final S/N table (and possibly S/N
     index) for this object or tile. */
  cltprm->sn        = &cltprm->clprm->sn[ cltprm->id ];
  cltprm->sn->ndim  = 1;                       /* Depends on `cltprm->sn' */
  cltprm->sn->type  = GAL_TYPE_FLOAT32;
  cltprm->sn->dsize = gal_data_malloc_array(GAL_TYPE_SIZE_T, 1);
  cltprm->sn->array = gal_data_malloc_array(cltprm->sn->type, tablen);
  cltprm->sn->size  = cltprm->sn->dsize[0] = tablen;      /* After dsize. */
  if(p->checkclumpsn)
    {
      cltprm->snind        = &cltprm->clprm->snind [ cltprm->id ];
      cltprm->snind->ndim  = 1;             /* Depends on `cltprm->snind' */
      cltprm->snind->type  = GAL_TYPE_UINT32;
      cltprm->snind->dsize = gal_data_malloc_array(GAL_TYPE_SIZE_T, 1);
      cltprm->snind->size  = cltprm->snind->dsize[0]=tablen;/* After dsize */
      cltprm->snind->array = gal_data_malloc_array(cltprm->snind->type,
                                                   tablen);
    }
  else cltprm->snind=NULL;


  /* Allocate the array to keep the raw information of each clump. Note the
     `+1' in `infodsize', this is because the labels begin with 1 and we
     want each label to have one row on the same label.*/
  cltprm->info=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, infodsize,
                              NULL, 1, p->cp.minmapsize, NULL, NULL, NULL);

  /* First get the raw information necessary for making the S/N table. */
  clumps_get_raw_info(cltprm);

  /* Calculate the signal to noise for successful detections */
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
      if( Ni>p->segsnminarea && I>O )   /* This is O, not 0 (zero). */
        {
          /* Here we have done sky subtraction once. However, if the sky
             was already subtracted (informed by the user), then the
             varience should be multiplied by 2.  */
          var = ( (p->skysubtracted ? 2.0f : 1.0f)
                  * row[INFO_INSTD] * row[INFO_INSTD] );

          /* Calculate the Signal to noise ratio, if we are on the
             noise regions, we don't care about the IDs of the clumps
             anymore, so store the Signal to noise ratios contiguously
             (for easy sorting and etc). Note that counter will always
             be smaller and equal to i. */
          ind = sky0_det1 ? i : counter++;
          if(cltprm->snind) indarr[ind]=i;
          snarr[ind]=( sqrt(Ni/p->cpscorr)*(I-O)
                       / sqrt( (I>0?I:-1*I) + (O>0?O:-1*O) + var ) );
        }
      else
        {
          ind = sky0_det1 ? i : counter++;
          if(cltprm->snind) indarr[ind]=i;
          snarr[ind]=NAN;
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
  size_t i, c, ind, tind, numsky, *indarr;
  size_t *scoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);
  size_t *icoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim);


  /* Initialize the parameters for this thread. */
  cltprm.clprm   = clprm;
  cltprm.topinds = clprm->sky0_det1 ? NULL : NULL;


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
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


      /* Find the number of detected pixels over this tile. Since this is
         the binary image, this is just the sum of all the pixels. */
      tmp=gal_statistics_sum(tile);
      numdet=*((double *)(tmp->array));
      gal_data_free(tmp);


      /* See if this tile should be used or not (has enough undetected
         pixels). */
      numsky=tile->size-numdet;
      if( (float)numsky/(float)tile->size > p->minbfrac )
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

          /* Add the index of every sky element to the array of indexs. */
          c=0;
          indarr=cltprm.indexs->array;
          GAL_TILE_PARSE_OPERATE({
              /* This pixel's index over all the image. */
              ind = (uint32_t *)i - (uint32_t *)(p->clabel->array);
              gal_dimension_index_to_coord(ind, ndim, dsize, icoord);

              /* If the pixel is on the edge, set it as river and
                 don't include it in the indexs. */
              if( icoord[0]==scoord[0]
                  || icoord[0]==scoord[0]+tile->dsize[0]-1
                  || icoord[1]==scoord[1]
                  || icoord[1]==scoord[1]+tile->dsize[1]-1 )
                *(uint32_t *)i=CLUMPS_RIVER;

              /* This pixel is not on the edge, check if it had a value
                 of `0' in the binary image (is not detected), and if
                 so, then add it to the list of indexs. */
              else if(binary[ind]==0)
                indarr[c++]=gal_data_ptr_dist(tile->block->array, i,
                                              p->clabel->type);
            }, tile, NULL, 0, 1);

          /* Correct the number of indexs. */
          cltprm.indexs->size=cltprm.indexs->dsize[0]=c;

          /* Generate the clumps over this region. */
          clumps_oversegment(&cltprm);

          /* Correct the river pixels */
          GAL_TILE_PARSE_OPERATE({if(*i==CLUMPS_RIVER) *i=0;}, tile,
                                 NULL, 0, 1);

          /* Make the clump S/N table. */
          clumps_make_sn_table(&cltprm);

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



















/**********************************************************************/
/*****************         High level functins        *****************/
/**********************************************************************/
/* The job of this function is to find the best signal to noise value to
   use as a threshold to detect real clumps.

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
clumps_on_undetected_sn(struct noisechiselparams *p)
{
  struct clumps_params clprm;

  /* Initialize/allocate the clump parameters structure,  */
  clprm.p=p;
  clprm.sky0_det1=0;
  clprm.sn=gal_data_array_calloc(p->ltl.tottiles);
  clprm.snind = ( p->checkclumpsn
                  ? gal_data_array_calloc(p->ltl.tottiles) : NULL );


  /* Spin off the threads to start the work. */
  gal_threads_spin_off(clumps_find_make_sn_table, &clprm, p->ltl.tottiles,
                       p->cp.numthreads);


  /* If the user wanted to see the steps. */
  if(p->segmentationname)
    gal_fits_img_write(p->clabel, p->segmentationname, NULL, PROGRAM_STRING);

  /* Clean up. */
  gal_data_array_free(clprm.sn, p->ltl.tottiles, 1);
  if(p->checkclumpsn) gal_data_array_free(clprm.snind, p->ltl.tottiles, 1);
}
