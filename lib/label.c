/*********************************************************************
label -- Work on labeled (integer valued) datasets.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/list.h>
#include <gnuastro/qsort.h>
#include <gnuastro/label.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>









/****************************************************************
 *****************         Internal          ********************
 ****************************************************************/
static void
label_check_type(gal_data_t *in, uint8_t needed_type, char *variable,
                 const char *func)
{
  if(in->type!=needed_type)
    error(EXIT_FAILURE, 0, "%s: the '%s' dataset has '%s' type, but it "
          "must have a '%s' type.\n\n"
          "You can use 'gal_data_copy_to_new_type' or "
          "'gal_data_copy_to_new_type_free' to convert your input dataset "
          "to this type before calling this function", func, variable,
          gal_type_name(in->type, 1), gal_type_name(needed_type, 1));
}




















/****************************************************************
 *****************          Indexs           ********************
 ****************************************************************/
/* Put the indexs of each labeled region into an array of 'gal_data_t's
   (where each element is a dataset containing the respective label's
   indexs). */
gal_data_t *
gal_label_indexs(gal_data_t *labels, size_t numlabs, size_t minmapsize,
                 int quietmmap)
{
  size_t i, *areas;
  int32_t *a, *l, *lf;
  gal_data_t *max, *labindexs;

  /* Sanity check. */
  label_check_type(labels, GAL_TYPE_INT32, "labels", __func__);

  /* If the user hasn't given the number of labels, find it (maximum
     label). */
  if(numlabs==0)
    {
      max=gal_statistics_maximum(labels);
      numlabs=*((int32_t *)(max->array));
      gal_data_free(max);
    }
  labindexs=gal_data_array_calloc(numlabs+1);

  /* Find the area in each detected object (to see how much space we need
     to allocate). If blank values are present, an extra check is
     necessary, so to get faster results when there aren't any blank
     values, we'll also do a check. */
  areas=gal_pointer_allocate(GAL_TYPE_SIZE_T, numlabs+1, 1, __func__,
                             "areas");
  lf=(l=labels->array)+labels->size;
  do
    if(*l>0)  /* Only labeled regions: *l==0 (undetected), *l<0 (blank). */
      ++areas[*l];
  while(++l<lf);

  /* For a check.
  for(i=0;i<numlabs+1;++i)
    printf("detection %zu: %zu\n", i, areas[i]);
  exit(0);
  */

  /* Allocate/Initialize the dataset containing the indexs of each
     object. We don't want the labels of the non-detected regions
     (areas[0]). So we'll set that to zero.*/
  for(i=1;i<numlabs+1;++i)
    gal_data_initialize(&labindexs[i], NULL, GAL_TYPE_SIZE_T, 1,
                        &areas[i], NULL, 0, minmapsize, quietmmap,
                        NULL, NULL, NULL);

  /* Put the indexs into each dataset. We will use the areas array again,
     but this time, use it as a counter. */
  memset(areas, 0, (numlabs+1)*sizeof *areas);
  lf=(a=l=labels->array)+labels->size;
  do
    if(*l>0)  /* No undetected regions (*l==0), or blank (<0) */
      ((size_t *)(labindexs[*l].array))[ areas[*l]++ ] = l-a;
  while(++l<lf);

  /* Clean up and return. */
  free(areas);
  return labindexs;
}




















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
   neighboring labeled region that pixel will be a 'river' pixel.

   DON'T FORGET: SET THE FLAGS FOR CONV EQUAL TO INPUT IN SEGMENT.

*/
size_t
gal_label_watershed(gal_data_t *values, gal_data_t *indexs,
                    gal_data_t *labels, size_t *topinds, int min0_max1)
{
  size_t ndim=values->ndim;

  int hasblank;
  float *arr=values->array;
  gal_list_sizet_t *Q=NULL, *cleanup=NULL;
  size_t *a, *af, ind, *dsize=values->dsize;
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  int32_t n1, nlab, rlab, curlab=1, *labs=labels->array;

  /* Sanity checks */
  label_check_type(values, GAL_TYPE_FLOAT32, "values", __func__);
  label_check_type(indexs, GAL_TYPE_SIZE_T,  "indexs", __func__);
  label_check_type(labels, GAL_TYPE_INT32,   "labels", __func__);
  if( gal_dimension_is_different(values, labels) )
    error(EXIT_FAILURE, 0, "%s: the 'values' and 'labels' arguments must "
          "have the same size", __func__);
  if(indexs->ndim!=1)
    error(EXIT_FAILURE, 0, "%s: 'indexs' has to be a 1D array, but it is "
          "%zuD", __func__, indexs->ndim);


  /* See if there are blank values in the input dataset. */
  hasblank=gal_blank_present(values, 0);


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
                                                         values->type),
                                  GAL_TYPE_INVALID, 2, checkdsize,
                                  NULL, 0, 0, NULL, NULL, NULL);
  tile->block=values;
  gal_checkset_writable_remove(filename, 0, 0);
  crop=gal_data_copy(tile);
  gal_fits_img_write(crop, filename, NULL, PROGRAM_NAME);
  gal_data_free(crop);
  printf("blank: %u\nriver: %u\ntmpcheck: %u\ninit: %u\n",
         (int32_t)GAL_BLANK_INT32, (int32_t)GAL_LABEL_RIVER,
         (int32_t)GAL_LABEL_TMPCHECK, (int32_t)GAL_LABEL_INIT);
  tile->array=gal_tile_block_relative_to_other(tile, labels);
  tile->block=labels;
  **********************************************/


  /* If the size of the indexs is zero, then this function is pointless. */
  if(indexs->size==0) return 0;


  /* If the indexs aren't already sorted (by the value they correspond to),
     sort them given indexs based on their flux ('gal_qsort_index_arr' is
     defined as static in 'gnuastro/qsort.h') */
  if( !( (indexs->flag & GAL_DATA_FLAG_SORT_CH)
        && ( indexs->flag
             & (GAL_DATA_FLAG_SORTED_I
                | GAL_DATA_FLAG_SORTED_D) ) ) )
    {
      gal_qsort_index_single=values->array;
      qsort(indexs->array, indexs->size, sizeof(size_t),
            ( min0_max1
              ? gal_qsort_index_single_float32_d
              : gal_qsort_index_single_float32_i) );
    }


  /* Initialize the region we want to over-segment. */
  af=(a=indexs->array)+indexs->size;
  do labs[*a]=GAL_LABEL_INIT; while(++a<af);


  /* Go over all the given indexs and pull out the clumps. */
  af=(a=indexs->array)+indexs->size;
  do
    /* When regions of a constant flux or masked regions exist, some later
       indexs (although they have same flux) will be filled before hand. If
       they are done, there is no need to do them again. */
    if(labs[*a]==GAL_LABEL_INIT)
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
                    "we can fix this problem. 'Q' and 'cleanup' should be "
                    "NULL but while checking the equal flux regions they "
                    "aren't", __func__, PACKAGE_BUGREPORT);

            /* Add this pixel to a queue. */
            gal_list_sizet_add(&Q, *a);
            gal_list_sizet_add(&cleanup, *a);
            labs[*a] = GAL_LABEL_TMPCHECK;

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
                         nlab=labs[ nind ];

                         /* This neighbor's label isn't zero. */
                         if(nlab)
                           {
                             /* If this neighbor has not been labeled yet
                                and has an equal flux, add it to the queue
                                to expand the studied region.*/
                             if( nlab==GAL_LABEL_INIT && arr[nind]==arr[*a] )
                               {
                                 labs[nind]=GAL_LABEL_TMPCHECK;
                                 gal_list_sizet_add(&Q, nind);
                                 gal_list_sizet_add(&cleanup, nind);
                               }
                             else
                               n1=( nlab>0

                                    /* If this neighbor has a positive
                                       nlab, it belongs to another object,
                                       so if 'n1' has not been set for the
                                       whole region (n1==0), put 'nlab'
                                       into 'n1'. If 'n1' has been set and
                                       is different from 'nlab' then this
                                       whole equal flux region should be a
                                       wide river because it is connecting
                                       two connected regions.*/
                                    ? ( n1
                                        ? (n1==nlab ? n1 : GAL_LABEL_RIVER)
                                        : nlab )

                                    /* If the data has blank pixels, see if
                                       the neighbor is blank. If so, set
                                       the label to a river. Checking for
                                       the presence of blank values in the
                                       dataset can be done outside this
                                       loop (or even outside this function
                                       if flags are set). So to help the
                                       compiler optimize the program, we'll
                                       first use the pre-checked value. */
                                    : ( ( hasblank && isnan(arr[nind]) )
                                        ? GAL_LABEL_RIVER
                                        : n1 ) );
                           }

                         /* If this neigbour has a label of zero, then we
                            are on the edge of the indexed region (the
                            neighbor is not in the initial list of pixels
                            to segment). When over-segmenting the noise and
                            the detections, 'label' is zero for the parts
                            of the image that we are not interested in
                            here. */
                         else labs[*a]=GAL_LABEL_RIVER;
                       }
                   } );
              }

            /* Set the label that is to be given to this equal flux
               region. If 'n1' was set to any value, then that label should
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
                if( labs[ ind ]==GAL_LABEL_TMPCHECK ) labs[ ind ]=rlab;
              }
          }

        /* The flux of this pixel is not the same as the next sorted
           flux, so simply find the label for this object. */
        else
          {
            /* 'n1' is the label of the first labeled neighbor found, so
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
                 /* When 'n1' has already been set as a river, there is no
                    point in looking at the other neighbors. */
                 if(n1!=GAL_LABEL_RIVER)
                   {
                     /* For easy reading. */
                     nlab=labs[ nind ];

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

                                /* If the data has blank pixels, see if the
                                   neighbor is blank. If so, set the label
                                   to a river. Checking for the presence of
                                   blank values in the dataset can be done
                                   outside this loop (or even outside this
                                   function if flags are set). So to help
                                   the compiler optimize the program, we'll
                                   first use the pre-checked value. */
                                : ( ( hasblank && isnan(arr[nind]) )
                                    ? GAL_LABEL_RIVER
                                    : n1 ) )

                            /* 'nlab==0' (the neighbor lies in the other
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
               and has been set to 'GAL_LABEL_RIVER' before) or all its
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
            labs[ *a ] = rlab;
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
                   labs[*a]);
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




















/**********************************************************************/
/*************             Clump significance             *************/
/**********************************************************************/
static int
label_clump_significance_sanity(gal_data_t *values, gal_data_t *std,
                                gal_data_t *label, gal_data_t *indexs,
                                struct gal_tile_two_layer_params *tl,
                                gal_data_t *sig, const char *func)
{
  size_t *a, *af;
  float first=NAN, second=NAN, *f=values->array;

  /* Type of values. */
  if( values->type!=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "%s: the values dataset must have a 'float' "
          "type, but it has a '%s' type", func,
          gal_type_name(values->type, 1));

  /* Type of standard deviation. */
  if( std->type!=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "%s: the standard deviation dataset must have a "
          "'float' ('float32') type, but it has a '%s' type", func,
          gal_type_name(std->type, 1));

  /* Type of labels image. */
  if( label->type!=GAL_TYPE_INT32 )
    error(EXIT_FAILURE, 0, "%s: the labels dataset must have an 'int32' "
          "type, but it has a '%s' type", func,
          gal_type_name(label->type, 1));

  /* Dimentionality of the values dataset. */
  if( values->ndim>3 )
    error(EXIT_FAILURE, 0, "%s: currently only supports 1, 2 or 3 "
          "dimensional datasets, but a %zu-dimensional dataset is given",
          func, values->ndim);

  /* Type of indexs image. */
  if( indexs->type!=GAL_TYPE_SIZE_T )
    error(EXIT_FAILURE, 0, "%s: the indexs dataset must have a 'size_t' "
          "type, but it has a '%s' type", func,
          gal_type_name(label->type, 1));

  /* Dimensionality of indexs (must be 1D). */
  if( indexs->ndim!=1 )
    error(EXIT_FAILURE, 0, "%s: the indexs dataset must be a 1D dataset, "
          "but it has %zu dimensions", func, indexs->ndim);

  /* Similar sizes between values and standard deviation. */
  if( gal_dimension_is_different(values, label) )
    error(EXIT_FAILURE, 0, "%s: the values and label arrays don't have the "
          "same size.", func);

  /* Size of the standard deviation. */
  if( !( std->size==1
         || std->size==values->size
         || (tl && std->size==tl->tottiles) ) )
    error(EXIT_FAILURE, 0, "%s: the standard deviation dataset has %zu "
          "elements. But it can only have one of these sizes: 1) a "
          "single value (used for the whole dataset), 2) The size of "
          "the values dataset (%zu elements, one value for each "
          "element), 3) The size of the number of tiles in the input "
          "tessellation (when a tessellation is given)",
          func, std->size, values->size);

  /* If the 'array' and 'dsize' elements of 'sig' have already been set. */
  if(sig->array)
    error(EXIT_FAILURE, 0, "%s: the dataset that will contain the "
          "significance values must have NULL pointers for its 'array' "
          "and 'dsize' pointers (they will be allocated here)", func);

  /* See if the clumps are to be built starting from local maxima or local
     minima. */
  af=(a=indexs->array)+indexs->size;
  do
    /* A label may have NAN values. */
    if( !isnan(f[*a]) )
      {
        if( isnan(first) )
          first=f[*a];
        else
          {
            if( isnan(second) )
              {
                /* Note that the elements may have equal values, so for
                   'second', we want the first non-blank AND different
                   value. */
                if( f[*a]!=first )
                  second=f[*a];
              }
            else
              break;
          }
      }
  while(++a<af);

  /* Note that if all the values are blank or there is only one value
     covered by all the indexs, then both (or one) of 'first' or 'second'
     will be NAN. In either case, the significance measure is not going to
     be meaningful if we assume the clumps start from the maxima or
     minima. So we won't check if they are NaN or not.*/
  return first>second ? 1 : 0;
}





/* In this function we want to find the general information for each clump
   in an over-segmented labeled array. The signal in each clump is the
   average signal inside it subtracted by the average signal in the river
   pixels around it. So this function will go over all the pixels in the
   object (already found in deblendclumps()) and add them appropriately.

   The output is an array of size cltprm->numinitial*INFO_NCOLS. as listed
   below.*/
enum infocols
  {
    INFO_STD,            /* Standard deviation.                           */
    INFO_INAREA,         /* Tatal area within clump.                      */
    INFO_RIVAREA,        /* Tatal area within rivers around clump.        */
    INFO_PEAK_RIVER,     /* Peak (min or max) river value around a clump. */
    INFO_PEAK_CENTER,    /* Peak (min or max) clump value.                */

    INFO_NCOLS,          /* Total number of columns in the 'info' table.  */
  };
static void
label_clump_significance_raw(gal_data_t *values_d, gal_data_t *std_d,
                             gal_data_t *label_d, gal_data_t *indexs,
                             struct gal_tile_two_layer_params *tl,
                             double *info)
{
  size_t ndim=values_d->ndim, *dsize=values_d->dsize;

  double *row;
  size_t i, *a, *af, ii, coord[3];
  size_t nngb=gal_dimension_num_neighbors(ndim);
  int32_t nlab, *ngblabs, *label=label_d->array;
  float *values=values_d->array, *std=std_d->array;
  size_t *dinc=gal_dimension_increment(ndim, dsize);

  /* Allocate the array to keep the neighbor labels of river pixels. */
  ngblabs=gal_pointer_allocate(GAL_TYPE_INT32, nngb, 0, __func__, "ngblabs");

  /* Go over all the pixels in this region. */
  af=(a=indexs->array)+indexs->size;
  do
    if( !isnan(values[ *a ]) )
      {
        /* This pixel belongs to a clump. */
        if( label[ *a ]>0 )
          {
            /* For easy reading. */
            row = &info [ label[*a] * INFO_NCOLS ];

            /* Add this pixel to this clump's area. */
            ++row[ INFO_INAREA ];

            /* In the loop 'INFO_INAREA' is just the pixel counter of this
               clump. The pixels are sorted by flux (decreasing for
               positive clumps and increasing for negative). So the second
               extremum value is just the second pixel of the clump. */
            if( row[ INFO_INAREA ]==1.0f )
              row[ INFO_PEAK_CENTER ] = values[*a];
          }

        /* This pixel belongs to a river (has a value of zero and isn't
           blank). */
        else
          {
            /* We are on a river pixel. So the value of this pixel has to
               be added to any of the clumps in touches. But since it might
               touch a labeled region more than once, we use 'ngblabs' to
               keep track of which label we have already added its value
               to. 'ii' is the number of different labels this river pixel
               has already been considered for. 'ngblabs' will keep the list
               labels. */
            ii=0;
            memset(ngblabs, 0, nngb*sizeof *ngblabs);

            /* Look into the 8-connected neighbors (recall that a
               connectivity of 'ndim' means all pixels touching it (even on
               one vertice). */
            GAL_DIMENSION_NEIGHBOR_OP(*a, ndim, dsize, ndim, dinc, {
                /* This neighbor's label. */
                nlab=label[ nind ];

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
                        row = &info[ nlab * INFO_NCOLS ];

                        ++row[INFO_RIVAREA];
                        if( row[INFO_RIVAREA]==1.0f )
                          {
                            /* Get the maximum river value. */
                            row[INFO_PEAK_RIVER] = values[*a];

                            /* Get the standard deviation. */
                            if(std_d->size==1 || std_d->size==values_d->size)
                              row[INFO_STD]=std_d->size==1?std[0]:std[*a];
                            else
                              {
                                gal_dimension_index_to_coord(*a, ndim, dsize,
                                                             coord);
                                row[INFO_STD]=
                                  std[gal_tile_full_id_from_coord(tl,coord)];
                              }
                          }
                      }
                  }
              } );
          }
      }
  while(++a<af);

  /* Clean up. */
  free(dinc);
  free(ngblabs);
}




/* Make an S/N table for the clumps in a given region. */
void
gal_label_clump_significance(gal_data_t *values, gal_data_t *std,
                             gal_data_t *label, gal_data_t *indexs,
                             struct gal_tile_two_layer_params *tl,
                             size_t numclumps, size_t minarea, int variance,
                             int keepsmall, gal_data_t *sig,
                             gal_data_t *sigind)
{
  double *info;
  int max1_min0;
  float *sigarr;
  double C, R, S, *row;
  int32_t *indarr=NULL;
  size_t i, ind, counter=0;
  size_t tablen=numclumps+1;

  /* If there were no initial clumps, then ignore this function. */
  if(numclumps==0) { sig->size=0; return; }

  /* Basic sanity checks. */
  max1_min0=label_clump_significance_sanity(values, std, label, indexs,
                                            tl, sig, __func__);

  /* Allocate the arrays to keep the final significance measure (and
     possibly the indexs). */
  sig->ndim  = 1;                        /* Depends on 'cltprm->sn' */
  sig->type  = GAL_TYPE_FLOAT32;
  if(sig->dsize==NULL)
    sig->dsize = gal_pointer_allocate(GAL_TYPE_SIZE_T, 1, 0, __func__,
                                      "sig->dsize");
  sig->array = gal_pointer_allocate(sig->type, tablen, 0, __func__,
                                    "sig->array");
  sig->size  = sig->dsize[0] = tablen;  /* MUST BE AFTER dsize. */
  info=gal_pointer_allocate(GAL_TYPE_FLOAT64, tablen*INFO_NCOLS, 1,
                            __func__, "info");
  if( sigind )
    {
      sigind->ndim  = 1;
      sigind->type  = GAL_TYPE_INT32;
      sigind->dsize = gal_pointer_allocate(GAL_TYPE_SIZE_T, 1, 0, __func__,
                                           "sigind->dsize");
      sigind->size  = sigind->dsize[0] = tablen;/* After dsize */
      sigind->array = gal_pointer_allocate(sigind->type, tablen, 0, __func__,
                                           "sigind->array");
    }


  /* First, get the raw information necessary for making the S/N table. */
  label_clump_significance_raw(values, std, label, indexs, tl, info);


  /* Calculate the signficance value for successful clumps */
  sigarr=sig->array;
  if(keepsmall) sigarr[0]=NAN;
  if(sigind) indarr=sigind->array;
  for(i=1;i<tablen;++i)
    {
      /* For readability. */
      row = &info[ i * INFO_NCOLS ];

      /* If we have a sufficient area and any rivers were actually found
         for this clump, then do the measurement. */
      if( row[ INFO_INAREA ]>minarea && row[ INFO_RIVAREA ])
        {
          /* Set the index to write the values. If 'keepsmall' is not
             called, we don't care about the IDs of the clumps anymore, so
             store the signal-to-noise ratios contiguously. Note that
             counter will always be smaller and equal to i. */
          ind = keepsmall ? i : counter++;

          /* For easy reading. */
          R   = row[ INFO_PEAK_RIVER  ];
          C   = row[ INFO_PEAK_CENTER ];
          S   = variance ? sqrt(row[ INFO_STD ]) : row[ INFO_STD ];

          /* NEGATIVE VALUES: Rivers are also defined on the edges of the
             image and on pixels touching blank pixels. In a strong
             gradient, such sitations can cause the river to be
             larger/smaller than the minimum/maximum within the clump. This
             can only happen in very strong gradients, so for now, I think
             it is safe to ignore that clump (its negative value will
             automatically discard it). Later, if we find a problem with
             this, we'll have to figure out a better solution. */
          if(sigind) indarr[ind]=i;
          sigarr[ind] = ( max1_min0 ? (C-R) : (R-C) ) / S;
        }
      else
        {
          /* Only over detections, we should put a NaN when the S/N isn't
             calculated.  */
          if(keepsmall)
            {
              sigarr[i]=NAN;
              if(sigind) indarr[i]=i;
            }
        }
    }


  /* If we don't want to keep the small clumps, the size of the S/N table
     has to be corrected. */
  if(keepsmall==0)
    {
      sig->dsize[0] = sig->size = counter;
      if(sigind) sigind->dsize[0] = sigind->size = counter;
    }


  /* Clean up. */
  free(info);
}




















/**********************************************************************/
/*************               Growing labels               *************/
/**********************************************************************/
/* Grow the given labels without creating new ones. */
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
  label_check_type(indexs, GAL_TYPE_SIZE_T, "indexs", __func__);
  label_check_type(labels, GAL_TYPE_INT32,  "labels", __func__);
  if(indexs->ndim!=1)
    error(EXIT_FAILURE, 0, "%s: 'indexs' has to be a 1D array, but it is "
          "%zuD", __func__, indexs->ndim);

  /* The basic idea is this: after growing, not all the blank pixels are
     necessarily filled, for example the pixels might belong to two regions
     above the growth threshold. So the pixels in between them (which are
     below the threshold will not ever be able to get a label, even if they
     are in the indexs list). Therefore, the safest way we can terminate
     the loop of growing the objects is to stop it when the number of
     pixels left to fill in this round (thisround) equals the number of
     blanks.

     To start the loop, we set 'thisround' to one more than the number of
     indexed pixels. Note that it will be corrected immediately after the
     loop has started, it is just important to pass the 'while'. */
  thisround=ninds+1;
  while( thisround > ninds )
    {
      /* 'thisround' will keep the number of pixels to be inspected in this
         round. 'ninds' will count the number of pixels left without a
         label by the end of this round. Since 'ninds' comes from the
         previous loop (or outside, for the first round) it has to be saved
         in 'thisround' to begin counting a fresh. */
      thisround=ninds;
      ninds=0;

      /* Go over all the available indexs. NOTE: while the 'indexs->array'
         pointer remains unchanged, 'indexs->size' can/will change (get
         smaller) in every loop. */
      sf = (s=indexs->array) + indexs->size;
      do
        {
          /* We'll begin by assuming the nearest neighbor of this pixel
             has no label (has a value of 0). */
          n1=0;

          /* Check the neighbors of this pixel. Note that since this
             macro has multiple loops within it, we can't use
             break. We'll use the 'searchngb' variable instead. */
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
                             ('withrivers==0'), then there is no point in
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
             we just need to write in the value of 'n1'. */
          if(n1)
            {
              /* Set the label. */
              olabel[*s]=n1;

              /* If this pixel is a river (can only happen when
                 'withrivers' is zero), keep it in the loop, because we
                 want the 'indexs' dataset to contain all non-positive
                 (non-labeled) pixels, including rivers. */
              if(n1==GAL_LABEL_RIVER)
                iarray[ ninds++ ] = *s;
            }
          else
            iarray[ ninds++ ] = *s;

          /* Correct the size of the 'indexs' dataset. */
          indexs->size = indexs->dsize[0] = ninds;
        }
      while(++s<sf);
    }

  /* Clean up. */
  free(dinc);
}
