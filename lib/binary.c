/*********************************************************************
binary -- Work on binary (0 and 1 valued) datasets.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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
#include <gnuastro/tile.h>
#include <gnuastro/blank.h>
#include <gnuastro/binary.h>
#include <gnuastro/dimension.h>
#include <gnuastro/linkedlist.h>









/*********************************************************************/
/*****************      Erosion and dilation      ********************/
/*********************************************************************/
static void
binary_erode_dilate_2d_4con(gal_data_t *input, int dilate0_erode1)
{
  uint8_t f, b, *pt, *fpt, *byt=input->array;
  size_t i, j, ind, nr=input->dsize[0], nc=input->dsize[1];

  /* Do a sanity check: */
  if(dilate0_erode1!=1 && dilate0_erode1!=0)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can "
          "fix this problem. The value to `dilate0_erode1' is %u while it "
          "should be 0 or 1", __func__, PACKAGE_BUGREPORT, dilate0_erode1);

  /* Set the foreground and background values. */
  if(dilate0_erode1==0) {f=1; b=0;}
  else                  {f=0; b=1;}

  /* Check the 4 corners: */
  if(byt[0]==b && (byt[1]==f || byt[nc]==f) )
    byt[0]=GAL_BINARY_TMP_VALUE;

  if(byt[nc-1]==b && (byt[nc-2]==f || byt[2*nc-1]==f))
    byt[nc-1]=GAL_BINARY_TMP_VALUE;

  if(byt[(nr-1)*nc]==b
     && (byt[(nr-2)*nc]==f || byt[(nr-1)*nc+1]==f) )
    byt[(nr-1)*nc]=GAL_BINARY_TMP_VALUE;

  if(byt[nr*nc-1]==b
     && (byt[nr*nc-2]==f || byt[nr*nc-1-nc]==f) )
    byt[nr*nc-1]=GAL_BINARY_TMP_VALUE;

  /* Check the 4 sides: */
  for(j=1;j<nc-1;++j)
    if(byt[j]==b
       && (byt[j+1]==f || byt[j-1]==f || byt[j+nc]==f) )
      byt[j]=GAL_BINARY_TMP_VALUE;

  for(j=1;j<nc-1;++j)
    {
      ind=(nr-1)*nc+j;
      if(byt[ind]==b
         && (byt[ind+1]==f || byt[ind-1]==f || byt[ind-nc]==f) )
        byt[ind]=GAL_BINARY_TMP_VALUE;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=i*nc;
      if(byt[ind]==b
         && (byt[ind+1]==f || byt[ind+nc]==f || byt[ind-nc]==f) )
        byt[ind]=GAL_BINARY_TMP_VALUE;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=(i+1)*nc-1;
      if(byt[ind]==b
         && (byt[ind-1]==f || byt[ind+nc]==f || byt[ind-nc]==f) )
        byt[ind]=GAL_BINARY_TMP_VALUE;
    }

  /* Check the body: */
  for(i=1;i<nr-1;++i)
    for(j=1;j<nc-1;++j)
      {
        ind=i*nc+j;
        if(byt[ind]==b
           && (byt[ind-1]==f     || byt[ind+1]==f
               || byt[ind+nc]==f || byt[ind-nc]==f) )
          byt[ind]=GAL_BINARY_TMP_VALUE;
      }

  /* Set all the changed pixels to the proper values: */
  fpt=(pt=byt)+nr*nc;
  do *pt = *pt==GAL_BINARY_TMP_VALUE ? f : *pt; while(++pt<fpt);
}





/* 8 connected dilation and erosion. b0_f1==0: Dilate the
   foreground. b0_f1==1: Erode the foreground. */
static void
binary_erode_dilate_2d_8con(gal_data_t *input, unsigned char dilate0_erode1)
{
  uint8_t f, b, *pt, *fpt, *byt=input->array;
  size_t i, j, ind, nr=input->dsize[0], nc=input->dsize[1];

  /* Do a sanity check: */
  if(dilate0_erode1!=1 && dilate0_erode1!=0)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can fix "
          "this problem. The value to dilate0_erode1 is %u while it should "
          "be 0 or 1", __func__, PACKAGE_BUGREPORT, dilate0_erode1);

  /* Set the foreground and background values: */
  if(dilate0_erode1==0) {f=1; b=0;}
  else         {f=0; b=1;}

  /* Check the 4 corners: */
  if(byt[0]==b && (byt[1]==f
                   || byt[nc]==f || byt[nc+1]==f) )
    byt[0]=GAL_BINARY_TMP_VALUE;

  if(byt[nc-1]==b && (byt[nc-2]==f
                      || byt[2*nc-1]==f
                      || byt[2*nc-2]==f) )
    byt[nc-1]=GAL_BINARY_TMP_VALUE;

  if(byt[(nr-1)*nc]==b
     && ( byt[(nr-2)*nc]==f || byt[(nr-1)*nc+1]==f
          || byt[(nr-2)*nc+1]==f) )
    byt[(nr-1)*nc]=GAL_BINARY_TMP_VALUE;

  if(byt[nr*nc-1]==b
     && ( byt[nr*nc-2]==f || byt[nr*nc-1-nc]==f
          || byt[nr*nc-2-nc]==f) )
    byt[nr*nc-1]=GAL_BINARY_TMP_VALUE;

  /* Check the 4 sides: */
  for(j=1;j<nc-1;++j)
    if(byt[j]==b
       && ( byt[j+1]==f || byt[j-1]==f || byt[j+nc]==f
            || byt[j-1+nc]==f || byt[j+1+nc]==f) )
      byt[j]=GAL_BINARY_TMP_VALUE;

  for(j=1;j<nc-1;++j)
    {
      ind=(nr-1)*nc+j;
      if(byt[ind]==b
         && ( byt[ind+1]==f || byt[ind-1]==f || byt[ind-nc]==f
              || byt[ind-1-nc]==f || byt[ind+1-nc]==f) )
        byt[ind]=GAL_BINARY_TMP_VALUE;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=i*nc;
      if(byt[ind]==b
         && ( byt[ind+1]==f || byt[ind+nc]==f || byt[ind-nc]==f
              || byt[ind+1-nc]==f || byt[ind+1+nc]==f) )
        byt[ind]=GAL_BINARY_TMP_VALUE;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=(i+1)*nc-1;
      if(byt[ind]==b
         && (byt[ind-1]==f || byt[ind+nc]==f || byt[ind-nc]==f
             || byt[ind-1-nc]==f || byt[ind-1+nc]==f) )
        byt[ind]=GAL_BINARY_TMP_VALUE;
    }

  /* Check the body: */
  for(i=1;i<nr-1;++i)
    for(j=1;j<nc-1;++j)
      {
        ind=i*nc+j;
        if(byt[ind]==b
           && (byt[ind-1]==f        || byt[ind+1]==f
               || byt[ind+nc]==f    || byt[ind-nc]==f
               || byt[ind-1-nc]==f  || byt[ind+1+nc]==f
               || byt[ind-1+nc]==f  || byt[ind+1-nc]==f) )
          byt[ind]=GAL_BINARY_TMP_VALUE;
      }

  /* Set all the changed pixels to the proper values: */
  fpt=(pt=byt)+nr*nc;
  do *pt = *pt==GAL_BINARY_TMP_VALUE ? f : *pt; while(++pt<fpt);
}





/* Erode a binary dataset any number of times. If `inplace' is given a
   value of `1', then do the erosion within the already allocated space,
   otherwise, allocate a new array and save the result into that.

   This function will only work on the elements with a value of 1 or 0. It
   will leave all the rest unchanged. Also note that it only works on
   `uint8_t' type datasets. So if the input doesn't have that type, it is
   going to copy it this type and return the newlyallocated dataset. So
   when the input's type isn't `uint8_t', `inplace' is irrelevant. */
static gal_data_t *
binary_erode_dilate(gal_data_t *input, size_t num, int connectivity,
                    int inplace, int d0e1)
{
  size_t counter;
  gal_data_t *binary;
  size_t *dinc=gal_dimension_increment(input->ndim, input->dsize);

  /* Currently this only works on blocks. */
  if(input->block)
    error(EXIT_FAILURE, 0, "%s: currently only works on a fully "
          "allocated block of memory, but the input is a tile (its `block' "
          "element is not NULL)", __func__);

  /* Set the dataset to work on. */
  binary = ( (inplace && input->type==GAL_TYPE_UINT8)
             ? input :
             gal_data_copy_to_new_type(input, GAL_TYPE_UINT8) );

  /* Go over every element and do the erosion. */
  switch(binary->ndim)
    {
    case 2:
      for(counter=0;counter<num;++counter)
        switch(connectivity)
          {
          case 4: binary_erode_dilate_2d_4con(binary, d0e1); break;
          case 8: binary_erode_dilate_2d_8con(binary, d0e1); break;
          default:
            error(EXIT_FAILURE, 0, "%s: %d not acceptable for connectivity "
                  "in a 2D dataset", __func__, connectivity);
          }
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: currently doesn't work on %zu "
            "dimensional datasets", __func__, binary->ndim);
    }

  /* Clean up and return. */
  free(dinc);
  return binary;
}





gal_data_t *
gal_binary_erode(gal_data_t *input, size_t num, int connectivity,
                 int inplace)
{
  return binary_erode_dilate(input, num, connectivity, inplace, 1);
}





gal_data_t *
gal_binary_dilate(gal_data_t *input, size_t num, int connectivity,
                  int inplace)
{
  return binary_erode_dilate(input, num, connectivity, inplace, 0);
}





gal_data_t *
gal_binary_open(gal_data_t *input, size_t num, int connectivity,
                int inplace)
{
  gal_data_t *out;

  /* First do the necessary number of erosions. */
  out=gal_binary_erode(input, num, connectivity, inplace);

  /* If `inplace' was called, then `out' is the same as `input', if it
     wasn't, then `out' is a newly allocated array. In any case, we should
     dilate in the same allocated space. */
  gal_binary_dilate(input, num, connectivity, 1);

  /* Return the output dataset. */
  return out;
}




















/*********************************************************************/
/*****************      Connected components      ********************/
/*********************************************************************/
/* Find the connected components in the input binary dataset `binary'
   through the breadth first search algorithm. `binary' has to have an
   `uint8' datatype and only zero and non-zero values in it will be
   distinguished. The output dataset (which will contain a label on each
   pixel) maybe already allocated (with type `int32'). If `*out!=NULL', the
   labels will be reset to zero before the start and the labels will be
   written into it. If `*out==NULL', the necessary dataset will be
   allocated here and put into it. */
size_t
gal_binary_connected_components(gal_data_t *binary, gal_data_t **out,
                                int connectivity)
{
  int32_t *l;
  uint8_t *b, *bf;
  gal_data_t *lab;
  size_t p, i, curlab=1;
  gal_list_sizet_t *Q=NULL;
  size_t *dinc=gal_dimension_increment(binary->ndim, binary->dsize);

  /* Two small sanity checks. */
  if(binary->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "%s: the input data set type must be `uint8'",
          __func__);
  if(binary->block)
    error(EXIT_FAILURE, 0, "%s: currently, the input data structure to "
          "must not be a tile", __func__);


  /* Prepare the dataset for the labels. */
  if(*out)
    {
      /* Use the given dataset.  */
      lab=*out;

      /* Make sure the given dataset has the same size as the input. */
      if( gal_data_dsize_is_different(binary, lab) )
        error(EXIT_FAILURE, 0, "%s: the `binary' and `out' datasets must have "
              "the same size", __func__);

      /* Make sure it has a `int32' type. */
      if( lab->type!=GAL_TYPE_INT32 )
        error(EXIT_FAILURE, 0, "%s: the `out' dataset must have `int32' type"
              "but the array you have given is `%s' type", __func__,
              gal_type_name(lab->type, 1));

      /* Reset all its values to zero. */
      memset(lab->array, 0, lab->size * gal_type_sizeof(lab->type));
    }
  else
    lab=*out=gal_data_alloc(NULL, GAL_TYPE_INT32, binary->ndim,
                            binary->dsize, binary->wcs, 1,
                            binary->minmapsize, NULL, "labels", NULL);


  /* Initialize the labels array. If we have blank pixels in the byt
     array, then give them the blank labeled array. Note that since
     their value will not be 0, they will also not be labeled. */
  l=lab->array;
  bf=(b=binary->array)+binary->size; /* Library must have no side effect,   */
  if( gal_blank_present(binary, 0) ) /* So blank flag should not be changed.*/
    do *l++ = *b==GAL_BLANK_UINT8 ? GAL_BLANK_INT32 : 0; while(++b<bf);


  /* Go over all the pixels and do a breadth-first: any pixel that is not
     labeled is used to label the full object by checking neighbors before
     going onto the next pixels. */
  l=lab->array;
  b=binary->array;
  for(i=0;i<binary->size;++i)
    /* Check if this pixel is already labeled. */
    if( b[i] && l[i]==0 )
      {
        /* This is the first pixel of this connected region that we have
           got to. */
        l[i]=curlab;

        /* Add this pixel to the queue of pixels to work with. */
        gal_list_sizet_add(&Q, i);

        /* While a pixel remains in the queue, continue labelling and
           searching for neighbors. */
        while(Q!=NULL)
          {
            /* Pop an element from the queue. */
            p=gal_list_sizet_pop(&Q);

            /* Go over all its neighbors and add them to the list if they
               haven't already been labeled. */
            GAL_DIMENSION_NEIGHBOR_OP(p, binary->ndim, binary->dsize,
                                      connectivity, dinc,
              {
                if( b[ nind ] && l[ nind ]==0 )
                  {
                    l[ nind ] = curlab;
                    gal_list_sizet_add(&Q, nind);
                  }
              } );
          }

        /* This object has been fully labeled, so increment the current
           label. */
        ++curlab;
      }


  /* Clean up and return the total number. */
  free(dinc);
  return curlab-1;
}





/* Given an adjacency matrix (which should be binary), find the number of
   connected objects and return an array of new labels for each old
   label. In other words, this function will find the objects that are
   connected (possibly through a third object) and in the output array, the
   respective elements for all these objects is going to have the same
   value. The total number of connected labels is put into the place
   pointed to by `numconnected'.

   Labels begin from 1 (0 is kept for non-labeled regions usually). So if
   you have 3 initial objects/labels, the input matrix to this function
   should have a size of 4x4. The first (label 0) row and column are not
   going to be parsed/checked.

   The adjacency matrix needs to be completely filled (on both sides of the
   diagonal) for this function.

   If the input adjacency matrix has a size of `amsize * amsize', the
   output will have a size of `amsize' with each index having a new label
   in its place. */
gal_data_t *
gal_binary_connected_adjacency_matrix(gal_data_t *adjacency,
                                      size_t *numconnected)
{
  gal_data_t *newlabs_d;
  gal_list_sizet_t *Q=NULL;
  int32_t *newlabs, curlab=1;
  uint8_t *adj=adjacency->array;
  size_t i, j, p, num=adjacency->dsize[0];

  /* Some small sanity checks. */
  if(adjacency->type != GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "%s: input must have type `uint8'. However, the "
          "input dataset has type of `%s'", __func__,
          gal_type_name(adjacency->type, 1));

  if(adjacency->ndim != 2)
    error(EXIT_FAILURE, 0, "%s: input must be 2-dimensional (a matrix)."
          "However, the input dataset has %zu dimensions", __func__,
          adjacency->ndim);

  if(adjacency->dsize[0] != adjacency->dsize[1])
    error(EXIT_FAILURE, 0, "%s: input must be square (same length in both "
          "dimensions). However, the input dataset has a size of %zu x %zu",
          __func__, adjacency->dsize[0], adjacency->dsize[1]);


  /* Allocate (and clear) the output datastructure. */
  newlabs_d=gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &num, NULL, 1,
                           adjacency->minmapsize, NULL, NULL, NULL);
  newlabs=newlabs_d->array;


  /* Go over the input matrix and apply the same principle as we used to
     identify connected components in an image: through a queue, find those
     elements that are connected. */
  for(i=1;i<num;++i)
    if(newlabs[i]==0)
      {
        /* Add this old label to the list that must be corrected. */
        gal_list_sizet_add(&Q, i);

        /* Continue while the list has elements. */
        while(Q!=NULL)
          {
            /* Pop the top old-label from the list. */
            p=gal_list_sizet_pop(&Q);

            /* If it has already been labeled then ignore it. */
            if( newlabs[p]!=curlab )
              {
                /* Give it the new label. */
                newlabs[p]=curlab;

                /* Go over the adjacecny matrix row for this touching
                   object and see if there are any not-yet-labeled objects
                   that are touching it. */
                for(j=1;j<num;++j)
                  if( adj[ p*num+j ] && newlabs[j]==0 )
                    gal_list_sizet_add(&Q, j);
              }
          }

        /* Increment the current label. */
        ++curlab;
      }


  /* For a check.
  printf("=== Old labels --> new labels ===\n");
  for(i=1;i<num;++i) printf("%zu: %u\n", i, newlabs[i]);
  */

  /* Return the output. */
  *numconnected = curlab-1;
  return newlabs_d;
}




















/*********************************************************************/
/*****************            Fill holes          ********************/
/*********************************************************************/
/* Make the array that is the inverse of the input byt of fill
   holes. The inverse array will also be 4 pixels larger in both
   dimensions. This is because we might also want to fill those holes
   that are touching the side of the image. One pixel for a pixel that
   is one pixel away from the image border. Another pixel for those
   objects that are touching the image border. */
static gal_data_t *
binary_make_padded_inverse(gal_data_t *input, gal_data_t **outtile)
{
  uint8_t *in;
  size_t i, startind;
  gal_data_t *inv, *tile;
  size_t *dsize=gal_data_malloc_array(GAL_TYPE_SIZE_T, input->ndim);
  size_t *startcoord=gal_data_malloc_array(GAL_TYPE_SIZE_T, input->ndim);


  /* Set the size of the padded inverse image and the coordinates of the
     start. We want the inverse to be padded on the edges of each dimension
     by 2 pixels, so each dimension should be padded by 4 pixels. */
  for(i=0;i<input->ndim;++i)
    {
      startcoord[i]=2;
      dsize[i]=input->dsize[i]+4;
    }


  /* Allocate the inverse dataset and initialize it to 1 (mainly for the
     edges, the inner region will be set afterwards).

     PADDING MUST BE INITIALIZED WITH 1: This is done so the connected body
     of 1 valued pixels (after inversion) gets a label of 1 after labeling
     the connected components and any hole, will get a value larger than
     1. */
  inv=gal_data_alloc(NULL, GAL_TYPE_UINT8, input->ndim, dsize, NULL, 0,
                     input->minmapsize, "INVERSE", "binary", NULL);
  memset(inv->array, 1, inv->size);


  /* Define a tile to fill the central regions of the inverse. */
  startind=gal_dimension_coord_to_index(input->ndim, inv->dsize, startcoord);
  tile=gal_data_alloc(gal_data_ptr_increment(inv->array, startind, inv->type),
                      inv->type, input->ndim, input->dsize, NULL, 0, 0, NULL,
                      NULL, NULL);
  *outtile=tile;
  tile->block=inv;


  /* Put the input's flags into the inverted array and the tile. */
  inv->flag = tile->flag = input->flag;


  /* Fill the central regions. */
  in=input->array;
  GAL_TILE_PARSE_OPERATE({*i = *in==GAL_BLANK_UINT8 ? *in : !*in; ++in;},
                         tile, NULL, 0, 0);


  /* Clean up and return. */
  free(dsize);
  free(startcoord);
  return inv;
}





/* Fill all the holes in an input unsigned char array that are bounded
   within a 4-connected region.

   The basic method is this:

   1. An inverse image is created:

        * For every pixel in the input that is 1, the inverse is 0.

        * The inverse image has two extra pixels on each edge to
          ensure that all the inv[i]==1 pixels around the image are
          touching each other and a diagonal object passing through
          the image does not cause the inv[i]==1 pixels on the edges
          of the image to get a different label.

   2. The 8 connected regions in this inverse image are found.

   3. Since we had a 2 pixel padding on the edges of the image, we
      know for sure that all labeled regions with a label of 1 are
      actually connected `holes' in the input image.

      Any pixel with a label larger than 1, is therefore a bounded
      hole that is not 8-connected to the rest of the holes.  */
void
gal_binary_fill_holes(gal_data_t *input)
{
  uint8_t *in;
  gal_data_t *inv, *tile, *holelabs=NULL;

  /* A small sanity check. */
  if( input->type != GAL_TYPE_UINT8 )
    error(EXIT_FAILURE, 0, "%s: input must have `uint8' type, but its "
          "input dataset has `%s' type", __func__,
          gal_type_name(input->type, 1));


  /* Make the inverse image. */
  inv=binary_make_padded_inverse(input, &tile);


  /* Label the 8-connected (connectivity==2) holes. */
  gal_binary_connected_components(inv, &holelabs, 2);


  /* Any pixel with a label larger than 1 is a hole in the input image and
     we should invert the respective pixel. To do it, we'll use the tile
     that was defined before, just change its block and array.*/
  in=input->array;
  tile->array=gal_tile_block_relative_to_other(tile, holelabs);
  tile->block=holelabs; /* has to be after correcting `tile->array'. */


  /* The type of the tile is already known (it is `int32_t') and we have no
     output, so we'll just put `int' as a place-holder. In this way we can
     avoid the switch statement of GAL_TILE_PARSE_OPERATE, and directly
     use the workhorse macro `GAL_TILE_PO_OISET'. */
  GAL_TILE_PO_OISET(int32_t, int, {
      *in = *i>1 && *i!=GAL_BLANK_INT32 ? 1 : *in;
      ++in;
    }, tile, NULL, 0, 0);


  /* Clean up and return. */
  tile->array=NULL;
  gal_data_free(inv);
  gal_data_free(tile);
  gal_data_free(holelabs);
}
