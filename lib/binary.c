/*********************************************************************
binary -- Work on binary (0 and 1 valued) datasets.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>









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
          "fix this problem. The value to 'dilate0_erode1' is %u while it "
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
  else                  {f=0; b=1;}

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





/* This is a general erosion and dilation function. It is less efficient
   than the more specialized cases above. */
static void
binary_erode_dilate_general(gal_data_t *input, unsigned char dilate0_erode1,
                            int connectivity)
{
  uint8_t f, b, *pt, *fpt, *byt=input->array;
  size_t i, *dinc=gal_dimension_increment(input->ndim, input->dsize);

  /* Do a sanity check: */
  if(dilate0_erode1!=1 && dilate0_erode1!=0)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can fix "
          "this problem. The value to dilate0_erode1 is %u while it should "
          "be 0 or 1", __func__, PACKAGE_BUGREPORT, dilate0_erode1);

  /* Set the foreground and background values: */
  if(dilate0_erode1==0) {f=1; b=0;}
  else                  {f=0; b=1;}

  /* Go over the neighbors of each pixel. */
  for(i=0;i<input->size;++i)
    if(byt[i]==b)
      GAL_DIMENSION_NEIGHBOR_OP(i, input->ndim, input->dsize, connectivity,
                                dinc,{
                                  if(byt[i]!=GAL_BINARY_TMP_VALUE
                                     && byt[nind]==f)
                                    byt[i]=GAL_BINARY_TMP_VALUE;
                                });

  /* Set all the changed pixels to the proper values: */
  fpt=(pt=byt)+input->size;
  do *pt = *pt==GAL_BINARY_TMP_VALUE ? f : *pt; while(++pt<fpt);
}





/* Erode a binary dataset any number of times. If 'inplace' is given a
   value of '1', then do the erosion within the already allocated space,
   otherwise, allocate a new array and save the result into that.

   This function will only work on the elements with a value of 1 or 0. It
   will leave all the rest unchanged. Also note that it only works on
   'uint8_t' type datasets. So if the input doesn't have that type, it is
   going to copy it this type and return the newlyallocated dataset. So
   when the input's type isn't 'uint8_t', 'inplace' is irrelevant. */
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
          "allocated block of memory, but the input is a tile (its 'block' "
          "element is not NULL)", __func__);

  /* Set the dataset to work on. */
  binary = ( (inplace && input->type==GAL_TYPE_UINT8)
             ? input
             : gal_data_copy_to_new_type(input, GAL_TYPE_UINT8) );

  /* Go over every element and do the erosion. */
  switch(binary->ndim)
    {
    case 2:
      for(counter=0;counter<num;++counter)
        switch(connectivity)
          {
          case 1: binary_erode_dilate_2d_4con(binary, d0e1); break;
          case 2: binary_erode_dilate_2d_8con(binary, d0e1); break;
          default:
            error(EXIT_FAILURE, 0, "%s: %d not acceptable for connectivity "
                  "in a 2D dataset", __func__, connectivity);
          }
      break;

    case 3:
      for(counter=0;counter<num;++counter)
        binary_erode_dilate_general(binary, d0e1, connectivity);
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

  /* If 'inplace' was called, then 'out' is the same as 'input', if it
     wasn't, then 'out' is a newly allocated array. In any case, we should
     dilate in the same allocated space. */
  gal_binary_dilate(input, num, connectivity, 1);

  /* Return the output dataset. */
  return out;
}




















/*********************************************************************/
/*****************      Connected components      ********************/
/*********************************************************************/
/* Find connected components in an intput dataset. */
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
    error(EXIT_FAILURE, 0, "%s: the input data set type must be 'uint8'",
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
      if( gal_dimension_is_different(binary, lab) )
        error(EXIT_FAILURE, 0, "%s: the 'binary' and 'out' datasets must "
              "have the same size", __func__);

      /* Make sure it has a 'int32' type. */
      if( lab->type!=GAL_TYPE_INT32 )
        error(EXIT_FAILURE, 0, "%s: the 'out' dataset must have 'int32' type"
              "but the array you have given is '%s' type", __func__,
              gal_type_name(lab->type, 1));

      /* Reset all its values to zero. */
      memset(lab->array, 0, lab->size * gal_type_sizeof(lab->type));
    }
  else
    lab=*out=gal_data_alloc(NULL, GAL_TYPE_INT32, binary->ndim,
                            binary->dsize, binary->wcs, 1,
                            binary->minmapsize, binary->quietmmap,
                            NULL, "labels", NULL);


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





/* Put the indexs of connected labels in a list of 'gal_data_t's, each with
   a one-dimensional array that has the indexs of that connected
   component.*/
#define BINARY_CONINDEX_VAL 2
gal_data_t *
gal_binary_connected_indexs(gal_data_t *binary, int connectivity)
{
  uint8_t *b, *bf;
  gal_data_t *lines=NULL;
  size_t p, i, onelabnum, *onelabarr;
  gal_list_sizet_t *Q=NULL, *onelab=NULL;
  size_t *dinc=gal_dimension_increment(binary->ndim, binary->dsize);

  /* Small sanity checks. */
  if(binary->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "%s: the input data set type must be 'uint8'",
          __func__);
  if(binary->block)
    error(EXIT_FAILURE, 0, "%s: currently, the input data structure to "
          "must not be a tile", __func__);

  /* Go over all the pixels and do a breadth-first search. */
  b=binary->array;
  for(i=0;i<binary->size;++i)
    /* A pixel that has already been recorded is given a value of
       'BINARY_CONINDEX_VAL'. */
    if( b[i]==1 )
      {
        /* Add this pixel to the queue of pixels to work with. */
	b[i]=BINARY_CONINDEX_VAL;
        gal_list_sizet_add(&Q, i);
        gal_list_sizet_add(&onelab, i);

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
                if( b[nind]==1 )
                  {
		    b[nind]=BINARY_CONINDEX_VAL;
                    gal_list_sizet_add(&Q, nind);
		    gal_list_sizet_add(&onelab, nind);
                  }
              } );
          }

	/* Parsing has finished, put all the indexs into an array. */
	onelabarr=gal_list_sizet_to_array(onelab, 1, &onelabnum);
	gal_list_data_add_alloc(&lines, onelabarr, GAL_TYPE_SIZE_T, 1,
				&onelabnum, NULL, 0, -1, 1, NULL, NULL, NULL);

	/* Clean up. */
	gal_list_sizet_free(onelab);
	onelab=NULL;
      }

  /* Reverse the order. */
  gal_list_data_reverse(&lines);

  /* For a check
  {
    gal_data_t *test=lines->next;
    size_t *b, *bf;
    bf=(b=test->array)+test->size;
    do printf("%zu\n", *b++);
    while(b<bf);
    exit(0);
  }
  */

  /* Set all the '2' values back to '1'. */
  bf=(b=binary->array)+binary->size;
  do if(*b==BINARY_CONINDEX_VAL) *b=1; while(++b<bf);

  /* Clean up and return the total number. */
  free(dinc);
  return lines;
}





/* Given an adjacency matrix (which should be binary), find the number of
   connected objects and return an array of new labels for each old
   label.  */
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
    error(EXIT_FAILURE, 0, "%s: input must have type 'uint8'. However, the "
          "input dataset has type of '%s'", __func__,
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
                           adjacency->minmapsize, adjacency->quietmmap,
                           NULL, NULL, NULL);
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

                /* Go over the adjacency matrix row for this touching
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
  size_t *startcoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim, 0,
                                          __func__, "startcoord");
  size_t *dsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, input->ndim, 0,
                                     __func__, "dsize");


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
                     input->minmapsize, input->quietmmap,
                     "INVERSE", "binary", NULL);
  memset(inv->array, 1, inv->size);


  /* Define a tile to fill the central regions of the inverse. */
  startind=gal_dimension_coord_to_index(input->ndim, inv->dsize, startcoord);
  tile=gal_data_alloc(gal_pointer_increment(inv->array, startind, inv->type),
                      inv->type, input->ndim, input->dsize, NULL, 0, 0, 0,
                      NULL, NULL, NULL);
  *outtile=tile;
  tile->block=inv;


  /* Put the input's flags into the inverted array and the tile. */
  inv->flag = tile->flag = input->flag;


  /* Fill the central regions. */
  in=input->array;
  GAL_TILE_PARSE_OPERATE( tile, NULL, 0, 0,
                          {*i = *in==GAL_BLANK_UINT8 ? *in : !*in; ++in;} );


  /* Clean up and return. */
  free(dsize);
  free(startcoord);
  return inv;
}





gal_data_t *
gal_binary_holes_label(gal_data_t *input, int connectivity,
                       size_t *numholes)
{
  size_t d;
  int32_t *lab;
  gal_data_t *inv, *tile, *holelabs=NULL;

  /* A small sanity check. */
  if( input->type != GAL_TYPE_UINT8 )
    error(EXIT_FAILURE, 0, "%s: input must have 'uint8' type, but its "
          "input dataset has '%s' type", __func__,
          gal_type_name(input->type, 1));


  /* Make the inverse image. */
  inv=binary_make_padded_inverse(input, &tile);


  /* Label the holes. Recall that the first label is just the undetected
     regions, so we should subtract that from the total number.*/
  *numholes=gal_binary_connected_components(inv, &holelabs, connectivity);
  *numholes -= 1;


  /* Any pixel with a label larger than 1 is a hole in the input image and
     we should invert the respective pixel. To do it, we'll use the tile
     that was defined before, just change its block and array.*/
  tile->array=gal_tile_block_relative_to_other(tile, holelabs);
  tile->block=holelabs; /* has to be after correcting 'tile->array'. */


  /* The type of the tile is already known (it is 'int32_t') and we have no
     output/other, so we'll just put 'int' as a place-holder. In this way
     we can avoid the switch statement of GAL_TILE_PARSE_OPERATE, and
     directly use the workhorse macro 'GAL_TILE_PO_OISET'. */
  lab=(holelabs)->array;
  GAL_TILE_PO_OISET(int32_t, int, tile, NULL, 0, 0, {
      *lab++ = ( *i
                 ? ( *i==1
                     ? 0        /* Originally, was background. */
                     : *i-1 )   /* Real label: -1 (background has 1). */
                 : -1 );        /* Originally, was foreground. */
    });


  /* Clean up */
  tile->array=NULL;
  gal_data_free(inv);
  gal_data_free(tile);


  /* Correct the sizes of the hole labels array. We have already filled the
     from the start, effectively removing the paddings. Therefore, ee will
     just correct the sizes and we won't bother actually re-allocating the
     array size in memory. According to the GNU C library's description
     after 'realloc': "In several allocation implementations, making a
     block smaller sometimes necessitates copying it, so it can fail if no
     other space is available.". The extra padding is only 2 pixels wide,
     thus, the extra space is negligible compared to the actual array. So
     it isn't worth possibly having to copy the whole array to another
     location. Later, when we free the space, the kernel knows how much the
     allocated size is. */
  for(d=0;d<input->ndim;++d)
    holelabs->dsize[d] = input->dsize[d];
  holelabs->size=input->size;


  /* Return the number of holes.  */
  return holelabs;
}





/* Fill all the holes in an input unsigned char array.

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
      actually connected 'holes' in the input image.

      Any pixel with a label larger than 1, is therefore a bounded
      hole that is not 8-connected to the rest of the holes.  */
void
gal_binary_holes_fill(gal_data_t *input, int connectivity, size_t maxsize)
{
  uint8_t *in;
  uint32_t *i, *fi;
  size_t numholes, *sizes;
  gal_data_t *inv, *tile, *holelabs=NULL;

  /* Small sanity checks. */
  if( input->type != GAL_TYPE_UINT8 )
    error(EXIT_FAILURE, 0, "%s: input must have 'uint8' type, but its "
          "input dataset has '%s' type", __func__,
          gal_type_name(input->type, 1));
  if(connectivity<1 || connectivity>input->ndim)
    error(EXIT_FAILURE, 0, "%s: connectivity value %d is not acceptable. "
          "It has to be between 1 and the number of input's dimensions "
          "(%zu)", __func__, connectivity, input->ndim);


  /* Make the inverse image. */
  inv=binary_make_padded_inverse(input, &tile);


  /* Label the holes */
  numholes=gal_binary_connected_components(inv, &holelabs, connectivity);


  /* Any pixel with a label larger than 1 is a hole in the input image and
     we should invert the respective pixel. To do it, we'll use the tile
     that was defined before, just change its block and array.*/
  in=input->array;
  tile->array=gal_tile_block_relative_to_other(tile, holelabs);
  tile->block=holelabs; /* has to be after correcting 'tile->array'. */

  /* If the user wants to only fill holes to a certain size, then remove
     those with a larger size. */
  if(maxsize<-1)
    {
      /* Allocate space to keep the size of each hole: */
      sizes=gal_pointer_allocate(GAL_TYPE_SIZE_T, numholes+1, 1, __func__,
                                 "sizes");
      fi=(i=holelabs->array)+holelabs->size; do ++sizes[*i]; while(++i<fi);

      /* Set those labels with a larger size to 1 (treat it as
         background). */
      fi=(i=holelabs->array)+holelabs->size;
      do
        if(*i!=GAL_BLANK_INT32) *i = sizes[*i]>maxsize ? 1 : *i;
      while(++i<fi);

      /* Clean up. */
      free(sizes);
    }

  /* The type of the tile is already known (it is 'int32_t') and we have no
     output, so we'll just put 'int' as a place-holder. In this way we can
     avoid the switch statement of GAL_TILE_PARSE_OPERATE, and directly use
     the workhorse macro 'GAL_TILE_PO_OISET'. */
  GAL_TILE_PO_OISET(int32_t, int, tile, NULL, 0, 0, {
      *in = *i>1 && *i!=GAL_BLANK_INT32 ? 1 : *in;
      ++in;
    });


  /* Clean up and return. */
  tile->array=NULL;
  gal_data_free(inv);
  gal_data_free(tile);
  gal_data_free(holelabs);
}
