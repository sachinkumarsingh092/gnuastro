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

#include <gnuastro/tile.h>
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
    error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can fix "
          "this problem. In binary_2d_4con, the value to `dilate0_erode1' "
          "is %u while it should be 0 or 1", PACKAGE_BUGREPORT,
          dilate0_erode1);

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
    error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can fix "
          "this problem. In dilate0_erode1_4con (binary.c), the value to "
          "dilate0_erode1 is %u while it should be 0 or 1", PACKAGE_BUGREPORT,
          dilate0_erode1);

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
    error(EXIT_FAILURE, 0, "`gal_binary_erode' currently only works on a "
          "fully allocated block of memory");

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
            error(EXIT_FAILURE, 0, "%d not acceptable for connectivity in "
                  "a 2D dataset", connectivity);
          }
      break;

    default:
      error(EXIT_FAILURE, 0, "`gal_binary_erode' currently doesn't work on "
            "%zu dimensional datasets", binary->ndim);
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
   pixel) maybe already allocated (with type `uint32'). If `*out!=NULL',
   the labels will be reset to zero before the start and the labels will be
   written into it. If `*out==NULL', the necessary dataset will be
   allocated here and put into it. */
size_t
gal_binary_connected_components(gal_data_t *binary, gal_data_t **out,
                                int connectivity)
{
  uint32_t *l;
  uint8_t *b, *bf;
  gal_data_t *lab;
  size_t p, i, curlab=1;
  struct gal_linkedlist_sll *Q=NULL;
  size_t *dinc=gal_dimension_increment(binary->ndim, binary->dsize);

  /* Two small sanity checks. */
  if(binary->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "the input data structure to "
          "`gal_binary_connected_components' must be `uint8' type");
  if(binary->block)
    error(EXIT_FAILURE, 0, "currently, the input data structure to "
          "`gal_binary_connected_components' must not be a tile");


  /* Prepare the dataset for the labels. */
  if(*out)
    {
      /* Use the given dataset.  */
      lab=*out;

      /* Make sure the given dataset has the same size as the input. */
      if( gal_data_dsize_is_different(binary, lab) )
        error(EXIT_FAILURE, 0, "the `binary' and `out' datasets must have "
              "the same size in `gal_binary_connected_components'");

      /* Make sure it has a `uint32' type. */
      if( lab->type!=GAL_TYPE_UINT32 )
        error(EXIT_FAILURE, 0, "the `out' dataset in "
              "`gal_binary_connected_components' must have `uint32' type"
              "but the array you have given is `%s' type",
              gal_type_to_string(lab->type, 1));

      /* Reset all its values to zero. */
      memset(lab->array, 0, lab->size * gal_type_sizeof(lab->type));
    }
  else
    lab=*out=gal_data_alloc(NULL, GAL_TYPE_UINT32, binary->ndim,
                            binary->dsize, binary->wcs, 1,
                            binary->minmapsize, NULL, "labels", NULL);


  /* Initialize the labels array. If we have blank pixels in the byt
     array, then give them the blank labeled array. Note that since
     their value will not be 0, they will also not be labeled. */
  l=lab->array;
  bf=(b=binary->array)+binary->size;
  if( gal_blank_present(binary) )
    do *l++ = *b==GAL_BLANK_UINT8 ? GAL_BLANK_UINT32 : 0; while(++b<bf);


  /* Go over all the pixels and do a breadth-first: any pixel that is not
     labeled is used to label the full object by checking neighbors before
     going onto the next pixels. */
  l=lab->array;
  b=binary->array;
  for(i=0;i<binary->size;++i)
    /* Check if this pixel is already labeled. */
    if( b[i] && !l[i] )
      {
        /* This is the first pixel of this connected region that we have
           got to. */
        l[i]=curlab;

        /* Add this pixel to the queue of pixels to work with. */
        gal_linkedlist_add_to_sll(&Q, i);

        /* While a pixel remains in the queue, continue labelling and
           searching for neighbors. */
        while(Q!=NULL)
          {
            /* Pop an element from the queue. */
            gal_linkedlist_pop_from_sll(&Q, &p);

            /* Go over all its neighbors and add them to the list if they
               haven't already been labeled. */
            GAL_DIMENSION_NEIGHBOR_OP(p, binary->ndim, binary->dsize,
                                      connectivity, dinc,
              {
                if( b[ nind ] && !l[ nind ] )
                  {
                    l[ nind ] = curlab;
                    gal_linkedlist_add_to_sll(&Q, nind);
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
