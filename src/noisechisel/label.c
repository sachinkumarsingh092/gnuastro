/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/neighbors.h>
#include <gnuastro/linkedlist.h>

#include "main.h"

#include "label.h"
#include "binary.h"










/**********************************************************************/
/**************     Connected component labeling    *******************/
/**********************************************************************/
/* Find the connected components in an image based on the breadth
   first algorithm. The output curlab is one more than the actual
   number of objects in the array. The input image is a unsigned char
   binary, either background (==0) or foreground (==1), image. This
   function will find and label the different components of all
   pixels labeled b0_f1.*/
size_t
BF_concmp(unsigned char *byt, long *lab, size_t s0, size_t s1,
          int anyblank, const size_t connectivity)
{
  struct gal_linkedlist_sll *Q=NULL;
  long *l=lab, curlab=1; /* Current label */
  size_t i, p, size=s0*s1, s0t1=s0-1, s1t1=s1-1;
  unsigned char *b=byt, *bf=byt+s0*s1, counter, bl, br, tl, tr;


  /* A simple sanity check. */
  if(connectivity!=4 && connectivity!=8)
    error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can fix "
          "the problem. For some reason, the value to connectivity in "
          "BF_concmp (label.c) is %lu which is not recognized",
          PACKAGE_BUGREPORT, connectivity);


  /* Initialize the labels array. If we have blank pixels in the byt
     array, then give them the blank labeled array. Note that since
     their value will not be 0, they will also not be labeled. */
  if(anyblank)
    do *l++ = *b==GAL_FITS_BYTE_BLANK ? GAL_FITS_LONG_BLANK
         : 0; while(++b<bf);
  else
    memset(lab, 0, size*sizeof *lab);

  /* Go over all the pixels: */
  if(connectivity==4)
    {
      for(i=0;i<size;++i)
        /* Check if it is not needed or already done: */
        if(byt[i] && !lab[i])
          {
            lab[i]=curlab;
            add_to_sll(&Q, i);
            while(Q!=NULL)
              {
                /* Pop from the queue */
                pop_from_sll(&Q, &p);

                /* Check the four connected neighbors: */
                if(p/s1>0    && byt[p-s1] && !lab[p-s1])
                  {add_to_sll(&Q, p-s1); lab[p-s1]=curlab;}
                if(p/s1<s0t1 && byt[p+s1] && !lab[p+s1])
                  {add_to_sll(&Q, p+s1); lab[p+s1]=curlab;}
                if(p%s1>0    && byt[p-1]  && !lab[p-1])
                  {add_to_sll(&Q, p-1); lab[p-1]=curlab;}
                if(p%s1<s1t1 && byt[p+1]  && !lab[p+1])
                  {add_to_sll(&Q, p+1); lab[p+1]=curlab;}
              }
            ++curlab;
          }
    }
  else
    {
      for(i=0;i<size;++i)
        /* Check if it is not needed or already done: */
        if(byt[i] && !lab[i])
          {
            lab[i]=curlab;
            add_to_sll(&Q, i);

            while(Q!=NULL)
              {

                pop_from_sll(&Q, &p);

                /* Set the counters for the corners in 8-connectivity */
                counter=bl=br=tl=tr=0;

                /* Check the four connected neighbors: */
                if(p/s1>0)
                  {
                    ++counter; ++bl; ++br;
                    if(byt[p-s1] && !lab[p-s1])
                      {add_to_sll(&Q, p-s1); lab[p-s1]=curlab;}
                  }
                if(p/s1<s0t1)
                  {
                    ++counter; ++tl; ++tr;
                    if(byt[p+s1] && !lab[p+s1])
                      {add_to_sll(&Q, p+s1); lab[p+s1]=curlab;}
                  }
                if(p%s1>0)
                  {
                    ++counter; ++bl; ++tl;
                    if(byt[p-1]  && !lab[p-1])
                      {add_to_sll(&Q, p-1); lab[p-1]=curlab;}
                  }
                if(p%s1<s1t1)
                  {
                    ++counter; ++tr; ++br;
                    if(byt[p+1]  && !lab[p+1])
                      {add_to_sll(&Q, p+1); lab[p+1]=curlab;}
                  }

                if(counter==4)  /* All four corners are in the image. */
                  {
                    if(byt[p-s1-1] && !lab[p-s1-1])
                      {add_to_sll(&Q, p-s1-1); lab[p-s1-1]=curlab;}
                    if(byt[p-s1+1] && !lab[p-s1+1])
                      {add_to_sll(&Q, p-s1+1); lab[p-s1+1]=curlab;}
                    if(byt[p+s1-1] && !lab[p+s1-1])
                      {add_to_sll(&Q, p+s1-1); lab[p+s1-1]=curlab;}
                    if(byt[p+s1+1] && !lab[p+s1+1])
                      {add_to_sll(&Q, p+s1+1); lab[p+s1+1]=curlab;}
                  }
                else            /* At least one corner isn't in the image. */
                  {
                    if(bl==2 && byt[p-s1-1] && !lab[p-s1-1])
                      {add_to_sll(&Q, p-s1-1); lab[p-s1-1]=curlab;}
                    if(br==2 && byt[p-s1+1] && !lab[p-s1+1])
                      {add_to_sll(&Q, p-s1+1); lab[p-s1+1]=curlab;}
                    if(tl==2 && byt[p+s1-1] && !lab[p+s1-1])
                      {add_to_sll(&Q, p+s1-1); lab[p+s1-1]=curlab;}
                    if(tr==2 && byt[p+s1+1] && !lab[p+s1+1])
                      {add_to_sll(&Q, p+s1+1); lab[p+s1+1]=curlab;}
                  }
              }
            ++curlab;
          }
    }

  return curlab;
}





/* This function applied the same principles of the above function but
   on an adjacency matrix. Its ouput is an array the size of one side
   of the adjacency matrix that will have the label each object should
   have.

   The adjacency matrix should be zero (for no connection) and
   non-zero for connection (note that it should be symmetric).
*/
size_t
BF_concomp_AdjMatrix(int *adj, size_t numside, long **outnewlabs)
{
  size_t i, j, p;
  struct gal_linkedlist_sll *Q=NULL;
  long *newlabs, curlab=1;

  errno=0;
  newlabs=calloc(numside, sizeof *newlabs);
  if(newlabs==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for newlabs in "
          "BF_concomp_AdjMatrix (label.c)", numside*sizeof *newlabs);

  for(i=1;i<numside;++i)
    if(newlabs[i]==0)
      {
        add_to_sll(&Q, i);
        while(Q!=NULL)
          {
            pop_from_sll(&Q, &p);
            if(newlabs[p]!=curlab)
              {
                newlabs[p]=curlab;
                for(j=1;j<numside;++j)
                  if( adj[p*numside+j] && newlabs[j]==0 )
                    add_to_sll(&Q, j);
              }
          }
        ++curlab;
      }

  /* For a check:
  for(i=1;i<numside;++i)
    printf("%lu: %ld\n", i, newlabs[i]);
  */

  *outnewlabs=newlabs;
  return curlab;
}




















/**********************************************************************/
/**************            Work on labels           *******************/
/**********************************************************************/
/* Find the areas of labeled regions in an array. `numlabs`, which is
   the number of labels in the array has to be one larger than the
   largest label in the array. */
void
labareas(long *lab, size_t size, size_t numlabs, size_t **outareas)
{
  size_t *areas;
  long *last=lab+size;

  /* Allocate the areas array: */
  errno=0;
  areas=*outareas=calloc(numlabs, sizeof *areas);
  if(areas==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for areas in labareas (label.c)",
          numlabs*sizeof *areas);

  /* Find the area of each label. */
  do if(ISINDEXABLELABEL) ++areas[*lab]; while(++lab<last);

  /* For a check:
  for(i=0;i<numlabs;++i)
    printf("%lu: %ld\n", i, a[i]);
  printf("\n\n\n\n\n");
  */
}





void
removesmallarea_relabel(long *in, unsigned char *byt, size_t size,
                        size_t *numlabs, size_t minarea)
{
  size_t i, *areas;
  long *newlabs, curlab=1;

  /* Allocate the space to keep the new labels. */
  errno=0; newlabs=calloc(*numlabs, sizeof *newlabs);
  if(newlabs==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for newlabs in "
          "removesmallarea_relabel (label.c)", *numlabs*sizeof *newlabs);

  /* Find the areas: */
  labareas(in, size, *numlabs, &areas);
  areas[0]=0;

  /* Set the new labels: */
  for(i=1;i<*numlabs;++i)
    if(areas[i]>minarea)
      newlabs[i]=curlab++;

  if(byt)
    {
      for(i=0;i<size;++i)
        if(in[i]!=GAL_FITS_LONG_BLANK)
          byt[i] = (in[i]=newlabs[in[i]]) > 0;
    }
  else
    {
      for(i=0;i<size;++i)
        if(in[i]!=GAL_FITS_LONG_BLANK)
          in[i]=newlabs[in[i]];
    }

  *numlabs=curlab;

  free(areas);
  free(newlabs);
}





/* Make an array of pointers to arrays that keep the indexes of
   different labeled regions in the image. Note that the label zero
   is not going to be considered. */
void
labindexs(long *inlab, size_t size, size_t numlabs, size_t **outareas,
          size_t ***outlabinds)
{
  long *lab, *lp;
  size_t i, *areas, *counters, **labinds;

  /* Allocate the pointers to the objects index array and set the
     pointer to label zero to NULL. */
  errno=0;
  labinds=*outlabinds=malloc(numlabs*sizeof *labinds);
  if(labinds==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for labinds in labindexs "
          "(label.c)", numlabs*sizeof *labinds);
  labinds[0]=NULL;

  /* Find the area of each label: */
  labareas(inlab, size, numlabs, outareas);
  areas=*outareas;
  areas[0]=0;

  /* For each label, counter will keep the position of the last filled
     element. */
  errno=0;
  counters=calloc(numlabs, sizeof *counters);
  if(counters==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for counters in labindexs "
          "(label.c)", numlabs*sizeof *counters);

  /* Allocate space for each label's indexs: */
  for(i=1;i<numlabs;++i)
    {
      errno=0;
      labinds[i]=malloc(areas[i]*sizeof **labinds);
      if(labinds[i]==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for labinds[%lu] in "
              "labindexs (label.c)", areas[i]*sizeof **labinds, i);
    }

  /* Fill in the indexs array. */
  lp=(lab=inlab)+size;
  do
    if(ISINDEXABLELABEL)
      labinds[*lab][counters[*lab]++]=lab-inlab;
  while(++lab<lp);
  labinds[0]=NULL;

  /* For a check:
  {
    size_t j;
    for(i=1;i<numlabs;++i)
      {
        printf("Lab: %lu\n", i);
        for(j=0;j<areas[i];++j)
          printf("  %lu\n", labinds[i][j]);
      }
    exit(0);
  }
  */
  free(counters);
}
