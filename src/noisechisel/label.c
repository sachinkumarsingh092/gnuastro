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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include "neighbors.h"
#include "linkedlist.h"
#include "fitsarrayvv.h"

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
          const size_t connectivity)
{
  struct sll *Q=NULL;
  long curlab=1; /* Current label */
  unsigned char counter, bl, br, tl, tr;
  size_t i, p, size=s0*s1, s0t1=s0-1, s1t1=s1-1;

  /* A simple sanity check. */
  if(connectivity!=4 && connectivity!=8)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix "
          "the problem. For some reason, the value to connectivity in "
          "BF_concmp (label.c) is %lu which is not recognized.",
          PACKAGE_BUGREPORT, connectivity);

  /* Set all the labels to zero. */
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
                    if(bl==2 && byt[p-s1-1] && !lab[p-s1-1])
                      {add_to_sll(&Q, p-s1-1); lab[p-s1-1]=curlab;}
                    if(br==2 && byt[p-s1+1] && !lab[p-s1+1])
                      {add_to_sll(&Q, p-s1+1); lab[p-s1+1]=curlab;}
                    if(tl==2 && byt[p+s1-1] && !lab[p+s1-1])
                      {add_to_sll(&Q, p+s1-1); lab[p+s1-1]=curlab;}
                    if(tr==2 && byt[p+s1+1] && !lab[p+s1+1])
                      {add_to_sll(&Q, p+s1+1); lab[p+s1+1]=curlab;}
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
  struct sll *Q=NULL;
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
labareas(long *in, size_t size, size_t numlabs, size_t **areas)
{
  size_t i, *a;

  errno=0; a=*areas=calloc(numlabs, sizeof *a);
  if(a==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for areas in labareas (label.c)",
          numlabs*sizeof *a);

  for(i=0;i<size;++i)
    ++a[in[i]];

  /*
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
    for(i=0;i<size;++i)
      byt[i] = (in[i]=newlabs[in[i]]) > 0;
  else
    for(i=0;i<size;++i)
      in[i]=newlabs[in[i]];

  *numlabs=curlab;

  free(areas);
  free(newlabs);
}





/* Calculate the Signal to noise ratio for all the labels in
   `labinmesh'. The sky subtracted input image and the standard
   deviations image for each pixel are over the full image: p->img and
   p->std. However, the array keeping the labels of each object is
   only over the mesh under consideration.
*/
void
detlabelsn(struct noisechiselparams *p, long *labinmesh, size_t *numlabs,
           size_t start, size_t s0, size_t s1, float **outsntable)
{
  float *img=p->img;
  double *fluxs, err;
  struct meshparams *smp=&p->smp;
  size_t i, ind, *xys, r=0, is1=p->smp.s1, counter=0;
  float *f, *ff, *s, ave, *sntable, cpscorr=p->cpscorr;
  long *areas, *l=labinmesh, detsnminarea=p->detsnminarea;


  /* Allocate the necessary arrays to keep the label information: */
  errno=0; sntable=*outsntable=calloc(*numlabs, sizeof *sntable);
  if(sntable==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for sntable in labelsn (label.c)",
          *numlabs*sizeof *sntable);
  errno=0; fluxs=calloc(*numlabs, sizeof *fluxs);
  if(fluxs==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for fluxs in labelsn (label.c)",
          *numlabs*sizeof *fluxs);
  errno=0; xys=calloc(*numlabs*2, sizeof *xys);
  if(xys==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for xys in labelsn (label.c)",
          *numlabs*sizeof *xys);
  errno=0; areas=calloc(*numlabs, sizeof *areas);
  if(areas==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for areas in labelsn (label.c)",
          *numlabs*sizeof *areas);


  /* For a check:
  arraytofitsimg("test.fits", "Labs", LONG_IMG, labinmesh,
                 s0, s1, 0, NULL, NULL, SPACK_STRING);
  */


  /* Go over the arrays and fill in the information. Note that since
     labinmesh is over the mesh and not the full image, setting it
     once in the definition was enough. */
  do
    {
      ff = ( f = img + start + r++ * is1 ) + s1;
      do
        {
          if(*l && !isnan(*f))     /* There is a label on this mesh. */
            {
              ++areas[*l];
              fluxs[*l]   += *f;
              xys[*l*2]   += (f-img)/is1;
              xys[*l*2+1] += (f-img)%is1;

            }
          ++l; ++s;
        }
      while(++f<ff);
    }
  while(r<s0);


  /* calculate the signal to noise for successful detections: */
  for(i=1;i<*numlabs;++i)
    if(areas[i]>detsnminarea && (ave=fluxs[i]/areas[i]) >0 )
      {
        /* Find the weighted center of this object and get the
           standard deviation at this point. Note that the standard
           deviation on the grid was stored in smp->garray2. The error
           should then be taken to the power of two and if the sky is
           subtracted, a 2 should be multiplied to it.*/
        ind=(xys[i*2]/areas[i]) * is1 + xys[i*2+1]/areas[i];
        err = smp->garray2[imgindextomeshid(smp, ind)];
        err *= p->skysubtracted ? err : 2.0f*err;

        /* Set the index in the sntable to store the Signal to noise
           ratio. When we are dealing with the noise, we only want the
           non-zero signal to noise values, so we will just use a
           counter. But for initial detections, it is very important
           that their Signal to noise ratio be placed in the same
           index as their label. */
        ind = p->b0f1 ? i : counter++;
        sntable[ind]=sqrt( (float)(areas[i])/cpscorr ) * ave / sqrt(ave+err);
      }
  if(p->b0f1==0) *numlabs=counter;


  /* For a check:
  for(i=0;i<*numlabs;++i)
    printf("%lu: %-10ld %-10.3f %-10.3f %-10.3f\n",
           i, areas[i], fluxs[i], stds[i]/areas[i], sntable[i]);
  */


  /* Clean up */
  free(xys);
  free(fluxs);
  free(areas);
}
