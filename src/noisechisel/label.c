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

#include "linkedlist.h"

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
  size_t i, p, size, s0t1, s1t1;
  long curlab=1; /* Current label */

  if(connectivity!=4)
    error(EXIT_FAILURE, 0, "Currently BF_concmp (label.c) is only set for "
          "4 connectivity.");

  size=s0*s1;
  s0t1=s0-1;
  s1t1=s1-1;
  memset(lab, 0, size*sizeof *lab);

  /* Go over all the pixels: */
  for(i=0;i<size;++i)

    /* Check if it is not needed or already done: */
    if(byt[i] && !lab[i])
      {
	add_to_sll(&Q, i);

	while(Q!=NULL)
	  {
	    pop_from_sll(&Q, &p);
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
