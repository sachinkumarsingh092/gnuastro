/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (gnuastro) package.

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

#include <math.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <gsl/gsl_sort.h>

#include "fitsarrayvv.h"

#include "main.h"
#include "imgwarp.h"



/***************************************************************/
/**************            MACROS             ******************/
/***************************************************************/
/* Multiply a 2 element vector with a transformation matrix and put
   the result in the 2 element output array. It is assumed that the
   input is from a flat coordinate systemy. */
#define mappoint(V, T, O)                                       \
  {                                                             \
    (O)[0]=( ( (T)[0]*(V)[0] + (T)[1]*(V)[1] + (T)[2] )         \
             / ( (T)[6]*(V)[0] + (T)[7]*(V)[1] + (T)[8] ) );    \
    (O)[1]=( ( (T)[3]*(V)[0] + (T)[4]*(V)[1] + (T)[5] )         \
             / ( (T)[6]*(V)[0] + (T)[7]*(V)[1] + (T)[8] ) );    \
  }                                                             \





/* The cross product of two points from the center. */
#define crossproduct(A, B) ( (A)[0]*(B)[1] - (B)[0]*(A)[1] )




/* Find the cross product (2*area) between three points. Each point is
   assumed to be a pointer that has atleast two values within it. */
#define tricrossproduct(A, B, C)                  \
  ( ( (B)[0]-(A)[0] ) * ( (C)[1]-(A)[1] ) -       \
    ( (C)[0]-(A)[0] ) * ( (B)[1]-(A)[1] ) )       \





/* We have the line A-B. We want to see if C is to the left of this
   line or to its right. This function will return 1 if it is to the
   left. It uses the basic property of vector multiplication: If the
   three points are anti-clockwise (the point is to the left), then
   the vector multiplication is positive, if it is negative, then it
   is clockwise (c is to the right).

   Ofcourse it is very important that A be below or equal to B in both
   the X and Y directions. The rounding error might give
   -0.0000000000001 (I didn't count the number of zeros!!) instead of
   zero for the area. Zero would indicate that they are on the same
   line in this case this should give a true result.
*/
#define pleftofline(A, B, C)                            \
  tricrossproduct((A), (B), (C))<-ROUNDERR ? 0 : 1



















/***************************************************************/
/**************       Basic operations        ******************/
/***************************************************************/

/* We have a simple polygon (that can result from projection, so its
   edges don't collide or it doesn't have holes) and we want to order
   its corners in a counter an anticlockwise fashion. This is
   necessary for finding its area later and for checking if a pixel is
   within it. Depending on the transformation, the corners can have
   practically any order even if before the transformation, they were
   ordered.

   The input is an array containing the coordinates (two values) of
   each corner. `n' is the number of corners. So the length of the
   input should be 2*n. After this function the sides are going to be
   sorted in an anti-clockwise fashion.

   This is very similar to the Graham scan in finding the Convex
   Hull. However, in projection we will never have the left condition
   below (where this algorithm will get to E before D), we will always
   have the right case or E won't exist!

                  D --------C          D------------- C
                    \      |         E /            |
                     \E    |           \            |
                     /     |            \           |
                   A--------B             A ---------B

   This is because we are always going to be calculating the area of
   the overlap between a quadrilateral and the pixel grid or the
   quadrilaterral its self.

   MAXPOLYGONCORNERS is defined so there will be no need to allocate
   these temporary arrays seprately.
*/
void
orderedpolygoncorners(double *in, size_t n)
{
  size_t ordinds[MAXPOLYGONCORNERS];
  double angles[MAXPOLYGONCORNERS], ordvalues[2*MAXPOLYGONCORNERS];
  size_t i, tmp, aindexs[MAXPOLYGONCORNERS], tindexs[MAXPOLYGONCORNERS];

  if(n>MAXPOLYGONCORNERS)
    error(EXIT_FAILURE, 0, "Most probably a bug! The number of corners "
          "given to `orderedpolygoncorners' is more than %d. This is an "
          "internal value and cannot be set from the outside. Most probably "
          "Some bug has caused this un-normal value. Please contact us at "
          PACKAGE_BUGREPORT" so we can solve this problem.",
          MAXPOLYGONCORNERS);

  /* For a check:
  printf("\n\nBefore sorting:\n");
  for(i=0;i<n;++i)
    printf("%lu: %.3f, %.3f\n", i, in[i*2], in[i*2+1]);
  printf("\n");
  */

  /* Find the point with the smallest Y (if there is two of them, the
     one with the smallest X too.). This is necessary because if the
     angles are not found relative to this point, the ordering of the
     corners might not be correct in non-trivial cases. The number of
     points is usually very small (less than 10), so this is
     insignificant in the grand scheme of processes that have to take
     place in warping the image. */
  gsl_sort_index(ordinds, in+1, 2, n);
  if( in[ ordinds[0]*2+1 ] == in[ ordinds[1]*2+1 ]
     && in[ ordinds[0]*2] > in[ ordinds[1]*2 ])
    {
      tmp=ordinds[0];
      ordinds[0]=ordinds[1];
      ordinds[1]=tmp;
    }


  /* We only have `n-1' more elements to sort, use the angle of the
     line between the three remaining points and the first point. */
  for(i=0;i<n-1;++i)
    angles[i]=atan2( in[ ordinds[i+1]*2+1 ] - in[ ordinds[0]*2+1 ],
                     in[ ordinds[i+1]*2 ]   - in[ ordinds[0]*2   ] );
  /* For a check:
  for(i=0;i<n-1;++i)
    printf("%lu: %.3f degrees\n", ordinds[i+1], angles[i]*180/M_PI);
  printf("\n");
  */

  /* Sort the angles into the correct order, we need an extra array to
     temporarily keep the newly angle-ordered indexs. Without it we
     are going to loose half of the ordinds indexs! */
  gsl_sort_index(aindexs, angles, 1, n-1);
  for(i=0;i<n-1;++i) tindexs[i]=ordinds[aindexs[i]+1];
  for(i=0;i<n-1;++i) ordinds[i+1]=tindexs[i];

  /* Put the values in the correct order: */
  for(i=0;i<n;++i)
    {
      ordvalues[i*2]=in[ordinds[i]*2];
      ordvalues[i*2+1]=in[ordinds[i]*2+1];
    }
  for(i=0;i<2*n;++i) in[i]=ordvalues[i];


  /* For a check:
  printf("\nAfter sorting:\n");
  for(i=0;i<n;++i)
    printf("%lu: %.3f, %.3f\n", ordinds[i], in[i*2], in[i*2+1]);
  printf("\n");
  */
}





/* The area of a polygon is the sum of the vector products of all the
   vertices in a counterclockwise order. See the Wikipedia page for
   Polygon for more information.

   `v' points to an array of doubles which keep the positions of the
   vertices such that v[0] and v[1] are the positions of the first
   corner to be considered.

   We will start from the edge connecting the last pixel to the first
   pixel for the first step of the loop, for the rest, j is always
   going to be one less than i.
 */
double
polygonarea(double *v, size_t n)
{
  double sum=0.0f;
  size_t i=0, j=n-1;

  while(i<n)
    {
      sum+=crossproduct(v+j*2, v+i*2);
      j=i++;
    }
  return sum/2.0f;
}





/* We have a quadrilateral (polygon with four sides) whose vertices
   are in the array `v'. Such that v[0], v[1] are the two coordinates
   of the first vertice. The vertices also have to be sorted in a
   counter clockwise fashion. We also have a point (with coordinates
   p[0], p[1]) and we want to see if it is inside the polygon or
   not.

   If the point is inside the polygon, it will always be to the left
   of the edge connecting the two vertices when the vertices are
   traversed in order. See the comments above `polygonarea' for an
   explanation about i and j and the loop.*/
int
pinquadrilateral(double *v, double *p)
{
  size_t i=0, j=3;

  while(i<4)
    {
      if( pleftofline(v+j*2, v+i*2, p) ) j=i++;
      else return 0;
    }
  return 1;
}




















/***************************************************************/
/**************      Processing function      ******************/
/***************************************************************/
void *
imgwarponthread(void *inparam)
{
  struct iwpparams *iwp=(struct iwpparams*)inparam;
  struct imgwarpparams *p=iwp->p;

  size_t is0=p->is0, is1=p->is1;
  double *outfpixval=p->outfpixval;
  double ocrn[8], icrn[8], *output=p->output;
  size_t i, j, ind, os1=p->os1, *extinds=p->extinds;

  for(i=0;(ind=iwp->indexs[i])!=NONTHRDINDEX;++i)
    {
      /* Set the corners of this output pixel. */
      ocrn[0]=ind%os1+outfpixval[0];     ocrn[1]=ind/os1+outfpixval[1];
      ocrn[2]=ind%os1+1+outfpixval[0];   ocrn[3]=ind/os1+outfpixval[1];
      ocrn[4]=ind%os1+outfpixval[0];     ocrn[5]=ind/os1+1+outfpixval[1];
      ocrn[6]=ind%os1+1+outfpixval[0];   ocrn[7]=ind/os1+1+outfpixval[1];

      /* Transform the four corners to the input image */
      for(j=0;j<4;++j) mappoint(&ocrn[j*2], p->inverse, &icrn[j*2]);

      /* For a check:
      printf("\n\n\n%lu, %lu:\n", ind%os1+1, ind/os1+1);
      for(j=0;j<4;++j) printf("%f, %f\n", icrn[j*2], icrn[j*2+1]);
      printf("%.16f, %.16f, %.16f, %.16f.\n", icrn[extinds[0]],
             icrn[extinds[1]], icrn[extinds[2]], icrn[extinds[3]]);
      */

      /* In case the four extremes of this output pixel are outside
         the range of the input image, then go onto the next pixel. To
         be completely outside the image all four corners have to be
         outside the image range.

         Rounding error might cause the maximum values to be
         0.00000000001 (I haven't counted the zeros!!!). So the
         maximum points are closer than ROUNDERR to zero and still
         postive, consider them to be negative.
      */
      if( icrn[extinds[0]]>is0 || icrn[extinds[1]]<=ROUNDERR
          || icrn[extinds[2]]>is1 || icrn[extinds[3]]<=ROUNDERR )
        continue;
      /* For a test: */
      else
        output[ind]=1.0f;



      /* Find the borders of the pixels which are to be used. */


    }



  /* Wait until all other threads finish. */
  if(p->cp.numthreads>1)
    pthread_barrier_wait(iwp->b);

  return NULL;
}


















/***************************************************************/
/**************          Preparations         ******************/
/***************************************************************/
/* Do all the preparations.

   Make the output array. We transform the four corners of the image
   into the output space. To find the four sides of the image.

   About fpixel and lpixel. The point is that we don't want to spend
   time, transforming any pixels which we know will not be in the
   input image.

   Find the proper order of transformed pixel corners from the output
   array to the input array. The order is fixed for all the pixels in
   the image altough the scale might change.
*/
void
imgwarppreparations(struct imgwarpparams *p)
{
  double *d, *df, output[8];
  size_t i, *extinds=p->extinds;
  double ocrn[8]={0,0,1,0,0,1,1,1}, icrn[8]={0,0,0,0,0,0,0,0};
  double xmin=DBL_MAX, xmax=-DBL_MAX, ymin=DBL_MAX, ymax=-DBL_MAX;
  double input[8]={0.0f, 0.0f, p->is1, 0.0f, 0.0f, p->is0, p->is1, p->is0};

  /* Find the range of pixels of the input image. */
  for(i=0;i<4;++i)
    {
      mappoint(&input[i*2], p->matrix, &output[i*2]);
      if(output[i*2]<xmin)     xmin = output[i*2];
      if(output[i*2]>xmax)     xmax = output[i*2];
      if(output[i*2+1]<ymin)   ymin = output[i*2+1];
      if(output[i*2+1]>ymax)   ymax = output[i*2+1];
    }
  /* For a check:
  for(i=0;i<4;++i)
      printf("(%.3f, %.3f) --> (%.3f, %.3f)\n",
             input[i*2], input[i*2+1], output[i*2], output[i*2+1]);
  printf("xmin: %.3f\nxmax: %.3f\nymin: %.3f\nymax: %.3f\n",
         xmin, xmax, ymin, ymax);
  */

  /* Set the final size of the image. The X axis is horizontal. we are
     using the bottom left corner, that is why we are adding a 1. */
  p->outfpixval[0]=floor(xmin);
  p->outfpixval[1]=floor(ymin);
  p->os0=(long)ymax-(long)ymin+1;
  p->os1=(long)xmax-(long)xmin+1;
  if(ymin*ymax<0.0f) ++p->os0;
  if(xmin*xmax<0.0f) ++p->os1;
  /* For a check:
  printf("Wrapped:\n");printf("os1: %lu\nos0: %lu\noutfpixval=(%.3f,%.3f)\n",
         p->os1, p->os0, p->outfpixval[0], p->outfpixval[1]);
  */

  /* In case the point (0.0f, 0.0f) has moved and the user has asked
     to incorporte that shift (by not calling the --wrap option), then
     change the relevant parameters.

     output[0], output[1] show the coordinates of the new origin.

     NOTE: To incorporate non-integer shifts, the borders of the
     output pixels have to have no fractional values.
  */
  if(p->wrap==0)
    for(i=0;i<2;++i)
    {
      if(output[i]>0.0f) p->outfpixval[i]-=(long)output[i];
      if(i==0) p->os1+=labs(output[i]);
      else     p->os0+=labs(output[i]);
    }
  /* For a check:
  printf("Corrected:\nos1: %lu\nos0: %lu\noutfpixval=(%.3f,%.3f)\n",
         p->os1, p->os0, p->outfpixval[0], p->outfpixval[1]);
  */

  /* We now know the size of the output and the starting and ending
     coordinates in the output image (bottom left corners of pixels)
     for the transformation. */
  errno=0;
  p->output=malloc(p->os0*p->os1*sizeof *p->output);
  if(p->output==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for the output array",
          p->os0*p->os1*sizeof *p->output);
  df=(d=p->output)+p->os0*p->os1; do *d++=NAN; while(d<df);


  /* Order the corners of the inverse-transformed pixel (from the
     output to the input) in an anti-clockwise transformation. In a
     general homographic transform, the scales of the output pixels
     may change, but the relative positions of the corners will
     not.  */
  for(i=0;i<4;++i) mappoint(&ocrn[i*2], p->inverse, &icrn[i*2]);
  orderedpolygoncorners(icrn, 4);


  /* Find which index after transformation will have the minimum and
     maximum positions along the two axises. We can't use the starting
     loop because that is based on the input image which can be
     not-a-square! So we do it here where pixels are squares. */
  xmin=DBL_MAX; xmax=-DBL_MAX; ymin=DBL_MAX; ymax=-DBL_MAX;
  for(i=0;i<4;++i)
    {
      if(icrn[i*2]<xmin)   { xmin=icrn[i*2];   extinds[0]=i*2;   }
      if(icrn[i*2]>xmax)   { xmax=icrn[i*2];   extinds[1]=i*2;   }
      if(icrn[i*2+1]<ymin) { ymin=icrn[i*2+1]; extinds[2]=i*2+1; }
      if(icrn[i*2+1]>ymax) { ymax=icrn[i*2+1]; extinds[3]=i*2+1; }
    }

  printf("\n%f\n", polygonarea(icrn, 4));
}




















/***************************************************************/
/**************       Outside function        ******************/
/***************************************************************/
void
imgwarp(struct imgwarpparams *p)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct iwpparams *iwp;
  size_t nt=p->cp.numthreads;
  size_t i, nb, *indexs, thrdcols;


  /* Array keeping thread parameters for each thread. */
  errno=0;
  iwp=malloc(nt*sizeof *iwp);
  if(iwp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in imgwarp "
          "(imgwarp.c) for iwp", nt*sizeof *iwp);


  /* Prepare the output array and all the necessary things: */
  imgwarppreparations(p);


  /* Distribute the output pixels into the threads: */
  distinthreads(p->os0*p->os1, nt, &indexs, &thrdcols);


  /* Start the convolution. */
  if(nt==1)
    {
      iwp[0].p=p;
      iwp[0].indexs=indexs;
      imgwarponthread(&iwp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
	 (that spinns off the nt threads) is also a thread, so the
	 number the barrier should be one more than the number of
	 threads spinned off. */
      if(p->os0*p->os1<nt) nb=p->os0*p->os1+1;
      else nb=nt+1;
      attrbarrierinit(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=NONTHRDINDEX)
          {
            iwp[i].p=p;
            iwp[i].b=&b;
            iwp[i].indexs=&indexs[i*thrdcols];
	    err=pthread_create(&t, &attr, imgwarponthread, &iwp[i]);
	    if(err)
	      error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }

  /* Save the output: */
  arraytofitsimg(p->cp.output, "Warped", DOUBLE_IMG, p->output, p->os0,
                 p->os1, NULL, SPACK_STRING);

  /* Free the allocated spaces: */
  free(iwp);
  free(indexs);
  free(p->output);
}
