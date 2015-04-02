/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <stdlib.h>
#include <gsl/gsl_sort.h>

#include "polygon.h"


/***************************************************************/
/**************       Basic operations        ******************/
/***************************************************************/

/* We have a simple polygon (that can result from projection, so its
   edges don't collide or it doesn't have holes) and we want to order
   its corners in an anticlockwise fashion. This is necessary for
   clipping it and finding its area later. Depending on the
   transformation, the corners can have practically any order even if
   before the transformation, they were ordered.

   The input is an array containing the coordinates (two values) of
   each corner. `n' is the number of corners. So the length of the
   input should be 2*n. The output is an array with `n' elements
   specifying the indexs in order. The reason the indexes are output
   is that for all the pixels in the image, in a homographic
   transform, the order is the same. So the input is unchanged, only
   `n' values will be put in the ordinds array. Such that calling the
   input coordinates in the following fashion will give an
   anti-clockwise order for 4 points for example:

   1st vertice: in[ordinds[0]*2], in[ordinds[0]*2+1]
   2nd vertice: in[ordinds[1]*2], in[ordinds[1]*2+1]
   3rd vertice: in[ordinds[2]*2], in[ordinds[2]*2+1]
   4th vertice: in[ordinds[3]*2], in[ordinds[3]*2+1]

   This is very similar to the Graham scan in finding the Convex
   Hull. However, in projection we will never have a concave polygon
   (the left condition below, where this algorithm will get to E
   before D), we will always have a convex polygon (right case) or E
   won't exist!

                 Convex Polygon        Concave Polygon

                  D --------C          D------------- C
                    \      |         E /            |
                     \E    |           \            |
                     /     |            \           |
                   A--------B             A ---------B

   This is because we are always going to be calculating the area of
   the overlap between a quadrilateral and the pixel grid or the
   quadrilaterral its self.

   MAXPOLYGONCORNERS is defined so there will be no need to allocate
   these temporary arrays seprately. Since we are dealing with pixels,
   the polygon can't really have too many vertices.
*/
void
orderedpolygoncorners(double *in, size_t n, size_t *ordinds)
{
  double angles[MAXPOLYGONCORNERS];
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
     corners might not be correct in non-trivial cases. */
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

  /* For a check:
  printf("\nAfter sorting:\n");
  for(i=0;i<n;++i)
    printf("%lu: %.3f, %.3f\n", ordinds[i], in[ordinds[i]*2],
           in[ordinds[i]*2+1]);
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
  return fabs(sum)/2.0f;
}





/* We have a polygon with `n' sides whose vertices are in the array
   `v' (with 2*n elements). Such that v[0], v[1] are the two coordinates
   of the first vertice. The vertices also have to be sorted in a
   counter clockwise fashion. We also have a point (with coordinates
   p[0], p[1]) and we want to see if it is inside the polygon or not.

   If the point is inside the polygon, it will always be to the left
   of the edge connecting the two vertices when the vertices are
   traversed in order. See the comments above `polygonarea' for an
   explanation about i and j and the loop.*/
int
pinpolygon(double *v, double *p, size_t n)
{
  size_t i=0, j=n-1;

  while(i<n)
    {
      if( pleftofline(&v[j*2], &v[i*2], p) )
        j=i++;
      else return 0;
    }
  return 1;
}





/* Similar to pinpolygon, except that if the point is on one of the
   edges of a polygon, this will return 0. */
int
ppropinpolygon(double *v, double *p, size_t n)
{
  size_t i=0, j=n-1;

  while(i<n)
    {
      /* For a check:
      printf("(%.2f, %.2f), (%.2f, %.2f), (%.2f, %.2f)=%f\n",
             v[j*2], v[j*2+1], v[i*2], v[i*2+1], p[0], p[1],
             tricrossproduct(&v[j*2], &v[i*2], p));
      */
      if( ppropleftofline(&v[j*2], &v[i*2], p) )
        j=i++;
      else
        return 0;
    }
  return 1;
}





/* Find the intersection of a line segment (Aa--Ab) and an infinite
   line (Ba--Bb) and put the intersection point in the output `o'. All
   the points are assumed to be two element double arrays already
   allocated outside this function.

   Return values:
    1: Intersection exists, value was put in the output.
    0: Intersection doesn't exist, no output written.
   -1: All four points are on the same line, no output written.
*/
int
seginfintersection(double *Aa, double *Ab, double *Ba, double *Bb,
                   double *o)
{
  int Aacollinear=pcollinearwithline(Ba, Bb, Aa);
  int Abcollinear=pcollinearwithline(Ba, Bb, Ab);

  /* If all four points lie on the same line, there is infinite
     intersections, so return -1. If three of the points are
     collinear, then the point on the line segment that is collinear
     with the infinite line is the intersection point.*/
  if( Aacollinear && Abcollinear)
    return -1;
  else if( Aacollinear || Abcollinear)
    {
      if(Aacollinear)  { o[0]=Aa[0]; o[1]=Aa[1]; }
      else             { o[0]=Ab[0]; o[1]=Ab[1]; }
      return 1;
    }

  /* None of Aa or Ab are collinear with the Ba--Bb line. They can
     only have an intersection if Aa and Ab are on opposite sides of
     the Ba--Bb. If they are on the same side of Ba--Bb, then there is
     no intersection (atleast within the line segment range
     (Aa--Ab).

     The intersection of two lines is calculated from the determinant
     form in the Wikipedia article of "line-line intersection" where

     x1=Ba[0]     x2=Bb[0]     x3=Aa[0]      x4=Ab[0]
     y1=Ba[1]     y2=Bb[1]     y3=Aa[1]      y4=Ab[1]

     Note that the denominators and the parenthesis with the
     subtraction of multiples are the same.
  */
  if ( ppropleftofline(Ba, Bb, Aa) ^ ppropleftofline(Ba, Bb, Ab) )
    {
      /* Find the intersection point of two infinite lines (we assume
         Aa-Ab is infinite in calculating this): */
      o[0]=( ( (Ba[0]*Bb[1]-Ba[1]*Bb[0])*(Aa[0]-Ab[0])
               - (Ba[0]-Bb[0])*(Aa[0]*Ab[1]-Aa[1]*Ab[0]) )
             / ( (Ba[0]-Bb[0])*(Aa[1]-Ab[1]) - (Ba[1]-Bb[1])*(Aa[0]-Ab[0]) )
             );
      o[1]=( ( (Ba[0]*Bb[1]-Ba[1]*Bb[0])*(Aa[1]-Ab[1])
               - (Ba[1]-Bb[1])*(Aa[0]*Ab[1]-Aa[1]*Ab[0]) )
             / ( (Ba[0]-Bb[0])*(Aa[1]-Ab[1]) - (Ba[1]-Bb[1])*(Aa[0]-Ab[0]) )
             );

      /* Now check if the output is in the same range as Aa--Ab. */
      if( o[0]>=minoftwo(Aa[0], Ab[0])-ROUNDERR
          && o[0]<=maxoftwo(Aa[0], Ab[0])+ROUNDERR
          && o[1]>=minoftwo(Aa[1], Ab[1])-ROUNDERR
          && o[1]<=maxoftwo(Aa[1], Ab[1])+ROUNDERR )
        return 1;
      else
        return 0;
    }
  else
    return 0;

}






/* Clip (find the overlap of) two polygons. This function uses the
   Sutherland-Hodgman polygon clipping psudocode from Wikipedia:

   List outputList = subjectPolygon;
   for (Edge clipEdge in clipPolygon) do
     List inputList = outputList;
     outputList.clear();
     Point S = inputList.last;
     for (Point E in inputList) do
        if (E inside clipEdge) then
           if (S not inside clipEdge) then
              outputList.add(ComputeIntersection(S,E,clipEdge));
           end if
           outputList.add(E);
        else if (S inside clipEdge) then
           outputList.add(ComputeIntersection(S,E,clipEdge));
        end if
        S = E;
     done
   done

   The difference is that we are not using lists, but arrays to keep
   polygon vertices. The two polygons are called Subject (`s') and
   Clip (`c') with `n' and `m' vertices respectively.

   The output is stored in `o' and the number of elements in the
   output are stored in what `*numcrn' (for number of corners) points
   to.*/
void
polygonclip(double *s, size_t n, double *c, size_t m,
            double *o, size_t *numcrn)
{
  double in[2*MAXPOLYGONCORNERS], *S, *E;
  size_t t, ii=m-1, i=0, jj, j, outnum, innum;

  /*
  if(n>MAXPOLYGONCORNERS || m>MAXPOLYGONCORNERS)
    error(EXIT_FAILURE, 0, "The two polygons given to the function "
          "polygonclip in polygon.c have %lu and %lu vertices. They cannot "
          "have any values larger than %lu.", n, m, MAXPOLYGONCORNERS);
  */

  /* 2*outnum because for each vertice, there are two elements. */
  outnum=n; for(t=0;t<2*outnum;++t) o[t]=s[t];

  while(i<m)                    /* clipEdge: c[ii*2] -- c[i*2]. */
    {
      /*
      printf("#################\nclipEdge %lu\n", i);
      printf("(%.3f, %.3f) -- (%.3f, %.3f)\n----------------\n",
             c[ii*2], c[ii*2+1], c[i*2], c[i*2+1]);
      */
      innum=outnum; for(t=0;t<2*innum;++t) in[t]=o[t];
      outnum=0;

      j=0;
      jj=innum-1;
      while(j<innum)
        {
          S=&in[jj*2];                      /* Starting point. */
          E=&in[j*2];                       /* Ending point.   */

          if( ppropleftofline(&c[ii*2], &c[i*2], E) )
            {
              if( ppropleftofline(&c[ii*2], &c[i*2], S)==0 )
                if( seginfintersection(S, E, &c[ii*2], &c[i*2], &o[2*outnum])
                    >0) ++outnum;
              o[2*outnum]=E[0]; o[2*outnum+1]=E[1]; ++outnum;
            }
          else if( ppropleftofline(&c[ii*2], &c[i*2], S) )
            if( seginfintersection(S, E, &c[ii*2], &c[i*2], &o[2*outnum])>0 )
              ++outnum;
          /*
          {
            size_t k;
            printf("(%.3f, %.3f) -- (%.3f, %.3f): %lu\n",
                   S[0], S[1], E[0], E[1], outnum);
            for(k=0;k<outnum;++k)
              printf("\t (%.3f, %.3f)\n", o[k*2], o[k*2+1]);
          }
          */

          jj=j++;
        }
      ii=i++;
      /*
      printf("outnum: %lu\n\n\n\n", outnum);
      */
    }
  *numcrn=outnum;
}
