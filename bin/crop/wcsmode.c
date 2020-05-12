/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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

#include <math.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/wcs.h>
#include <gnuastro/pointer.h>

#include "main.h"

#include "onecrop.h"










/*******************************************************************/
/****************        Check for ui.c        *********************/
/*******************************************************************/
/* This function is called from ui.c. Its job is to check the WCS
   values of this  */
void
wcsmode_check_prepare(struct cropparams *p, struct inputimgs *img)
{
  double *pixscale;
  struct wcsprm *wcs=img->wcs;
  size_t *dsize=img->dsize, ndim=img->ndim;

  /* For two dimensions, we have four corners (2 numbers for each) and for
     three dimensions we have 8 corners (3 numbers for each), so we'll just
     assume the largest. */
  int i, status[8]={0,0,0,0,0,0,0,0}, ncorners=0;
  double imgcrd[24], phi[8], theta[8], pixcrd[24];


  /* Check if the image is aligned with the WCS coordinates. Note that
     because of small floating point errors, some programs might still keep
     very small values in the off-diagonal matrix elements. */
  if( fabs(wcs->pc[1]/wcs->pc[3])>1e-6 || fabs(wcs->pc[2]/wcs->pc[3])>1e-6 )
    error(EXIT_FAILURE, 0, "%s: HDU %s: is not aligned to the "
          "celestial coordinates. The first FITS axis should be "
          "along the Right Ascension and the second FITS axis "
          "should be along the declination.\n\n"
          "Gnuastro's Warp program can align it with the following "
          "command:\n\n"
          "    $ astwarp %s --hdu=%s --align\n",
          img->name, p->cp.hdu, img->name, p->cp.hdu);
  if(wcs->pc[0]>0)
    error(EXIT_FAILURE, 0, "%s: HDU %s: An increase in the first "
          "FITS axis pixel coordinates should be a decrease in the "
          "RA. You have to flip the image along the second axis "
          "before running Crop", img->name, p->cp.hdu);
  if(wcs->pc[3]<0)
    error(EXIT_FAILURE, 0, "%s: HDU %s: An increase in the second "
          "FITS axis pixel coordinates should translate to an "
          "increase in the declination. You have to flip the "
          "image along the first axis before running Crop",
          img->name, p->cp.hdu);


  /* Check the nature of the coordinates, currently we can only support RA
     and Dec, other modes haven't been checked. */
  if( strcmp(wcs->ctype[0], "RA---TAN")
      || strcmp(wcs->ctype[1], "DEC--TAN") )
    error(EXIT_FAILURE, 0, "currently the only WCS types usable are "
          "'RA---TAN' and 'DEC--TAN' for the first and second axises "
          "respectively. The WCS types of '%s' (hdu %s) are '%s' and '%s' "
          "respectively", img->name, p->cp.hdu, wcs->ctype[0], wcs->ctype[1]);


  /* Check if the pixels are actually a square, then compare the resolution
     with the other input images. Due to floating point errors, some very
     small differences might exist in the pixel scale, so break out with an
     error only if the pixel scales are more different than 1e-6. */
  pixscale=gal_wcs_pixel_scale(wcs);
  if(pixscale==NULL)
    error(EXIT_FAILURE, 0, "the pixel scale couldn't be deduced from the "
          "WCS");
  if( fabs(pixscale[0]-pixscale[1])/pixscale[0] > 1e-6 )
    error(EXIT_FAILURE, 0, "%s: HDU %s: The pixel scale along "
          "the two image axises is not the same. The first axis "
          "is %.15g deg/pixel, while the second is %.15g",
          img->name, p->cp.hdu, pixscale[1], pixscale[0]);
  if(p->pixscale)
    {
      for(i=0;i<ndim;++i)
        if(p->pixscale[i] != pixscale[i])
          error(EXIT_FAILURE, 0, "%s (hdu %s): has resolution of %g along "
                "dimension %d. However, previously checked input(s) had "
                "a resolution of %g in this dimension", img->name, p->cp.hdu,
                pixscale[i], i+1, p->pixscale[i]);
      free(pixscale);
    }
  else
    p->pixscale=pixscale;


  /* Set the coordinates of the dataset's corners. Note that 'dsize' is in
     C order, while pixcrd is in FITS order.*/
  switch(ndim)
    {
    case 2:
      ncorners=4;
      pixcrd[0] = 1;          pixcrd[1] = 1;
      pixcrd[2] = dsize[1];   pixcrd[3] = 1;
      pixcrd[4] = 1;          pixcrd[5] = dsize[0];
      pixcrd[6] = dsize[1];   pixcrd[7] = dsize[0];
      break;
    case 3:
      ncorners=8;
      pixcrd[0]  = 1;         pixcrd[1]  = 1;         pixcrd[2]  = 1;
      pixcrd[3]  = dsize[2];  pixcrd[4]  = 1;         pixcrd[5]  = 1;
      pixcrd[6]  = 1;         pixcrd[7]  = dsize[1];  pixcrd[8]  = 1;
      pixcrd[9]  = dsize[2];  pixcrd[10] = dsize[1];  pixcrd[11] = 1;
      pixcrd[12] = 1;         pixcrd[13] = 1;         pixcrd[14] = dsize[0];
      pixcrd[15] = dsize[2];  pixcrd[16] = 1;         pixcrd[17] = dsize[0];
      pixcrd[18] = 1;         pixcrd[19] = dsize[1];  pixcrd[20] = dsize[0];
      pixcrd[21] = dsize[2];  pixcrd[22] = dsize[1];  pixcrd[23] = dsize[0];
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: %zu dimensional datasets not supported",
            __func__, ndim);
    }


  /* Get the coordinates of the corners of the dataset in WCS.  */
  wcsp2s(wcs, ncorners, ndim, pixcrd, imgcrd, phi, theta,
         img->corners, status);


  /* Check if there was no error in the conversion. */
  for(i=0;i<ncorners;++i)
    if(status[i])
      error(EXIT_FAILURE, 0, "wcsp2s ERROR %d in row %d of pixcrd: %s",
            i, status[i], wcs_errmsg[status[i]]);


  /* Fill in the size of the dataset in WCS from the first pixel in the
     image. Note that 'dsize' is in C axises, while the 'pixscale',
     'corners' and 'sized' are in FITS axises. */
  if(ndim==2)
    {
      img->sized[0] = ( img->dsize[1] * p->pixscale[0]
                        / cos( img->corners[1] * M_PI / 180 ) );
      img->sized[1] = img->dsize[0] * p->pixscale[1];
    }
  else /* 3D */
    {
      img->sized[0] = ( img->dsize[2] * p->pixscale[0]             /* RA  */
                        / cos( img->corners[1] * M_PI / 180 ) );
      img->sized[1] = img->dsize[1] * p->pixscale[1];              /* Dec */
      img->sized[2] = img->dsize[0] * p->pixscale[2];              /* 3D  */
    }


  /* In case the image crosses the equator, we will calculate these values
     here so later on, we don't have to calculate them on every check. See
     the explanation above 'point_in_dataset'.

     Note that in both 2D and 3D data, the declination is in the second
     coordinate (index 1). */
  if( img->corners[1] * (img->corners[1]+img->sized[1]) < 0 )
    {
      /* re in the comments above 'point_in_dataset'. */
      img->equatorcorr[0]=img->corners[0]
        -0.5*img->sized[0]*(1-cos(img->corners[1]*M_PI/180));

      /* sre in the comments above 'point_in_dataset'. */
      img->equatorcorr[1]=img->sized[0]*cos(img->corners[1]*M_PI/180);
    }


  /* Just to check:
  printf("\n\n%s:\n", img->name);
  if(ndim==2)
    printf("(%.10f, %.10f)\n"
           "(%.10f, %.10f)\n"
           "(%.10f, %.10f)\n"
           "(%.10f, %.10f)\n\n",
           img->corners[0], img->corners[1],
           img->corners[2], img->corners[3],
           img->corners[4], img->corners[5],
           img->corners[6], img->corners[7]);
  else
    printf("(%.10f, %.10f, %.10f)\n"
           "(%.10f, %.10f, %.10f)\n"
           "(%.10f, %.10f, %.10f)\n"
           "(%.10f, %.10f, %.10f)\n"
           "(%.10f, %.10f, %.10f)\n"
           "(%.10f, %.10f, %.10f)\n"
           "(%.10f, %.10f, %.10f)\n"
           "(%.10f, %.10f, %.10f)\n\n",
           img->corners[0],  img->corners[1],  img->corners[2],
           img->corners[3],  img->corners[4],  img->corners[5],
           img->corners[6],  img->corners[7],  img->corners[8],
           img->corners[9],  img->corners[10], img->corners[11],
           img->corners[12], img->corners[13], img->corners[14],
           img->corners[15], img->corners[16], img->corners[17],
           img->corners[18], img->corners[19], img->corners[20],
           img->corners[21], img->corners[22], img->corners[23] );
  exit(0);
  */
}




















/*******************************************************************/
/************        Check if WCS is in image         **************/
/*******************************************************************/
/* Set the four sides around the point of interest in RA and Dec.

   NOTE: When the image is aligned with the celestial coordinates (current
   working paradigm), the declination is measured on a great circle, while
   the right ascension is not. So we have to consider this in calculating
   the difference in RA.
*/
void
wcsmode_crop_corners(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;

  size_t i, ndim=p->imgs->ndim;
  double minra=FLT_MAX, mindec=FLT_MAX;
  double maxra=-FLT_MAX, maxdec=-FLT_MAX;
  double r, d, l, dr, h[MAXDIM], hr[MAXDIM];
  size_t rmini=-1, rmaxi=-1, dmini=-1, dmaxi=-1;

  /* Set the four corners of the WCS region. */
  if(p->polygon)
    {
      /* A small sanity check. */
      if(ndim!=2)
        error(EXIT_FAILURE, 0, "%s: polygon crops are currently only "
              "supported on 2D datasets, the input dataset is %zuD",
              __func__, ndim);

      /* Find their minimum and maximum values. */
      for(i=0;i<p->nvertices;++i)
        {
          if(p->wpolygon[i*2]>maxra) maxra=p->wpolygon[i*2];
          if(p->wpolygon[i*2]<minra) minra=p->wpolygon[i*2];
          if(p->wpolygon[i*2+1]>maxdec) maxdec=p->wpolygon[i*2+1];
          if(p->wpolygon[i*2+1]<mindec) mindec=p->wpolygon[i*2+1];
        }

      /* Set the corners: */
      crp->corners[0] = maxra;  crp->corners[1] = mindec; /* Bottom Left  */
      crp->corners[2] = minra;  crp->corners[3] = mindec; /* Bottom Right */
      crp->corners[4] = maxra;  crp->corners[5] = maxdec; /* Top Left     */
      crp->corners[6] = minra;  crp->corners[7] = maxdec; /* Top Right    */
    }
  else
    {
      /* Set the RA and Dec to use as center. */
      r=crp->world[0]=p->centercoords[0][crp->out_ind];
      d=crp->world[1]=p->centercoords[1][crp->out_ind];
      if(ndim==3) l=crp->world[2]=p->centercoords[2][crp->out_ind];


      /* Calculate the declination in radians for easy readability. */
      dr=d*M_PI/180;

      /* Set the half width in each dimension. For the angular dimensions,
         also calculate it in radians. */
      hr[0] = ( h[0] = ((double *)(p->width->array))[0] / 2 ) * M_PI / 180;
      hr[1] = ( h[1] = ((double *)(p->width->array))[1] / 2 ) * M_PI / 180;
      if(ndim==3)
        h[2] = ((double *)(p->width->array))[2] / 2;

      /* Set the corners of this crop. */
      switch(ndim)
        {
        case 2:
          crp->corners[0] = r+h[0]/cos(dr-hr[1]);
          crp->corners[1] = d-h[1];                   /* Bottom left.  */

          crp->corners[2] = r-h[0]/cos(dr-hr[1]);
          crp->corners[3] = d-h[1];                   /* Bottom Right. */

          crp->corners[4] = r+h[0]/cos(dr+hr[1]);
          crp->corners[5] = d+h[1];                   /* Top Left.     */

          crp->corners[6] = r-h[0]/cos(dr+hr[1]);
          crp->corners[7] = d+h[1];                   /* Top Right.    */
          break;

        case 3:
          /* Note that the third dimension is assumed to be independent of
             the first two. So the first two coordinates of its corners in
             the front and back (on the two faces in the third dimension),
             are equal.  */
          crp->corners[0]  = crp->corners[12] = r+h[0]/cos(dr-hr[1]);
          crp->corners[1]  = crp->corners[13] = d-h[1];
          crp->corners[2]  = l-h[2];                /* Bottom left front. */

          crp->corners[3]  = crp->corners[15] = r-h[0]/cos(dr-hr[1]);
          crp->corners[4]  = crp->corners[16] = d-h[1];
          crp->corners[5]  = l-h[2];                /* Bottom right front.*/

          crp->corners[6]  = crp->corners[18] = r+h[0]/cos(dr+hr[1]);
          crp->corners[7]  = crp->corners[19] = d+h[1];
          crp->corners[8]  = l-h[2];                /* Top Left front.    */

          crp->corners[9]  = crp->corners[21] = r-h[0]/cos(dr+hr[1]);
          crp->corners[10] = crp->corners[22] = d+h[1];
          crp->corners[11] = l-h[2];                /* Top right front.   */

          crp->corners[14] = l+h[2];                /* Bottom left back.  */
          crp->corners[17] = l+h[2];                /* Bottom right back. */
          crp->corners[20] = l+h[2];                /* Top left back.     */
          crp->corners[23] = l+h[2];                /* Top right back.     */
          break;
        }
    }


  /* Set the bottom width and height of the crop in degrees. Note that the
     width changes as the height changes, so here we want the height and
     the lowest declination. Note that in 2D on the bottom edge, corners[0]
     is the maximum RA and corners[2] is the minimum RA.  For all the 2D
     region, corners[5] is one of the maximum declinations and corners[1]
     is one of the the minimum declinations.

     North and south hemispheres are no problem: When using the center,
     they are set properly (in any hemisphere) and for a polygon, the
     minimums and maximums are automatically found. */
  rmini = ndim;                 /* First element in second corner. */
  rmaxi = 0;                    /* First element.                  */
  dmini = 1;                    /* Second element.                 */
  dmaxi = ndim==2 ? 5 : 7;      /* Second element in third corner. */
  crp->sized[0]=( (crp->corners[rmaxi]-crp->corners[rmini])
                  / cos(crp->corners[dmini]*M_PI/180) );
  crp->sized[1]=crp->corners[dmaxi]-crp->corners[dmini];
  if(ndim==3)
    crp->sized[2] = crp->corners[14] - crp->corners[2];


  /* In case the crop crosses the equator, then we need these two
     corrections. See the complete explanations above 'point_in_dataset'. */
  if(crp->corners[1]*(crp->corners[1]+crp->sized[1]) < 0 )
    {
      /* re in the explanations above 'point_in_dataset'. */
      crp->equatorcorr[0]=crp->corners[0]
        -0.5*crp->sized[0]*(1-cos(crp->corners[1]*M_PI/180));

      /* sre in the explanations above 'point_in_dataset'. */
      crp->equatorcorr[1]=crp->sized[0]*cos(crp->corners[1]*M_PI/180);
    }


  /* Just to check:
  if(ndim==2)
    {
      printf("\n\n%g, %g:\n", r, d);
      printf("\t(%.10f, %.10f)\n"
             "\t(%.10f, %.10f)\n"
             "\t(%.10f, %.10f)\n"
             "\t(%.10f, %.10f)\n\n",
             crp->corners[0], crp->corners[1],
             crp->corners[2], crp->corners[3],
             crp->corners[4], crp->corners[5],
             crp->corners[6], crp->corners[7]);
    }
  else
    {
      printf("\n\n%g, %g, %g:\n", r, d, l);
      printf("\t(%.10f, %.10f, %g)\n"
             "\t(%.10f, %.10f, %g)\n"
             "\t(%.10f, %.10f, %g)\n"
             "\t(%.10f, %.10f, %g)\n"
             "\t(%.10f, %.10f, %g)\n"
             "\t(%.10f, %.10f, %g)\n"
             "\t(%.10f, %.10f, %g)\n"
             "\t(%.10f, %.10f, %g)\n\n",
             crp->corners[0],  crp->corners[1],  crp->corners[2],
             crp->corners[3],  crp->corners[4],  crp->corners[5],
             crp->corners[6],  crp->corners[7],  crp->corners[8],
             crp->corners[9],  crp->corners[10], crp->corners[11],
             crp->corners[12], crp->corners[13], crp->corners[14],
             crp->corners[15], crp->corners[16], crp->corners[17],
             crp->corners[18], crp->corners[19], crp->corners[20],
             crp->corners[21], crp->corners[22], crp->corners[23] );
    }
  exit(0);
  */
}





/* We have the polygon vertices in WCS coordinates and need to change them
   to one input image's pixel coordinates. */
void
fillcrpipolygon(struct onecropparams *crp)
{
  struct cropparams *p=crp->p;
  gal_data_t *tmp, *coords=NULL;
  size_t i, d, ndim=p->imgs->ndim;

  /* Allocate the necessary arrays for each column. */
  for(d=0;d<ndim;++d)
    gal_list_data_add_alloc(&coords, NULL, GAL_TYPE_FLOAT64, 1, &p->nvertices,
                            NULL, 0, -1, 1, NULL, NULL, NULL);


  /* Fill in the world coordinate columns. */
  for(i=0;i<p->nvertices;++i)
    {
      d=0;
      for(tmp=coords;tmp!=NULL;tmp=tmp->next)
        ((double *)(tmp->array))[i] = p->wpolygon[ i * ndim + d++ ];
    }


  /* Convert them to image coordinates. */
  gal_wcs_world_to_img(coords, p->imgs[crp->in_ind].wcs, 1);


  /* Allocate the image polygon array, and put the image polygon vertice
     values into it. */
  crp->ipolygon=gal_pointer_allocate(GAL_TYPE_FLOAT64, ndim*p->nvertices, 0,
                                     __func__, "crp->ipolygon");
  for(i=0;i<p->nvertices;++i)
    {
      d=0;
      for(tmp=coords;tmp!=NULL;tmp=tmp->next)
        crp->ipolygon[ i * ndim + d++ ] = ((double *)(tmp->array))[i];
    }


  /* Clean up. */
  gal_list_data_free(coords);
}





/* BASICS:
   =======

   An image is a rectangle, but the sky is on a globe. When the images
   are aligned to the celstial coordinates, (as we have required here
   in wcscheckprepare) the first FITS axis shows change in RA, while
   the second axis shows change in Dec. The declination always changes
   along a great circle, so there is no problem. But unless
   declination is constrained to zero, RA changes on small circles.

   See the rectangle below, assume it is an image. To check if a given
   point is within the same declination as this rectangle is very
   simple, since d3==d4 and d1==d2. Note that r1>r2 and r3>r4 (because
   right ascension increases to the east).

       (r3,d3)    ------------------ (r4,d4)
                  |                |
                  |                |
                  |                |
                  |                |
       (r1,d1)    ------------------ (r2,d2)

   But for RA, the same number of pixels on each declination,
   corresponds to different ranges in right ascention. As the
   declination gets higher in the northern hemisphere (where the
   declination rises towards the top of the image) r1-r2 becomes
   smaller than r3-r4. So, in terms of coverage in RA and Dec, this
   box should rather be shown like this trapezoid (exaggerated):

                 --------------------
                 |                  |
                  |                |         (Northern hemisphere)
                   |              |
                    |            |
                    --------------

   On the southern hemisphere it should be shown like this:

                   ----------------
                   |              |
                  |                |         (Southern hemisphere)
                 |                  |
                |                    |
                ----------------------

   The functional form of the change is the inverse cosine, so:

           (r3-r4)=(r1-r2)/cos(d3-d1)            (North)

           (r1-r2)=(r3-r4)/cos(d1-d3)            (South)

   QUESTION:
   ========
   Is a given point at the RA and Dec of (rp,dp) inside our
   rectangular image?



   IMAGE IS FULLY WITHIN ONE HEMISPHERE
   ------------------------------------

   Our reference point for the image is the first pixel in the image,
   which by convention is the (r1,d1) point in the rectangle above. We
   also have the angular size of the rectangular image as 'sr', 'sd'
   (for "size in RA" and "size in dec"). NOTE: This has to be We also
   assume r1+sr and d1+sd are the distances to the last pixels in our
   rectangular image.

   As explained above, to check the declination range, everything is
   very easy:
                                        dp>=d1     &&      dp<=d1+sd

   For RA, things become a little more complicated (recall that
   r1>r3). 'n' is defined as half of the extra space between the top
   and bottom lines of the two trapezoids.

   (North) n=0.5*sr*(1/cos(dp-d1)-1) ==> rp<=r1+n   &&   rp>=r1-sr-n

   (South) n=0.5*sr*(1-cos(dp-d1))   ==> rp<=r1-n   &&   rp>=r1-sr+n



   IMAGE CROSSES THE EQUATOR
   -------------------------

   When d1*(d1+sd)<0, the image crosses the equator (d1 is negative
   and d1+sd is positive). In this case, we define 're' and 'sre' as
   an equivalent of r1 and sr but on the equator:

       re=r1-0.5*sr*(1-cos(d1))   &&   sre=sr*cos(d1)

   then  for all  the  points  with negative  declination  we use  the
   (South) equations of above as before and for those points that have
   a positive declination,  we use the North formula  but replacing r1
   with re, d1 with 0 and sr with sre.


   INPUTS
   ------
   p[]: Point coordinates (rp and dp above).
   i[]: Coordinates of first pixel in image. (r1 and d1 above).
   s[]: Size/width of box (sr and sd above).
   c[]: Corrections if equator is passed, (se and sre above).


   IMPORTANT: It is assumed that the dimensions are ordered with:

      0: RA
      1: Dec
      2: Third dimension (independent of RA and Dec).          */
static int
point_in_dataset(double *p, double *i, double *s, double *c, size_t ndim)
{
  double n;

  /* If there is a third dimension, then first check that. Note that the
     third dimension is assumed to be indendent of the first two. */
  if(ndim==3 && ( p[2]<i[2] || p[2]>i[2]+s[2] ) )
    return 0;

  /* In the RA and Dec checks, first check the declination. If it is not in
     range, you can safely return 0. */
  if(p[1]>=i[1] && p[1]<=i[1]+s[1])
    {
      if(p[1]<=0)            /* Point is in southern hemisphere, it     */
        {                    /* doesn't matter if image passes equator! */
          n=0.5f*s[0]*( 1 - cos((p[1]-i[1])*M_PI/180) );
          if(p[0]<=i[0]-n   &&   p[0]>=i[0]-s[0]+n)
            return 1;
        }
      else                      /* Point is in the northern hemisphere. */
        {
          if( i[1] * (s[1]+i[1]) > 0 ) /* Image does not cross equator. */
            {
              n=0.5f*s[0]*( 1/cos((p[1]-i[1])*M_PI/180) - 1);
              if(p[0]<=i[0]+n   &&   p[0]>=i[0]-s[0]-n)
                return 1;
            }
          else                             /* Image crosses the equator.*/
            {
              n=0.5f*c[1]*( 1/cos((p[1]     )*M_PI/180) - 1);
              if(p[0]<=c[0]+n   &&   p[0]>=c[0]-c[1]-n)
                return 1;
            }
        }
    }

  return 0;
}






/* Is there an overlap between this crop box and the survey image?
   This function will return 0 if there isn't and 1 if there is.

   We don't want to necessarily assume that the crop box is smaller
   than the survey images. If we made that assumption, we only had to
   check if the corners of the crop are in the image. When we allow
   the input survey images to be smaller than the crop box, it might
   happen that none of the corners of the crop are in the image but
   there is an overlap (the survey image is completely within the crop
   box). So we have to check both.  */
int
wcsmode_overlap(struct onecropparams *crp)
{
  double *d, *fd;
  double *i, *s, *c;                /* for clear viewing. */
  struct cropparams *p=crp->p;
  size_t ndim=crp->p->imgs->ndim;

  /* First check if the corners of the crop are in the image.*/
  s=p->imgs[crp->in_ind].sized;
  i=p->imgs[crp->in_ind].corners;
  c=p->imgs[crp->in_ind].equatorcorr;
  fd=(d=crp->corners) + (ndim==2 ? 8 : 24);
  do
    {
      /* As long as one of the crop corners are in the image, we know there
         is overlap and can return a true value. We don't need to check all
         corners. */
      if( point_in_dataset(d, i, s, c, ndim) ) return 1;
      d+=ndim;
    }
  while(d<fd);

  /* None of the crop box corners where within the image. Now check if
     the image corners are within the crop.*/
  s=crp->sized;
  i=crp->corners;
  c=crp->equatorcorr;
  fd=(d=p->imgs[crp->in_ind].corners) + (ndim==2 ? 8 : 24);
  do
    {
      if( point_in_dataset(d, i, s, c, ndim) ) return 1;
      d+=ndim;
    }
  while(d<fd);

  /* If control reaches here, there was no overlap. */
  return 0;
}
