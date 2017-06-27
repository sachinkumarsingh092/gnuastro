/*********************************************************************
Functions to that only use WCSLIB functionality.
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

#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include <gsl/gsl_linalg.h>

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/dimension.h>








/*************************************************************
 ***********               Read WCS                ***********
 *************************************************************/
/* Read the WCS information from the header. Unfortunately, WCS lib is
   not thread safe, so it needs a mutex. In case you are not using
   multiple threads, just pass a NULL pointer as the mutex.

   After you finish with this WCS, you should free the space with:

   status = wcsvfree(&nwcs,&wcs);

   If the WCS structure is not recognized, then this function will
   return a NULL pointer for the wcsprm structure and a zero for
   nwcs. It will also report the fact to the user in stderr.

   ===================================
   WARNING: wcspih IS NOT THREAD SAFE!
   ===================================
   Don't call this function within a thread or use a mutex.
*/
struct wcsprm *
gal_wcs_read_fitsptr(fitsfile *fptr, size_t hstartwcs, size_t hendwcs,
                     int *nwcs)
{
  /* Declaratins: */
  int nkeys=0, status=0;
  struct wcsprm *wcs=NULL;
  char *fullheader, *to, *from;
  int relax    = WCSHDR_all; /* Macro: use all informal WCS extensions. */
  int ctrl     = 0;          /* Don't report why a keyword wasn't used. */
  int nreject  = 0;          /* Number of keywords rejected for syntax. */

  /* CFITSIO function: */
  if( fits_hdr2str(fptr, 1, NULL, 0, &fullheader, &nkeys, &status) )
    gal_fits_io_error(status, NULL);

  /* Only consider the header keywords in the current range: */
  if(hendwcs>hstartwcs)
    {
      /* Mark the last character in the desired region. */
      fullheader[hendwcs*(FLEN_CARD-1)]='\0';
      /*******************************************************/
      /******************************************************
      printf("%s\n", fullheader);
      ******************************************************/
      /*******************************************************/

      /* Shift all the characters to the start of the string. */
      if(hstartwcs)                /* hstartwcs!=0 */
        {
          to=fullheader;
          from=&fullheader[hstartwcs*(FLEN_CARD-1)-1];
          while(*from++!='\0') *to++=*from;
        }

      nkeys=hendwcs-hstartwcs;

      /*******************************************************/
      /******************************************************
      printf("\n\n\n###############\n\n\n\n\n\n");
      printf("%s\n", &fullheader[1*(FLEN_CARD-1)]);
      exit(0);
      ******************************************************/
      /*******************************************************/
    }


  /* WCSlib function to parse the FITS headers. */
  status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, nwcs, &wcs);
  if(status)
    {
      fprintf(stderr, "\n##################\n"
              "WCSLIB Warning: wcspih ERROR %d: %s.\n"
              "##################\n",
              status, wcs_errmsg[status]);
      wcs=NULL; *nwcs=0;
    }
  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in fitsarrayvv.c for freeing "
                           "the memory used to keep all the headers");


  /* Set the internal structure: */
  if(wcs)
    {
      /* CTYPE is a mandatory WCS keyword, so if it hasn't been given (its
         '\0'), then the headers didn't have a WCS structure. However,
         WCSLIB still fills in the basic information (for example the
         dimensionality of the dataset). */
      if(wcs->ctype[0][0]=='\0')
        {
          wcsfree(wcs);
          wcs=NULL;
          *nwcs=0;
        }
      else
        {
          /* Set the WCS structure. */
          status=wcsset(wcs);
          if(status)
            {
              fprintf(stderr, "\n##################\n"
                      "WCSLIB Warning: wcsset ERROR %d: %s.\n"
                      "##################\n",
                      status, wcs_errmsg[status]);
              wcsfree(wcs);
              wcs=NULL;
              *nwcs=0;
            }
          else
            /* A correctly useful WCS is present. When no PC matrix
               elements were present in the header, the default PC matrix
               (a unity matrix) is used. In this case WCSLIB doesn't set
               `altlin' (and gives it a value of 0). In Gnuastro, later on,
               we might need to know the type of the matrix used, so in
               such a case, we will set `altlin' to 1. */
            if(wcs->altlin==0) wcs->altlin=1;
        }
    }


  /* Return the WCS structure. */
  return wcs;
}





struct wcsprm *
gal_wcs_read(char *filename, char *hdu, size_t hstartwcs,
             size_t hendwcs, int *nwcs)
{
  int status=0;
  fitsfile *fptr;
  struct wcsprm *wcs;

  /* Check HDU for realistic conditions: */
  fptr=gal_fits_hdu_open_format(filename, hdu, 0);

  /* Read the WCS information: */
  wcs=gal_wcs_read_fitsptr(fptr, hstartwcs, hendwcs, nwcs);

  /* Close the FITS file and return. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
  return wcs;
}




















/**************************************************************/
/**********              Utilities                 ************/
/**************************************************************/
/* Copy a given WSC structure into another one. */
struct wcsprm *
gal_wcs_copy(struct wcsprm *wcs)
{
  struct wcsprm *out;

  /* If the input WCS is NULL, return a NULL WCS. */
  if(wcs)
    {
      /* Allocate the output WCS structure. */
      errno=0;
      out=malloc(sizeof *out);
      if(out==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `out'",
              __func__, sizeof *out);

      /* Initialize the allocated WCS structure. The WCSLIB manual says "On
         the first invokation, and only the first invokation, wcsprm::flag
         must be set to -1 to initialize memory management"*/
      out->flag=-1;
      wcsini(1, wcs->naxis, out);

      /* Copy the input WCS to the output WSC structure. */
      wcscopy(1, wcs, out);
    }
  else
    out=NULL;

  /* Return the final output. */
  return out;
}





/* Using the block data structure of the tile, add a WCS structure for
   it. In many cases, tiles are created for internal processing, so there
   is no need to keep their WCS. Hence for preformance reasons, when
   creating the tiles they don't have any WCS structure. When needed, this
   function can be used to add a WCS structure to the tile by copying the
   WCS structure of its block and correcting its starting points. If the
   tile already has a WCS structure, this function won't do anything.*/
void
gal_wcs_on_tile(gal_data_t *tile)
{
  size_t i, start_ind, ndim=tile->ndim;
  gal_data_t *block=gal_tile_block(tile);
  size_t *coord=gal_data_malloc_array(GAL_TYPE_SIZE_T, ndim, __func__,
                                      "coord");

  /* If the tile already has a WCS structure, don't do anything. */
  if(tile->wcs) return;
  else
    {
      /* Copy the block's WCS into the tile. */
      tile->wcs=gal_wcs_copy(block->wcs);

      /* Find the coordinates of the tile's starting index. */
      start_ind=gal_data_ptr_dist(block->array, tile->array, block->type);
      gal_dimension_index_to_coord(start_ind, ndim, block->dsize, coord);

      /* Correct the copied WCS structure. Note that crpix is indexed in
         the FITS/Fortran order while coord is ordered in C, it also starts
         counting from 1, not zero. */
      for(i=0;i<ndim;++i)
        tile->wcs->crpix[i] -= coord[ndim-1-i];
      /*
      printf("start_ind: %zu\n", start_ind);
      printf("coord: %zu, %zu\n", coord[1]+1, coord[0]+1);
      printf("CRPIX: %f, %f\n", tile->wcs->crpix[0], tile->wcs->crpix[1]);
      */
    }

  /* Clean up. */
  free(coord);
}





/* Return the Warping matrix of the given WCS structure. This will be the
   final matrix irrespective of the type of storage in the WCS
   structure. Recall that the FITS standard has several methods to store
   the matrix, which is up to this function to account for and return the
   final matrix. The output is an allocated DxD matrix where `D' is the
   number of dimensions. */
double *
gal_wcs_warp_matrix(struct wcsprm *wcs)
{
  double *out;
  size_t i, j, size=wcs->naxis*wcs->naxis;

  /* Allocate the necessary array. */
  errno=0;
  out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `out'",
          __func__, size*sizeof *out);

  /* Fill in the array. */
  if(wcs->altlin & 0x1)          /* Has a PCi_j array. */
    {
      for(i=0;i<wcs->naxis;++i)
        for(j=0;j<wcs->naxis;++j)
          out[i*wcs->naxis+j] = wcs->cdelt[i] * wcs->pc[i*wcs->naxis+j];
    }
  else if(wcs->altlin & 0x2)     /* Has CDi_j array.   */
    {
      for(i=0;i<size;++i)
        out[i]=wcs->cd[i];
    }
  else
    error(EXIT_FAILURE, 0, "%s: currently only PCi_ja and CDi_ja keywords "
          "are recognized", __func__);

  /* Return the result */
  return out;
}





/* According to the FITS standard, in the `PCi_j' WCS formalism, the matrix
   elements m_{ij} are encoded in the `PCi_j' keywords and the scale
   factors are encoded in the `CDELTi' keywords. There is also another
   formalism (the `CDi_j' formalism) which merges the two into one
   matrix.

   However, WCSLIB's internal operations are apparently done in the `PCi_j'
   formalism. So its outputs are also all in that format by default. When
   the input is a `CDi_j', WCSLIB will still read the image into the
   `PCi_j' formalism and the `CDELTi's are set to 1. This function will
   decompose the two matrices to give a reasonable `CDELTi' and `PCi_j' in
   such cases. */
void
gal_wcs_decompose_pc_cdelt(struct wcsprm *wcs)
{
  double *ps;
  size_t i, j;

  /* The correction is only needed when the PC matrix is filled. */
  if(wcs->pc)
    {
      /* Get the pixel scale. */
      ps=gal_wcs_pixel_scale(wcs);

      /* The PC matrix and the CDELT elements might both contain scale
         elements (during processing the scalings might be added only to PC
         elements for example). So to be safe, we first multiply them into
         one matrix. */
      for(i=0;i<wcs->naxis;++i)
        for(j=0;j<wcs->naxis;++j)
          wcs->pc[i*wcs->naxis+j] *= wcs->cdelt[i];

      /* Set the CDELTs. */
      for(i=0; i<wcs->naxis; ++i)
        wcs->cdelt[i] = ps[i];

      /* Correct the PCi_js */
      for(i=0;i<wcs->naxis;++i)
        for(j=0;j<wcs->naxis;++j)
          wcs->pc[i*wcs->naxis+j] /= ps[i];

      /* Clean up. */
      free(ps);

      /* According to the `wcslib/wcs.h' header: "In particular, wcsset()
         resets wcsprm::cdelt to unity if CDi_ja is present (and no
         PCi_ja).". So apparently, when the input is a `CDi_j', it might
         expect the `CDELTi' elements to be 1.0. But we have changed that
         here, so we will correct the `altlin' element of the WCS structure
         to make sure that WCSLIB only looks into the `PCi_j' and `CDELTi'
         and makes no assumptioins about `CDELTi'. */
      wcs->altlin=1;
    }
}





/* The distance (along a great circle) on a sphere between two points
   is calculated here. Since the pixel sides are usually very small,
   we won't be using the direct formula:

   cos(distance)=sin(d1)*sin(d2)+cos(d1)*cos(d2)*cos(r1-r2)

   We will be using the haversine formula which better considering
   floating point errors (from Wikipedia:)

   sin^2(distance)/2=sin^2( (d1-d2)/2 )+cos(d1)*cos(d2)*sin^2( (r1-r2)/2 )

   Inputs and outputs are all in degrees.
*/
double
gal_wcs_angular_distance_deg(double r1, double d1, double r2, double d2)
{
  /* Convert degrees to radians. */
  double r1r=r1*M_PI/180, d1r=d1*M_PI/180;
  double r2r=r2*M_PI/180, d2r=d2*M_PI/180;

  /* To make things easier to read: */
  double a=sin( (d1r-d2r)/2 );
  double b=sin( (r1r-r2r)/2 );

  /* Return the result: */
  return 2*asin( sqrt( a*a + cos(d1r)*cos(d2r)*b*b) ) * 180/M_PI;
}




/* Return the pixel scale of the dataset in units of the WCS. */
double *
gal_wcs_pixel_scale(struct wcsprm *wcs)
{
  gsl_vector S;
  gsl_matrix A, V;
  size_t n=wcs->naxis;
  double *a, *v, *pixscale;

  /* Allocate space for the `v' and `pixscale' arrays. */
  errno=0; v=malloc(n*n*sizeof *v);
  if(v==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `v'",
          __func__, n*n*sizeof *v);
  errno=0; pixscale=malloc(n*sizeof *pixscale);
  if(pixscale==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `pixscale'",
          __func__, n*sizeof *pixscale);

  /* Write the full matrix into an array, irrespective of what type it was
     stored in the wcsprm structure (`PCi_j' style or `CDi_j' style). */
  a=gal_wcs_warp_matrix(wcs);

  /* Fill in the necessary GSL vector and Matrix structures. */
  S.size=n;     S.stride=1;                S.data=pixscale;
  V.size1=n;    V.size2=n;    V.tda=n;     V.data=v;
  A.size1=n;    A.size2=n;    A.tda=n;     A.data=a;

  /* Run GSL's Singular Value Decomposition, using one-sided Jacobi
     orthogonalization which computes the singular (scale) values to a
     higher relative accuracy.*/
  gsl_linalg_SV_decomp_jacobi(&A, &V, &S);

  /* Clean up and return. */
  free(a);
  free(v);
  return pixscale;
}





/* Report the arcsec^2 area of the pixels in the image based on the
   WCS information in that image. */
double
gal_wcs_pixel_area_arcsec2(struct wcsprm *wcs)
{
  double out;
  double *pixscale;

  /* A small sanity check. Later, when higher dimensions are necessary, we
     can find which ones correlate to RA and Dec and use them to find the
     pixel area in arcsec^2. */
  if(wcs->naxis!=2)
    error(EXIT_FAILURE, 0, "%s: currently only 2D datasets supported. "
          "The input WCS has %d dimensions", __func__, wcs->naxis);

  /* Get the pixel scales along each axis in degrees, then multiply. */
  pixscale=gal_wcs_pixel_scale(wcs);

  /* Clean up and return the result. */
  out = pixscale[0] * pixscale[1] * 3600.0f * 3600.0f;
  free(pixscale);
  return out;
}



















/**************************************************************/
/**********            Array conversion            ************/
/**************************************************************/
/* Convert an array of world coordinates to image coordinates. Note that in
   Gnuastro, each column is treated independently, so the inputs are
   separate. If `*x==NULL', or `*y==NULL', then space will be allocated for
   them, otherwise, it is assumed that space has already been
   allocated. Note that they must each be a 1 dimensional array.

   You can do the conversion in place: just pass the same array as you give
   to RA and Dec to X and Y. */
void
gal_wcs_world_to_img(struct wcsprm *wcs, double *ra, double *dec,
                     double **x, double **y, size_t size)
{
  size_t i;
  int status, *stat, ncoord=size, nelem=2;
  double *phi, *theta, *world, *pixcrd, *imgcrd;

  /* Allocate all the necessary arrays. */
  phi    = gal_data_malloc_array( GAL_TYPE_FLOAT64, size, __func__, "phi");
  stat   = gal_data_calloc_array( GAL_TYPE_INT32,   size, __func__, "stat");
  theta  = gal_data_malloc_array( GAL_TYPE_FLOAT64, size, __func__, "theta");
  world  = gal_data_malloc_array( GAL_TYPE_FLOAT64, 2*size, __func__,
                                  "world");
  imgcrd = gal_data_malloc_array( GAL_TYPE_FLOAT64, 2*size, __func__,
                                  "imgcrd");
  pixcrd = gal_data_malloc_array( GAL_TYPE_FLOAT64, 2*size, __func__,
                                  "pixcrd");

  /* Write the values into the allocated contiguous array. */
  for(i=0;i<size;++i) { world[i*2]=ra[i]; world[i*2+1]=dec[i]; }

  /* Use WCSLIB's `wcss2p'. */
  status=wcss2p(wcs, ncoord, nelem, world, phi, theta, imgcrd, pixcrd, stat);
  if(status)
    error(EXIT_FAILURE, 0, "%s: wcss2p ERROR %d: %s", __func__, status,
          wcs_errmsg[status]);

  /* For a sanity check:
  printf("\n\ngal_wcs_world_to_img sanity check:\n");
  for(i=0;i<size;++i)
    printf("world (%f, %f) --> pix (%f, %f), [stat: %d]\n",
           world[i*2], world[i*2+1], pixcrd[i*2], pixcrd[i*2+1], stat[i]);
  */

  /* Allocate the output arrays if they were not already allocated. */
  if(*x==NULL) *x=gal_data_malloc_array(GAL_TYPE_FLOAT64, size, __func__,"x");
  if(*y==NULL) *y=gal_data_malloc_array(GAL_TYPE_FLOAT64, size, __func__,"y");

  /* Put the values into the output arrays. */
  for(i=0;i<size;++i)
    {
      (*x)[i] = stat[i] ? NAN : pixcrd[i*2];
      (*y)[i] = stat[i] ? NAN : pixcrd[i*2+1];
    }

  /* Clean up. */
  free(phi);
  free(stat);
  free(theta);
  free(world);
  free(pixcrd);
}





/* Similar to `gal_wcs_world_to_img' but converts image coordinates into
   world coordinates. */
void
gal_wcs_img_to_world(struct wcsprm *wcs, double *x, double *y,
                     double **ra, double **dec, size_t size)
{
  size_t i;
  int status, *stat, ncoord=size, nelem=2;
  double *phi, *theta, *world, *pixcrd, *imgcrd;

  /* Allocate all the necessary arrays. */
  phi    = gal_data_malloc_array( GAL_TYPE_FLOAT64, size, __func__, "phi");
  stat   = gal_data_calloc_array( GAL_TYPE_INT32,   size, __func__, "stat");
  theta  = gal_data_malloc_array( GAL_TYPE_FLOAT64, size, __func__, "theta");
  world  = gal_data_malloc_array( GAL_TYPE_FLOAT64, 2*size, __func__,
                                  "world");
  imgcrd = gal_data_malloc_array( GAL_TYPE_FLOAT64, 2*size, __func__,
                                  "imgcrd");
  pixcrd = gal_data_malloc_array( GAL_TYPE_FLOAT64, 2*size, __func__,
                                  "pixcrd");

  /* Write the values into the allocated contiguous array. */
  for(i=0;i<size;++i) { pixcrd[i*2]=x[i]; pixcrd[i*2+1]=y[i]; }

  /* Use WCSLIB's wcsp2s for the conversion. */
  status=wcsp2s(wcs, ncoord, nelem, pixcrd, imgcrd, phi, theta, world, stat);
  if(status)
    error(EXIT_FAILURE, 0, "%s: wcsp2s ERROR %d: %s", __func__, status,
          wcs_errmsg[status]);

  /* For a sanity check:
  printf("\n\ngal_wcs_img_to_world sanity check:\n");
  for(i=0;i<size;++i)
    printf("img (%f, %f) --> world (%f, %f), [stat: %d]\n",
           pixcrd[i*2], pixcrd[i*2+1], world[i*2], world[i*2+1], stat[i]);
  */

  /* Allocate the output arrays if they were not already allocated. */
  if(*ra==NULL)
    *ra  = gal_data_malloc_array(GAL_TYPE_FLOAT64, size, __func__, "ra");
  if(*dec==NULL)
    *dec = gal_data_malloc_array(GAL_TYPE_FLOAT64, size, __func__, "dec");

  /* Put the values into the output arrays. */
  for(i=0;i<size;++i)
    {
      (*ra)[i]  = stat[i] ? NAN : world[i*2];
      (*dec)[i] = stat[i] ? NAN : world[i*2+1];
    }

  /* Clean up. */
  free(phi);
  free(stat);
  free(theta);
  free(world);
  free(pixcrd);
}
