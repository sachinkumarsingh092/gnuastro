/*********************************************************************
Functions to that only use WCSLIB functionality.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/permutation.h>








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
  int sumcheck;
  size_t i, fulllen;
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

  /* Set the internal structure: */
  if(wcs)
    {
      /* It may happen that the WCS-related keyword values are stored as
         strings (they have single-quotes around them). In this case,
         WCSLIB will read the CRPIX and CRVAL values as zero. When this
         happens do a small check and abort, while informing the user about
         the problem. */
      sumcheck=0;
      for(i=0;i<wcs->naxis;++i)
        {sumcheck += (wcs->crval[i]==0.0f) + (wcs->crpix[i]==0.0f);}
      if(sumcheck==wcs->naxis*2)
        {
          /* We only care about the first set of characters in each
             80-character row, so we don't need to parse the last few
             characters anyway. */
          fulllen=strlen(fullheader)-12;
          for(i=0;i<fulllen;++i)
            if( strncmp(fullheader+i, "CRVAL1  = '", 11) == 0 )
              fprintf(stderr, "WARNING: WCS Keyword values are not "
                      "numbers.\n\n"
                      "WARNING: The values to the WCS-related keywords are "
                      "enclosed in single-quotes. In the FITS standard "
                      "this is how string values are stored, therefore "
                      "WCSLIB is unable to read them AND WILL PUT ZERO IN "
                      "THEIR PLACE (creating a wrong WCS in the output). "
                      "Please update the respective keywords of the input "
                      "to be numbers (see next line).\n\n"
                      "WARNING: You can do this with Gnuastro's 'astfits' "
                      "program and the '--update' option. The minimal WCS "
                      "keywords that need a numerical value are: 'CRVAL1', "
                      "'CRVAL2', 'CRPIX1', 'CRPIX2', 'EQUINOX' and "
                      "'CD%%_%%' (or 'PC%%_%%', where the %% are integers), "
                      "please see the FITS standard, and inspect your FITS "
                      "file to identify the full set of keywords that you "
                      "need correct (for example PV%%_%% keywords).\n\n");
        }

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
          /* For a check.
          printf("flag: %d\n", wcs->flag);
          printf("naxis: %d\n", wcs->naxis);
          printf("crpix: %f, %f\n", wcs->crpix[0], wcs->crpix[1]);
          printf("pc: %f, %f, %f, %f\n", wcs->pc[0], wcs->pc[1], wcs->pc[2],
                 wcs->pc[3]);
          printf("cdelt: %f, %f\n", wcs->cdelt[0], wcs->cdelt[1]);
          printf("crval: %f, %f\n", wcs->crval[0], wcs->crval[1]);
          printf("cunit: %s, %s\n", wcs->cunit[0], wcs->cunit[1]);
          printf("ctype: %s, %s\n", wcs->ctype[0], wcs->ctype[1]);
          printf("lonpole: %f\n", wcs->lonpole);
          printf("latpole: %f\n", wcs->latpole);
          */

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
               'altlin' (and gives it a value of 0). In Gnuastro, later on,
               we might need to know the type of the matrix used, so in
               such a case, we will set 'altlin' to 1. */
            if(wcs->altlin==0) wcs->altlin=1;
        }
    }

  /* Clean up and return. */
  status=0;
  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in freeing the memory used to "
                      "keep all the headers");
  return wcs;
}





struct wcsprm *
gal_wcs_read(char *filename, char *hdu, size_t hstartwcs,
             size_t hendwcs, int *nwcs)
{
  int status=0;
  fitsfile *fptr;
  struct wcsprm *wcs;

  /* Make sure we are dealing with a FITS file. */
  if( gal_fits_name_is_fits(filename) == 0 )
    return NULL;

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
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'out'",
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





/* Remove the algorithm part of CTYPE (anything after, and including, a
   '-') if necessary. */
static void
wcs_ctype_noalgorithm(char *str)
{
  size_t i, len=strlen(str);
  for(i=0;i<len;++i) if(str[i]=='-') { str[i]='\0'; break; }
}




/* See if the CTYPE string ends with TAN. */
static int
wcs_ctype_has_tan(char *str)
{
  size_t len=strlen(str);

  return !strcmp(&str[len-3], "TAN");
}





/* Remove dimension. */
#define WCS_REMOVE_DIM_CHECK 0
void
gal_wcs_remove_dimension(struct wcsprm *wcs, size_t fitsdim)
{
  size_t c, i, j, naxis;

  /* If the WCS structure is NULL, just return. */
  if(wcs==NULL) return;

  /* Sanity check. */
  naxis=wcs->naxis;
  if(fitsdim==0 || fitsdim>naxis)
    error(EXIT_FAILURE, 0, "%s: requested dimension (fitsdim=%zu) must be "
          "larger than zero and smaller than the number of dimensions in "
          "the given WCS structure (%zu)", __func__, fitsdim, naxis);

  /**************************************************/
#if WCS_REMOVE_DIM_CHECK
  printf("\n\nfitsdim: %zu\n", fitsdim);
  printf("\n##################\n");
  /*
  wcs->pc[0]=0;   wcs->pc[1]=1;   wcs->pc[2]=2;
  wcs->pc[3]=3;   wcs->pc[4]=4;   wcs->pc[5]=5;
  wcs->pc[6]=6;   wcs->pc[7]=7;   wcs->pc[8]=8;
  */
  for(i=0;i<wcs->naxis;++i)
    {
      for(j=0;j<wcs->naxis;++j)
        printf("%-5g", wcs->pc[i*wcs->naxis+j]);
      printf("\n");
    }
#endif
  /**************************************************/


  /* First loop over the arrays. */
  for(i=0;i<naxis;++i)
    {
      /* The dimensions are in FITS order, but counting starts from 0, so
         we'll have to subtract 1 from 'fitsdim'. */
      if(i>fitsdim-1)
        {
          /* 1-D arrays. */
          if(wcs->crpix) wcs->crpix[i-1] = wcs->crpix[i];
          if(wcs->cdelt) wcs->cdelt[i-1] = wcs->cdelt[i];
          if(wcs->crval) wcs->crval[i-1] = wcs->crval[i];
          if(wcs->crota) wcs->crota[i-1] = wcs->crota[i];
          if(wcs->crder) wcs->crder[i-1] = wcs->crder[i];
          if(wcs->csyer) wcs->csyer[i-1] = wcs->csyer[i];

          /* The strings are all statically allocated, so we don't need to
             check. */
          memcpy(wcs->cunit[i-1], wcs->cunit[i], 72);
          memcpy(wcs->ctype[i-1], wcs->ctype[i], 72);
          memcpy(wcs->cname[i-1], wcs->cname[i], 72);

          /* For 2-D arrays, just bring up all the rows. We'll fix the
             columns in a second loop. */
          for(j=0;j<naxis;++j)
            {
              if(wcs->pc) wcs->pc[ (i-1)*naxis+j ] = wcs->pc[ i*naxis+j ];
              if(wcs->cd) wcs->cd[ (i-1)*naxis+j ] = wcs->cd[ i*naxis+j ];
            }
        }
    }


  /**************************************************/
#if WCS_REMOVE_DIM_CHECK
  printf("\n###### Respective row removed (replaced).\n");
  for(i=0;i<wcs->naxis;++i)
    {
      for(j=0;j<wcs->naxis;++j)
        printf("%-5g", wcs->pc[i*wcs->naxis+j]);
      printf("\n");
    }
#endif
  /**************************************************/


  /* Second loop for 2D arrays. */
  c=0;
  for(i=0;i<naxis;++i)
    for(j=0;j<naxis;++j)
      if(j!=fitsdim-1)
        {
          if(wcs->pc) wcs->pc[ c ] = wcs->pc[ i*naxis+j ];
          if(wcs->cd) wcs->cd[ c ] = wcs->cd[ i*naxis+j ];
          ++c;
        }


  /* Correct the total number of dimensions in the WCS structure. */
  naxis = wcs->naxis -= 1;


  /* The 'TAN' algorithm needs two dimensions. So we need to remove it when
     it can cause confusion. */
  switch(naxis)
    {
    /* The 'TAN' algorithm cannot be used for any single-dimensional
       dataset. So we'll have to remove it if it exists. */
    case 1:
      wcs_ctype_noalgorithm(wcs->ctype[0]);
      break;

    /* For any other dimensionality, 'TAN' should be kept only when exactly
       two dimensions have it. */
    default:

      c=0;
      for(i=0;i<naxis;++i)
        if( wcs_ctype_has_tan(wcs->ctype[i]) )
          ++c;

      if(c!=2)
        for(i=0;i<naxis;++i)
          if( wcs_ctype_has_tan(wcs->ctype[i]) )
            wcs_ctype_noalgorithm(wcs->ctype[i]);
      break;
    }



  /**************************************************/
#if WCS_REMOVE_DIM_CHECK
  printf("\n###### Respective column removed.\n");
  for(i=0;i<naxis;++i)
    {
      for(j=0;j<naxis;++j)
        printf("%-5g", wcs->pc[i*naxis+j]);
      printf("\n");
    }
  printf("\n###### One final string\n");
  for(i=0;i<naxis;++i)
    printf("%s\n", wcs->ctype[i]);
  exit(0);
#endif
  /**************************************************/
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
  size_t *coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                     "coord");

  /* If the tile already has a WCS structure, don't do anything. */
  if(tile->wcs) return;
  else
    {
      /* Copy the block's WCS into the tile. */
      tile->wcs=gal_wcs_copy(block->wcs);

      /* Find the coordinates of the tile's starting index. */
      start_ind=gal_pointer_num_between(block->array, tile->array,
                                        block->type);
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
   final matrix. The output is an allocated DxD matrix where 'D' is the
   number of dimensions. */
double *
gal_wcs_warp_matrix(struct wcsprm *wcs)
{
  double *out, crota2;
  size_t i, j, size=wcs->naxis*wcs->naxis;

  /* Allocate the necessary array. */
  errno=0;
  out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'out'",
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
  else if(wcs->altlin & 0x4)     /* Has CROTAi array.   */
    {
      /* Basic sanity checks. */
      if(wcs->naxis!=2)
        error(EXIT_FAILURE, 0, "%s: CROTAi currently on works in 2 "
              "dimensions.", __func__);
      if(wcs->crota[0]!=0.0)
        error(EXIT_FAILURE, 0, "%s: CROTA1 is not zero", __func__);

      /* CROTAi keywords are depreciated in the FITS standard. However, old
         files may still use them. For a full description of CROTAi
         keywords and their history (along with the conversion equations
         here), please see the link below:

         https://fits.gsfc.nasa.gov/users_guide/users_guide/node57.html

         Just note that the equations of the link above convert CROTAi to
         PC. But here we want the "final" matrix (after multiplication by
         the 'CDELT' values). So to speed things up, we won't bother
         dividing and then multiplying by the same CDELT values in the
         off-diagonal elements. */
      crota2=wcs->crota[1];
      out[0] = wcs->cdelt[0] * cos(crota2);
      out[1] = -1 * wcs->cdelt[1] *sin(crota2);
      out[2] = wcs->cdelt[0] * sin(crota2);
      out[3] = wcs->cdelt[1] * cos(crota2);
    }
  else
    error(EXIT_FAILURE, 0, "%s: currently only PCi_ja and CDi_ja keywords "
          "are recognized", __func__);

  /* Return the result */
  return out;
}





/* According to the FITS standard, in the 'PCi_j' WCS formalism, the matrix
   elements m_{ij} are encoded in the 'PCi_j' keywords and the scale
   factors are encoded in the 'CDELTi' keywords. There is also another
   formalism (the 'CDi_j' formalism) which merges the two into one
   matrix.

   However, WCSLIB's internal operations are apparently done in the 'PCi_j'
   formalism. So its outputs are also all in that format by default. When
   the input is a 'CDi_j', WCSLIB will still read the image into the
   'PCi_j' formalism and the 'CDELTi's are set to 1. This function will
   decompose the two matrices to give a reasonable 'CDELTi' and 'PCi_j' in
   such cases. */
void
gal_wcs_decompose_pc_cdelt(struct wcsprm *wcs)
{
  double *ps;
  size_t i, j;

  /* If there is on WCS, then don't do anything. */
  if(wcs==NULL) return;

  /* The correction is only needed when the PC matrix is filled. */
  if(wcs->pc)
    {
      /* Get the pixel scale. */
      ps=gal_wcs_pixel_scale(wcs);
      if(ps==NULL) return;

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

      /* According to the 'wcslib/wcs.h' header: "In particular, wcsset()
         resets wcsprm::cdelt to unity if CDi_ja is present (and no
         PCi_ja).". So apparently, when the input is a 'CDi_j', it might
         expect the 'CDELTi' elements to be 1.0. But we have changed that
         here, so we will correct the 'altlin' element of the WCS structure
         to make sure that WCSLIB only looks into the 'PCi_j' and 'CDELTi'
         and makes no assumptioins about 'CDELTi'. */
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
  int warning_printed;
  gal_data_t *pixscale;
  size_t i, j, n, maxj, *permutation;
  double jvmax, *a, *out, *v, maxrow, minrow;

  /* Only continue if a WCS exists. */
  if(wcs==NULL) return NULL;


  /* Write the full WCS rotation matrix into an array, irrespective of what
     style it was stored in the wcsprm structure ('PCi_j' style or 'CDi_j'
     style). */
  a=gal_wcs_warp_matrix(wcs);


  /* Now that everything is good, we can allocate the necessary memory. */
  n=wcs->naxis;
  v=gal_pointer_allocate(GAL_TYPE_FLOAT64, n*n, 0, __func__, "v");
  permutation=gal_pointer_allocate(GAL_TYPE_SIZE_T, n, 0, __func__,
                                   "permutation");
  pixscale=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &n, NULL,
                          0, -1, 1, NULL, NULL, NULL);


  /* To avoid confusing issues with floating point errors being written in
     the non-diagonal elements of the FITS header PC or CD matrices, we
     need to check if the minimum and maximum values in each row are not
     several orders of magnitude apart.

     Note that in some cases (for example a spectrum), one axis might be in
     degrees (value around 1e-5) and the other in angestroms (value around
     1e-10). So we can't look at the minimum and maximum of the whole
     matrix. However, in such cases, people will probably not warp/rotate
     the image to mix the coordinates. So the important thing to check is
     the minimum and maximum (non-zero) values in each row. */
  warning_printed=0;
  for(i=0;i<n;++i)
    {
      /* Find the minimum and maximum values in each row. */
      minrow=FLT_MAX;
      maxrow=-FLT_MAX;
      for(j=0;j<n;++j)
        if(a[i*n+j]!=0.0) /* We aren't concerned with 0 valued elements. */
          {
            /* We have to use absolutes because in cases like RA, the
               diagonal values in different rows may have different signs*/
            if(fabs(a[i*n+j])<minrow) minrow=fabs(a[i*n+j]);
            if(fabs(a[i*n+j])>maxrow) maxrow=fabs(a[i*n+j]);
          }

      /* Do the check, print warning and make correction. */
      if(maxrow!=minrow && maxrow/minrow>1e5 && warning_printed==0)
        {
          fprintf(stderr, "\nWARNING: The input WCS matrix (possibly taken "
                  "from the FITS header keywords starting with 'CD' or 'PC') "
                  "contains values with very different scales (more than "
                  "10^5 different). This is probably due to floating point "
                  "errors. These values might bias the pixel scale (and "
                  "subsequent) calculations.\n\n"
                  "You can see the respective matrix with one of the "
                  "following two commands (depending on how the FITS file "
                  "was written). Recall that if the desired extension/HDU "
                  "isn't the default, you can choose it with the '--hdu' "
                  "(or '-h') option before the '|' sign in these commands."
                  "\n\n"
                  "    $ astfits file.fits -p | grep 'PC._.'\n"
                  "    $ astfits file.fits -p | grep 'CD._.'\n\n"
                  "You can delete the ones with obvious floating point "
                  "error values using the following command (assuming you "
                  "want to delete 'CD1_2' and 'CD2_1'). Afterwards, you can "
                  "re-run your original command to remove this warning "
                  "message and possibly correct errors that it might have "
                  "caused.\n\n"
                  "    $ astfits file.fits --delete=CD1_2 --delete=CD2_1\n\n"
                  );
          warning_printed=1;
        }
    }


  /* Fill in the necessary GSL vector and Matrix structures. */
  S.size=n;     S.stride=1;                S.data=pixscale->array;
  V.size1=n;    V.size2=n;    V.tda=n;     V.data=v;
  A.size1=n;    A.size2=n;    A.tda=n;     A.data=a;


  /* Run GSL's Singular Value Decomposition, using one-sided Jacobi
     orthogonalization which computes the singular (scale) values to a
     higher relative accuracy. */
  gsl_linalg_SV_decomp_jacobi(&A, &V, &S);


  /* The raw pixel scale array produced from the singular value
     decomposition above is ordered based on values, not the input. So when
     the pixel scales in all the dimensions aren't the same (the units of
     the dimensions differ), the order of the values in 'pixelscale' will
     not necessarily correspond to the input's dimensions.

     To correct the order, we can use the 'V' matrix to find the original
     position of the pixel scale values and then use permutation to
     re-order it correspondingly. The column with the largest (absolute)
     value will be taken as the one to be used for each row. */
  for(i=0;i<n;++i)
    {
      /* Find the column with the maximum value. */
      maxj=-1;
      jvmax=-FLT_MAX;
      for(j=0;j<n;++j)
        if(fabs(v[i*n+j])>jvmax)
          {
            maxj=j;
            jvmax=fabs(v[i*n+j]);
          }

      /* Use the column with the maximum value for this dimension. */
      permutation[i]=maxj;
    }


  /* Apply the permutation described above. */
  gal_permutation_apply(pixscale, permutation);


  /* Clean up and return. */
  free(a);
  free(v);
  free(permutation);
  out=pixscale->array;
  pixscale->array=NULL;
  gal_data_free(pixscale);
  return out;
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
  if(wcs->naxis!=2) return NAN;

  /* Check if the units of the axis are degrees or not. Currently all FITS
     images I have worked with use 'deg' for degrees. If other alternatives
     exist, we can add corrections later. */
  if( strcmp("deg", wcs->cunit[0]) || strcmp("deg", wcs->cunit[1]) )
    return NAN;

  /* Get the pixel scales along each axis in degrees, then multiply. */
  pixscale=gal_wcs_pixel_scale(wcs);
  if(pixscale==NULL) return NAN;

  /* Clean up and return the result. */
  out = pixscale[0] * pixscale[1] * 3600.0f * 3600.0f;
  free(pixscale);
  return out;
}




















/**************************************************************/
/**********            Array conversion            ************/
/**************************************************************/
/* Some sanity checks for the WCS conversion functions. */
static void
wcs_convert_sanity_check_alloc(gal_data_t *coords, struct wcsprm *wcs,
                               const char *func, int **stat, double **phi,
                               double **theta, double **world,
                               double **pixcrd, double **imgcrd)
{
  gal_data_t *tmp;
  size_t ndim=0, firstsize=0, size=coords->size;

  /* Make sure a WCS structure is actually given. */
  if(wcs==NULL)
    error(EXIT_FAILURE, 0, "%s: input WCS structure is NULL", func);

  for(tmp=coords; tmp!=NULL; tmp=tmp->next)
    {
      /* Count how many coordinates are given. */
      ++ndim;

      /* Check the type of the input. */
      if(tmp->type!=GAL_TYPE_FLOAT64)
        error(EXIT_FAILURE, 0, "%s: input coordinates must have 'float64' "
              "type", func);

      /* Make sure it has a single dimension. */
      if(tmp->ndim!=1)
        error(EXIT_FAILURE, 0, "%s: input coordinates for each dimension "
              "must each be one dimensional. Coordinate dataset %zu of the "
              "inputs has %zu dimensions", func, ndim, tmp->ndim);

      /* See if all inputs have the same size. */
      if(ndim==1) firstsize=tmp->size;
      else
        if(firstsize!=tmp->size)
          error(EXIT_FAILURE, 0, "%s: all input coordinates must have the "
                "same number of elements. Coordinate dataset %zu has %zu "
                "elements while the first coordinate has %zu", func, ndim,
                tmp->size, firstsize);
    }

  /* See if the number of coordinates given corresponds to the dimensions
     of the WCS structure. */
  if(ndim!=wcs->naxis)
    error(EXIT_FAILURE, 0, "%s: the number of input coordinates (%zu) does "
          "not match the dimensions of the input WCS structure (%d)", func,
          ndim, wcs->naxis);

  /* Allocate all the necessary arrays. */
  *phi    = gal_pointer_allocate( GAL_TYPE_FLOAT64, size,      0, __func__,
                                  "phi");
  *stat   = gal_pointer_allocate( GAL_TYPE_INT32,   size,      1, __func__,
                                  "stat");
  *theta  = gal_pointer_allocate( GAL_TYPE_FLOAT64, size,      0, __func__,
                                  "theta");
  *world  = gal_pointer_allocate( GAL_TYPE_FLOAT64, ndim*size, 0, __func__,
                                  "world");
  *imgcrd = gal_pointer_allocate( GAL_TYPE_FLOAT64, ndim*size, 0, __func__,
                                  "imgcrd");
  *pixcrd = gal_pointer_allocate( GAL_TYPE_FLOAT64, ndim*size, 0, __func__,
                                  "pixcrd");
}





/* In Gnuastro, each column (coordinate for WCS conversion) is treated as a
   separate array in a 'gal_data_t' that are linked through a linked
   list. But in WCSLIB, the input is a single array (with multiple
   columns). This function will convert between the two. */
static void
wcs_convert_list_to_array(gal_data_t *list, double *array, int *stat,
                          size_t ndim, int listtoarray)
{
  size_t i, d=0;
  gal_data_t *tmp;

  for(tmp=list; tmp!=NULL; tmp=tmp->next)
    {
      /* Put all this coordinate's values into the single array that is
         input into or output from WCSLIB. */
      for(i=0;i<list->size;++i)
        {
          if(listtoarray)
            array[i*ndim+d] = ((double *)(tmp->array))[i];
          else
            ((double *)(tmp->array))[i] = stat[i] ? NAN : array[i*ndim+d];
        }

      /* Increment the dimension. */
      ++d;
    }
}





/* Prepare the output of the WCS conversion functions. */
static gal_data_t *
wcs_convert_prepare_out(gal_data_t *coords, struct wcsprm *wcs, int inplace)
{
  size_t i;
  gal_data_t *out=NULL;
  if(inplace)
    out=coords;
  else
    for(i=0;i<wcs->naxis;++i)
      gal_list_data_add_alloc(&out, NULL, GAL_TYPE_FLOAT64, 1,
                              &coords->size, NULL, 0, coords->minmapsize,
                              coords->quietmmap, wcs->ctype[i], wcs->cunit[i],
                              NULL);
  return out;
}





/* Convert world coordinates to image coordinates given the input WCS
   structure. The input must be a linked list of data structures of float64
   ('double') type. The top element of the linked list must be the first
   coordinate and etc. If 'inplace' is non-zero, then the output will be
   written into the input's allocated space. */
gal_data_t *
gal_wcs_world_to_img(gal_data_t *coords, struct wcsprm *wcs, int inplace)
{
  gal_data_t *out;
  int status, *stat=NULL, ncoord=coords->size, nelem;
  double *phi=NULL, *theta=NULL, *world=NULL, *pixcrd=NULL, *imgcrd=NULL;

  /* Some sanity checks. */
  wcs_convert_sanity_check_alloc(coords, wcs, __func__, &stat, &phi, &theta,
                                 &world, &pixcrd, &imgcrd);
  nelem=wcs->naxis; /* We have to make sure a WCS is given first. */


  /* Write the values from the input list of separate columns into a single
     array (WCSLIB input). */
  wcs_convert_list_to_array(coords, world, stat, wcs->naxis, 1);


  /* Use WCSLIB's wcss2p for the conversion. */
  status=wcss2p(wcs, ncoord, nelem, world, phi, theta, imgcrd, pixcrd, stat);
  if(status)
    error(EXIT_FAILURE, 0, "%s: wcss2p ERROR %d: %s", __func__, status,
          wcs_errmsg[status]);


  /* For a sanity check.
  {
    size_t i;
    printf("\n\n%s sanity check:\n", __func__);
    for(i=0;i<coords->size;++i)
      printf("(%g, %g, %g) --> (%g, %g, %g), [stat: %d]\n",
              world[i*3],  world[i*3+1 ], world[i*3+2],
             pixcrd[i*3], pixcrd[i*3+1], pixcrd[i*3+2], stat[i]);
  }
  */


  /* Allocate the output arrays if they were not already allocated. */
  out=wcs_convert_prepare_out(coords, wcs, inplace);


  /* Write the output from a single array (WCSLIB output) into the output
     list of this function. */
  wcs_convert_list_to_array(out, pixcrd, stat, wcs->naxis, 0);


  /* Clean up. */
  free(phi);
  free(stat);
  free(theta);
  free(world);
  free(pixcrd);

  /* Return the output list of coordinates. */
  return out;
}





/* Similar to 'gal_wcs_world_to_img'. */
gal_data_t *
gal_wcs_img_to_world(gal_data_t *coords, struct wcsprm *wcs, int inplace)
{
  gal_data_t *out;
  int status, *stat=NULL, ncoord=coords->size, nelem;
  double *phi=NULL, *theta=NULL, *world=NULL, *pixcrd=NULL, *imgcrd=NULL;

  /* Some sanity checks. */
  wcs_convert_sanity_check_alloc(coords, wcs, __func__, &stat, &phi, &theta,
                                 &world, &pixcrd, &imgcrd);
  nelem=wcs->naxis; /* We have to make sure a WCS is given first. */


  /* Write the values from the input list of separate columns into a single
     array (WCSLIB input). */
  wcs_convert_list_to_array(coords, pixcrd, stat, wcs->naxis, 1);


  /* Use WCSLIB's wcsp2s for the conversion. */
  status=wcsp2s(wcs, ncoord, nelem, pixcrd, imgcrd, phi, theta, world, stat);
  if(status)
    error(EXIT_FAILURE, 0, "%s: wcsp2s ERROR %d: %s", __func__, status,
          wcs_errmsg[status]);


  /* For a sanity check.
  {
    size_t i;
    printf("\n\n%s sanity check:\n", __func__);
    for(i=0;i<coords->size;++i)
      printf("(%g, %g, %g) --> (%g, %g, %g), [stat: %d]\n",
             pixcrd[i*3], pixcrd[i*3+1], pixcrd[i*3+2],
             world[i*3],  world[i*3+1],  world[i*3+2], stat[i]);
  }
  */


  /* Allocate the output arrays if they were not already allocated. */
  out=wcs_convert_prepare_out(coords, wcs, inplace);


  /* Write the output from a single array (WCSLIB output) into the output
     list of this function. */
  wcs_convert_list_to_array(out, world, stat, wcs->naxis, 0);


  /* Clean up. */
  free(phi);
  free(stat);
  free(theta);
  free(world);
  free(pixcrd);


  /* Return the output list of coordinates. */
  return out;
}
