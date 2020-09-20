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
#include <wcslib/wcsmath.h>

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>
#include <gnuastro/permutation.h>

#include <gnuastro-internal/checkset.h>

#if GAL_CONFIG_HAVE_WCSLIB_DIS_H
#include <wcslib/dis.h>
#include <gnuastro-internal/wcsdistortion.h>
#endif





/* Static functions on for this file. */
static void
gal_wcs_to_cd(struct wcsprm *wcs);





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
  int fixstatus[NWCSFIX]={0};/* For the various wcsfix checks.          */
  int relax    = WCSHDR_all; /* Macro: use all informal WCS extensions. */
  int ctrl     = 0;          /* Don't report why a keyword wasn't used. */
  int nreject  = 0;          /* Number of keywords rejected for syntax. */
  int fixctrl  = 1;          /* Correct non-standard units in wcsfix.   */
  void *fixnaxis = NULL;     /* For now disable cylfix() with this      */
                             /* (because it depends on image size).     */

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

          /* Fix non-standard WCS features. */
          if( wcsfix(fixctrl, fixnaxis, wcs, fixstatus) )
            {
              if(fixstatus[CDFIX])
                error(0, 0, "%s: (warning) wcsfix status for cdfix: %d",
                      __func__, fixstatus[CDFIX]);
              if(fixstatus[DATFIX])
                error(0, 0, "%s: (warning) wcsfix status for datfix: %d",
                      __func__, fixstatus[DATFIX]);
#if GAL_CONFIG_HAVE_WCSLIB_OBSFIX
              if(fixstatus[OBSFIX])
                error(0, 0, "%s: (warning) wcsfix status for obsfix: %d",
                      __func__, fixstatus[OBSFIX]);
#endif
              if(fixstatus[UNITFIX])
                error(0, 0, "%s: (warning) wcsfix status for unitfix: %d",
                      __func__, fixstatus[UNITFIX]);
              if(fixstatus[SPCFIX])
                error(0, 0, "%s: (warning) wcsfix status for spcfix: %d",
                      __func__, fixstatus[SPCFIX]);
              if(fixstatus[CELFIX])
                error(0, 0, "%s: (warning) wcsfix status for celfix: %d",
                      __func__, fixstatus[CELFIX]);
              if(fixstatus[CYLFIX])
                error(0, 0, "%s: (warning) wcsfix status for cylfix: %d",
                      __func__, fixstatus[CYLFIX]);
            }

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
            {
              /* A correctly useful WCS is present. When no PC matrix
                 elements were present in the header, the default PC matrix
                 (a unity matrix) is used. In this case WCSLIB doesn't set
                 'altlin' (and gives it a value of 0). In Gnuastro, later on,
                 we might need to know the type of the matrix used, so in
                 such a case, we will set 'altlin' to 1. */
              if(wcs->altlin==0) wcs->altlin=1;
            }
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





/* Extract the dimension name from CTYPE. */
char *
gal_wcs_dimension_name(struct wcsprm *wcs, size_t dimension)
{
  size_t i;
  char *out;

  /* Make sure a WCS pointer actually exists. */
  if(wcs==NULL) return NULL;

  /* Make sure the requested dimension is not larger than the number of
     dimensions in the WCS. */
  if(dimension >= wcs->naxis) return NULL;

  /* Make a copy of the CTYPE value and set the first occurance of '-' to
     '\0', to avoid the projection type. */
  gal_checkset_allocate_copy(wcs->ctype[dimension], &out);
  for(i=0;i<strlen(out);++i) if(out[i]=='-') out[i]='\0';

  /* Return the output array. */
  return out;
}



















/*************************************************************
 ***********               Write WCS               ***********
 *************************************************************/
void
gal_wcs_write_in_fitsptr(fitsfile *fptr, struct wcsprm *wcs)
{
  char *wcsstr;
  int tpvdist, status=0, nkeyrec;

  /* Prepare the main rotation matrix. Note that for TPV distortion, WCSLIB
     versions 7.3 and before couldn't deal with the CDELT keys, so to be
     safe, in such cases, we'll remove the effect of CDELT in the
     'gal_wcs_to_cd' function. */
  tpvdist=wcs->lin.disseq && !strcmp(wcs->lin.disseq->dtype[1], "TPV");
  if( tpvdist ) gal_wcs_to_cd(wcs);
  else          gal_wcs_decompose_pc_cdelt(wcs);

  /* Clean up small errors in the PC matrix and CDELT values. */
  gal_wcs_clean_errors(wcs);

  /* Convert the WCS information to text. */
  status=wcshdo(WCSHDO_safe, wcs, &nkeyrec, &wcsstr);
  if(status)
    error(0, 0, "%s: WARNING: WCSLIB error, no WCS in output.\n"
          "wcshdu ERROR %d: %s", __func__, status,
          wcs_errmsg[status]);
  else
    gal_fits_key_write_wcsstr(fptr, wcsstr, nkeyrec);
  status=0;

   /* WCSLIB is going to write PC+CDELT keywords in any case. But when we
      have a TPV distortion, it is cleaner to use a CD matrix. Also,
      including and before version 7.3, WCSLIB wouldn't convert coordinates
      properly if the PC matrix is used with the TPV distortion. So to help
      users with WCSLIB 7.3 or earlier, we need to conver the PC matrix to
      CD. 'gal_wcs_to_cd' function already made sure that CDELT=1, so
      effectively the CD matrix and PC matrix are equivalent, we just need
      to convert the keyword names and delete the CDELT keywords. Note that
      zero-valued PC/CD elements may not be present, so we'll manually set
      'status' to zero and not let CFITSIO crash.*/
  if(tpvdist)
    {
      status=0; fits_modify_name(fptr, "PC1_1", "CD1_1", &status);
      status=0; fits_modify_name(fptr, "PC1_2", "CD1_2", &status);
      status=0; fits_modify_name(fptr, "PC2_1", "CD2_1", &status);
      status=0; fits_modify_name(fptr, "PC2_2", "CD2_2", &status);
      status=0; fits_delete_str(fptr, "CDELT1", &status);
      status=0; fits_delete_str(fptr, "CDELT2", &status);
      status=0;
      fits_write_comment(fptr, "The CD matrix is used instead of the "
                         "PC+CDELT due to conflicts with TPV distortion "
                         "in WCSLIB 7.3 (released on 2020/06/03) and "
                         "ealier. By default Gnuastro will write "
                         "PC+CDELT matrices because the rotation (PC) and "
                         "pixel-scale (CDELT) are separate; providing "
                         "more physically relevant metadata for human "
                         "readers (PC+CDELT is also the default format "
                         "of WCSLIB).", &status);
    }
}





void
gal_wcs_write(struct wcsprm *wcs, char *filename,
              char *extname, gal_fits_list_key_t *headers,
              char *program_string)
{
  int status=0;
  size_t ndim=0;
  fitsfile *fptr;
  long *naxes=NULL;

  /* Small sanity checks */
  if(wcs==NULL)
    error(EXIT_FAILURE, 0, "%s: input WCS is NULL", __func__);
  if( gal_fits_name_is_fits(filename)==0 )
    error(EXIT_FAILURE, 0, "%s: not a FITS suffix", filename);

  /* Open the file for writing */
  fptr=gal_fits_open_to_write(filename);

  /* Create the FITS file. */
  fits_create_img(fptr, gal_fits_type_to_bitpix(GAL_TYPE_UINT8),
                  ndim, naxes, &status);
  gal_fits_io_error(status, NULL);

  /* Remove the two comment lines put by CFITSIO. Note that in some cases,
     it might not exist. When this happens, the status value will be
     non-zero. We don't care about this error, so to be safe, we will just
     reset the status variable after these calls. */
  fits_delete_key(fptr, "COMMENT", &status);
  fits_delete_key(fptr, "COMMENT", &status);
  status=0;

  /* If an extension name was requested, add it. */
  if(extname)
    fits_write_key(fptr, TSTRING, "EXTNAME", extname, "", &status);

  /* Write the WCS structure. */
  gal_wcs_write_in_fitsptr(fptr, wcs);

  /* Write all the headers and the version information. */
  gal_fits_key_write_version_in_ptr(&headers, program_string, fptr);

  /* Close the FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}



















/*************************************************************
 ***********              Distortions              ***********
 *************************************************************/
int
gal_wcs_distortion_from_string(char *distortion)
{
  if(      !strcmp(distortion,"TPD") ) return GAL_WCS_DISTORTION_TPD;
  else if( !strcmp(distortion,"SIP") ) return GAL_WCS_DISTORTION_SIP;
  else if( !strcmp(distortion,"TPV") ) return GAL_WCS_DISTORTION_TPV;
  else if( !strcmp(distortion,"DSS") ) return GAL_WCS_DISTORTION_DSS;
  else if( !strcmp(distortion,"WAT") ) return GAL_WCS_DISTORTION_WAT;
  else
    error(EXIT_FAILURE, 0, "WCS distortion name '%s' not recognized, "
          "currently recognized names are 'TPD', 'SIP', 'TPV', 'DSS' and "
          "'WAT'", distortion);

  /* Control should not reach here. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return GAL_WCS_DISTORTION_INVALID;
}





char *
gal_wcs_distortion_to_string(int distortion)
{
  /* Return the proper literal string. */
  switch(distortion)
    {
    case GAL_WCS_DISTORTION_TPD: return "TPD";
    case GAL_WCS_DISTORTION_SIP: return "SIP";
    case GAL_WCS_DISTORTION_TPV: return "TPV";
    case GAL_WCS_DISTORTION_DSS: return "DSS";
    case GAL_WCS_DISTORTION_WAT: return "WAT";
    default:
      error(EXIT_FAILURE, 0, "WCS distortion id '%d' isn't recognized",
            distortion);
    }

  /* Control should not reach here. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NULL;
}





/* Check the type of distortion present and return the appropriate
   integer based on `enum gal_wcs_distortion`.

   Parameters:
    struct wcsprm *wcs - The wcs parameters of the fits file.

   Return:
    int out_distortion - The type of distortion present. */
int
gal_wcs_distortion_identify(struct wcsprm *wcs)
{
#if GAL_CONFIG_HAVE_WCSLIB_DIS_H
  struct disprm *dispre=NULL;
  struct disprm *disseq=NULL;

  /* Small sanity check. */
  if(wcs==NULL) return GAL_WCS_DISTORTION_INVALID;

  /* To help in reading. */
  disseq=wcs->lin.disseq;
  dispre=wcs->lin.dispre;

  /* Check if distortion present. */
  if( disseq==NULL && dispre==NULL ) return GAL_WCS_DISTORTION_INVALID;

  /* Check the type of distortion.

     As mentioned in the WCS paper IV section 2.4.2 available at
     https://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf, the
     DPja and DQia keywords are used to record the parameters required by
     the prior and sequent distortion functions respectively.

     Now, as mentioned in dis.h file reference section in WCSLIB manual
     given here
     https://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/dis_8h.html, TPV,
     DSS, and WAT are sequent polynomial distortions, while SIP is prior
     polynomial distortion. TPD is a superset of all these distortions and
     hence can be used both as a prior and sequent distortion polynomial.

     References and citations:
     "Representations of distortions in FITS world coordinate systems",
     Calabretta, M.R. et al. (WCS Paper IV, draft dated 2004/04/22)
      */

  if( dispre != NULL )
    {
      if(      !strcmp(*dispre->dtype, "SIP") ) return GAL_WCS_DISTORTION_SIP;
      else if( !strcmp(*dispre->dtype, "TPD") ) return GAL_WCS_DISTORTION_TPD;
      else
        error(EXIT_FAILURE, 0, "%s: distortion '%s' isn't recognized in "
              "the 'dispre' structure of the given 'wcsprm'", __func__,
              *dispre->dtype);
    }
  else if( disseq != NULL )
    {
      if(      !strcmp(*disseq->dtype, "TPV") ) return GAL_WCS_DISTORTION_TPV;
      else if( !strcmp(*disseq->dtype, "TPD") ) return GAL_WCS_DISTORTION_TPD;
      else if( !strcmp(*disseq->dtype, "DSS") ) return GAL_WCS_DISTORTION_DSS;
      else if( !strcmp(*disseq->dtype, "WAT") ) return GAL_WCS_DISTORTION_WAT;
      else
        error(EXIT_FAILURE, 0, "%s: distortion '%s' isn't recognized in "
              "the 'disseq' structure of the given 'wcsprm'", __func__,
              *dispre->dtype);
    }

  /* Control should not reach here. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to fix "
        "the problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
#else
  /* The 'wcslib/dis.h' isn't present. */
  error(EXIT_FAILURE, 0, "%s: the installed version of WCSLIB on this "
        "system doesn't have the 'dis.h' header! Thus Gnuastro can't do "
        "distortion-related operations on the world coordinate system "
        "(WCS). To use these features, please upgrade your version "
        "of WCSLIB and re-build Gnuastro", __func__);
#endif
  return GAL_WCS_DISTORTION_INVALID;
}










/* Convert a given distrotion type to other.

  Parameters:
    struct wcsprm *wcs - The wcs parameters of the fits file.
    int out_distortion - The desired output distortion.
    size_t* fitsize    - The size of the array along each dimension.

  Return:
    struct wcsprm *outwcs - The transformed wcs parameters in the
                            required distortion type. */
struct wcsprm *
gal_wcs_distortion_convert(struct wcsprm *inwcs, int outdisptype,
                           size_t *fitsize)
{
#if GAL_CONFIG_HAVE_WCSLIB_DIS_H
  struct wcsprm *outwcs=NULL;
  int indisptype=gal_wcs_distortion_identify(inwcs);

  /* Make sure we have a PC+CDELT structure in the input WCS. */
  gal_wcs_decompose_pc_cdelt(inwcs);

  /* If the input and output types are the same, just copy the input,
     otherwise, do the conversion. */
  if(indisptype==outdisptype) outwcs=gal_wcs_copy(inwcs);
  else
    switch(indisptype)
      {
      /* If there is no distortion in the input, just return a
         newly-allocated copy. */
      case GAL_WCS_DISTORTION_INVALID: outwcs=gal_wcs_copy(inwcs); break;

      /* Input's distortion is SIP. */
      case GAL_WCS_DISTORTION_SIP:
        switch(outdisptype)
          {
          case GAL_WCS_DISTORTION_TPV:
            outwcs=gal_wcsdistortion_sip_to_tpv(inwcs);
            break;
          default:
            error(EXIT_FAILURE, 0, "%s: conversion from %s to %s is not yet "
                  "supported. Please contact us at '%s'", __func__,
                  gal_wcs_distortion_to_string(indisptype),
                  gal_wcs_distortion_to_string(outdisptype),
                  PACKAGE_BUGREPORT);
              }
        break;

      /* Input's distortion is TPV. */
      case GAL_WCS_DISTORTION_TPV:
        switch(outdisptype)
          {
          case GAL_WCS_DISTORTION_SIP:
            if(fitsize==NULL)
              error(EXIT_FAILURE, 0, "%s: the size array is necessary "
                    "for this conversion", __func__);
            outwcs=gal_wcsdistortion_tpv_to_sip(inwcs, fitsize);
            break;
          default:
            error(EXIT_FAILURE, 0, "%s: conversion from %s to %s is not yet "
                  "supported. Please contact us at '%s'", __func__,
                  gal_wcs_distortion_to_string(indisptype),
                  gal_wcs_distortion_to_string(outdisptype),
                  PACKAGE_BUGREPORT);
          }
        break;

      /* Input's distortion is not yet supported.. */
      case GAL_WCS_DISTORTION_TPD:
      case GAL_WCS_DISTORTION_DSS:
      case GAL_WCS_DISTORTION_WAT:
        error(EXIT_FAILURE, 0, "%s: input %s distortion is not yet "
              "supported. Please contact us at '%s'", __func__,
              gal_wcs_distortion_to_string(indisptype),
              PACKAGE_BUGREPORT);

      /* A bug! This distortion is not yet recognized. */
      default:
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
              "the problem. The identifier '%d' is not recognized as a "
              "distortion", __func__, PACKAGE_BUGREPORT, indisptype);
      }

  /* Return the converted WCS. */
  return outwcs;
#else
  /* The 'wcslib/dis.h' isn't present. */
  error(EXIT_FAILURE, 0, "%s: the installed version of WCSLIB on this "
        "system doesn't have the 'dis.h' header! Thus Gnuastro can't do "
        "distortion-related operations on the world coordinate system "
        "(WCS). To use these features, please upgrade your version "
        "of WCSLIB and re-build Gnuastro", __func__);
  return NULL;
#endif
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




/* Clean up small/negligible errros that are clearly caused by measurement
   errors in the PC and CDELT elements. */
void
gal_wcs_clean_errors(struct wcsprm *wcs)
{
  double crdcheck=NAN;
  size_t i, crdnum=0, ndim=wcs->naxis;
  double mean, crdsum=0, sum=0, min=FLT_MAX, max=0;
  double *pc=wcs->pc, *cdelt=wcs->cdelt, *crder=wcs->crder;

  /* First clean up CDELT: if the CRDER keyword is set, then we'll use that
     as a reference, if not, we'll use the absolute floating point error
     defined in 'GAL_WCS_FLTERR'. */
  for(i=0; i<ndim; ++i)
    {
      sum+=cdelt[i];
      if(cdelt[i]>max) max=cdelt[i];
      if(cdelt[i]<min) min=cdelt[i];
      if(crder[i]!=UNDEFINED) {++crdnum; crdsum=crder[i];}
    }
  mean=sum/ndim;
  crdcheck = crdnum ? crdsum/crdnum : GAL_WCS_FLTERROR;
  if( (max-min)/mean < crdcheck )
    for(i=0; i<ndim; ++i)
      cdelt[i]=mean;

  /* Now clean up the PC elements. If the diagonal elements are too close
     to 0, 1, or -1, set them to 0 or 1 or -1. */
  if(pc)
    for(i=0;i<ndim*ndim;++i)
      {
        if(      fabs(pc[i] -  0 ) < GAL_WCS_FLTERROR ) pc[i]=0;
        else if( fabs(pc[i] -  1 ) < GAL_WCS_FLTERROR ) pc[i]=1;
        else if( fabs(pc[i] - -1 ) < GAL_WCS_FLTERROR ) pc[i]=-1;
      }
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

  /* The correction is only needed when the PC matrix is filled. Note that
     internally, WCSLIB always uses the PC matrix, even when the input has
     a CD matrix.*/
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





/* Set the WCS structure to use the CD matrix. */
static void
gal_wcs_to_cd(struct wcsprm *wcs)
{
  size_t i, j, naxis;

  /* If there is on WCS, then don't do anything. */
  if(wcs==NULL) return;

  /* 'wcs->altlin' identifies which rotation element is being used (PCi_j,
     CDi_J or CROTAi). For PCi_j, the first bit will be 1 (==1), for CDi_j,
     the second bit is 1 (==2) and for CROTAi, the third bit is 1 (==4). */
  naxis=wcs->naxis;
  switch(wcs->altlin)
    {
   /* PCi_j: Convert it to CDi_j. */
    case 1:

      /* Fill in the CD matrix and correct the PC and CDELT arrays. We have
         to do this because ultimately, WCSLIB will be writing the PC and
         CDELT keywords, even when 'altlin' is 2. So effectively we have to
         multiply the PC and CDELT matrices, then set cdelt=1 in all
         dimensions. This is actually how WCSLIB reads a FITS header with
         only a CD matrix. */
      for(i=0;i<naxis;++i)
        {
          for(j=0;j<naxis;++j)
            wcs->cd[i*naxis+j] = wcs->pc[i*naxis+j] *= wcs->cdelt[i];
          wcs->cdelt[i]=1;
        }

      /* Set the altlin to be the CD matrix and free the PC matrix. */
      wcs->altlin=2;
      break;

    /* CDi_j: No need to do any conversion. */
    case 2: return; break;

    /* CROTAi: not yet supported. */
    case 4:
      error(0, 0, "%s: WARNING: Conversion of 'CROTAi' keywords to the CD "
            "matrix is not yet supported (for lack of time!), please "
            "contact us at %s to add this feature. But this may not cause a "
            "problem at all, so please check if the output's WCS is "
            "reasonable", __func__, PACKAGE_BUGREPORT);
      break;

    /* The value isn't supported! */
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
            "problem. The value %d for wcs->altlin isn't recognized",
            __func__, PACKAGE_BUGREPORT, wcs->altlin);
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
      if(maxrow!=minrow
         && maxrow/minrow>1e5    /* The difference between elements is large */
         && maxrow/minrow<GAL_WCS_FLTERROR
         && warning_printed==0)
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





int
gal_wcs_coverage(char *filename, char *hdu, size_t *ondim,
                 double **ocenter, double **owidth, double **omin,
                 double **omax)
{
  fitsfile *fptr;
  struct wcsprm *wcs;
  int nwcs=0, type, status=0;
  char *name=NULL, *unit=NULL;
  gal_data_t *tmp, *coords=NULL;
  size_t i, ndim, *dsize=NULL, numrows;
  double *x, *y, *z, *min, *max, *center, *width;

  /* Read the desired WCS. */
  wcs=gal_wcs_read(filename, hdu, 0, 0, &nwcs);

  /* If a WCS doesn't exist, return NULL. */
  if(wcs==NULL) return 0;

  /* Make sure the input HDU is an image. */
  if( gal_fits_hdu_format(filename, hdu) != IMAGE_HDU )
    error(EXIT_FAILURE, 0, "%s (hdu %s): is not an image HDU, the "
          "'--skycoverage' option only applies to image extensions",
          filename, hdu);

  /* Get the array information of the image. */
  fptr=gal_fits_hdu_open(filename, hdu, READONLY);
  gal_fits_img_info(fptr, &type, ondim, &dsize, &name, &unit);
  fits_close_file(fptr, &status);
  ndim=*ondim;

  /* Abort if we have more than 3 dimensions. */
  if(ndim==1 || ndim>3) return 0;

  /* Allocate the output datasets. */
  center=*ocenter=gal_pointer_allocate(GAL_TYPE_FLOAT64, ndim, 0, __func__,
                                       "ocenter");
  width=*owidth=gal_pointer_allocate(GAL_TYPE_FLOAT64, ndim, 0, __func__,
                                     "owidth");
  min=*omin=gal_pointer_allocate(GAL_TYPE_FLOAT64, ndim, 0, __func__,
                                 "omin");
  max=*omax=gal_pointer_allocate(GAL_TYPE_FLOAT64, ndim, 0, __func__,
                                 "omax");

  /* Now that we have the number of dimensions in the image, allocate the
     space needed for the coordinates. */
  numrows = (ndim==2) ? 5 : 9;
  for(i=0;i<ndim;++i)
    {
      tmp=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &numrows, NULL, 0,
                         -1, 1, NULL, NULL, NULL);
      tmp->next=coords;
      coords=tmp;
    }

  /* Fill in the coordinate arrays, Note that 'dsize' is ordered in C
     dimensions, for the WCS conversion, we need to have the dimensions
     ordered in FITS/Fortran order. */
  switch(ndim)
    {
    case 2:
      x=coords->array;  y=coords->next->array;
      x[0] = 1;         y[0] = 1;
      x[1] = dsize[1];  y[1] = 1;
      x[2] = 1;         y[2] = dsize[0];
      x[3] = dsize[1];  y[3] = dsize[0];
      x[4] = dsize[1]/2 + (dsize[1]%2 ? 1 : 0.5f);
      y[4] = dsize[0]/2 + (dsize[0]%2 ? 1 : 0.5f);
      break;
    case 3:
      x=coords->array; y=coords->next->array; z=coords->next->next->array;
      x[0] = 1;        y[0] = 1;              z[0]=1;
      x[1] = dsize[2]; y[1] = 1;              z[1]=1;
      x[2] = 1;        y[2] = dsize[1];       z[2]=1;
      x[3] = dsize[2]; y[3] = dsize[1];       z[3]=1;
      x[4] = 1;        y[4] = 1;              z[4]=dsize[0];
      x[5] = dsize[2]; y[5] = 1;              z[5]=dsize[0];
      x[6] = 1;        y[6] = dsize[1];       z[6]=dsize[0];
      x[7] = dsize[2]; y[7] = dsize[1];       z[7]=dsize[0];
      x[8] = dsize[2]/2 + (dsize[2]%2 ? 1 : 0.5f);
      y[8] = dsize[1]/2 + (dsize[1]%2 ? 1 : 0.5f);
      z[8] = dsize[0]/2 + (dsize[0]%2 ? 1 : 0.5f);
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
            "fix the problem. 'ndim' of %zu is not recognized.",
            __func__, PACKAGE_BUGREPORT, ndim);
    }

  /* For a check:
  printf("IMAGE COORDINATES:\n");
  for(i=0;i<numrows;++i)
    if(ndim==2)
      printf("%-15g%-15g\n", x[i], y[i]);
    else
      printf("%-15g%-15g%-15g\n", x[i], y[i], z[i]);
  */

  /* Convert to the world coordinate system. */
  gal_wcs_img_to_world(coords, wcs, 1);

  /* For a check:
  printf("\nWORLD COORDINATES:\n");
  for(i=0;i<numrows;++i)
    if(ndim==2)
      printf("%-15g%-15g\n", x[i], y[i]);
    else
      printf("%-15g%-15g%-15g\n", x[i], y[i], z[i]);
  */

  /* Get the minimum and maximum values in each dimension. */
  tmp=gal_statistics_minimum(coords);
  min[0] = ((double *)(tmp->array))[0];      gal_data_free(tmp);
  tmp=gal_statistics_maximum(coords);
  max[0] = ((double *)(tmp->array))[0];      gal_data_free(tmp);
  tmp=gal_statistics_minimum(coords->next);
  min[1] = ((double *)(tmp->array))[0];      gal_data_free(tmp);
  tmp=gal_statistics_maximum(coords->next);
  max[1] = ((double *)(tmp->array))[0];      gal_data_free(tmp);
  if(ndim>2)
    {
      tmp=gal_statistics_minimum(coords->next->next);
      min[2] = ((double *)(tmp->array))[0];      gal_data_free(tmp);
      tmp=gal_statistics_maximum(coords->next->next);
      max[2] = ((double *)(tmp->array))[0];      gal_data_free(tmp);
    }

  /* Write the center and width. */
  switch(ndim)
    {
    case 2:
      center[0]=x[4];         center[1]=y[4];
      width[0]=max[0]-min[0]; width[1]=max[1]-min[1];
      break;
    case 3:
      center[0]=x[8];         center[1]=y[8];         center[2]=z[8];
      width[0]=max[0]-min[0]; width[1]=max[1]-min[1]; width[2]=max[2]-min[2];
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to solve the "
            "problem. The value %zu is not a recognized dimension", __func__,
            PACKAGE_BUGREPORT, ndim);
    }

  /* Clean up and return success. */
  free(dsize);
  wcsfree(wcs);
  return 1;
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


  /* For a check.
  {
    size_t i;
    printf("\n\n%s sanity check (%d dimensions):\n", __func__, nelem);
    for(i=0;i<coords->size;++i)
      switch(nelem)
        {
        case 2:
          printf("(%-10g %-10g) --> (%-10g %-10g), [stat: %d]\n",
                 pixcrd[i*2], pixcrd[i*2+1],
                 world[i*2],  world[i*2+1], stat[i]);
          break;
        case 3:
          printf("(%g, %g, %g) --> (%g, %g, %g), [stat: %d]\n",
                 pixcrd[i*3], pixcrd[i*3+1], pixcrd[i*3+2],
                 world[i*3],  world[i*3+1],  world[i*3+2], stat[i]);
          break;
        }
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
