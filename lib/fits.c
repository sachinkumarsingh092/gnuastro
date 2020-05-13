/*********************************************************************
Functions to work with FITS image data.
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
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_version.h>

#include <gnuastro/git.h>
#include <gnuastro/wcs.h>
#include <gnuastro/list.h>
#include <gnuastro/fits.h>
#include <gnuastro/tile.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>
#include <gnuastro-internal/fixedstringmacros.h>










/*************************************************************
 **************        Reporting errors:       ***************
 *************************************************************/
void
gal_fits_io_error(int status, char *message)
{
  char defmessage[]="Error in CFITSIO, see above.";
  if(status)
    {
      fits_report_error(stderr, status);
      if(message)
        error(EXIT_FAILURE, 0, "%s", message);
      else
        error(EXIT_FAILURE, 0, "%s", defmessage);
    }
}




















/*************************************************************
 **************           FITS names           ***************
 *************************************************************/
/* IMPORTANT NOTE: if other compression suffixes are add to this function,
   include them in 'gal_checkset_automatic_output', so the compression
   suffix can be skipped when the user doesn't specify an output
   filename.*/
int
gal_fits_name_is_fits(char *name)
{
  size_t len;

  if(name)
    {
      len=strlen(name);
      if ( (    len>=3 && strcmp(&name[len-3], "fit"     ) == 0 )
           || ( len>=4 && strcmp(&name[len-4], "fits"    ) == 0 )
           || ( len>=7 && strcmp(&name[len-7], "fits.gz" ) == 0 )
           || ( len>=6 && strcmp(&name[len-6], "fits.Z"  ) == 0 )
           || ( len>=3 && strcmp(&name[len-3], "imh"     ) == 0 )
           || ( len>=7 && strcmp(&name[len-7], "fits.fz" ) == 0 ) )
        return 1;
      else
        return 0;
    }
  else return 0;
}





/* IMPORTANT NOTE: if other compression suffixes are add to this function,
   include them in 'gal_checkset_automatic_output', so the compression
   suffix can be skipped when the user doesn't specify an output
   filename.*/
int
gal_fits_suffix_is_fits(char *suffix)
{
  char *nodot;

  if(suffix)
    {
      nodot=suffix[0]=='.' ? (suffix+1) : suffix;
      if ( strcmp(   nodot, "fit"     ) == 0
           || strcmp(nodot, "fits"    ) == 0
           || strcmp(nodot, "fits.gz" ) == 0
           || strcmp(nodot, "fits.Z"  ) == 0
           || strcmp(nodot, "imh"     ) == 0
           || strcmp(nodot, "fits.fz" ) == 0 )
        return 1;
      else
        return 0;
    }
  else return 0;
}





/* If the name is a FITS name, then put a '(hdu: ...)' after it and return
   the string. If it isn't a FITS file, just print the name. Note that the
   space is allocated. */
char *
gal_fits_name_save_as_string(char *filename, char *hdu)
{
  char *name;

  /* Small sanity check. */
  if(filename==NULL)
    gal_checkset_allocate_copy("stdin", &name);
  else
    {
      if( gal_fits_name_is_fits(filename) )
        {
          if( asprintf(&name, "%s (hdu: %s)", filename, hdu)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else gal_checkset_allocate_copy(filename, &name);
    }
  return name;
}




















/*************************************************************
 **************           Type codes           ***************
 *************************************************************/
uint8_t
gal_fits_bitpix_to_type(int bitpix)
{
  switch(bitpix)
    {
    case BYTE_IMG:                  return GAL_TYPE_UINT8;
    case SBYTE_IMG:                 return GAL_TYPE_INT8;
    case USHORT_IMG:                return GAL_TYPE_UINT16;
    case SHORT_IMG:                 return GAL_TYPE_INT16;
    case ULONG_IMG:                 return GAL_TYPE_UINT32;
    case LONG_IMG:                  return GAL_TYPE_INT32;
    case LONGLONG_IMG:              return GAL_TYPE_INT64;
    case FLOAT_IMG:                 return GAL_TYPE_FLOAT32;
    case DOUBLE_IMG:                return GAL_TYPE_FLOAT64;
    default:
      error(EXIT_FAILURE, 0, "%s: bitpix value of %d not recognized",
            __func__, bitpix);
    }
  return 0;
}





int
gal_fits_type_to_bitpix(uint8_t type)
{
  switch(type)
    {
    case GAL_TYPE_UINT8:       return BYTE_IMG;
    case GAL_TYPE_INT8:        return SBYTE_IMG;
    case GAL_TYPE_UINT16:      return USHORT_IMG;
    case GAL_TYPE_INT16:       return SHORT_IMG;
    case GAL_TYPE_UINT32:      return ULONG_IMG;
    case GAL_TYPE_INT32:       return LONG_IMG;
    case GAL_TYPE_INT64:       return LONGLONG_IMG;
    case GAL_TYPE_FLOAT32:     return FLOAT_IMG;
    case GAL_TYPE_FLOAT64:     return DOUBLE_IMG;

    case GAL_TYPE_BIT:
    case GAL_TYPE_STRLL:
    case GAL_TYPE_STRING:
    case GAL_TYPE_UINT64:
    case GAL_TYPE_COMPLEX32:
    case GAL_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, "%s: type %s not recognized for FITS image "
            "BITPIX", __func__, gal_type_name(type, 1));

    default:
      error(EXIT_FAILURE, 0, "%s: type value of %d not recognized",
            __func__, type);
    }
  return 0;
}





/* The values to the TFORM header keywords of FITS binary tables are single
   letter capital letters, but that is useless in identifying the data type
   of the column. So this function will do the conversion based on the
   CFITSIO manual.*/
char
gal_fits_type_to_bin_tform(uint8_t type)
{
  switch(type)
    {
    /* Recognized by CFITSIO. */
    case GAL_TYPE_STRING:      return 'A';
    case GAL_TYPE_BIT:         return 'X';
    case GAL_TYPE_UINT8:       return 'B';
    case GAL_TYPE_INT8:        return 'S';
    case GAL_TYPE_UINT16:      return 'U';
    case GAL_TYPE_INT16:       return 'I';
    case GAL_TYPE_UINT32:      return 'V';
    case GAL_TYPE_INT32:       return 'J';
    case GAL_TYPE_INT64:       return 'K';
    case GAL_TYPE_FLOAT32:     return 'E';
    case GAL_TYPE_FLOAT64:     return 'D';
    case GAL_TYPE_COMPLEX32:   return 'C';
    case GAL_TYPE_COMPLEX64:   return 'M';

    /* Not recognized by CFITSIO. */
    case GAL_TYPE_UINT64:
      error(EXIT_FAILURE, 0, "%s: type %s not recognized for FITS binary "
            "table TFORM", __func__, gal_type_name(type, 1));
      break;

    /* Wrong type code. */
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized", __func__, type);
    }

  error(EXIT_FAILURE, 0, "%s: s bug! Please contact us so we can fix this. "
        "Control must not reach the end of this function", __func__);
  return '\0';
}





int
gal_fits_type_to_datatype(uint8_t type)
{
  int w=0;

  switch(type)
    {
    /* Recognized CFITSIO types. */
    case GAL_TYPE_BIT:              return TBIT;
    case GAL_TYPE_UINT8:            return TBYTE;
    case GAL_TYPE_INT8:             return TSBYTE;
    case GAL_TYPE_FLOAT32:          return TFLOAT;
    case GAL_TYPE_FLOAT64:          return TDOUBLE;
    case GAL_TYPE_COMPLEX32:        return TCOMPLEX;
    case GAL_TYPE_COMPLEX64:        return TDBLCOMPLEX;
    case GAL_TYPE_STRING:           return TSTRING;

    /* Types that depend on the host system. The C standard says that the
       'short', 'int' and 'long' types are ATLEAST 2, 2, 4 bytes, so to be
       safe, we will check all of them for the 32-bit types.*/
    case GAL_TYPE_UINT16:
      w=2;
      if     ( sizeof(short)    == w )   return TUSHORT;
      else if( sizeof(int)      == w )   return TUINT;
      break;

    case GAL_TYPE_INT16:
      w=2;
      if     ( sizeof(short)    == w )   return TSHORT;
      else if( sizeof(int)      == w )   return TINT;
      break;

    /* On 32-bit systems, the length of 'int' and 'long' are both
       32-bits. But CFITSIO's LONG type is preferred because it is designed
       to be 32-bit. Its 'INT' type is not clearly defined and caused
       problems when reading keywords.*/
    case GAL_TYPE_UINT32:
      w=4;
      if     ( sizeof(long)     == w )   return TULONG;
      else if( sizeof(int)      == w )   return TUINT;
      else if( sizeof(short)    == w )   return TUSHORT;
      break;

    /* Similar to UINT32 above. */
    case GAL_TYPE_INT32:
      w=4;
      if     ( sizeof(long)     == w )   return TLONG;
      else if( sizeof(int)      == w )   return TINT;
      else if( sizeof(short)    == w )   return TSHORT;
      break;

    case GAL_TYPE_UINT64:
      w=8;
      if     ( sizeof(long)     == w )   return TULONG;
      break;

    case GAL_TYPE_INT64:
      w=8;
      if     ( sizeof(long)     == w )   return TLONG;
      else if( sizeof(LONGLONG) == w )   return TLONGLONG;
      break;

    /* Wrong type. */
    default:
      error(EXIT_FAILURE, 0, "%s: type code %d is not a recognized",
            __func__, type);
    }

  /* If control reaches, here, there was a problem with the host types. */
  if(w)
    error(EXIT_FAILURE, 0, "%s: this system doesn't have a %d byte integer "
          "type, so type '%s' cannot be written to FITS", __func__, w,
          gal_type_name(type, 1));
  else
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can "
          "fix the problem. Control must not have reached the end for the "
          "given type '%s'", __func__, PACKAGE_BUGREPORT,
          gal_type_name(type, 1));
  return -1;
}





uint8_t
gal_fits_datatype_to_type(int datatype, int is_table_column)
{
  int inttype;

  switch(datatype)
    {
    case TBIT:            return GAL_TYPE_BIT;
    case TBYTE:           return GAL_TYPE_UINT8;
    case TSBYTE:          return GAL_TYPE_INT8;
    case TFLOAT:          return GAL_TYPE_FLOAT32;
    case TDOUBLE:         return GAL_TYPE_FLOAT64;
    case TCOMPLEX:        return GAL_TYPE_COMPLEX32;
    case TDBLCOMPLEX:     return GAL_TYPE_COMPLEX64;
    case TSTRING:         return GAL_TYPE_STRING;

    /* Sizes that depend on the host system. */
    case TUSHORT:
      switch( sizeof(short) )
        {
        case 2:           return GAL_TYPE_UINT16;
        case 4:           return GAL_TYPE_UINT32;
        case 8:           return GAL_TYPE_UINT64;
        }
      break;

    case TSHORT:
      switch( sizeof(short) )
        {
        case 2:           return GAL_TYPE_INT16;
        case 4:           return GAL_TYPE_INT32;
        case 8:           return GAL_TYPE_INT64;
        }
      break;

    case TUINT:
      switch( sizeof(int) )
        {
        case 2:           return GAL_TYPE_UINT16;
        case 4:           return GAL_TYPE_UINT32;
        case 8:           return GAL_TYPE_UINT64;
        }
      break;

    case TINT:
      switch( sizeof(int) )
        {
        case 2:           return GAL_TYPE_INT16;
        case 4:           return GAL_TYPE_INT32;
        case 8:           return GAL_TYPE_INT64;
        }
      break;

    case TULONG:
      switch( sizeof(long) )
        {
        case 4:           return GAL_TYPE_UINT32;
        case 8:           return GAL_TYPE_UINT64;
        }
      break;

    case TLONG: /* ==TINT32BIT when in a table column. */
      if(is_table_column) return GAL_TYPE_INT32;
      else
        switch( sizeof(long) )
          {
          case 4:         return GAL_TYPE_INT32;
          case 8:         return GAL_TYPE_INT64;
          }
      break;

    case TLONGLONG:
      return GAL_TYPE_INT64;
      break;

    /* The TLOGICAL depends on the context: for keywords, it is int32, for
       table columns, it is int8. */
    case TLOGICAL:
      switch( sizeof(int) )
        {
        case 2: inttype=GAL_TYPE_INT16;  break;
        case 4: inttype=GAL_TYPE_INT32;  break;
        case 8: inttype=GAL_TYPE_INT64;  break;
        }
      return is_table_column ? GAL_TYPE_INT8 : inttype;
      break;

    /* A bug! */
    default:
      error(EXIT_FAILURE, 0, "%s: %d is not a recognized CFITSIO datatype",
            __func__, datatype);
    }

  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can fix "
        "this. Control must not have reached the end of this function.",
        __func__, PACKAGE_BUGREPORT);
  return GAL_TYPE_INVALID;
}





/* When there is a BZERO or TZERO and BSCALE or TSCALE keywords, then the
   type that must be used to store the actual values of the pixels may be
   different from the type from BITPIX. This function does the necessary
   corrections.*/
static void
fits_type_correct(int *type, double bscale, char *bzero_str)
{
  double bzero;
  int tofloat=1;
  char *tailptr, *bzero_u64="9223372036854775808";

  /* If bzero_str is given and 'bscale=1.0' the case might be that we are
     dealing with an integer dataset that just needs a different sign. */
  if(bzero_str && bscale==1.0f)
    {
      /* Read the 'bzero' string as a 'double' number. */
      bzero  = strtod(bzero_str, &tailptr);
      if(tailptr==bzero_str)
        error(EXIT_FAILURE, 0, "%s: BZERO value '%s' couldn't be "
              "parsed as a number", __func__, bzero_str);

      /* Work based on type. For the default conversions defined by the
         FITS standard to change the signs of integers, make the proper
         correction, otherwise set the type to float. */
      switch(*type)
        {
        case GAL_TYPE_UINT8:
          if(bzero == -128.0f)      { *type = GAL_TYPE_INT8;   tofloat=0; }
          break;

        case GAL_TYPE_INT16:
          if(bzero == 32768)        { *type = GAL_TYPE_UINT16; tofloat=0; }
          break;

        case GAL_TYPE_INT32:
          if(bzero == 2147483648LU) { *type = GAL_TYPE_UINT32; tofloat=0; }
          break;

        /* The 'bzero' value for unsigned 64-bit integers has 19 decimal
           digits, but a 64-bit floating point ('double' type) can only
           safely store 15 decimal digits. As a result, the safest way to
           check the 'bzero' value for this type is to compare it as a
           string. But all integers nearby (for example
           '9223372036854775807') get rounded to this same value (when
           stored as 'double'). So we will also check the parsed number and
           if it equals this number, a warning will be printed. */
        case GAL_TYPE_INT64:
          if( !strcmp(bzero_str, bzero_u64) )
            {*type = GAL_TYPE_UINT64; tofloat=0;}
          else
            if( bzero == 9223372036854775808LLU )
              {
                fprintf(stderr, "\nWARNING in %s: the BZERO header keyword "
                        "value ('%s') is very close (but not exactly equal) "
                        "to '%s'. The latter value in the FITS standard is "
                        "used to identify that the dataset should be read as "
                        "unsigned 64-bit integers instead of signed 64-bit "
                        "integers. Depending on your version of CFITSIO, "
                        "it may be read as a signed or unsigned 64-bit "
                        "integer array\n\n", __func__, bzero_str, bzero_u64);
                tofloat=0;
              }
          break;

          /* For the other types (when 'BSCALE=1.0f'), currently no
             correction is necessary, maybe later we can check if the
             scales are integers and set the integer output type to the
             smallest type that can allow the scaled values. */
        default: tofloat=0;
        }
    }

  /* If the type must be a float, then do the conversion. */
  if(tofloat) *type=GAL_TYPE_FLOAT32;
}




















/**************************************************************/
/**********                  HDU                   ************/
/**************************************************************/
fitsfile *
gal_fits_open_to_write(char *filename)
{
  int status=0;
  long naxes=0;
  fitsfile *fptr;

  /* When the file exists just open it. Otherwise, create the file. But we
     want to leave the first extension as a blank extension and put the
     image in the next extension to be consistent between tables and
     images. */
  if(access(filename,F_OK) == -1 )
    {
      /* Create the file. */
      if( fits_create_file(&fptr, filename, &status) )
        gal_fits_io_error(status, NULL);

      /* Create blank extension. */
      if( fits_create_img(fptr, BYTE_IMG, 0, &naxes, &status) )
        gal_fits_io_error(status, NULL);

      /* Close the blank extension. */
      if( fits_close_file(fptr, &status) )
        gal_fits_io_error(status, NULL);
    }

  /* Open the file, ready for later steps. */
  if( fits_open_file(&fptr, filename, READWRITE, &status) )
    gal_fits_io_error(status, NULL);

  /* Return the pointer. */
  return fptr;
}





size_t
gal_fits_hdu_num(char *filename)
{
  fitsfile *fptr;
  int num, status=0;

  /* We don't need to check for an error everytime, because we don't
     make any non CFITSIO usage of the output. It is necessary to
     check every CFITSIO call, only when you will need to use the
     outputs. */
  fits_open_file(&fptr, filename, READONLY, &status);

  fits_get_num_hdus(fptr, &num, &status);

  fits_close_file(fptr, &status);

  gal_fits_io_error(status, NULL);

  return num;
}





/* Given the filename and HDU, this function will return the CFITSIO code
   for the type of data it contains (table, or image). The CFITSIO codes
   are:

       IMAGE_HDU:    An image HDU.
       ASCII_TBL:    An ASCII table HDU.
       BINARY_TBL:   BINARY TABLE HDU.       */
int
gal_fits_hdu_format(char *filename, char *hdu)
{
  fitsfile *fptr;
  int hdutype, status=0;

  /* Open the HDU. */
  fptr=gal_fits_hdu_open(filename, hdu, READONLY);

  /* Check the type of the given HDU: */
  if (fits_get_hdu_type(fptr, &hdutype, &status) )
    gal_fits_io_error(status, NULL);

  /* Clean up and return.. */
  if( fits_close_file(fptr, &status) )
    gal_fits_io_error(status, NULL);
  return hdutype;
}





/* Open a given HDU and return the FITS pointer. 'iomode' determines how
   the FITS file will be opened: only to read or to read and write. You
   should use the macros given by the CFITSIO header:

     READONLY:   read-only.
     READWRITE:  read and write.         */
fitsfile *
gal_fits_hdu_open(char *filename, char *hdu, int iomode)
{
  int status=0;
  char *ffname;
  fitsfile *fptr;

  /* Add hdu to filename: */
  if( asprintf(&ffname, "%s[%s#]", filename, hdu)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);

  /* Open the FITS file: */
  if( fits_open_file(&fptr, ffname, iomode, &status) )
    {
      switch(status)
        {
        /* Since the default HDU is '1', when the file only has one
           extension, this error is common, so we will put a special
           notice. */
        case END_OF_FILE:
          if( hdu[0]=='1' && hdu[1]=='\0' )
            error(EXIT_FAILURE, 0, "%s: only has one HDU.\n\n"
                  "You should tell Gnuastro's command-line "
                  "programs to look for data in the primary HDU with the "
                  "'--hdu=0' option (or '-h0'). For library users, you can "
                  "put \"0\" (a string literal) for the function's HDU "
                  "argument. For more, see the FOOTNOTE below.\n\n"
                  "Pro TIP: if your desired HDU has a name (value to "
                  "'EXTNAME' keyword), it is best to just use that name "
                  "with '--hdu' instead of relying on a counter. You can "
                  "see the list of HDUs in a FITS file (with their data "
                  "format, type, size and possibly HDU name) using "
                  "Gnuastro's 'astfits' program, for example:\n\n"
                  "    astfits %s\n\n"
                  "FOOTNOTE -- When writing a new FITS file, Gnuastro leaves "
                  "the pimary HDU only for metadata. The output datasets "
                  "(tables, images or cubes) are written after the primary "
                  "HDU. In this way the keywords of the the first HDU can be "
                  "used as metadata of the whole file (which may contain many "
                  "extensions, this is stipulated in the FITS standard). "
                  "Usually the primary HDU keywords contains the option names "
                  "and values that the program was run with. Because of this, "
                  "Gnuastro's default HDU to read data in a FITS file is the "
                  "second (or '--hdu=1'). This error is commonly caused when "
                  "the FITS file wasn't created by Gnuastro or by a program "
                  "respecting this convention.", filename, filename);
          break;

        /* Generic error below is fine for this case */
        case BAD_HDU_NUM:
          break;

        /* In case an un-expected error occurs, use the general CFITSIO
           reporting that we have already prepared. */
        default:
          gal_fits_io_error(status, "opening the given extension/HDU in "
                            "the given file");
        }

      error(EXIT_FAILURE, 0, "%s: extension/HDU '%s' doesn't exist. Please "
            "run the following command to see the extensions/HDUs in "
            "'%s':\n\n"
            "    $ astfits %s\n\n"
            "The respective HDU number (or name, when present) may be used "
            "with the '--hdu' option in Gnuastro's programs (or the 'hdu' "
            "argument in Gnuastro's libraries) to open the respective HDU.",
            filename, hdu, filename, filename);
    }

  /* Clean up and the pointer. */
  free(ffname);
  return fptr;
}





/* Check the desired HDU in a FITS image and also if it has the
   desired type. */
fitsfile *
gal_fits_hdu_open_format(char *filename, char *hdu, int img0_tab1)
{
  fitsfile *fptr;
  int status=0, hdutype;

  /* A small sanity check. */
  if(hdu==NULL)
    error(EXIT_FAILURE, 0, "no HDU specified for %s", filename);

  /* Open the HDU. */
  fptr=gal_fits_hdu_open(filename, hdu, READONLY);

  /* Check the type of the given HDU: */
  if (fits_get_hdu_type(fptr, &hdutype, &status) )
    gal_fits_io_error(status, NULL);

  /* Check if the type of the HDU is the expected type. We could have
     written these as && conditions, but this is easier to read, it makes
     no meaningful difference to the compiler. */
  if(img0_tab1)
    {
      if(hdutype==IMAGE_HDU)
        error(EXIT_FAILURE, 0, "%s (hdu: %s): is not a table",
              filename, hdu);
    }
  else
    {
      if(hdutype!=IMAGE_HDU)
        error(EXIT_FAILURE, 0, "%s (hdu: %s): not an image",
              filename, hdu);
    }

  /* Clean up and return. */
  return fptr;
}





















/**************************************************************/
/**********            Header keywords             ************/
/**************************************************************/
/* Return allocated pointer to the blank value to use in a FITS file header
   keyword. */
void *
gal_fits_key_img_blank(uint8_t type)
{
  void *out=NULL, *tocopy=NULL;
  uint8_t u8=0;          /* Equivalent of minimum in signed   with BZERO. */
  int16_t s16=INT16_MAX; /* Equivalend of maximum in unsigned with BZERO. */
  int32_t s32=INT32_MAX; /* Equivalend of maximum in unsigned with BZERO. */
  int64_t s64=INT64_MAX; /* Equivalend of maximum in unsigned with BZERO. */

  switch(type)
    {
    /* Types with no special treatment: */
    case GAL_TYPE_BIT:
    case GAL_TYPE_UINT8:
    case GAL_TYPE_INT16:
    case GAL_TYPE_INT32:
    case GAL_TYPE_INT64:
    case GAL_TYPE_FLOAT32:
    case GAL_TYPE_FLOAT64:
    case GAL_TYPE_COMPLEX32:
    case GAL_TYPE_COMPLEX64:
    case GAL_TYPE_STRING:
    case GAL_TYPE_STRLL:
      out = gal_blank_alloc_write(type);
      break;

    /* Types that need values from their opposite-signed types. */
    case GAL_TYPE_INT8:      tocopy=&u8;      break;
    case GAL_TYPE_UINT16:    tocopy=&s16;     break;
    case GAL_TYPE_UINT32:    tocopy=&s32;     break;
    case GAL_TYPE_UINT64:    tocopy=&s64;     break;

    default:
      error(EXIT_FAILURE, 0, "%s: type code %u not recognized as a Gnuastro "
            "data type", __func__, type);
    }

  /* If 'gal_blank_alloc_write' wasn't used (copy!=NULL), then allocate the
     necessary space and fill it in. Note that the width of the signed and
     unsigned values doesn't differ, so we can use the actual input
     type. */
  if(tocopy)
    {
      out = gal_pointer_allocate(type, 1, 0, __func__, "out");
      memcpy(out, tocopy, gal_type_sizeof(type));
    }

  /* Return. */
  return out;
}





/* CFITSIO doesn't remove the two single quotes around the string value, so
   the strings it reads are like: 'value ', or 'some_very_long_value'. To
   use the value, it is commonly necessary to remove the single quotes (and
   possible extra spaces). This function will modify the string in its own
   allocated space. You can use this to later free the original string (if
   it was allocated). */
void
gal_fits_key_clean_str_value(char *string)
{
  int end;       /* Has to be int because size_t is always >=0. */
  char *c, *cf;

  /* Start from the second last character (the last is a single quote) and
     go down until you hit a non-space character. This will also work when
     there is no space characters between the last character of the value
     and ending single-quote: it will be set to '\0' after this loop. */
  for(end=strlen(string)-2;end>=0;--end)
    if(string[end]!=' ')
      break;

  /* Shift all the characters after the first one (which is a ''' back by
     one and put the string ending characters on the 'end'th element. */
  cf=(c=string)+end; do *c=*(c+1); while(++c<cf);
  *cf='\0';
}




/* Fill the 'tm' structure (defined in 'time.h') with the values derived
   from a FITS format date-string and return the (optional) sub-second
   information as a double.

   The basic FITS string is defined under the 'DATE' keyword in the FITS
   standard. For the more complete format which includes timezones, see the
   W3 standard: https://www.w3.org/TR/NOTE-datetime */
char *
gal_fits_key_date_to_struct_tm(char *fitsdate, struct tm *tp)
{
  int hasT=0, hassq=0, usesdash=0, usesslash=0, hasZ=0;
  char *C, *cc, *c=NULL, *cf, *subsec=NULL, *nosubsec=fitsdate;

  /* Initialize the 'tm' structure to all-zero elements. In particular, The
     FITS standard times are written in UTC, so, the time zone ('tm_zone'
     element, which specifies number of seconds to shift for the time zone)
     has to be zero. The day-light saving flag ('isdst' element) also has
     to be set to zero. */
  tp->tm_sec=tp->tm_min=tp->tm_hour=tp->tm_mday=tp->tm_mon=tp->tm_year=0;
  tp->tm_wday=tp->tm_yday=tp->tm_isdst=tp->tm_gmtoff=0;
  tp->tm_zone=NULL;

  /* According to the FITS standard the 'T' in the middle of the date and
     time of day is optional (the time is not mandatory). */
  cf=(c=fitsdate)+strlen(fitsdate);
  do
    switch(*c)
      {
      case 'T':  hasT=1;      break; /* With 'T' HH:MM:SS are defined.    */
      case '-':  usesdash=1;  break; /* Day definition: YYYY-MM-DD.       */
      case '/':  usesslash=1; break; /* Day definition(old): DD/MM/YY.    */
      case '\'': hassq=1;     break; /* Wholly Wrapped in a single-quote. */
      case 'Z':  hasZ=1;      break; /* When ends in 'Z', means UTC. See  */
                                   /* https://www.w3.org/TR/NOTE-datetime */

      /* In case we have sub-seconds in the string, we need to remove it
         because 'strptime' doesn't recognize sub-second resolution.*/
      case '.':
        /* Allocate space (by copying the remaining full string and adding
           a '\0' where necessary. */
        gal_checkset_allocate_copy(c, &subsec);
        gal_checkset_allocate_copy(fitsdate, &nosubsec);

        /* Parse the sub-second part and remove it from the copy. */
        cc=nosubsec+(c-fitsdate);
        for(C=subsec+1;*C!='\0';C++)
          if(!isdigit(*C)) {*cc++=*C; *C='\0';}
        *cc='\0';
      }
  while(++c<cf);

  /* Convert this date into seconds since 1970/01/01, 00:00:00. */
  c = ( (usesdash==0 && usesslash==0)
        ? NULL
        : ( usesdash
            ? ( hasT
                ? ( hasZ
                    ? strptime(nosubsec, hassq?"'%FT%TZ'":"%FT%TZ", tp)
                    : strptime(nosubsec, hassq?"'%FT%T'":"%FT%T", tp) )
                : strptime(nosubsec, hassq?"'%F'"   :"%F"   , tp))
            : ( hasT
                ? ( hasZ
                    ? strptime(nosubsec, hassq?"'%d/%m/%yT%TZ'":"%d/%m/%yT%TZ", tp)
                    : strptime(nosubsec, hassq?"'%d/%m/%yT%T'":"%d/%m/%yT%T", tp))
                : strptime(nosubsec, hassq?"'%d/%m/%y'"   :"%d/%m/%y"   , tp)
                )
            )
        );

  /* The value might have sub-seconds. In that case, 'c' will point to a
     '.' and we'll have to parse it as double. */
  if( c==NULL || (*c!='.' && *c!='\0') )
    error(EXIT_FAILURE, 0, "'%s' isn't in the FITS date format.\n\n"
          "According to the FITS standard, the date must be in one of "
          "these formats:\n"
          "   YYYY-MM-DD\n"
          "   YYYY-MM-DDThh:mm:ss\n"
          "   YYYY-MM-DDThh:mm:ssZ   (Note the 'Z',  see *) \n"
          "   DD/MM/YY               (Note the 'YY', see ^)\n"
          "   DD/MM/YYThh:mm:ss\n"
          "   DD/MM/YYThh:mm:ssZ\n\n"
          "[*]: The 'Z' is interpreted as being in the UTC Timezone.\n"
          "[^]: Gnuastro's FITS library (this program), interprets the "
          "older (two character for year) format, year values 68 to 99 as "
          "the years 1969 to 1999 and values 0 to 68 as the years 2000 to "
          "2068.", fitsdate);

  /* If the subseconds were removed (and a new string was allocated), free
     that extra new string. */
  if(nosubsec!=fitsdate) free(nosubsec);

  /* Return the subsecond value. */
  return subsec;
}





/* Convert the FITS standard date format (as a string, already read from
   the keywords) into number of seconds since 1970/01/01, 00:00:00. Very
   useful to avoid calendar issues like number of days in a different
   months or leap years and etc. The remainder of the format string
   (sub-seconds) will be put into the two pointer arguments: 'subsec' is in
   double-precision floating point format and  */
size_t
gal_fits_key_date_to_seconds(char *fitsdate, char **subsecstr,
                             double *subsec)
{
  time_t t;
  char *tmp;
  struct tm tp;
  void *outptr=subsec;

  /* Fill in the 'tp' elements with values read from the string. */
  tmp=gal_fits_key_date_to_struct_tm(fitsdate, &tp);

  /* If the user cared about the remainder (sub-second string), then set it
     and convert it to a double type. */
  if(subsecstr && tmp)
    {
      /* Set the output pointer. */
      *subsecstr=tmp;

      /* Convert the remainder string to double-precision floating point
         (if the given pointer isn't NULL). */
      if(subsec)
        if( gal_type_from_string(&outptr, tmp, GAL_TYPE_FLOAT64) )
          error(EXIT_FAILURE, 0, "%s: the sub-second portion of '%s' (or "
                "'%s') couldn't be read as a number", __func__, fitsdate,
                tmp);
    }

  /* Convert the 'tm' structure to 'time_t'. Note that the system's
     timezone and daylight saving need to be subtracted from the output of
     'mktime'. Otherwise the result will be different on different
     host-system timezones (which is not what we want here: bug #57995). */
  t=mktime(&tp)-timezone-daylight;

  /* Return the value and set the output pointer. */
  return (size_t)t;
}





/* Read the keyword values from a FITS pointer. The input should be a
   linked list of 'gal_data_t'. Before calling this function, you just have
   to set the 'name' and desired 'type' values of each element in the list
   to the keyword you want it to keep the value of. The given 'name' value
   will be directly passed to CFITSIO to read the desired keyword. This
   function will allocate space to keep the value. Here is one example of
   using this function:

      gal_data_t *keysll=gal_data_array_calloc(N);

      for(i=0;i<N-2;++i) keysll[i]->next=keysll[i+1];

      \\ Put a name and type for each element.

      gal_fits_key_read_from_ptr(fptr, keysll, 0, 0);

      \\ use the values as you like.

      gal_data_array_free(keysll, N, 1);

   If the 'array' pointer of each keyword's dataset is not NULL, then it is
   assumed that the space has already been allocated. If it is NULL, then
   space will be allocated internally here.

   Strings need special consideration: the reason is that generally,
   'gal_data_t' needs to also allow for array of strings (as it supports
   arrays of integers for example). Hence two allocations will be done here
   (one if 'array!=NULL') and 'keysll[i].array' must be interpretted as
   'char **': one allocation for the pointer, one for the actual
   characters. You don't have to worry about the freeing,
   'gal_data_array_free' will free both allocations. So to read a string,
   one easy way would be the following:

      char *str, **strarray;
      strarr = keysll[i].array;
      str    = strarray[0];

   If CFITSIO is unable to read a keyword for any reason the 'status'
   element of the respective 'gal_data_t' will be non-zero. You can check
   the successful reading of the keyword from the 'status' value in each
   keyword's 'gal_data_t'. If it is zero, then the keyword was found and
   succesfully read. Otherwise, it a CFITSIO status value. You can use
   CFITSIO's error reporting tools or 'gal_fits_io_error' for reporting the
   reason. A tip: when the keyword doesn't exist, then CFITSIO's status
   value will be 'KEY_NO_EXIST'.

   CFITSIO will start searching for the keywords from the last place in the
   header that it searched for a keyword. So it is much more efficient if
   the order that you ask for keywords is based on the order they are
   stored in the header.
 */
void
gal_fits_key_read_from_ptr(fitsfile *fptr, gal_data_t *keysll,
                           int readcomment, int readunit)
{
  void *valueptr;
  char **strarray;
  gal_data_t *tmp;

  /* Get the desired keywords. */
  for(tmp=keysll;tmp!=NULL;tmp=tmp->next)
    if(tmp->name)
      {
        /* Initialize the status: */
        tmp->status=0;

        /* For each keyword, this function stores one value currently. So
           set the size and ndim to 1. But first allocate dsize if it
           wasn't already allocated. */
        if(tmp->dsize==NULL)
          tmp->dsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, 1, 0, __func__,
                                          "tmp->dsize");
        tmp->ndim=tmp->size=tmp->dsize[0]=1;

        /* When the type is a string, 'tmp->array' is an array of pointers
           to a separately allocated piece of memory. So we have to
           allocate that space here. If its not a string, then the
           allocated space above is enough to keep the value.*/
        switch(tmp->type)
          {
          case GAL_TYPE_STRING:
            tmp->array=strarray=( tmp->array
                                  ? tmp->array
                                  : gal_pointer_allocate(tmp->type, 1, 0,
                                                         __func__,
                                                         "tmp->array") );
            errno=0;
            valueptr=strarray[0]=malloc(FLEN_VALUE * sizeof *strarray[0]);
            if(strarray[0]==NULL)
              error(EXIT_FAILURE, errno, "%s: %zu bytes for strarray[0]",
                    __func__, FLEN_VALUE * sizeof *strarray[0]);
            break;

          default:
            tmp->array=valueptr=( tmp->array
                                  ? tmp->array
                                  : gal_pointer_allocate(tmp->type, 1, 0,
                                                         __func__,
                                                         "tmp->array") );
          }

        /* Allocate space for the keyword comment if necessary.*/
        if(readcomment)
          {
            errno=0;
            tmp->comment=calloc(FLEN_COMMENT, sizeof *tmp->comment);
            if(tmp->comment==NULL)
              error(EXIT_FAILURE, errno, "%s: %zu bytes for tmp->comment",
                    __func__, FLEN_COMMENT * sizeof *tmp->comment);
          }
        else
          tmp->comment=NULL;

        /* Allocate space for the keyword unit if necessary. Note that
           since there is no precise CFITSIO length for units, we will use
           the 'FLEN_COMMENT' length for units too (theoretically, the unit
           might take the full remaining area in the keyword). Also note
           that the unit is only optional, so it needs a separate CFITSIO
           function call which is done here.*/
        if(readunit)
          {
            /* Allocate space for the unit and read it in. */
            errno=0;
            tmp->unit=calloc(FLEN_COMMENT, sizeof *tmp->unit);
            if(tmp->unit==NULL)
              error(EXIT_FAILURE, errno, "%s: %zu bytes for tmp->unit",
                    __func__, FLEN_COMMENT * sizeof *tmp->unit);
            fits_read_key_unit(fptr, tmp->name, tmp->unit, &tmp->status);

            /* If the string is empty, free the space and set it to NULL. */
            if(tmp->unit[0]=='\0') {free(tmp->unit); tmp->unit=NULL;}
          }
        else
          tmp->unit=NULL;

        /* Read the keyword and place its value in the pointer. */
        fits_read_key(fptr, gal_fits_type_to_datatype(tmp->type),
                      tmp->name, valueptr, tmp->comment, &tmp->status);

        /* If the comment was empty, free the space and set it to zero. */
        if(tmp->comment && tmp->comment[0]=='\0')
          {free(tmp->comment); tmp->comment=NULL;}

        /* Strings need to be cleaned (CFITSIO puts ''' around them with
           some (possiblly) extra space on the two ends of the string. */
      }
}





/* Same as 'gal_fits_read_keywords_fptr', but accepts the filename and HDU
   as input instead of an already opened CFITSIO 'fitsfile' pointer. */
void
gal_fits_key_read(char *filename, char *hdu, gal_data_t *keysll,
                  int readcomment, int readunit)
{
  size_t len;
  int status=0;
  char *ffname;
  fitsfile *fptr;

  /* Add hdu to filename: */
  errno=0;
  len=strlen(filename)+strlen(hdu)+4;
  ffname=malloc(len*sizeof *ffname);
  if(ffname==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu characters", __func__, len);
  sprintf(ffname, "%s[%s#]", filename, hdu);

  /* Open the FITS file: */
  if( fits_open_file(&fptr, ffname, READONLY, &status) )
    gal_fits_io_error(status, "reading this FITS file");

  /* Read the keywords. */
  gal_fits_key_read_from_ptr(fptr, keysll, readcomment, readunit);

  /* Close the FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);

  /* Clean up. */
  free(ffname);
}





/* Add on keyword to the list of header keywords that need to be added
   to a FITS file. In the end, the keywords will have to be freed, so
   it is important to know before hand if they were allocated or
   not. If not, they don't need to be freed.

   NOTE FOR STRINGS: the value should be the pointer to the string its-self
   (char *), not a pointer to a pointer (char **). */
void
gal_fits_key_list_add(gal_fits_list_key_t **list, uint8_t type,
                      char *keyname, int kfree, void *value, int vfree,
                      char *comment, int cfree, char *unit)
{
  gal_fits_list_key_t *newnode;

  /* Allocate space for the new node and fill it in. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);
  newnode->type=type;
  newnode->keyname=keyname;
  newnode->value=value;
  newnode->comment=comment;
  newnode->unit=unit;
  newnode->kfree=kfree;                /* Free pointers after using them. */
  newnode->vfree=vfree;
  newnode->cfree=cfree;

  newnode->next=*list;
  *list=newnode;
}





void
gal_fits_key_list_add_end(gal_fits_list_key_t **list, uint8_t type,
                          char *keyname, int kfree, void *value, int vfree,
                          char *comment, int cfree, char *unit)
{
  gal_fits_list_key_t *newnode, *tmp;

  /* Allocate space for the new node and fill it in. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocation of new node", __func__);
  newnode->type=type;
  newnode->keyname=keyname;
  newnode->value=value;
  newnode->comment=comment;
  newnode->unit=unit;
  newnode->kfree=kfree;            /* Free pointers after using them. */
  newnode->vfree=vfree;
  newnode->cfree=cfree;

  if(*list)         /* List is already full, add this node to the end */
    {
      /* After this line, tmp points to the last node. */
      tmp=*list; while(tmp->next!=NULL) tmp=tmp->next;
      tmp->next=newnode;
      newnode->next=NULL;
    }
  else                 /* List is empty */
    {
      newnode->next=*list;
      *list=newnode;
    }
}





void
gal_fits_key_list_reverse(gal_fits_list_key_t **list)
{
  gal_fits_list_key_t *in=*list, *tmp=in, *reversed=NULL;

  /* Only do the reversal if there is more than one element. */
  if(in && in->next)
    {
      /* Fill the 'reversed' list. */
      while(in!=NULL)
        {
          tmp=in->next;
          in->next=reversed;
          reversed=in;
          in=tmp;
        }

      /* Write the reversed list into the input pointer. */
      *list=reversed;
    }
}





/* Write a blank keyword field and a title under it in the specified FITS
   file pointer. */
void
gal_fits_key_write_title_in_ptr(char *title, fitsfile *fptr)
{
  size_t i;
  int status=0;
  char *cp, *cpf, blankrec[80], titlerec[80];

  /* Just in case title is 'NULL'. */
  if(title)
    {
      /* A small sanity check. */
      if( strlen(title) + strlen(GAL_FITS_KEY_TITLE_START) > 78 )
        fprintf(stderr, "%s: FITS keyword title '%s' is too long to be fully "
                "included in the keyword record (80 characters, where the "
                "title is prefixed with %zu space characters)",
                __func__, title, strlen(GAL_FITS_KEY_TITLE_START));

      /* Set the last element of the blank array. */
      cpf=blankrec+79;
      *cpf='\0';
      titlerec[79]='\0';
      cp=blankrec; do *cp=' '; while(++cp<cpf);

      /* Print the blank lines before the title */
      if(fits_write_record(fptr, blankrec, &status))
        gal_fits_io_error(status, NULL);

      /* Print the title */
      sprintf(titlerec, "%s%s", GAL_FITS_KEY_TITLE_START, title);
      for(i=strlen(titlerec);i<79;++i)
        titlerec[i]=' ';
      if(fits_write_record(fptr, titlerec, &status))
        gal_fits_io_error(status, NULL);
    }
}





/* Write a filename into keyword values. */
void
gal_fits_key_write_filename(char *keynamebase, char *filename,
                            gal_fits_list_key_t **list, int top1end0)
{
  char *keyname, *value;
  size_t numkey=1, maxlength;
  size_t i, j, len=strlen(filename), thislen;

  /* When you give string arguments, CFITSIO puts them within two ''s,
     so the actual length available is two less. It seems this length
     also didn't include the null character, so, ultimately you have
     to take three from it.*/
  maxlength=FLEN_VALUE-3;

  /* Parse the filename and see when it is necessary to break it. */
  i=0;
  while(i<len)
    {
      /* Set the keyname: */
      errno=0;
      keyname=malloc(FLEN_KEYWORD);
      if(keyname==NULL)
        error(EXIT_FAILURE, errno, "%s: %d bytes for 'keyname'", __func__,
              FLEN_KEYWORD);
      if(len<maxlength)
        strcpy(keyname, keynamebase);
      else
        sprintf(keyname, "%s_%zu", keynamebase, numkey++);

      /* Set the keyword value: */
      errno=0;
      thislen=strlen(&filename[i]);
      value=malloc(maxlength);
      if(value==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes", __func__,
              thislen);
      strncpy(value, &filename[i], maxlength);

      /* If the FROM string (=&filename[i]) in strncpy is shorter than
         SIZE (=maxlength), then the rest of the space will be filled
         with null characters. So we can use this to check if the full
         length was copied. */
      if(value[maxlength-1]=='\0')
        {
          if(top1end0)
            gal_fits_key_list_add(list, GAL_TYPE_STRING, keyname, 1,
                                  value, 1, NULL, 0, NULL);
          else
            gal_fits_key_list_add_end(list, GAL_TYPE_STRING, keyname, 1,
                                      value, 1, NULL, 0, NULL);
          break;
        }
      else
        {
          /* Find the last place in the copied array that contains a
             '/' and put j on the next character (so it can be turned
             into a null character.*/
          for(j=maxlength-2;j>0;--j)
            if(value[j]=='/')
              {
                value[j+1]='\0';
                break;
              }
          if(j==0)
            {
              /* So the loop doesn't continue after this. */
              i=len+1;

              /* There will only be one keyword, so we'll just use the
                 basename. */
              strcpy(keyname, keynamebase);

              /* Let the user know that  */
              error(0,0, "%s: WARNING: '%s' is too long to fit "
                    "into a FITS keyword value (max of %zu characters), "
                    "it will be truncated", __func__, filename,
                    maxlength);
            }

          /* Convert the last useful character and save the file name.*/
          if(top1end0)
            gal_fits_key_list_add(list, GAL_TYPE_STRING, keyname, 1,
                                  value, 1, NULL, 0, NULL);
          else
            gal_fits_key_list_add_end(list, GAL_TYPE_STRING, keyname, 1,
                                      value, 1, NULL, 0, NULL);
          i+=j+1;
        }
    }
}





/* Write the WCS header string into a FITS files*/
void
gal_fits_key_write_wcsstr(fitsfile *fptr, char *wcsstr, int nkeyrec)
{
  size_t i;
  int status=0;
  char *keystart;

  /* Write the title. */
  gal_fits_key_write_title_in_ptr("World Coordinate System (WCS)", fptr);

  /* Write the keywords one by one: */
  for(i=0;i<nkeyrec;++i)
    {
      /* Set the start of this header. */
      keystart=&wcsstr[i*80];

      /* Write it if it isn't blank (first character is a space), or not a
         comment (first 7 characters equal to 'COMMENT'). The reason is
         that WCSLIB adds a blank line and a 'COMMENT' keyword saying its
         own version. But Gnuastro writes the version of WCSLIB as a
         separate keyword along with all other important software, so it is
         redundant and just makes the keywrods hard to read by eye.*/
      if( keystart[0]!=' ' && strncmp(keystart, "COMMENT", 7) )
        fits_write_record(fptr, keystart, &status);
    }
  gal_fits_io_error(status, NULL);
}





/* Write the given list of header keywords into the specified HDU of the
   specified FITS file. */
void
gal_fits_key_write(gal_fits_list_key_t **keylist, char *title,
                   char *filename, char *hdu)
{
  int status=0;
  fitsfile *fptr=gal_fits_hdu_open(filename, hdu, READWRITE);

  /* Write the title */
  gal_fits_key_write_title_in_ptr(title, fptr);

  /* Write the keywords into the  */
  gal_fits_key_write_in_ptr(keylist, fptr);

  /* Close the input FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}





/* Write the keywords in the gal_fits_list_key_t linked list to the FITS
   file. Every keyword that is written is freed, that is why we need the
   pointer to the linked list (to correct it after we finish). */
void
gal_fits_key_write_in_ptr(gal_fits_list_key_t **keylist, fitsfile *fptr)
{
  int status=0;
  gal_fits_list_key_t *tmp, *ttmp;

  tmp=*keylist;
  while(tmp!=NULL)
    {
      /* Write the basic key value and comments. */
      if(tmp->value)
        {
          if( fits_update_key(fptr, gal_fits_type_to_datatype(tmp->type),
                              tmp->keyname, tmp->value, tmp->comment,
                              &status) )
            gal_fits_io_error(status, NULL);
        }
      else
        {
          if(fits_update_key_null(fptr, tmp->keyname, tmp->comment, &status))
            gal_fits_io_error(status, NULL);
        }

      /* Write the units if it was given. */
      if( tmp->unit
          && fits_write_key_unit(fptr, tmp->keyname, tmp->unit, &status) )
        gal_fits_io_error(status, NULL);

      /* Free the value pointer if desired: */
      if(tmp->kfree) free(tmp->keyname);
      if(tmp->vfree) free(tmp->value);
      if(tmp->cfree) free(tmp->comment);

      /* Keep the pointer to the next keyword and free the allocated
         space for this keyword.*/
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }

  /* Set it to NULL so it isn't mistakenly used later. */
  *keylist=NULL;
}





/* Write the given list of header keywords and version info into the give
   file. */
void
gal_fits_key_write_version(gal_fits_list_key_t **keylist, char *title,
                           char *filename, char *hdu)
{
  int status=0;
  fitsfile *fptr=gal_fits_hdu_open(filename, hdu, READWRITE);

  /* Write the given keys followed by the versions. */
  gal_fits_key_write_version_in_ptr(keylist, title, fptr);

  /* Close the input FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}





void
gal_fits_key_write_version_in_ptr(gal_fits_list_key_t **keylist, char *title,
                                  fitsfile *fptr)
{
  int status=0;
  char *gitdescribe;
  char cfitsioversion[20];

  /* Before WCSLIB 5.0, the wcslib_version function was not
     defined. Sometime in the future were everyone has moved to more
     recent versions of WCSLIB, we can remove this macro and its check
     in configure.ac.*/
#if GAL_CONFIG_HAVE_WCSLIB_VERSION == 1
  int wcslibvers[3];
  char wcslibversion[20];
  const char *wcslibversion_const;
#endif

  /* Small sanity check. */
  if(fptr==NULL)
    error(EXIT_FAILURE, 0, "%s: input FITS pointer is NULL", __func__);

  /* If any header keywords are specified, add them: */
  if(keylist && *keylist)
    {
      gal_fits_key_write_title_in_ptr(title?title:"Specific keys", fptr);
      gal_fits_key_write_in_ptr(keylist, fptr);
    }

  /* Print 'Versions and date' title. */
  gal_fits_key_write_title_in_ptr("Versions and date", fptr);

  /* Set the version of CFITSIO as a string. */
  sprintf(cfitsioversion, "%-.2f", CFITSIO_VERSION);

  /* Write all the information: */
  fits_write_date(fptr, &status);

  /* Write the version of CFITSIO */
  fits_update_key(fptr, TSTRING, "CFITSIO", cfitsioversion,
                  "CFITSIO version.", &status);

  /* Write the WCSLIB version. */
#if GAL_CONFIG_HAVE_WCSLIB_VERSION == 1
  wcslibversion_const=wcslib_version(wcslibvers);
  strcpy(wcslibversion, wcslibversion_const);
  fits_update_key(fptr, TSTRING, "WCSLIB", wcslibversion,
                  "WCSLIB version.", &status);
#endif

  /* Write the GSL version. */
  fits_update_key(fptr, TSTRING, "GSL", GSL_VERSION,
                  "GNU Scientific Library version.", &status);

  /* Write the Gnuastro version. */
  fits_update_key(fptr, TSTRING, "GNUASTRO", PACKAGE_VERSION,
                  "GNU Astronomy Utilities version.", &status);

  /* If we are in a version controlled directory and have libgit2
     installed, write the commit description into the FITS file. */
  gitdescribe=gal_git_describe();
  if(gitdescribe)
    {
      fits_update_key(fptr, TSTRING, "COMMIT", gitdescribe,
                      "Git's commit description in running dir.", &status);
      free(gitdescribe);
    }

  /* Report any error if a problem came up */
  gal_fits_io_error(status, NULL);
}






/* Write the given keywords into the given extension of the given file,
   ending it with version information. This is primarily intended for
   writing configuration settings of a program into the zero-th extension
   of the FITS file (which is empty when the FITS file is created by
   Gnuastro's program and this library). */
void
gal_fits_key_write_config(gal_fits_list_key_t **keylist, char *title,
                          char *extname, char *filename, char *hdu)
{
  int status=0;
  fitsfile *fptr=gal_fits_hdu_open(filename, hdu, READWRITE);

  /* Delete the two extra comment lines describing the FITS standard that
     CFITSIO puts in when it creates a new extension. We'll set status to 0
     afterwards so even if they don't exist, the program continues
     normally. */
  fits_delete_key(fptr, "COMMENT", &status);
  fits_delete_key(fptr, "COMMENT", &status);
  status=0;

  /* Put a name for the zero-th extension. */
  if(fits_write_key(fptr, TSTRING, "EXTNAME", extname, "", &status))
    gal_fits_io_error(status, NULL);

  /* Write all the given keywords. */
  gal_fits_key_write_version_in_ptr(keylist, title, fptr);

  /* Close the FITS file. */
  if( fits_close_file(fptr, &status) )
    gal_fits_io_error(status, NULL);
}




















/*************************************************************
 ***********            Array functions            ***********
 *************************************************************/

/* Note that the FITS standard defines any array as an 'image',
   irrespective of how many dimensions it has. This function will return
   the Gnuastro-type, the number of dimensions and size along each
   dimension of the image along with its name and units if necessary (not
   NULL). Note that '*dsize' will be allocated here, so it must not point
   to any already allocated space. */
void
gal_fits_img_info(fitsfile *fptr, int *type, size_t *ndim, size_t **dsize,
                  char **name, char **unit)
{
  double bscale=NAN;
  size_t i, dsize_key=1;
  char **str, *bzero_str=NULL;
  int bitpix, status=0, naxis;
  gal_data_t *key, *keysll=NULL;
  long naxes[GAL_FITS_MAX_NDIM];

  /* Get the BITPIX, number of dimensions and size of each dimension. */
  if( fits_get_img_param(fptr, GAL_FITS_MAX_NDIM, &bitpix, &naxis,
                         naxes, &status) )
    gal_fits_io_error(status, NULL);
  *ndim=naxis;


  /* Convert bitpix to Gnuastro's known types. */
  *type=gal_fits_bitpix_to_type(bitpix);


  /* Define the names of the possibly existing important keywords about the
     dataset. We are defining these in the opposite order to be read by
     CFITSIO. The way Gnuastro writes the FITS keywords, the output will
     first have 'BZERO', then 'BSCALE', then 'EXTNAME', then, 'BUNIT'.*/
  gal_list_data_add_alloc(&keysll, NULL, GAL_TYPE_STRING, 1, &dsize_key,
                          NULL, 0, -1, 1, "BUNIT", NULL, NULL);
  gal_list_data_add_alloc(&keysll, NULL, GAL_TYPE_STRING, 1, &dsize_key,
                          NULL, 0, -1, 1, "EXTNAME", NULL, NULL);
  gal_list_data_add_alloc(&keysll, NULL, GAL_TYPE_FLOAT64, 1, &dsize_key,
                          NULL, 0, -1, 1, "BSCALE", NULL, NULL);
  gal_list_data_add_alloc(&keysll, NULL, GAL_TYPE_STRING, 1, &dsize_key,
                          NULL, 0, -1, 1, "BZERO", NULL, NULL);
  gal_fits_key_read_from_ptr(fptr, keysll, 0, 0);


  /* Read the special keywords. */
  i=1;
  for(key=keysll;key!=NULL;key=key->next)
    {
      /* Recall that the order is the opposite (this is a last-in-first-out
         list. */
      if(key->status==0)
        {
        switch(i)
          {
          case 4: if(unit) {str = key->array; *unit = *str; *str=NULL;} break;
          case 3: if(name) {str = key->array; *name = *str; *str=NULL;} break;
          case 2: bscale = *(double *)(key->array);                     break;
          case 1: str = key->array; bzero_str = *str;                   break;
          default:
            error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                  "fix the problem. For some reason, there are more "
                  "keywords requested ", __func__, PACKAGE_BUGREPORT);
          }
        }
      ++i;
    }


  /* If a type correction is necessary, then do it. */
  if( !isnan(bscale) || bzero_str )
    fits_type_correct(type, bscale, bzero_str);


  /* Allocate the array to keep the dimension size and fill it in, note
     that its order is the opposite of naxes. */
  *dsize=gal_pointer_allocate(GAL_TYPE_INT64, *ndim, 0, __func__, "dsize");
  for(i=0; i<*ndim; ++i)
    (*dsize)[i]=naxes[*ndim-1-i];


  /* Clean up. Note that bzero_str, gets freed by 'gal_data_free' (which is
     called by 'gal_list_data_free'. */
  gal_list_data_free(keysll);
}





/* Get the basic array info to remove extra dimensions if necessary. */
size_t *
gal_fits_img_info_dim(char *filename, char *hdu, size_t *ndim)
{
  fitsfile *fptr;
  size_t *dsize=NULL;
  int status=0, type;

  /* Open the given header, read the basic image information and close it
     again. */
  fptr=gal_fits_hdu_open(filename, hdu, READONLY);
  gal_fits_img_info(fptr, &type, ndim, &dsize, NULL, NULL);
  if( fits_close_file(fptr, &status) ) gal_fits_io_error(status, NULL);

  return dsize;
}





/* Read a FITS image HDU into a Gnuastro data structure. */
gal_data_t *
gal_fits_img_read(char *filename, char *hdu, size_t minmapsize,
                  int quietmmap)
{
  void *blank;
  long *fpixel;
  fitsfile *fptr;
  gal_data_t *img;
  size_t i, ndim, *dsize;
  char *name=NULL, *unit=NULL;
  int status=0, type, anyblank;


  /* Check HDU for realistic conditions: */
  fptr=gal_fits_hdu_open_format(filename, hdu, 0);


  /* Get the info and allocate the data structure. */
  gal_fits_img_info(fptr, &type, &ndim, &dsize, &name, &unit);


  /* Check if there is any dimensions (the first header can sometimes have
     no images). */
  if(ndim==0)
    error(EXIT_FAILURE, 0, "%s (hdu: %s) has 0 dimensions! The most common "
          "cause for this is a wrongly specified HDU. In some FITS images, "
          "the first HDU doesn't have any data, the data is in subsequent "
          "extensions. So probably reading the second HDU (with '--hdu=1' "
          "or '-h1') will solve the problem (following CFITSIO's "
          "convention, currently HDU counting starts from 0)." , filename,
          hdu);


  /* Set the fpixel array (first pixel in all dimensions). Note that the
     'long' type will not be larger than 64-bits, so, we'll just assume it
     is 64-bits for space allocation. On 32-bit systems, this won't be a
     problem, the space will be written/read as 32-bit 'long' any way,
     we'll just have a few empty bytes that will be freed anyway at the end
     of this function. */
  fpixel=gal_pointer_allocate(GAL_TYPE_INT64, ndim, 0, __func__, "fpixel");
  for(i=0;i<ndim;++i) fpixel[i]=1;


  /* Allocate the space for the array and for the blank values. */
  img=gal_data_alloc(NULL, type, ndim, dsize, NULL, 0, minmapsize,
                     quietmmap, name, unit, NULL);
  blank=gal_blank_alloc_write(type);
  if(name) free(name);
  if(unit) free(unit);
  free(dsize);


  /* Read the image into the allocated array: */
  fits_read_pix(fptr, gal_fits_type_to_datatype(type), fpixel,
                img->size, blank, img->array, &anyblank, &status);
  if(status) gal_fits_io_error(status, NULL);
  free(fpixel);
  free(blank);


  /* Close the input FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);


  /* Return the filled data structure */
  return img;
}





/* The user has specified an input file + extension, and your program needs
   this input to be a special type. For such cases, this function can be
   used to convert the input file to the desired type. */
gal_data_t *
gal_fits_img_read_to_type(char *inputname, char *hdu, uint8_t type,
                          size_t minmapsize, int quietmmap)
{
  gal_data_t *in, *converted;

  /* Read the specified input image HDU. */
  in=gal_fits_img_read(inputname, hdu, minmapsize, quietmmap);

  /* If the input had another type, convert it to float. */
  if(in->type!=type)
    {
      converted=gal_data_copy_to_new_type(in, type);
      gal_data_free(in);
      in=converted;
    }

  /* Return the final structure. */
  return in;
}





gal_data_t *
gal_fits_img_read_kernel(char *filename, char *hdu, size_t minmapsize,
                         int quietmmap)
{
  size_t i;
  int check=0;
  double sum=0;
  gal_data_t *kernel;
  float *f, *fp, tmp;

  /* Read the image as a float and if it has a WCS structure, free it. */
  kernel=gal_fits_img_read_to_type(filename, hdu, GAL_TYPE_FLOAT32,
                                   minmapsize, quietmmap);
  if(kernel->wcs) { wcsfree(kernel->wcs); kernel->wcs=NULL; }

  /* Check if the size along each dimension of the kernel is an odd
     number. If they are all an odd number, then the for each dimension,
     check will be incremented once. */
  for(i=0;i<kernel->ndim;++i)
    check += kernel->dsize[i]%2;
  if(check!=kernel->ndim)
    error(EXIT_FAILURE, 0, "%s: the kernel image has to have an odd number "
          "of pixels in all dimensions (there has to be one element/pixel "
          "in the center). At least one of the dimensions of %s (hdu: %s) "
          "doesn't have an odd number of pixels", __func__, filename, hdu);

  /* If there are any NaN pixels, set them to zero and normalize it. A
     blank pixel in a kernel is going to make a completely blank output.*/
  fp=(f=kernel->array)+kernel->size;
  do
    {
      if(isnan(*f)) *f=0.0f;
      else          sum+=*f;
    }
  while(++f<fp);
  kernel->flag |= GAL_DATA_FLAG_BLANK_CH;
  kernel->flag &= ~GAL_DATA_FLAG_HASBLANK;
  f=kernel->array; do *f++ *= 1/sum; while(f<fp);

  /* Flip the kernel about the center (necessary for non-symmetric
     kernels). */
  f=kernel->array;
  for(i=0;i<kernel->size/2;++i)
    {
      tmp=f[i];
      f[i]=f[ kernel->size - i - 1 ];
      f[ kernel->size - i - 1 ]=tmp;
    }

  /* Return the kernel*/
  return kernel;
}





/* This function will write all the data array information (including its
   WCS information) into a FITS file, but will not close it. Instead it
   will pass along the FITS pointer for further modification. */
fitsfile *
gal_fits_img_write_to_ptr(gal_data_t *input, char *filename)
{
  void *blank;
  int64_t *i64;
  fitsfile *fptr;
  uint64_t *u64, *u64f;
  long fpixel=1, *naxes;
  char *wcsstr, *u64key;
  size_t i, ndim=input->ndim;
  int nkeyrec, hasblank, status=0, datatype=0;
  gal_data_t *i64data, *towrite, *block=gal_tile_block(input);

  /* Small sanity check. */
  if( gal_fits_name_is_fits(filename)==0 )
    error(EXIT_FAILURE, 0, "%s: not a FITS suffix", filename);

  /* If the input is a tile (isn't a contiguous region of memory), then
     copy it into a contiguous region. */
  towrite = input==block ? input : gal_data_copy(input);
  hasblank=gal_blank_present(towrite, 0);

  /* Allocate the naxis area. */
  naxes=gal_pointer_allocate( ( sizeof(long)==8
                                ? GAL_TYPE_INT64
                                : GAL_TYPE_INT32 ), ndim, 0, __func__,
                              "naxes");


  /* Open the file for writing */
  fptr=gal_fits_open_to_write(filename);


  /* Fill the 'naxes' array (in opposite order, and 'long' type): */
  for(i=0;i<ndim;++i) naxes[ndim-1-i]=towrite->dsize[i];


  /* Create the FITS file. Unfortunately CFITSIO doesn't have a macro for
     UINT64, TLONGLONG is only for (signed) INT64. So if the dataset has
     that type, we'll have to convert it to 'INT64' and in the mean-time
     shift its zero, we will then have to write the BZERO and BSCALE
     keywords accordingly. */
  if(block->type==GAL_TYPE_UINT64)
    {
      /* Allocate the necessary space. */
      i64data=gal_data_alloc(NULL, GAL_TYPE_INT64, ndim, towrite->dsize,
                             NULL, 0, block->minmapsize, block->quietmmap,
                             NULL, NULL, NULL);

      /* Copy the values while making the conversion. */
      i64=i64data->array;
      u64f=(u64=towrite->array)+towrite->size;
      if(hasblank)
        {
          do *i64++ = ( *u64==GAL_BLANK_UINT64
                        ? GAL_BLANK_INT64
                        : (*u64 + INT64_MIN) );
          while(++u64<u64f);
        }
      else
        do *i64++ = (*u64 + INT64_MIN); while(++u64<u64f);

      /* We can now use CFITSIO's signed-int64 type macros. */
      datatype=TLONGLONG;
      fits_create_img(fptr, LONGLONG_IMG, ndim, naxes, &status);
      gal_fits_io_error(status, NULL);


      /* Write the image into the file. */
      fits_write_img(fptr, datatype, fpixel, i64data->size, i64data->array,
                     &status);
      gal_fits_io_error(status, NULL);


      /* We need to write the BZERO and BSCALE keywords manually. VERY
         IMPORTANT: this has to be done after writing the array. We cannot
         write this huge integer as a variable, so we'll simply write the
         full record/card. It is just important that the string be larger
         than 80 characters, CFITSIO will trim the rest of the string. */
      u64key="BZERO   =  9223372036854775808 / Offset of data                                         ";
      fits_write_record(fptr, u64key, &status);
      u64key="BSCALE  =                    1 / Default scaling factor                                 ";
      fits_write_record(fptr, u64key, &status);
      gal_fits_io_error(status, NULL);
    }
  else
    {
      /* Set the datatype */
      datatype=gal_fits_type_to_datatype(block->type);

      /* Create the FITS file. */
      fits_create_img(fptr, gal_fits_type_to_bitpix(towrite->type),
                      ndim, naxes, &status);
      gal_fits_io_error(status, NULL);

      /* Write the image into the file. */
      fits_write_img(fptr, datatype, fpixel, towrite->size, towrite->array,
                     &status);
      gal_fits_io_error(status, NULL);
    }


  /* Remove the two comment lines put by CFITSIO. Note that in some cases,
     it might not exist. When this happens, the status value will be
     non-zero. We don't care about this error, so to be safe, we will just
     reset the status variable after these calls. */
  fits_delete_key(fptr, "COMMENT", &status);
  fits_delete_key(fptr, "COMMENT", &status);
  status=0;


  /* If we have blank pixels, we need to define a BLANK keyword when we are
     dealing with integer types. */
  if(hasblank)
    switch(towrite->type)
      {
      case GAL_TYPE_FLOAT32:
      case GAL_TYPE_FLOAT64:
        /* Do nothing! Since there are much fewer floating point types
           (that don't need any BLANK keyword), we are checking them.*/
        break;

      default:
        blank=gal_fits_key_img_blank(towrite->type);
        if(fits_write_key(fptr, datatype, "BLANK", blank,
                          "Pixels with no data.", &status) )
          gal_fits_io_error(status, "adding the BLANK keyword");
        free(blank);
      }


  /* Write the extension name to the header. */
  if(towrite->name)
    fits_write_key(fptr, TSTRING, "EXTNAME", towrite->name, "", &status);


  /* Write the units to the header. */
  if(towrite->unit)
    fits_write_key(fptr, TSTRING, "BUNIT", towrite->unit, "", &status);


  /* Write comments if they exist. */
  if(towrite->comment)
    fits_write_comment(fptr, towrite->comment, &status);

  /* If a WCS structure is present, write it in */
  if(towrite->wcs)
    {
      /* Decompose the 'PCi_j' matrix and 'CDELTi' vector. */
      gal_wcs_decompose_pc_cdelt(towrite->wcs);

      /* Convert the WCS information to text. */
      status=wcshdo(WCSHDO_safe, towrite->wcs, &nkeyrec, &wcsstr);
      if(status)
        error(0, 0, "%s: WARNING: WCSLIB error, no WCS in output.\n"
              "wcshdu ERROR %d: %s", __func__, status,
              wcs_errmsg[status]);
      else
        gal_fits_key_write_wcsstr(fptr, wcsstr, nkeyrec);
      status=0;
    }

  /* Report any errors if we had any */
  free(naxes);
  gal_fits_io_error(status, NULL);
  if(towrite!=input) gal_data_free(towrite);
  return fptr;
}





void
gal_fits_img_write(gal_data_t *data, char *filename,
                   gal_fits_list_key_t *headers, char *program_string)
{
  int status=0;
  fitsfile *fptr;

  /* Write the data array into a FITS file and keep it open: */
  fptr=gal_fits_img_write_to_ptr(data, filename);

  /* Write all the headers and the version information. */
  gal_fits_key_write_version_in_ptr(&headers, program_string, fptr);

  /* Close the FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}





void
gal_fits_img_write_to_type(gal_data_t *data, char *filename,
                           gal_fits_list_key_t *headers, char *program_string,
                           int type)
{
  /* If the input dataset is not the correct type, then convert it,
     otherwise, use the input data structure. */
  gal_data_t *towrite = (data->type==type
                         ? data
                         : gal_data_copy_to_new_type(data, type));

  /* Write the converted dataset into an image. */
  gal_fits_img_write(towrite, filename, headers, program_string);

  /* Free the dataset if it was allocated. */
  if(towrite!=data) gal_data_free(towrite);
}





/* This function is mainly useful when you want to make FITS files in
   parallel (from one main WCS structure, with just differing CRPIX) for
   two reasons:

      - When a large number of FITS images (with WCS) need to be created in
        parallel, it can be much more efficient to write the header's WCS
        keywords once at first, write them in the FITS file, then just
        correct the CRPIX values.

      - WCSLIB's header writing function is not thread safe. So when
        writing FITS images in parallel, we can't write the header keywords
        in each thread.   */
void
gal_fits_img_write_corr_wcs_str(gal_data_t *input, char *filename,
                                char *wcsstr, int nkeyrec, double *crpix,
                                gal_fits_list_key_t *headers,
                                char *program_string)
{
  int status=0;
  fitsfile *fptr;

  /* The data should not have any WCS structure for this function. */
  if(input->wcs)
    error(EXIT_FAILURE, 0, "%s: input must not have WCS meta-data",
          __func__);

  /* Write the data array into a FITS file and keep it open. */
  fptr=gal_fits_img_write_to_ptr(input, filename);

  /* Write the WCS headers into the FITS file. */
  gal_fits_key_write_wcsstr(fptr, wcsstr, nkeyrec);

  /* Update the CRPIX keywords. Note that we don't want to change the
     values in the WCS information of gal_data_t. Because, it often happens
     that things are done in parallel, so we don't want to touch the
     original version, we just want to change the copied version. */
  if(crpix)
    {
      fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix[0], NULL, &status);
      fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix[1], NULL, &status);
      if(input->ndim==3)
        fits_update_key(fptr, TDOUBLE, "CRPIX3", &crpix[2], NULL, &status);
      gal_fits_io_error(status, NULL);
    }

  /* Write all the headers and the version information. */
  gal_fits_key_write_version_in_ptr(&headers, program_string, fptr);

  /* Close the file and return. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}




















/**************************************************************/
/**********                 Table                  ************/
/**************************************************************/
/* Get the size of a table HDU. CFITSIO doesn't use size_t, also we want to
   check status here.*/
void
gal_fits_tab_size(fitsfile *fitsptr, size_t *nrows, size_t *ncols)
{
  long lnrows;
  int incols, status=0;

  /* Read the sizes and put them in. */
  fits_get_num_rows(fitsptr, &lnrows, &status);
  fits_get_num_cols(fitsptr, &incols, &status);
  *ncols=incols;
  *nrows=lnrows;

  /* Report an error if any was issued. */
  gal_fits_io_error(status, NULL);
}





int
gal_fits_tab_format(fitsfile *fitsptr)
{
  int status=0;
  char value[FLEN_VALUE];

  fits_read_key(fitsptr, TSTRING, "XTENSION", value, NULL, &status);

  if(status==0)
    {
      if(!strcmp(value, "TABLE"))
        return GAL_TABLE_FORMAT_AFITS;
      else if(!strcmp(value, "BINTABLE"))
        return GAL_TABLE_FORMAT_BFITS;
      else
        error(EXIT_FAILURE, 0, "%s: the 'XTENSION' keyword of this FITS "
              "table ('%s') doesn't have a standard value", __func__, value);
    }
  else
    {
      if(status==KEY_NO_EXIST)
        error(EXIT_FAILURE, 0, "%s: input fitsfile pointer isn't a table",
              __func__);
      else
        gal_fits_io_error(status, NULL);
    }

  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can fix it. "
        "Control should not have reached the end of this function", __func__,
        PACKAGE_BUGREPORT);
  return -1;
}





/* The general format of the TDISPn keywords in FITS is like this: 'Tw.p',
   where 'T' specifies the general format, 'w' is the width to be given to
   this column and 'p' is the precision. For integer types, percision is
   actually the minimum number of integers, for floats, it is the number of
   decimal digits beyond the decimal point. */
static void
set_display_format(char *tdisp, gal_data_t *data, char *filename, char *hdu,
                   char *keyname)
{
  int isanint=0;
  char *tailptr;

  /* First, set the general display format */
  switch(tdisp[0])
    {
    case 'A':
      data->disp_fmt=GAL_TABLE_DISPLAY_FMT_STRING;
      break;

    case 'I':
      isanint=1;
      data->disp_fmt=GAL_TABLE_DISPLAY_FMT_DECIMAL;
      break;

    case 'O':
      isanint=1;
      data->disp_fmt=GAL_TABLE_DISPLAY_FMT_OCTAL;
      break;

    case 'Z':
      isanint=1;
      data->disp_fmt=GAL_TABLE_DISPLAY_FMT_HEX;
      break;

    case 'F':
      data->disp_fmt=GAL_TABLE_DISPLAY_FMT_FLOAT;
      break;

    case 'E':
    case 'D':
      data->disp_fmt=GAL_TABLE_DISPLAY_FMT_EXP;
      break;

    case 'G':
      data->disp_fmt=GAL_TABLE_DISPLAY_FMT_GENERAL;
      break;

    default:
      error(EXIT_FAILURE, 0, "%s (hdu: %s): Format character '%c' in the "
            "value (%s) of the keyword %s not recognized in %s", filename, hdu,
            tdisp[0], tdisp, keyname, __func__);
    }

  /* Parse the rest of the string to see if a width and precision are given
     or not. */
  data->disp_width=strtol(&tdisp[1], &tailptr, 0);
  switch(*tailptr)
    {
    case '.':      /* Width is set, go onto finding the precision. */
      data->disp_precision = strtol(&tailptr[1], &tailptr, 0);
      if(*tailptr!='\0')
        error(EXIT_FAILURE, 0, "%s (hdu: %s): The value '%s' of the "
              "'%s' keyword could not recognized (it doesn't finish after "
              "the precision) in %s", filename, hdu, tdisp, keyname, __func__);
      break;

    case '\0':     /* No precision given, use a default value.     */
      data->disp_precision = ( isanint
                               ? GAL_TABLE_DEF_PRECISION_INT
                               : GAL_TABLE_DEF_PRECISION_FLT );
      break;

    default:
      error(EXIT_FAILURE, 0, "%s (hdu: %s): The value '%s' of the "
            "'%s' keyword could not recognized (it doesn't have a '.', or "
            "finish, after the width) in %s", filename, hdu, tdisp,
            keyname, __func__);
    }


}




/* The FITS standard for binary tables (not ASCII tables) does not allow
   unsigned types for short, int and long types, or signed char! So it has
   'TSCALn' and 'TZEROn' to scale the signed types to an unsigned type. It
   does this internally, but since we need to define our data type and
   allocate space for it before actually reading the array, it is necessary
   to do this setting here.  */
static void
fits_correct_bin_table_int_types(gal_data_t *allcols, int tfields,
                                 int *tscal, long long *tzero)
{
  size_t i;

  for(i=0;i<tfields;++i)
    {
      /* If TSCALn is not 1, the reason for it isn't to use a different
         signed/unsigned type, so don't change anything. */
      if(tscal[i]!=1) continue;

      /* For a check
      printf("Column %zu initial type: %s (s: %d, z: %lld)\n", i+1,
             gal_data_type_as_string(allcols[i].type, 1), tscal[i], tzero[i]);
      */

      /* Correct the type based on the initial read type and the value to
         tzero. If tzero is any other value, then again, its not a type
         conversion, so just ignore it. */
      if(allcols[i].type==GAL_TYPE_UINT8 && tzero[i]==INT8_MIN)
        allcols[i].type = GAL_TYPE_INT8;

      else if ( allcols[i].type==GAL_TYPE_INT16
                && tzero[i] == -(long long)INT16_MIN )
        allcols[i].type = GAL_TYPE_UINT16;

      else if (allcols[i].type==GAL_TYPE_INT32
               && tzero[i] ==  -(long long)INT32_MIN)
        allcols[i].type = GAL_TYPE_UINT32;

      /* For a check
      printf("Column %zu corrected type: %s\n", i+1,
             gal_data_type_as_string(allcols[i].type, 1));
      */
    }
}





/* See the descriptions of 'gal_table_info'. */
gal_data_t *
gal_fits_tab_info(char *filename, char *hdu, size_t *numcols,
                  size_t *numrows, int *tableformat)
{
  long repeat;
  int tfields;        /* The maximum number of fields in FITS is 999 */
  char *tailptr;
  fitsfile *fptr;
  size_t i, index;
  long long *tzero;
  gal_data_t *allcols;
  int status=0, datatype, *tscal;
  char keyname[FLEN_KEYWORD]="XXXXXXXXXXXXX", value[FLEN_VALUE];


  /* Open the FITS file and get the basic information. */
  fptr=gal_fits_hdu_open_format(filename, hdu, 1);
  *tableformat=gal_fits_tab_format(fptr);
  gal_fits_tab_size(fptr, numrows, numcols);


  /* Read the total number of fields, then allocate space for the data
     structure array and store the information within it. */
  fits_read_key(fptr, TINT, "TFIELDS", &tfields, NULL, &status);
  allcols=gal_data_array_calloc(tfields);


  /* See comments of 'fits_correct_bin_table_int_types'. Here we are
     allocating the space to keep these values. */
  errno=0;
  tscal=calloc(tfields, sizeof *tscal);
  if(tscal==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for tscal", __func__,
          tfields*sizeof *tscal);
  errno=0;
  tzero=calloc(tfields, sizeof *tzero);
  if(tzero==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for tzero", __func__,
          tfields*sizeof *tzero);


  /* Read all the keywords one by one and if they match, then put them in
     the correct value. Note that we are starting from keyword 9 because
     according to the FITS standard, the first 8 keys in a FITS table are
     reserved. */
  for(i=9; strcmp(keyname, "END"); ++i)
    {
      /* Read the next keyword. */
      fits_read_keyn(fptr, i, keyname, value, NULL, &status);

      /* For string valued keywords, CFITSIO's function above, keeps the
         single quotes around the value string, one before and one
         after. 'gal_fits_key_clean_str_value' will remove these single
         quotes and any possible trailing space within the allocated
         space.*/
      if(value[0]=='\'') gal_fits_key_clean_str_value(value);

      /* COLUMN DATA TYPE. According the the FITS standard, the value of
         TFORM is most generally in this format: 'rTa'. 'T' is actually a
         code of the datatype. 'r' is the 'repeat' counter and 'a' is
         depreciated. Currently we can only read repeat==1 cases. When no
         number exists before the defined capital letter, it defaults to 1,
         but if a number exists (for example '5D'), then the repeat is 5
         (there are actually five values in each column). Note that
         value[0] is a single quote.*/
      if(strncmp(keyname, "TFORM", 5)==0)
        {
          /* See which column this information was for and add it, if the
             index is larger than the number of columns, then ignore
             the . The FITS standard says there should be no extra TFORM
             keywords beyond the number of columns, but we don't want to be
             that strict here.*/
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)     /* Counting from zero was corrected above. */
            {
              /* The FITS standard's value to this option for FITS ASCII
                 and binary files differ. */
              if(*tableformat==GAL_TABLE_FORMAT_AFITS)
                fits_ascii_tform(value, &datatype, NULL, NULL, &status);
              else
                fits_binary_tform(value, &datatype, &repeat, NULL, &status);

              /* Write the type into the data structure. */
              allcols[index].type=gal_fits_datatype_to_type(datatype, 1);

              /* If we are dealing with a string type, we need to know the
                 number of bytes in both cases for printing later. */
              if( allcols[index].type==GAL_TYPE_STRING )
                {
                  if(*tableformat==GAL_TABLE_FORMAT_AFITS)
                    {
                      repeat=strtol(value+1, &tailptr, 0);
                      if(*tailptr!='\0')
                        error(EXIT_FAILURE, 0, "%s (hdu: %s): the value to "
                              "keyword '%s' ('%s') is not in 'Aw' format "
                              "(for strings) as required by the FITS "
                              "standard in %s", filename, hdu, keyname, value,
                              __func__);
                    }
                  allcols[index].disp_width=repeat;
                }
            }
        }

      /* COLUMN SCALE FACTOR. */
      else if(strncmp(keyname, "TSCAL", 5)==0)
        {
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            {
              tscal[index]=strtol(value, &tailptr, 0);
              if(*tailptr!='\0')
                error(EXIT_FAILURE, 0, "%s (hdu: %s): value to %s keyword "
                      "('%s') couldn't be read as a number in %s", filename,
                      hdu, keyname, value, __func__);
            }
        }

      /* COLUMN ZERO VALUE (for signed/unsigned types). */
      else if(strncmp(keyname, "TZERO", 5)==0)
        {
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            {
              tzero[index]=strtoll(value, &tailptr, 0);
              if(*tailptr!='\0')
                error(EXIT_FAILURE, 0, "%s (hdu: %s): value to %s keyword "
                      "('%s') couldn't be read as a number in %s", filename,
                      hdu, keyname, value, __func__);
            }
        }

      /* COLUMN NAME. All strings in CFITSIO start and finish with single
         quotation marks, CFITSIO puts them in itsself, so if we don't
         remove them here, we might have duplicates later, its easier to
         just remove them to have a simple string that might be used else
         where too (without the single quotes).*/
      else if(strncmp(keyname, "TTYPE", 5)==0)
        {
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            gal_checkset_allocate_copy(value, &allcols[index].name);
        }

      /* COLUMN UNITS. */
      else if(strncmp(keyname, "TUNIT", 5)==0)
        {
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            gal_checkset_allocate_copy(value, &allcols[index].unit);
        }

      /* COLUMN COMMENTS */
      else if(strncmp(keyname, "TCOMM", 5)==0)
        {
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            gal_checkset_allocate_copy(value, &allcols[index].comment);
        }

      /* COLUMN BLANK VALUE. Note that to interpret the blank value the
         type of the column must already have been defined for this column
         in previous keywords. Otherwise, there will be a warning and it
         won't be used. */
      else if(strncmp(keyname, "TNULL", 5)==0)
        {
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields )
            {
              if(allcols[index].type<0)
                fprintf(stderr, "%s (hdu: %s): %s is located before "
                        "TFORM%zu, so the proper type to read/store the "
                        "blank value cannot be deduced", filename, hdu,
                        keyname, index+1);
              else
                gal_tableintern_read_blank(&allcols[index], value);
            }
        }

      /* COLUMN DISPLAY FORMAT */
      else if(strncmp(keyname, "TDISP", 5)==0)
        {
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            set_display_format(value, &allcols[index], filename, hdu,
                               keyname);
        }

      /* Column zero. */
    }

  /* Correct integer types, then free the allocated arrays. */
  fits_correct_bin_table_int_types(allcols, tfields, tscal, tzero);
  free(tscal);
  free(tzero);

  /* Close the FITS file and report an error if we had any. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
  return allcols;
}





/* Read CFITSIO un-readable (INF, -INF or NAN) floating point values in
   FITS ASCII tables. */
static void
fits_tab_read_ascii_float_special(char *filename, char *hdu, fitsfile *fptr,
                                  gal_data_t *out, size_t colnum,
                                  size_t numrows, size_t minmapsize,
                                  int quietmmap)
{
  double tmp;
  char **strarr;
  char *tailptr;
  gal_data_t *strrows;
  size_t i, colwidth=50;
  int anynul=0, status=0;

  /* Allocate the dataset to keep the string values. */
  strrows=gal_data_alloc(NULL, GAL_TYPE_STRING, 1, &numrows, NULL, 0,
                         minmapsize, quietmmap, NULL, NULL, NULL);

  /* Allocate the space to keep the string values. */
  strarr=strrows->array;
  for(i=0;i<numrows;++i)
    {
      errno=0;
      strarr[i]=calloc(colwidth, sizeof *strarr[i]);
      if(strarr[i]==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for "
              "strarr[%zu]", __func__, colwidth * sizeof *strarr[i], i);
    }

  /* Read the column as a string. */
  fits_read_col(fptr, TSTRING, colnum, 1, 1, out->size, NULL, strrows->array,
                &anynul, &status);
  gal_fits_io_error(status, NULL);

  /* Convert the strings to float. */
  for(i=0;i<numrows;++i)
    {
      /* Parse the string, if its not readable as a special number (like
         'inf' or 'nan', then just read it as a NaN. */
      tmp=strtod(strarr[i], &tailptr);
      if(tailptr==strarr[i]) tmp=NAN;

      /* Write it into the output dataset. */
      if(out->type==GAL_TYPE_FLOAT32)
        ((float *)(out->array))[i]=tmp;
      else
        ((double *)(out->array))[i]=tmp;
    }

  /* Clean up. */
  gal_data_free(strrows);
}





/* Read the column indexs into a dataset. */
gal_data_t *
gal_fits_tab_read(char *filename, char *hdu, size_t numrows,
                  gal_data_t *allcols, gal_list_sizet_t *indexll,
                  size_t minmapsize, int quietmmap)
{
  size_t i=0;
  void *blank;
  char **strarr;
  fitsfile *fptr;
  gal_data_t *out=NULL;
  gal_list_sizet_t *ind;
  int isfloat, status=0, anynul=0, hdutype;

  /* We actually do have columns to read. */
  if(numrows)
    {
      /* Open the FITS file */
      fptr=gal_fits_hdu_open_format(filename, hdu, 1);

      /* See if its a Binary or ASCII table (necessary for floating point
         blank values). */
      if (fits_get_hdu_type(fptr, &hdutype, &status) )
        gal_fits_io_error(status, NULL);

      /* Pop each index and read/store the array. */
      for(ind=indexll; ind!=NULL; ind=ind->next)
        {
          /* Allocate the necessary data structure (including the array)
             for this column. */
          gal_list_data_add_alloc(&out, NULL, allcols[ind->v].type, 1,
                                  &numrows, NULL, 0, minmapsize, quietmmap,
                                  allcols[ind->v].name, allcols[ind->v].unit,
                                  allcols[ind->v].comment);

          /* For a string column, we need an allocated array for each element,
             even in binary values. This value should be stored in the
             disp_width element of the data structure, which is done
             automatically in 'gal_fits_table_info'. */
          if(out->type==GAL_TYPE_STRING)
            {
              strarr=out->array;
              for(i=0;i<numrows;++i)
                {
                  errno=0;
                  strarr[i]=calloc(allcols[ind->v].disp_width+1,
                                   sizeof *strarr[i]);
                  if(strarr[i]==NULL)
                    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for "
                          "strarr[%zu]", __func__,
                          (allcols[ind->v].disp_width+1) * sizeof *strarr[i],
                          i);
                }
            }

          /* Allocate a blank value for the given type and read/store the
             column using CFITSIO. Note that for binary tables, we only
             need blank values for integer types. For binary floating point
             types, the FITS standard defines blanks as NaN (same as almost
             any other software like Gnuastro). However if a blank value is
             specified, CFITSIO will convert other special numbers like
             'inf' to NaN also. We want to be able to distringuish 'inf'
             and NaN here, so for floating point types in binary tables, we
             won't define any blank value. In ASCII tables, CFITSIO doesn't
             read the 'NAN' values (that it has written itself) unless we
             specify a blank pointer/value. */
          isfloat = ( out->type==GAL_TYPE_FLOAT32
                      || out->type==GAL_TYPE_FLOAT64 );
          blank = ( ( hdutype==BINARY_TBL && isfloat )
                    ? NULL
                    : gal_blank_alloc_write(out->type) );
          fits_read_col(fptr, gal_fits_type_to_datatype(out->type), ind->v+1,
                        1, 1, out->size, blank, out->array, &anynul, &status);

          /* In the ASCII table format, CFITSIO might not be able to read
             'INF' or '-INF'. In this case, it will set status to 'BAD_C2D'
             or 'BAD_C2F'. So, we'll use our own parser for the column
             values. */
          if( hdutype==ASCII_TBL
              && isfloat
              && (status==BAD_C2D || status==BAD_C2F) )
            {
              fits_tab_read_ascii_float_special(filename, hdu, fptr, out,
                                                ind->v+1, numrows,
                                                minmapsize, quietmmap);
              status=0;
            }

          /* Clean up and sanity check. */
          if(blank) free(blank);
          gal_fits_io_error(status, NULL);
        }

      /* Close the FITS file */
      fits_close_file(fptr, &status);
      gal_fits_io_error(status, NULL);
    }

  /* There are no rows to read ('numrows==NULL'). Make an empty-sized
     array. */
  else
    {
      /* We are setting a 1-element array to avoid any allocation
         errors. Then we are freeing the allocated spaces and correcting
         the sizes. */
      numrows=1;
      for(ind=indexll; ind!=NULL; ind=ind->next)
        {
          /* Do the allocation. */
          gal_list_data_add_alloc(&out, NULL, allcols[ind->v].type, 1,
                                  &numrows, NULL, 0, minmapsize, quietmmap,
                                  allcols[ind->v].name, allcols[ind->v].unit,
                                  allcols[ind->v].comment);

          /* Correct the array and sizes. */
          out->size=0;
          free(out->array);
          free(out->dsize);
          out->dsize=out->array=NULL;
        }
    }

  /* Return the output. */
  return out;
}





/* This function will allocate new copies for all elements to have the same
   length as the maximum length and set all trailing elements to '\0' for
   those that are shorter than the length. The return value is the
   allocated space. If the dataset is not a string, the returned value will
   be -1 (largest number of 'size_t'). */
static size_t
fits_string_fixed_alloc_size(gal_data_t *data)
{
  size_t i, j, maxlen=0;
  char *tmp, **strarr=data->array;

  /* Return 0 if the dataset is not a string. */
  if(data->type!=GAL_TYPE_STRING)
    return -1;

  /* Get the maximum length. */
  for(i=0;i<data->size;++i)
    maxlen = strlen(strarr[i])>maxlen ? strlen(strarr[i]) : maxlen;

  /* For all elements, check the length and if they aren't equal to maxlen,
     then allocate a maxlen sized array and put the values in. */
  for(i=0;i<data->size;++i)
    {
      /* Allocate (and clear) the space for the new string. We want it to
         be cleared, so when the strings are smaller, the rest of the space
         is filled with '\0' (ASCII for 0) values.*/
      errno=0;
      tmp=calloc(maxlen+1, sizeof *strarr[i]);
      if(tmp==NULL)
        error(EXIT_FAILURE, 0, "%s: %zu bytes for tmp", __func__,
              (maxlen+1)*sizeof *strarr[i]);

      /* Put the old array into the newly allocated space. 'tmp' was
         cleared (all values set to '\0', so we don't need to set the final
         one explicity after the copy.*/
      for(j=0;strarr[i][j]!='\0';++j)
        tmp[j]=strarr[i][j];

      /* Free the old array and put in the new one. */
      free(strarr[i]);
      strarr[i]=tmp;
    }

  /* Return the allocated space. */
  return maxlen+1;
}





static void
fits_table_prepare_arrays(gal_data_t *cols, size_t numcols, int tableformat,
                          char ***outtform, char ***outttype,
                          char ***outtunit)
{
  size_t i=0;
  gal_data_t *col;
  char fmt[2], lng[3];
  char *blank, **tform, **ttype, **tunit;


  /* Allocate the arrays to keep the 'tform' values */
  errno=0;
  tform=*outtform=malloc(numcols*sizeof *tform);
  if(tform==NULL)
    error(EXIT_FAILURE, 0, "%s: %zu bytes for tform", __func__,
          numcols*sizeof *tform);
  errno=0;
  ttype=*outttype=malloc(numcols*sizeof *ttype);
  if(ttype==NULL)
    error(EXIT_FAILURE, 0, "%s: %zu bytes for ttype", __func__,
          numcols*sizeof *ttype);
  errno=0;
  tunit=*outtunit=malloc(numcols*sizeof *tunit);
  if(tunit==NULL)
    error(EXIT_FAILURE, 0, "%s: %zu bytes for tunit", __func__,
          numcols*sizeof *tunit);


  /* Go over each column and fill in these arrays. */
  for(col=cols; col!=NULL; col=col->next)
    {
      /* Set the 'ttype' and 'tunit' values: */
      if( asprintf(&ttype[i], "%s", col->name ? col->name : "")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      if( asprintf(&tunit[i], "%s", col->unit ? col->unit : "")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);


      /* FITS's TFORM depends on the type of FITS table, so work
         differently. */
      switch(tableformat)
        {
        /* FITS ASCII table. */
        case GAL_TABLE_FORMAT_AFITS:

            /* Fill the printing format. */
            gal_tableintern_col_print_info(col, GAL_TABLE_FORMAT_AFITS,
                                           fmt, lng);

            /* We need to check if the blank value needs is larger than the
               expected width or not. Its initial width is set the output
               of the function above, but if the value is larger,
               'asprintf' (which is used) will make it wider. */
            blank = ( gal_blank_present(col, 0)
                      ? gal_blank_as_string(col->type, col->disp_width)
                      : NULL );

            /* Adjust the width. */
            if(blank)
              {
                col->disp_width = ( strlen(blank) > col->disp_width
                                    ? strlen(blank) : col->disp_width );
                free(blank);
              }

            /* Print the value to be used as TFORMn:  */
            switch(col->type)
              {
              case GAL_TYPE_STRING:
              case GAL_TYPE_UINT8:
              case GAL_TYPE_INT8:
              case GAL_TYPE_UINT16:
              case GAL_TYPE_INT16:
              case GAL_TYPE_UINT32:
              case GAL_TYPE_INT32:
              case GAL_TYPE_UINT64:
              case GAL_TYPE_INT64:
                if( asprintf(&tform[i], "%c%d", fmt[0], col->disp_width)<0 )
                  error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
                break;

              case GAL_TYPE_FLOAT32:
              case GAL_TYPE_FLOAT64:
                if( asprintf(&tform[i], "%c%d.%d", fmt[0], col->disp_width,
                             col->disp_precision)<0 )
                  error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
                break;

              default:
                error(EXIT_FAILURE, 0, "%s: col->type code %d not recognized",
                      __func__, col->type);
              }
          break;


        /* FITS binary table. */
        case GAL_TABLE_FORMAT_BFITS:

          /* If this is a string column, set all the strings to same size,
             then write the value of tform depending on the type. */
          col->disp_width=fits_string_fixed_alloc_size(col);
          fmt[0]=gal_fits_type_to_bin_tform(col->type);
          if( col->type==GAL_TYPE_STRING )
            {
              if( asprintf(&tform[i], "%d%c", col->disp_width, fmt[0])<0 )
                error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
            }
          else
            {
              if( asprintf(&tform[i], "%c", fmt[0])<0 )
                error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
            }
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: tableformat code %d not recognized",
                __func__, tableformat);
        }


      /* Increment the column index. */
      ++i;
    }
}





/* Write the TNULLn keywords into the FITS file. Note that this depends on
   the type of the table: for an ASCII table, all the columns need it. For
   a binary table, only the non-floating point ones (even if they don't
   have NULL values) must have it. */
static void
fits_write_tnull_tcomm(fitsfile *fptr, gal_data_t *col, int tableformat,
                       size_t colnum, char *tform)
{
  void *blank;
  int status=0;
  char *c, *keyname, *bcomment;

  /* Write the NULL value */
  switch(tableformat)
    {
    case GAL_TABLE_FORMAT_AFITS:

      /* Print the keyword and value. */
      if( asprintf(&keyname, "TNULL%zu", colnum)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      blank=gal_blank_as_string(col->type, col->disp_width);

      /* When in exponential form ('tform' starting with 'E'), CFITSIO
         writes a NaN value as 'NAN', but when in floating point form
         ('tform' starting with 'F'), it writes it as 'nan'. So in the
         former case, we need to convert the string to upper case. */
      if(tform[0]=='E' || tform[0]=='e')
        for(c=blank; *c!='\0'; ++c) *c=toupper(*c);

      /* Write in the header. */
      fits_write_key(fptr, TSTRING, keyname, blank,
                     "blank value for this column", &status);

      /* Clean up. */
      free(keyname);
      free(blank);
      break;

    case GAL_TABLE_FORMAT_BFITS:
      /* FITS binary tables don't accept NULL values for floating point or
         string columns. For floating point is must be NaN and for strings
         it is a blank string. */
      if( col->type!=GAL_TYPE_FLOAT32
          && col->type!=GAL_TYPE_FLOAT64
          && col->type!=GAL_TYPE_STRING )
        {
          blank=gal_blank_alloc_write(col->type);
          if( asprintf(&keyname, "TNULL%zu", colnum)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          fits_write_key(fptr, gal_fits_type_to_datatype(col->type),
                         keyname, blank, "blank value for this column",
                         &status);
          gal_fits_io_error(status, NULL);
          free(keyname);
          free(blank);
        }
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: tableformat code %d not recognized",
            __func__, tableformat);
    }

  /* Write the comments if there is any. */
  if(col->comment)
    {
      if( asprintf(&keyname, "TCOMM%zu", colnum)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      if( asprintf(&bcomment, "comment for field %zu", colnum)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      fits_write_key(fptr, TSTRING, keyname, col->comment,
                     bcomment, &status);
      gal_fits_io_error(status, NULL);
      free(keyname);
      free(bcomment);
    }
}





/* Write the given columns (a linked list of 'gal_data_t') into a FITS
   table.*/
void
gal_fits_tab_write(gal_data_t *cols, gal_list_str_t *comments,
                   int tableformat, char *filename, char *extname)
{
  void *blank;
  fitsfile *fptr;
  gal_data_t *col;
  size_t i, numrows=-1;
  gal_list_str_t *strt;
  char **ttype, **tform, **tunit;
  int tbltype, numcols=0, status=0;


  /* Make sure all the input columns have the same number of elements */
  for(col=cols; col!=NULL; col=col->next)
    {
      if(numrows==-1) numrows=col->size;
      else if(col->size!=numrows)
        error(EXIT_FAILURE, 0, "%s: the number of records/rows in the input "
              "columns are not equal", __func__);
      ++numcols;
    }


  /* Open the FITS file for writing. */
  fptr=gal_fits_open_to_write(filename);


  /* Prepare necessary arrays and if integer type columns have blank
     values, write the TNULLn keywords into the FITS file. */
  fits_table_prepare_arrays(cols, numcols, tableformat,
                            &tform, &ttype, &tunit);


  /* Make the FITS file pointer. Note that tableformat was checked in
     'fits_table_prepare_arrays'. */
  tbltype = tableformat==GAL_TABLE_FORMAT_AFITS ? ASCII_TBL : BINARY_TBL;
  fits_create_tbl(fptr, tbltype, numrows, numcols, ttype, tform, tunit,
                  extname, &status);
  gal_fits_io_error(status, NULL);


  /* Write the columns into the file and also write the blank values into
     the header when necessary. */
  i=0;
  for(col=cols; col!=NULL; col=col->next)
    {
      /* Write the blank value into the header and return a pointer to
         it. Otherwise, */
      fits_write_tnull_tcomm(fptr, col, tableformat, i+1, tform[i]);

      /* Set the blank pointer if its necessary, note that strings don't
         need a blank pointer in a FITS ASCII table.*/
      blank = ( gal_blank_present(col, 0)
                ? gal_blank_alloc_write(col->type) : NULL );
      if(tableformat==GAL_TABLE_FORMAT_AFITS && col->type==GAL_TYPE_STRING)
        { if(blank) free(blank); blank=NULL; }

      /* Write the full column into the table. */
      fits_write_colnull(fptr, gal_fits_type_to_datatype(col->type),
                         i+1, 1, 1, col->size, col->array, blank, &status);
      gal_fits_io_error(status, NULL);

      /* Clean up and Increment the column counter. */
      if(blank) free(blank);
      ++i;
    }


  /* Write the comments if there were any. */
  for(strt=comments; strt!=NULL; strt=strt->next)
    fits_write_comment(fptr, strt->v, &status);


  /* Write all the headers and the version information. */
  gal_fits_key_write_version_in_ptr(NULL, NULL, fptr);


  /* Clean up and close the FITS file. Note that each element in the
     'ttype' and 'tunit' arrays just points to the respective string in the
     column data structure, the space for each element of the array wasn't
     allocated.*/
  for(i=0;i<numcols;++i)
    {
      free(tform[i]);
      free(ttype[i]);
      free(tunit[i]);
    }
  free(tform);
  free(ttype);
  free(tunit);
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}
