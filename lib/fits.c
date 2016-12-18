/*********************************************************************
Functions to work with FITS image data.
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
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gnuastro/git.h>
#include <gnuastro/fits.h>

#include "checkset.h"
#include "fixedstringmacros.h"










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
   include them in `gal_checkset_automatic_output', so the compression
   suffix can be skipped when the user doesn't specify an output
   filename.*/
int
gal_fits_name_is_fits(char *name)
{
  size_t len;
  len=strlen(name);
  if ( ( len>=4 && strcmp(&name[len-4], "fits") == 0 )
       || ( len>=7 && strcmp(&name[len-7], "fits.gz") == 0 )
       || ( len>=6 && strcmp(&name[len-6], "fits.Z") == 0 )
       || ( len>=3 && strcmp(&name[len-3], "imh") == 0 )
       || ( len>=7 && strcmp(&name[len-7], "fits.fz") == 0 ) )
    return 1;
  else
    return 0;
}





/* IMPORTANT NOTE: if other compression suffixes are add to this function,
   include them in `gal_checkset_automatic_output', so the compression
   suffix can be skipped when the user doesn't specify an output
   filename.*/
int
gal_fits_suffix_is_fits(char *suffix)
{
  if (strcmp(suffix, "fits") == 0 || strcmp(suffix, ".fits") == 0
      || strcmp(suffix, "fits.gz") == 0 || strcmp(suffix, ".fits.gz") == 0
      || strcmp(suffix, "fits.Z") == 0 || strcmp(suffix, ".fits.Z") == 0
      || strcmp(suffix, "imh") == 0 || strcmp(suffix, ".imh") == 0
      || strcmp(suffix, "fits.fz") == 0 || strcmp(suffix, ".fits.fz") == 0)
   return 1;
 else
   return 0;
}





/* We have the name of the input file. But in most cases, the files
   that should be used (for example a mask image) are other extensions
   in the same file. So the user only has to give the HDU. The job of
   this function is to determine which is the case and set othername
   to the appropriate value. */
void
gal_fits_file_or_ext_name(char *inputname, char *inhdu, int othernameset,
                          char **othername, char *ohdu, int ohduset,
                          char *type)
{
  if(othernameset)
    {
      /* In some cases, for example a mask image, both the name and
         HDU are optional. So just to be safe, we will check this all
         the time. */
      if(ohduset==0)
        error(EXIT_FAILURE, 0, "a %s image was specified (%s). However, "
              "no HDU is given for it. Please add a HDU. If you regularly "
              "use the same HDU as %s, you may consider adding it to "
              "the configuration file. For more information, please see "
              "the `Configuration files' section of the %s manual by "
              "running ` info gnuastro ' on the command-line", type,
              *othername, type, PACKAGE_NAME);
      if(strcmp(inputname, *othername)==0)
        {
          if(strcmp(ohdu, inhdu)==0)
            error(EXIT_FAILURE, 0, "the specified %s name and "
                  "input image name (%s) are the same while the input "
                  "image hdu name and mask hdu are also identical (%s)",
                  type, inputname, inhdu);
        }
    }
    else if(ohduset && strcmp(ohdu, inhdu))
      *othername=inputname;
    else
      *othername=NULL;
}




















/*************************************************************
 **************           Type codes           ***************
 *************************************************************/
int
gal_fits_bitpix_to_type(int bitpix)
{
  switch(bitpix)
    {
    case BYTE_IMG:
      return GAL_DATA_TYPE_UCHAR;
    case SBYTE_IMG:
      return GAL_DATA_TYPE_CHAR;
    case USHORT_IMG:
      return GAL_DATA_TYPE_USHORT;
    case SHORT_IMG:
      return GAL_DATA_TYPE_SHORT;
    case ULONG_IMG:
      return GAL_DATA_TYPE_ULONG;
    case LONG_IMG:
      return GAL_DATA_TYPE_LONG;
    case LONGLONG_IMG:
      return GAL_DATA_TYPE_LONGLONG;
    case FLOAT_IMG:
      return GAL_DATA_TYPE_FLOAT;
    case DOUBLE_IMG:
      return GAL_DATA_TYPE_DOUBLE;
    default:
      error(EXIT_FAILURE, 0, "bitpix value of %d not recognized in "
            "gal_fits_bitpix_to_type", bitpix);
    }
  return 0;
}





int
gal_fits_type_to_bitpix(int type)
{
  switch(type)
    {
    case GAL_DATA_TYPE_UCHAR:
      return BYTE_IMG;
    case GAL_DATA_TYPE_CHAR:
      return SBYTE_IMG;
    case GAL_DATA_TYPE_USHORT:
      return USHORT_IMG;
    case GAL_DATA_TYPE_SHORT:
      return SHORT_IMG;
    case GAL_DATA_TYPE_ULONG:
      return ULONG_IMG;
    case GAL_DATA_TYPE_LONG:
      return LONG_IMG;
    case GAL_DATA_TYPE_LONGLONG:
      return LONGLONG_IMG;
    case GAL_DATA_TYPE_FLOAT:
      return FLOAT_IMG;
    case GAL_DATA_TYPE_DOUBLE:
      return DOUBLE_IMG;
    default:
      error(EXIT_FAILURE, 0, "type value of %d not recognized in "
            "gal_fits_type_to_bitpix", type);
    }
  return 0;
}





/* The values to the TFORM header keyword are single letter capital
   letters, but that is useless in identifying the data type of the
   column. So this function will do the conversion based on the CFITSIO
   manual.

   Note that the characters are the same for ASCII or binary tables.*/
int
gal_fits_tform_to_type(char tform)
{
  switch(tform)
    {
    case 'X':
      return GAL_DATA_TYPE_BIT;
    case 'B':
      return GAL_DATA_TYPE_UCHAR;
    case 'S': case 'L':
      return GAL_DATA_TYPE_CHAR;
    case 'A':
      return GAL_DATA_TYPE_STRING;
    case 'V':
      return GAL_DATA_TYPE_UINT;
    case 'U':
      return GAL_DATA_TYPE_USHORT;
    case 'I':
      return GAL_DATA_TYPE_SHORT;
    case 'J':
      return GAL_DATA_TYPE_LONG;
    case 'K':
      return GAL_DATA_TYPE_LONGLONG;
    case 'E':
      return GAL_DATA_TYPE_FLOAT;
    case 'D':
      return GAL_DATA_TYPE_DOUBLE;
    case 'C':
      return GAL_DATA_TYPE_COMPLEX;
    case 'M':
      return GAL_DATA_TYPE_DCOMPLEX;
    default:
      error(EXIT_FAILURE, 0, "'%c' is not a recognized CFITSIO value for "
            "the TFORMn header keyword(s).", tform);
    }

  error(EXIT_FAILURE, 0, "A bug! Please contact us so we can fix this. "
        "For some reason, control has reached to the end of the "
        "gal_fits_tform_to_datatype function in fits.c.");
  return -1;
}





int
gal_fits_type_to_datatype(int type)
{
  switch(type)
    {
    case GAL_DATA_TYPE_BIT:
      return TBIT;

    case GAL_DATA_TYPE_UCHAR:
      return TBYTE;

    case GAL_DATA_TYPE_CHAR:
      return TSBYTE;

    case GAL_DATA_TYPE_STRING:
      return TSTRING;

    case GAL_DATA_TYPE_USHORT:
      return TUSHORT;

    case GAL_DATA_TYPE_SHORT:
      return TSHORT;

    case GAL_DATA_TYPE_UINT:
      return TUINT;

    case GAL_DATA_TYPE_INT:
      return TINT;

    case GAL_DATA_TYPE_ULONG:
      return TULONG;

    case GAL_DATA_TYPE_LONG:
      return TLONG;

    case GAL_DATA_TYPE_LONGLONG:
      return TLONGLONG;

    case GAL_DATA_TYPE_FLOAT:
      return TFLOAT;

    case GAL_DATA_TYPE_DOUBLE:
      return TDOUBLE;

    case GAL_DATA_TYPE_COMPLEX:
      return TCOMPLEX;

    case GAL_DATA_TYPE_DCOMPLEX:
      return TDBLCOMPLEX;

    default:
      error(EXIT_FAILURE, 0, "'%d' is not a recognized Gnuastro type. "
            "It was given to `gal_fits_type_to_datatype'.", type);
    }

  error(EXIT_FAILURE, 0, "A bug! Please contact us so we can fix this. "
        "For some reason, control has reached to the end of the "
        "gal_fits_type_to_datatype function in fits.c.");
  return -1;
}





int
gal_fits_datatype_to_type(int datatype)
{
  switch(datatype)
    {
    case TBIT:
      return GAL_DATA_TYPE_BIT;

    case TBYTE:
      return GAL_DATA_TYPE_UCHAR;

    case TSBYTE:
      return GAL_DATA_TYPE_CHAR;

    case TSTRING:
      return GAL_DATA_TYPE_STRING;

    case TUSHORT:
      return GAL_DATA_TYPE_USHORT;

    case TSHORT:
      return GAL_DATA_TYPE_SHORT;

    case TUINT:
      return GAL_DATA_TYPE_UINT;

    case TINT:
      return GAL_DATA_TYPE_INT;

    case TULONG:
      return GAL_DATA_TYPE_ULONG;

    case TLONG:
      return GAL_DATA_TYPE_LONG;

    case TLONGLONG:
      return GAL_DATA_TYPE_LONGLONG;

    case TFLOAT:
      return GAL_DATA_TYPE_FLOAT;

    case TDOUBLE:
      return GAL_DATA_TYPE_DOUBLE;

    case TCOMPLEX:
      return GAL_DATA_TYPE_COMPLEX;

    case TDBLCOMPLEX:
      return GAL_DATA_TYPE_DCOMPLEX;

    default:
      error(EXIT_FAILURE, 0, "'%d' is not a recognized CFITSIO datatype. "
            "It was given to `gal_fits_datatype_to_type'.", datatype);
    }

  error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix "
        "this. For some reason, control has reached to the end of the "
        "gal_fits_datatype_to_type function in fits.c.", PACKAGE_BUGREPORT);
  return -1;
}




















/*************************************************************
 **************        Get information         ***************
 *************************************************************/
void
gal_fits_num_hdus(char *filename, int *numhdu)
{
  int status=0;
  fitsfile *fptr;

  /* We don't need to check for an error everytime, because we don't
     make any non CFITSIO usage of the output. It is necessary to
     check every CFITSIO call, only when you will need to use the
     outputs. */
  fits_open_file(&fptr, filename, READONLY, &status);

  fits_get_num_hdus(fptr, numhdu, &status);

  fits_close_file(fptr, &status);

  gal_fits_io_error(status, NULL);
}





/* Note that the FITS standard defines any array as an `image',
   irrespective of how many dimensions it has. */
void
gal_fits_img_info(fitsfile *fptr, int *type, size_t *ndim, long **dsize)
{
  size_t i;
  int bitpix, status=0, naxis;
  long naxes[GAL_DATA_MAXDIM];

  /* Get the BITPIX, number of dimensions and size of each dimension. */
  if( fits_get_img_param(fptr, GAL_DATA_MAXDIM, &bitpix, &naxis,
                         naxes, &status) )
    gal_fits_io_error(status, NULL);
  *ndim=naxis;

  /* Convert bitpix to Gnuastro's known types. */
  *type=gal_fits_bitpix_to_type(bitpix);

  /* Allocate space for the size along each dimension. */
  errno=0;
  *dsize=malloc( *ndim * sizeof **dsize );
  if(*dsize==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for dsize in gal_fits_img_info",
          *ndim * sizeof **dsize);

  /* Put the size of each dimention into the output array. */
  for(i=0;i<*ndim;++i)
    (*dsize)[i]=naxes[i];
}




















/**************************************************************/
/**********                  HDU                   ************/
/**************************************************************/
/* Check the desired HDU in a FITS image and also if it has the
   desired type. */
fitsfile *
gal_fits_read_hdu(char *filename, char *hdu, unsigned char img0_tab1)
{
  size_t len;
  char *ffname;
  fitsfile *fptr;
  int status=0, hdutype;

  /* Add hdu to filename: */
  errno=0;
  len=strlen(filename)+strlen(hdu)+4;
  ffname=malloc(len*sizeof *ffname);
  if(ffname==NULL)
    error(EXIT_FAILURE, errno, "%zu characters", len);
  sprintf(ffname, "%s[%s#]", filename, hdu);

  /* Open the FITS file: */
  if( fits_open_file(&fptr, ffname, READONLY, &status) )
    gal_fits_io_error(status, "reading this FITS file");

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
  free(ffname);
  return fptr;
}




















/**************************************************************/
/**********            Header keywords             ************/
/**************************************************************/

/* Each keyword name, value, comment, type is kept within a Gnuastro data
   structure. The input should be an array of such data structures,
   either statically allocated like:

        gal_data_t keys[2];

   or dynamically allocated like:

        gal_data_t *keys;
        keys=malloc(2*sizeof *keys);

   Before calling this function, you just have to set the `title' and
   `type' values of the data structure. The given title value will be
   directly passed to CFITSIO to read the desired keyword.

   CFITSIO will start searching for the keywords from the last place in the
   header that it searched for a keyword. So it is much more efficient if
   the order that you ask for keywords is based on the order they are
   stored in the header.

   ABOUT THE STRING VALUES:

   The space for a string value is dynamically allocated within the
   `gal_fits_key' structure (to be `FLEN_VALUE' characters, `FLEN_VALUE' is
   defined by CFITSIO and is the largest possible number of characters).

   This function will not abort if CFITSIO is unable to read the keyword
   due to any reason. You can check the successful reading of the keyword
   from the `status' value in each keyword's data structure. If its zero,
   then the keyword was found and read as expected. Otherwise, this a
   CFITSIO status, so you use its error reporting tools or
   `gal_fits_io_error' for reporting the reason. If you have an alternative
   to the keyword, or its not mandatory, and the keyword doesn't exist,
   then the status value will be KEY_NO_EXIST (from CFITSIO).
 */
void
gal_fits_read_keywords_fptr(fitsfile *fptr, gal_data_t *keysll,
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

        /* When the type is a string, `tmp->array' is an array of pointers
           to a separately allocated piece of memory. So we have to
           allocate that space here. If its not a string, then the
           allocated space above is enough to keep the value.*/
        switch(tmp->type)
          {
          case GAL_DATA_TYPE_STRING:
            errno=0;
            strarray=tmp->array;
            valueptr=strarray[0]=malloc(FLEN_VALUE * sizeof *strarray[0]);
            if(strarray[0]==NULL)
              error(EXIT_FAILURE, errno, "%zu bytes for strarray[0] in "
                    "`gal_fits_read_keywords_fprt'",
                    FLEN_VALUE * sizeof *strarray[0]);
            break;

          default:
            valueptr=tmp->array;
          }

        /* Allocate space for the keyword comment if necessary.*/
        if(readcomment)
          {
            errno=0;
            tmp->comment=malloc(FLEN_COMMENT * sizeof *tmp->comment);
            if(tmp->comment==NULL)
              error(EXIT_FAILURE, errno, "%zu bytes for tmp->comment in "
                    "`gal_fits_read_keywords_fprt'",
                    FLEN_COMMENT * sizeof *tmp->comment);
          }
        else
          tmp->comment=NULL;

        /* Allocate space for the keyword unit if necessary. Note that
           since there is no precise CFITSIO length for units, we will use
           the `FLEN_COMMENT' length for units too (theoretically, the unit
           might take the full remaining area in the keyword). Also note
           that the unit is only optional, so it needs a separate CFITSIO
           function call which is done here.*/
        if(readunit)
          {
            /* Allocate space for the unit and read it in. */
            errno=0;
            tmp->unit=malloc(FLEN_COMMENT * sizeof *tmp->unit);
            if(tmp->unit==NULL)
              error(EXIT_FAILURE, errno, "%zu bytes for tmp->unit in "
                    "`gal_fits_read_keywords_fprt'",
                    FLEN_COMMENT * sizeof *tmp->unit);
            fits_read_key_unit(fptr, tmp->name, tmp->unit, &tmp->status);

            /* If the string is empty, free the space and set it to NULL. */
            if(tmp->unit[0]=='\0') {free(tmp->unit); tmp->unit=NULL;}
          }
        else
          tmp->unit=NULL;

        /* Read the keyword and place its value in the poitner. */
        fits_read_key(fptr, gal_fits_type_to_datatype(tmp->type),
                      tmp->name, valueptr, tmp->comment, &tmp->status);

        /* If the comment was empty, free the space and set it to zero. */
        if(tmp->comment && tmp->comment[0]=='\0')
          {free(tmp->comment); tmp->comment=NULL;}
      }
}





/* Same as `gal_fits_read_keywords_fptr', but accepts the filename and HDU
   as input instead of an already opened CFITSIO `fitsfile' pointer. */
void
gal_fits_read_keywords(char *filename, char *hdu, gal_data_t *keysll,
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
    error(EXIT_FAILURE, errno, "%zu characters", len);
  sprintf(ffname, "%s[%s#]", filename, hdu);

  /* Open the FITS file: */
  if( fits_open_file(&fptr, ffname, READONLY, &status) )
    gal_fits_io_error(status, "reading this FITS file");

  /* Read the keywords. */
  gal_fits_read_keywords_fptr(fptr, keysll, readcomment, readunit);

  /* Close the FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);

  /* Clean up. */
  free(ffname);
}





/* Add on keyword to the list of header keywords that need to be added
   to a FITS file. In the end, the keywords will have to be freed, so
   it is important to know before hand if they were allocated or
   not. If not, they don't need to be freed. */
void
gal_fits_add_to_key_ll(struct gal_fits_key_ll **list, int type,
                       char *keyname, int kfree, void *value, int vfree,
                       char *comment, int cfree, char *unit)
{
  struct gal_fits_key_ll *newnode;

  /* Allocate space for the new node and fill it in. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno,
          "linkedlist: new element in gal_fits_key_ll");
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
gal_fits_add_to_key_ll_end(struct gal_fits_key_ll **list, int type,
                           char *keyname, int kfree, void *value, int vfree,
                           char *comment, int cfree, char *unit)
{
  struct gal_fits_key_ll *newnode, *tmp;

  /* Allocate space for the new node and fill it in. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno,
          "linkedlist: new element in gal_fits_key_ll");
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
gal_fits_file_name_in_keywords(char *keynamebase, char *filename,
                               struct gal_fits_key_ll **list)
{
  char *keyname, *value;
  size_t numkey=1, maxlength;
  size_t i, j, len=strlen(filename), thislen;

  /* When you give string arguments, CFITSIO puts them within two ''s,
     so the actual length available is two less. It seems this length
     also didn't include the null character, so, ultimately you have
     to take three from it.*/
  maxlength=FLEN_VALUE-3;

  i=0;
  while(i<len)
    {
      /* Set the keyname: */
      errno=0;
      keyname=malloc(FLEN_KEYWORD);
      if(keyname==NULL)
        error(EXIT_FAILURE, errno, "%d bytes", FLEN_KEYWORD);
      sprintf(keyname, "%s_%zu", keynamebase, numkey++);

      /* Set the keyword value: */
      errno=0;
      thislen=strlen(&filename[i]);
      value=malloc(maxlength);
      if(value==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes", thislen);
      strncpy(value, &filename[i], maxlength);

      /* If the FROM string (=&filename[i]) in strncpy is shorter than
         SIZE (=maxlength), then the rest of the space will be filled
         with null characters. So we can use this to check if the full
         length was copied. */
      if(value[maxlength-1]=='\0')
        {
          gal_fits_add_to_key_ll_end(list, TSTRING, keyname, 1,
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
            error(EXIT_FAILURE, 0, "the filename `%sP has at least one "
                  "span of %zu characters without a `/`. It cannot be "
                  "written to the header of the output fits file",
                  filename, maxlength);

          /* Convert the last useful character and save the file name.*/
          gal_fits_add_to_key_ll_end(list, TSTRING, keyname, 1,
                                     value, 1, NULL, 0, NULL);
          i+=j+1;
        }
    }
}





/* Write the WCS and begin the part on this particular program's
   key words. */
void
gal_fits_add_wcs_to_header(fitsfile *fptr, char *wcsheader, int nkeyrec)
{
  size_t i;
  int h, status=0;
  char startblank[]="                      / ";
  char *cp, *cpf, blankrec[80], titlerec[80];

  /* Set the last element of the blank array. */
  cpf=blankrec+79;
  *cpf='\0';
  titlerec[79]='\0';
  cp=blankrec; do *cp=' '; while(++cp<cpf);

  /* Print the first two lines before the WCS header information. */
  if(fits_write_record(fptr, blankrec, &status))
    gal_fits_io_error(status, NULL);
  sprintf(titlerec, "%sWCS information", startblank);
  for(i=strlen(titlerec);i<79;++i)
    titlerec[i]=' ';
  if(fits_write_record(fptr, titlerec, &status))
    gal_fits_io_error(status, NULL);

  /* Write the keywords one by one: */
  for(h=0;h<nkeyrec-1;++h)
    fits_write_record(fptr, &wcsheader[h*80], &status);
  gal_fits_io_error(status, NULL);
}





/* Write the keywords in the gal_fits_key_ll linked list to the FITS
   file. Every keyword that is written is freed, that is why we need
   the pointer to the linked list (to correct it after we finish). */
void
gal_fits_update_keys(fitsfile *fptr, struct gal_fits_key_ll **keylist)
{
  int status=0;
  struct gal_fits_key_ll *tmp, *ttmp;

  tmp=*keylist;
  while(tmp!=NULL)
    {
      /* Write the information: */
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
      if(tmp->unit && fits_write_key_unit(fptr, tmp->keyname,
                                          tmp->unit, &status) )
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
  *keylist=NULL;
}





void
gal_fits_write_keys_version(fitsfile *fptr, struct gal_fits_key_ll *headers,
                            char *spack_string)
{
  size_t i;
  int status=0;
  char *gitdescribe;
  char cfitsioversion[20];
  char startblank[]="              / ";
  char *cp, *cpf, blankrec[80], titlerec[80];

  /* Before WCSLIB 5.0, the wcslib_version function was not
     defined. Sometime in the future were everyone has moved to more
     recent versions of WCSLIB, we can remove this macro and its check
     in configure.ac.*/
#ifdef GAL_CONFIG_HAVE_WCSLIB_VERSION
  int wcslibvers[3];
  char wcslibversion[20];
  const char *wcslibversion_const;
#endif

  /* Set the last element of the blank array. */
  cpf=blankrec+79;
  *cpf='\0';
  titlerec[79]='\0';
  cp=blankrec; do *cp=' '; while(++cp<cpf);

  /* If any header keywords are specified add them: */
  if(headers)
    {
      fits_write_record(fptr, blankrec, &status);
      sprintf(titlerec, "%s%s", startblank, spack_string);
      for(i=strlen(titlerec);i<79;++i) titlerec[i]=' ';
      fits_write_record(fptr, titlerec, &status);
      gal_fits_update_keys(fptr, &headers);
    }


  /* Start printing the version information */
  fits_write_record(fptr, blankrec, &status);
  sprintf(titlerec, "%sVersions and date", startblank);
  for(i=strlen(titlerec);i<79;++i) titlerec[i]=' ';
  fits_write_record(fptr, titlerec, &status);
  gal_fits_io_error(status, NULL);

  /* Set the version of CFITSIO as a string. */
  sprintf(cfitsioversion, "%-.2f", CFITSIO_VERSION);

  /* Write all the information: */
  fits_write_date(fptr, &status);

  /* Write the version of CFITSIO */
  fits_update_key(fptr, TSTRING, "CFITSIO", cfitsioversion,
                  "CFITSIO version.", &status);

  /* Write the WCSLIB version */
#ifdef GAL_CONFIG_HAVE_WCSLIB_VERSION
  wcslibversion_const=wcslib_version(wcslibvers);
  strcpy(wcslibversion, wcslibversion_const);
  fits_update_key(fptr, TSTRING, "WCSLIB", wcslibversion,
                  "WCSLIB version.", &status);
#endif

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




















/*************************************************************
 ***********       Read WCS from FITS pointer      ***********
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
void
gal_fits_read_wcs_from_pointer(fitsfile *fptr, int *nwcs, struct wcsprm **wcs,
                               size_t hstartwcs, size_t hendwcs)
{
  /* Declaratins: */
  int nkeys=0, status=0;
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

  /* WCSlib function */
  status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, nwcs, wcs);
  if(status)
    {
      fprintf(stderr, "\n##################\n"
              "WCSLIB Warning: wcspih ERROR %d: %s.\n"
              "##################\n",
              status, wcs_errmsg[status]);
      *wcs=NULL; *nwcs=0;
    }
  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in fitsarrayvv.c for freeing "
                           "the memory used to keep all the headers");

  /* Set the internal structure: */
  status=wcsset(*wcs);
  if(status)
    {
      fprintf(stderr, "\n##################\n"
              "WCSLIB Warning: wcsset ERROR %d: %s.\n"
              "##################\n",
            status, wcs_errmsg[status]);
      *wcs=NULL; *nwcs=0;
    }

  /* Initialize the wcsprm struct
  if ((status = wcsset(*wcs)))
    error(EXIT_FAILURE, 0, "wcsset ERROR %d: %s.\n", status,
          wcs_errmsg[status]);
  */
}





void
gal_fits_read_wcs(char *filename, char *hdu, size_t hstartwcs,
                  size_t hendwcs, int *nwcs, struct wcsprm **wcs)
{
  int status=0;
  fitsfile *fptr;

  /* Check HDU for realistic conditions: */
  fptr=gal_fits_read_hdu(filename, hdu, 0);

  /* Read the WCS information: */
  gal_fits_read_wcs_from_pointer(fptr, nwcs, wcs, hstartwcs, hendwcs);

  /* Close the FITS file: */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}

















/*************************************************************
 ***********            Array functions            ***********
 *************************************************************/

/* Read a FITS image into an array corresponding to fitstype and also
   save the size of the array.

   If the image has any null pixels, their number is returned by this
   function. The value that is placed for those pixels is defined by
   the macros in fitsarrayvv.h and depends on the type of the data.*/
gal_data_t *
gal_fits_read_img_hdu(char *filename, char *hdu, char *maskname,
                      char *mhdu, size_t minmapsize)
{
  void *blank;
  int anyblank;
  size_t i, ndim;
  fitsfile *fptr;
  int status=0, type;
  long *fpixel, *dsize, dsize_key=1;
  char **str, *name=NULL, *unit=NULL;
  gal_data_t *img, *mask, *keysll=NULL;


  /* Check HDU for realistic conditions: */
  fptr=gal_fits_read_hdu(filename, hdu, 0);


  /* Get the info and allocate the data structure. */
  gal_fits_img_info(fptr, &type, &ndim, &dsize);


  /* Check if there is any dimensions (the first header can sometimes have
     no images). */
  if(ndim==0)
    error(EXIT_FAILURE, 0, "%s (hdu: %s) has 0 dimensions! The most common "
          "cause for this is a wrongly specified HDU: in some FITS images, "
          "the first HDU doesn't have any data. So probably reading the "
          "second HDU (with `--hdu=1' or `-h1') will solve the problem. Note "
          "that currently HDU counting starts from 0." , filename, hdu);


  /* Set the fpixel array (first pixel in all dimensions): */
  errno=0;
  fpixel=malloc(ndim*sizeof *fpixel);
  if(fpixel==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for fpixel in "
          "`gal_fits_read_img_hdu'", ndim*sizeof *fpixel);
  for(i=0;i<ndim;++i) fpixel[i]=1;


  /* Read the possibly existing useful keywords. Note that the values are
     in allocated strings in the keys[i] data structures. Note that we need
     the linked list of keys to keep the `name' and `unit' pointers. We can
     free the linked list after `gal_data_alloc' has read/copied the
     values.*/
  gal_data_add_to_ll(&keysll, NULL, GAL_DATA_TYPE_STRING, 1, &dsize_key,
                        NULL, 0, -1, "EXTNAME", NULL, NULL);
  gal_data_add_to_ll(&keysll, NULL, GAL_DATA_TYPE_STRING, 1, &dsize_key,
                        NULL, 0, -1, "BUNIT", NULL, NULL);
  gal_fits_read_keywords_fptr(fptr, keysll, 0, 0);
  if(keysll->status==0)       {str=keysll->array;       unit=*str; }
  if(keysll->next->status==0) {str=keysll->next->array; name=*str; }


  /* Allocate the space for the array and for the blank values. */
  img=gal_data_alloc(NULL, type, (long)ndim, dsize, NULL, 0, minmapsize,
                     name, unit, NULL);
  blank=gal_data_alloc_blank(type);
  gal_data_free_ll(keysll);
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


  /* If a mask was specified, read it into the mask data structure, then
     set all the corresponding pixels of the input image to NaN. */
  if(maskname)
    {
      /* Read the mask HDU. */
      mask=gal_fits_read_img_hdu(maskname, mhdu, NULL, NULL, minmapsize);

      /* Apply the mask on the input. */
      gal_data_apply_mask(img, mask);

      /* Free the mask space. */
      gal_data_free(mask, 0);
    }

  /* Return the filled data structure */
  return img;
}





/* The user has specified an input file and a mask file. In the
   processing, all masked pixels will be converted to NaN pixels in
   the input image so we only have to deal with one array. Also since
   all processing is done on floating point arrays, the input is
   converted to floating point, irrespective of its input type. The
   original input bitpix will be stored so if you want to, you can
   return it back to the input type if you please. */
gal_data_t *
gal_fits_read_to_type(char *inputname, char *inhdu, char *maskname,
                      char *mhdu, int type, size_t minmapsize)
{
  gal_data_t *in, *converted;

  /* Read the specified input image HDU. */
  in=gal_fits_read_img_hdu(inputname, inhdu, maskname, mhdu, minmapsize);

  /* If the input had another type, convert it to float. */
  if(in->type!=type)
    {
      converted=gal_data_copy_to_new_type(in, type);
      gal_data_free(in, 0);
      in=converted;
    }

  /* Return the final structure. */
  return in;
}





gal_data_t *
gal_fits_read_float_kernel(char *inputname, char *inhdu, float **outkernel,
                           size_t *ins0, size_t *ins1)
{
  size_t i;
  int check=0;
  double sum=0;
  gal_data_t *kernel;
  float *f, *fp, tmp;

  /* Read the image as a float */
  kernel=gal_fits_read_to_type(inputname, inhdu, NULL, NULL,
                               GAL_DATA_TYPE_FLOAT, -1);

  /* Check if the size along each dimension of the kernel is an odd
     number. If they are all an odd number, then the for each dimension,
     check will be incremented once. */
  for(i=0;i<kernel->ndim;++i)
    check=kernel->dsize[i]%2;
  if(check!=kernel->ndim)
    error(EXIT_FAILURE, 0, "the kernel image has to have an odd number "
          "of pixels in all dimensions (there has to be one element/pixel "
          "in the center). At least one of the dimensions of %s (hdu: %s) "
          "doesn't have an odd number of pixels", inputname, inhdu);

  /* If there are any NaN pixels, set them to zero and normalize it. There
     are no blank pixels any more.*/
  fp=(f=kernel->array)+kernel->size;
  do
    {
      if(isnan(*f)) *f=0.0f;
      else          sum+=*f;
    }
  while(++f<fp);
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
gal_fits_write_img_fitsptr(gal_data_t *data, char *filename, char *extname)
{
  void *blank;
  long fpixel=1;
  fitsfile *fptr;
  char *wcsheader;
  int nkeyrec, status=0, datatype=gal_fits_type_to_datatype(data->type);

  /* Check if the file already exists. If it does, we want to add the array
     as a new extension. */
  if(access(filename,F_OK) != -1 )
    fits_open_file(&fptr,filename, READWRITE, &status);
  else
    fits_create_file(&fptr, filename, &status);

  /* Create the FITS file and put the image into it. */
  fits_create_img(fptr, gal_fits_type_to_bitpix(data->type),
                  data->ndim, data->dsize, &status);
  fits_write_img(fptr, datatype, fpixel, data->size, data->array, &status);

  /* If we have blank pixels, we need to define a BLANK keyword when we are
     dealing with integer types. */
  if(gal_data_has_blank(data))
    switch(data->type)
      {
      case GAL_DATA_TYPE_FLOAT:
      case GAL_DATA_TYPE_DOUBLE:
        /* Do nothing! Since there are much fewer floating point types
           (that don't need any BLANK keyword), we are checking them.*/
        break;

      default:
        blank=gal_data_alloc_blank(data->type);
        if(fits_write_key(fptr, datatype, "BLANK", blank,
                          "Pixels with no data.", &status) )
          gal_fits_io_error(status, "adding the BLANK keyword");
        free(blank);
      }

  /* Write the extension name to the header. */
  fits_write_key(fptr, TSTRING, "EXTNAME", extname, "", &status);
  gal_fits_io_error(status, NULL);

  /* If a WCS structure is present, write it in */
  if(data->wcs)
    {
      /* Convert the WCS information to text. */
      status=wcshdo(WCSHDO_safe, data->wcs, &nkeyrec, &wcsheader);
      if(status)
        error(EXIT_FAILURE, 0, "wcshdo ERROR %d: %s", status,
              wcs_errmsg[status]);
      gal_fits_add_wcs_to_header(fptr, wcsheader, nkeyrec);
    }

  /* Report any errors if we had any */
  gal_fits_io_error(status, NULL);
  return fptr;
}





void
gal_fits_write_img(gal_data_t *data, char *filename, char *extname,
                   struct gal_fits_key_ll *headers, char *spack_string)
{
  int status=0;
  fitsfile *fptr;

  /* Write the data array into a FITS file and keep it open: */
  fptr=gal_fits_write_img_fitsptr(data, filename, extname);

  /* Write all the headers and the version information. */
  gal_fits_write_keys_version(fptr, headers, spack_string);

  /* Close the FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}





void
gal_fits_write_img_update_crpix(gal_data_t *data, char *filename,
                                char *extname,
                                struct gal_fits_key_ll *headers,
                                double *crpix, char *spack_string)
{
  int status=0;
  fitsfile *fptr;

  /* Write the data array into a FITS file and keep it open: */
  fptr=gal_fits_write_img_fitsptr(data, filename, extname);

  /* Update the CRPIX keywords. Note that we don't want to change the
     values in the WCS information of gal_data_t. Because, it often happens
     that things are done in parallel, so we don't want to touch the
     original version, we just want to change the copied version. */
  if(crpix)
    {
      fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix[0],
                      NULL, &status);
      fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix[1],
                      NULL, &status);
      gal_fits_io_error(status, NULL);
    }

  /* Write all the headers and the version information. */
  gal_fits_write_keys_version(fptr, NULL, spack_string);

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
gal_fits_table_size(fitsfile *fitsptr, size_t *nrows, size_t *ncols)
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
gal_fits_table_type(fitsfile *fptr)
{
  int status=0;
  char value[FLEN_VALUE];

  fits_read_key(fptr, TSTRING, "XTENSION", value, NULL, &status);

  if(status==0)
    {
      if(!strcmp(value, "TABLE   "))
        return GAL_TABLE_TYPE_AFITS;
      else if(!strcmp(value, "BINTABLE"))
        return GAL_TABLE_TYPE_BFITS;
      else
        error(EXIT_FAILURE, 0, "The `XTENSION' keyword of this FITS file "
              "doesn't have a standard value (`%s')", value);
    }
  else
    {
      if(status==KEY_NO_EXIST)
        error(EXIT_FAILURE, 0, "the `gal_fits_table_type' function was "
              "called on a FITS extension which is not a table.");
      else
        gal_fits_io_error(status, NULL);
    }

  error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix it. "
        "for some reason, the control of `gal_fits_table_type' has reached "
        "the end of the function! This must not happen", PACKAGE_BUGREPORT);
  return -1;
}




static void
remove_trailing_space(char *str)
{
  size_t i;

  /* Start from the second last character (the last is a single quote) and
     go down until you hit a non-space character. */
  for(i=strlen(str)-2;i>0;--i)
    if(str[i]!=' ')
      break;

  /* If the string is empty, it will stop at the first character that is a
     single quote. So no need to check further. */
  str[i+1]='\0';
}





/* The general format of the TDISPn keywords in FITS is like this: `Tw.p',
   where `T' specifies the general format, `w' is the width to be given to
   this column and `p' is the precision. For integer types, percision is
   actually the minimum number of integers, for floats, it is the number of
   decimal digits beyond the decimal point.

 */
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
      error(EXIT_FAILURE, 0, "%s (hdu: %s): Format character `%c' in the "
            "value (%s) of the keywork %s not recognized", filename, hdu,
            tdisp[0], tdisp, keyname);
    }

  /* Parse the rest of the string to see if a width and precision are given
     or not. */
  data->disp_width=strtol(&tdisp[1], &tailptr, 0);
  switch(*tailptr)
    {
    case '.':      /* Width is set, go onto finding the precision. */
      data->disp_precision = strtol(&tailptr[1], &tailptr, 0);
      if(*tailptr!='\0')
        error(EXIT_FAILURE, 0, "%s (hdu: %s): The value `%s' of the "
              "`%s' keyword could not recognized (it doesn't finish after "
              "the precision)", filename, hdu, tdisp, keyname);
      break;

    case '\0':     /* No precision given, use a default value.     */
      data->disp_precision = ( isanint
                               ? GAL_TABLE_DEF_INT_PRECISION
                               : GAL_TABLE_DEF_FLT_PRECISION );
      break;

    default:
      error(EXIT_FAILURE, 0, "%s (hdu: %s): The value `%s' of the "
            "`%s' keyword could not recognized (it doesn't have a `.', or "
            "finish, after the width)", filename, hdu, tdisp,
            keyname);
    }


}





/* See the descriptions of `gal_table_info'. */
gal_data_t *
gal_fits_table_info(char *filename, char *hdu, size_t *numcols,
                    int *tabletype)
{
  long repeat;
  int tfields;        /* The maximum number of fields in FITS is 999 */
  size_t index;
  fitsfile *fptr;
  size_t i, numrows;
  gal_data_t *cols=NULL;
  int status=0, datatype;
  char *tailptr, keyname[FLEN_KEYWORD]="XXXXXXXXXXXXX", value[FLEN_VALUE];


  /* Open the FITS file and get the basic information. */
  fptr=gal_fits_read_hdu(filename, hdu, 1);
  *tabletype=gal_fits_table_type(fptr);
  gal_fits_table_size(fptr, &numrows, numcols);


  /* Read the total number of fields, then allocate space for the data
     structure and store the information within it. */
  fits_read_key(fptr, TINT, "TFIELDS", &tfields, NULL, &status);
  errno=0;
  cols=malloc(tfields*sizeof *cols);
  if(cols==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for cols in `gal_fits_table_info'",
          tfields*sizeof *cols);


  /* Save the number of rows as the data structure size and also length
     along the first (and only) dimension. */
  for(i=0;i<*numcols;++i) cols[i].size=numrows;


  /* Read all the keywords one by one and if they match, then put them in
     the correct value. Note that we are starting from keyword 9 because
     according to the FITS standard, the first 8 keys in a FITS table are
     reserved. */
  for(i=9; strcmp(keyname, "END"); ++i)
    {
      /* Read the next keyword. */
      fits_read_keyn(fptr, i, keyname, value, NULL, &status);

      /* COLUMN DATA TYPE. According the the FITS standard, the value of
         TFORM is most generally in this format: `rTa'. `T' is actually a
         code of the datatype. `r' is the `repeat' counter and `a' is
         depreciated. Currently we can only read repeat==1 cases. When no
         number exists before the defined capital letter, it defaults to 1,
         but if a number exists (for example `5D'), then the repeat is 5
         (there are actually five values in each column). Note that
         value[0] is a single quote.*/
      if(strncmp(keyname, "TFORM", 5)==0)
        {
          /* Remove the ending trailing space and quotation sign. */
          remove_trailing_space(value);
          if(*tabletype==GAL_TABLE_TYPE_AFITS)
            fits_ascii_tform(&value[1], &datatype, NULL, NULL, &status);
          else
            fits_binary_tform(&value[1], &datatype, &repeat, NULL, &status);

          /* Small sanity check. */
          if(repeat>1)
             error(EXIT_FAILURE, 0, "The repeat value of %s is %ld, "
                   "currently we can only use columns with a repeat "
                   "of 1. Please get in touch with us at %s to add this "
                   "feature", keyname, repeat, PACKAGE_BUGREPORT);

          /* See which column this information was for and add it. In the
             meantime, also do a sanity check. */
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)     /* Counting from zero was corrected above. */
            cols[index].type=gal_fits_datatype_to_type(datatype);
        }

      /* COLUMN NAME. All strings in CFITSIO start and finish with single
         quotation marks, CFITSIO puts them in itsself, so if we don't
         remove them here, we might have duplicates later, its easier to
         just remove them to have a simple string that might be used else
         where too (without the single quotes).*/
      else if(strncmp(keyname, "TTYPE", 5)==0)
        {
          remove_trailing_space(value);
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            gal_checkset_allocate_copy(&value[1], &cols[index].name);
        }

      /* COLUMN UNITS. */
      else if(strncmp(keyname, "TUNIT", 5)==0)
        {
          /* similar to tname, see above.*/
          remove_trailing_space(value);
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            gal_checkset_allocate_copy(&value[1], &cols[index].unit);
        }

      /* COLUMN COMMENTS */
      else if(strncmp(keyname, "TCOMM", 5)==0)
        {
          /* similar to tname, see above.*/
          remove_trailing_space(value);
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            gal_checkset_allocate_copy(&value[1], &cols[index].comment);
        }

      /* COLUMN DISPLAY FORMAT */
      else if(strncmp(keyname, "TDISP", 5)==0)
        {
          /* similar to tname, see above.*/
          remove_trailing_space(value);
          index = strtoul(&keyname[5], &tailptr, 10) - 1;
          if(index<tfields)
            set_display_format(&value[1], &cols[index], filename, hdu,
                               keyname);
        }
    }

  /* Close the FITS file and report an error if we had any. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
  return cols;
}





/* Read the column indexs given in the `indexll' linked list from a FITS
   table into a linked list of data structures, note that this is a
   low-level function, so the output data linked list is the inverse of the
   input indexs linked list. You can use */
gal_data_t *
gal_fits_table_read(char *filename, char *hdu, gal_data_t *colinfo,
                    struct gal_linkedlist_sll *indexll, int minmapsize)
{
  size_t ind;
  void *blank;
  long dsize[1];
  fitsfile *fptr;
  int status=0, anynul;
  gal_data_t *out=NULL, *col;

  /* Open the FITS file */
  fptr=gal_fits_read_hdu(filename, hdu, 1);

  /* Pop each index and read/store the array. */
  while(indexll!=NULL)
    {
      /* Pop the index. */
      gal_linkedlist_pop_from_sll(&indexll, &ind);

      /* Allocate the necessary data structure (including the array) for
         this column. */
      dsize[0]=colinfo[ind].size;
      col=gal_data_alloc(NULL, colinfo[ind].type, 1, dsize, NULL, 0,
                         minmapsize, colinfo[ind].name, colinfo[ind].unit,
                         colinfo[ind].comment);

      /* Allocate a blank value for the give type and read/store the
         column using CFITSIO. Afterwards, free the blank value. */
      blank=gal_data_alloc_blank(col->type);
      fits_read_col(fptr, gal_fits_type_to_datatype(col->type), ind+1, 1, 1,
                    col->size, blank, col->array, &anynul, &status);
      free(blank);

      /* Add the column to the final list of data structures. */
      col->next=out;
      out=col;
    }

  /* Close the FITS file */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
  return out;
}





/* Write the given columns (a linked list of `gal_data_t') into a FITS
   table.*/
void
gal_fits_table_write(gal_data_t *cols, char *comments, int tabletype,
                     char *filename, int dontdelete)
{

}
