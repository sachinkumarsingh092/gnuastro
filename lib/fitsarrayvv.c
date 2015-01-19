/*********************************************************************
Functions to convert a FITS array to a C array and vice versa.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include "config.h"
#include "checkset.h"
#include "fitsarrayvv.h"
#include "fixedstringmacros.h"


/*************************************************************
 **************        Reporting errors:       ***************
 *************************************************************/
void
fitsioerror(int status, char *message)
{
  char defmessage[]="Error in CFITSIO, see above.";
  if(status)
    {
      fits_report_error(stderr, status);
      if(message)
	error(EXIT_FAILURE, 0, message);
      else
	error(EXIT_FAILURE, 0, defmessage);
    }
}




















/*************************************************************
 **************      Acceptable FITS names     ***************
 *************************************************************/
int
nameisfits(char *name)
{
  size_t len;
  len=strlen(name);
  if (strcmp(&name[len-5], ".fits") == 0
      || strcmp(&name[len-8], ".fits.gz") == 0
      || strcmp(&name[len-7], ".fits.Z") == 0
      || strcmp(&name[len-4], ".imh") == 0)
    return 1;
  else
    return 0;
}




















/*************************************************************
 **************      BITPIX Dependancies       ***************
 *************************************************************/
void
imgbitpixsize(fitsfile *fptr, int *bitpix, long *naxes)
{
  int status=0, maxdim=10, naxis;

  if( fits_get_img_param(fptr, maxdim, bitpix, &naxis, naxes, &status) )
    fitsioerror(status, NULL);

  if(naxis!=2)
    error(EXIT_FAILURE, 0, "Currently only a 2 dimentional image array "
	  "is supported. Your array is %d dimention(s). %s", naxis,
	  naxis ? "Please contact us to add this feature." : "");
}





/* Set datatype (in CFITSIO) based on BITPIX. */
int
bitpixtodtype(int bitpix)
{
  switch(bitpix)
    {
    case BYTE_IMG:
      return TBYTE;
    case SHORT_IMG:
      return TSHORT;
    case LONG_IMG:
      return TLONG;
    case LONGLONG_IMG:
      return TLONGLONG;
    case FLOAT_IMG:
      return TFLOAT;
    case DOUBLE_IMG:
      return TDOUBLE;
    default:
      error(EXIT_FAILURE, 0, "Bitpix value of %d not recognized.",
	    bitpix);
    }
  return 0;
}





void *
bitpixnull(int bitpix)
{
  uint8_t *b;
  int16_t *s;
  int32_t *l;
  int64_t *L;
  float *f;
  double *d;

  errno=0;

  switch(bitpix)
    {
    case BYTE_IMG:
      b=malloc(sizeof(uint8_t));
      if(b==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(uint8_t));
      *b=FITSBYTENUL;
      return b;

    case SHORT_IMG:
      s=malloc(sizeof(uint16_t));
      if(s==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(int16_t));
      *s=FITSSHORTNUL;
      return s;

    case LONG_IMG:
      l=malloc(sizeof(uint32_t));
      if(l==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(int32_t));
      *l=FITSLONGNUL;
      return l;

    case LONGLONG_IMG:
      L=malloc(sizeof(uint64_t));
      if(L==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(int64_t));
      *L=FITSLLONGNUL;
      return L;

    case FLOAT_IMG:
      f=malloc(sizeof(float));
      if(f==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(float));
      *f=FITSFLOATNUL;
      return f;

    case DOUBLE_IMG:
      d=malloc(sizeof(double));
      if(d==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(double));
      *d=FITSFLOATNUL;
      return d;

    default:
      error(EXIT_FAILURE, 0, "Bitpix value of %d not recognized.",
	    bitpix);
    }

  return NULL;
}





/* Allocate an array based on the value of bitpix. */
void *
bitpixalloc(size_t size, int bitpix)
{
  void *array;

  /* Allocate space for the array to keep the image. */
  switch(bitpix)
    {
    case BYTE_IMG:
      size*=sizeof(uint8_t);
      break;

    case SHORT_IMG:
      size*=sizeof(int16_t);
      break;

    case LONG_IMG:
      size*=sizeof(int32_t);
      break;

    case LONGLONG_IMG:
      size*=sizeof(int64_t);
      break;

    case FLOAT_IMG:
      if(sizeof(float)!=4)
	error(EXIT_FAILURE, 0, "`float` is not 32bits on this machine.");
      size*=sizeof(float);
      break;

    case DOUBLE_IMG:
      if(sizeof(float)!=8)
	error(EXIT_FAILURE, 0, "`double` is not 64bits on this machine.");
      size*=sizeof(double);
      break;

    default:
      error(EXIT_FAILURE, 0, "Bitpix value of %d not recognized.", bitpix);
    }

  errno=0;
  array=malloc(size);
  if(array==NULL)
    error(EXIT_FAILURE, errno, "Array of %lu bytes", size);

  return array;
}





void
nultovalue(void *array, int bitpix, size_t size, void *value)
{
  uint8_t *b, *bf, bv=*(uint8_t *) value; /* Value will only be read from */
  int16_t *s, *sf, sv=*(int16_t *) value; /* one of these based on bitpix.*/
  int32_t *l, *lf, lv=*(int32_t *) value; /* Which the caller assigned.   */
  int64_t *L, *Lf, Lv=*(int64_t *) value; /* If there is any problem, it  */
  float   *f, *ff, fv=*(float   *) value; /* is their responsability, not */
  double  *d, *df, dv=*(double  *) value; /* this functions :-D.          */

  switch(bitpix)
    {
    case BYTE_IMG:
      bf=(b=array)+size;
      do if(*b==FITSBYTENUL) *b=bv; while(++b<bf);
      break;

    case SHORT_IMG:
      sf=(s=array)+size;
      do if(*s==FITSSHORTNUL) *s=sv; while(++s<sf);
      break;

    case LONG_IMG:
      lf=(l=array)+size;
      do if(*l==FITSLONGNUL) *l=lv; while(++l<lf);
      break;

    case LONGLONG_IMG:
      Lf=(L=array)+size;
      do if(*L==FITSLLONGNUL) *L=Lv; while(++L<Lf);
      break;

    case FLOAT_IMG:
      ff=(f=array)+size;
      do if(*f==FITSFLOATNUL) *f=fv; while(++f<ff);
      break;

    case DOUBLE_IMG:
      df=(d=array)+size;
      do if(*d==FITSFLOATNUL) *d=dv; while(++d<df);
      break;

    default:
      error(EXIT_FAILURE, 0, "A bug! Bitpix value of %d not recognized. "
	    "This should not happen here (convertnul in fitsarrayvv.c). "
	    "Please contact us to see how this happened.", bitpix);
    }
}





void
changetype(void *in, int inbitpix, size_t size, void **out, int outbitpix)
{
  uint8_t *b, *bf, *ob=in;
  int16_t *s, *sf, *os=in;
  int32_t *l, *lf, *ol=in;
  int64_t *L, *Lf, *oL=in;
  float *f, *ff, *of=in;
  double *d, *df, *od=in;

  /* Allocate space for the output and start filling it. */
  *out=bitpixalloc(size, outbitpix);
  switch(outbitpix)
    {
    case BYTE_IMG:
      switch(inbitpix)
	{
	case BYTE_IMG:
	  bf=(b=*out)+size; do *b=*ob++; while(++b<bf); return;
	case SHORT_IMG:
	  sf=(s=*out)+size; do *s=*ob++; while(++s<sf); return;
	case LONG_IMG:
	  lf=(l=*out)+size; do *l=*ob++; while(++l<lf); return;
	case LONGLONG_IMG:
	  Lf=(L=*out)+size; do *L=*ob++; while(++L<Lf); return;
	case FLOAT_IMG:
	  ff=(f=*out)+size; do *f=*ob++; while(++f<ff); return;
	case DOUBLE_IMG:
	  df=(d=*out)+size; do *d=*ob++; while(++d<df); return;
	default:
	  error(EXIT_FAILURE, 0, "A bug! In changetype (fitsarrayvv.c). "
		"BITPIX=%d of input not recognized. Please contact us so "
		"we can fix it.", inbitpix);
	}
      break;

    case SHORT_IMG:
      switch(inbitpix)
	{
	case BYTE_IMG:
	  bf=(b=*out)+size; do *b=*os++; while(++b<bf); return;
	case SHORT_IMG:
	  sf=(s=*out)+size; do *s=*os++; while(++s<sf); return;
	case LONG_IMG:
	  lf=(l=*out)+size; do *l=*os++; while(++l<lf); return;
	case LONGLONG_IMG:
	  Lf=(L=*out)+size; do *L=*os++; while(++L<Lf); return;
	case FLOAT_IMG:
	  ff=(f=*out)+size; do *f=*os++; while(++f<ff); return;
	case DOUBLE_IMG:
	  df=(d=*out)+size; do *d=*os++; while(++d<df); return;
	default:
	  error(EXIT_FAILURE, 0, "A bug! In changetype (fitsarrayvv.c). "
		"BITPIX=%d of input not recognized. Please contact us so "
		"we can fix it.", inbitpix);
	}
      break;

    case LONG_IMG:
      switch(inbitpix)
	{
	case BYTE_IMG:
	  bf=(b=*out)+size; do *b=*ol++; while(++b<bf); return;
	case SHORT_IMG:
	  sf=(s=*out)+size; do *s=*ol++; while(++s<sf); return;
	case LONG_IMG:
	  lf=(l=*out)+size; do *l=*ol++; while(++l<lf); return;
	case LONGLONG_IMG:
	  Lf=(L=*out)+size; do *L=*ol++; while(++L<Lf); return;
	case FLOAT_IMG:
	  ff=(f=*out)+size; do *f=*ol++; while(++f<ff); return;
	case DOUBLE_IMG:
	  df=(d=*out)+size; do *d=*ol++; while(++d<df); return;
	default:
	  error(EXIT_FAILURE, 0, "A bug! In changetype (fitsarrayvv.c). "
		"BITPIX=%d of input not recognized. Please contact us so "
		"we can fix it.", inbitpix);
	}
      break;

    case LONGLONG_IMG:
      switch(inbitpix)
	{
	case BYTE_IMG:
	  bf=(b=*out)+size; do *b=*oL++; while(++b<bf); return;
	case SHORT_IMG:
	  sf=(s=*out)+size; do *s=*oL++; while(++s<sf); return;
	case LONG_IMG:
	  lf=(l=*out)+size; do *l=*oL++; while(++l<lf); return;
	case LONGLONG_IMG:
	  Lf=(L=*out)+size; do *L=*oL++; while(++L<Lf); return;
	case FLOAT_IMG:
	  ff=(f=*out)+size; do *f=*oL++; while(++f<ff); return;
	case DOUBLE_IMG:
	  df=(d=*out)+size; do *d=*oL++; while(++d<df); return;
	default:
	  error(EXIT_FAILURE, 0, "A bug! In changetype (fitsarrayvv.c). "
		"BITPIX=%d of input not recognized. Please contact us so "
		"we can fix it.", inbitpix);
	}
      break;

    case FLOAT_IMG:
      switch(inbitpix)
	{
	case BYTE_IMG:
	  bf=(b=*out)+size; do *b=*of++; while(++b<bf); return;
	case SHORT_IMG:
	  sf=(s=*out)+size; do *s=*of++; while(++s<sf); return;
	case LONG_IMG:
	  lf=(l=*out)+size; do *l=*of++; while(++l<lf); return;
	case LONGLONG_IMG:
	  Lf=(L=*out)+size; do *L=*of++; while(++L<Lf); return;
	case FLOAT_IMG:
	  ff=(f=*out)+size; do *f=*of++; while(++f<ff); return;
	case DOUBLE_IMG:
	  df=(d=*out)+size; do *d=*of++; while(++d<df); return;
	default:
	  error(EXIT_FAILURE, 0, "A bug! In changetype (fitsarrayvv.c). "
		"BITPIX=%d of input not recognized. Please contact us so "
		"we can fix it.", inbitpix);
	}
      break;

    case DOUBLE_IMG:
      switch(inbitpix)
	{
	case BYTE_IMG:
	  bf=(b=*out)+size; do *b=*od++; while(++b<bf); return;
	case SHORT_IMG:
	  sf=(s=*out)+size; do *s=*od++; while(++s<sf); return;
	case LONG_IMG:
	  lf=(l=*out)+size; do *l=*od++; while(++l<lf); return;
	case LONGLONG_IMG:
	  Lf=(L=*out)+size; do *L=*od++; while(++L<Lf); return;
	case FLOAT_IMG:
	  ff=(f=*out)+size; do *f=*od++; while(++f<ff); return;
	case DOUBLE_IMG:
	  df=(d=*out)+size; do *d=*od++; while(++d<df); return;
	default:
	  error(EXIT_FAILURE, 0, "A bug! In changetype (fitsarrayvv.c). "
		"BITPIX=%d of input not recognized. Please contact us so "
		"we can fix it.", inbitpix);
	}
      break;


    default:
      error(EXIT_FAILURE, 0, "A bug! Output Bitpix value of %d is not "
	    "recognized. This should not happen here (changetype in "
	    "fitsarrayvv.c). Please contact us to see how this happened.",
	    outbitpix);
    }
}
















/*************************************************************
 **************      Number of extensions:     ***************
 *************************************************************/
void
numhduensions(char *filename, int *numhdu)
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

  fitsioerror(status, NULL);
}




















/**************************************************************/
/**********       Check FITS image HDUs:      ************/
/**************************************************************/
char *
hdutypestring(int hdutype)
{
  switch(hdutype)
    {
    case IMAGE_HDU:
      return "an Image";
    case ASCII_TBL:
      return "an ASCII table";
    case BINARY_TBL:
      return "a binary table";
      break;
    default:
      error(EXIT_FAILURE, 0, "HDU code %d in CFITSIO not recognized.",
	    hdutype);
    }
  return NULL;
}





/* Check the desired HDU in a FITS image and also if it has the
   desired type. */
void
readfitshdu(char *filename, char *hdu, int desiredtype, fitsfile **outfptr)
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
    error(EXIT_FAILURE, errno, "%lu characters", len);
  sprintf(ffname, "%s[%s#]", filename, hdu);

  /* Open the FITS file: */
  if( fits_open_file(outfptr, ffname, READONLY, &status) )
    fitsioerror(status, "Reading this FITS file.");
  fptr=*outfptr;

  /* Check the Type of the given HDU: */
  if (fits_get_hdu_type(fptr, &hdutype, &status) )
    fitsioerror(status, NULL);

  if(hdutype!=desiredtype)
    error(EXIT_FAILURE, 0, "%s: HDU %s is %s, not %s.",
	  filename, hdu, hdutypestring(hdutype),
	  hdutypestring(desiredtype));

  free(ffname);
}




















/*************************************************************
 ******************         Header          ******************
 *************************************************************/
void
add_to_fitsheaderllend(struct fitsheaderll **list, int datatype,
		       char *keyname, int kfree, void *value, int vfree,
		       char *comment, int cfree, char *unit)
{
  struct fitsheaderll *newnode, *tmp;

  /* Allocate space for the new node and fill it in. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno,
	  "linkedlist: New element in fitsheaderll");
  newnode->datatype=datatype;
  newnode->keyname=keyname;
  newnode->value=value;
  newnode->comment=comment;
  newnode->unit=unit;
  newnode->kfree=kfree;		/* Free pointers after using them. */
  newnode->vfree=vfree;
  newnode->cfree=cfree;

  if(*list)	 /* List is already full, add this node to the end */
    {
      /* After this line, tmp points to the last node. */
      tmp=*list; while(tmp->next!=NULL) tmp=tmp->next;
      tmp->next=newnode;
      newnode->next=NULL;
    }
  else		 /* List is empty */
    {
      newnode->next=*list;
      *list=newnode;
    }
}





void
filenameinkeywords(char *keynamebase, char *filename,
		  struct fitsheaderll **list)
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
      sprintf(keyname, "%s_%lu", keynamebase, numkey++);

      /* Set the keyword value: */
      errno=0;
      thislen=strlen(&filename[i]);
      value=malloc(maxlength);
      if(value==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", thislen);
      strncpy(value, &filename[i], maxlength);

      /* If the FROM string (=&filename[i]) in strncpy is shorter than
	 SIZE (=maxlength), then the rest of the space will be filled
	 with null characters. So we can use this to check if the full
	 length was copied. */
      if(value[maxlength-1]=='\0')
	{
	  add_to_fitsheaderllend(list, TSTRING, keyname, 1, value, 1,
				 NULL, 0, NULL);
	  break;
	}
      else
	{
	  /* Find the last place in the copied array that contains a
	     '/' and put j on the next character (so it can be turned
	     into a null character.*/
	  for(j=maxlength-1;j>0;--j)
	    if(value[j]=='/')
	      {
		value[j+1]='\0';
		break;
	      }
	  if(j==0)
	    error(EXIT_FAILURE, 0, "The filename `%sP has at least one span "
		  "of %lu characters without a `/`. It cannot be written "
		  "to the header of the output fits file.", filename,
		  maxlength);

	  /* Convert the last useful character and save the file name.*/
	  add_to_fitsheaderllend(list, TSTRING, keyname, 1, value, 1,
				 NULL, 0, NULL);
	  i+=j+1;
	}
    }
}





/* Write the WCS and begin the part on this particular program's
   key words. */
void
addwcstoheader(fitsfile *fptr, char *wcsheader, int nkeyrec)
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
    fitsioerror(status, NULL);
  sprintf(titlerec, "%sWCS information", startblank);
  for(i=strlen(titlerec);i<79;++i)
    titlerec[i]=' ';
  if(fits_write_record(fptr, titlerec, &status))
    fitsioerror(status, NULL);

  /* Write the keywords one by one: */
  for(h=0;h<nkeyrec-1;++h)
    fits_write_record(fptr, &wcsheader[h*80], &status);
  fitsioerror(status, NULL);
}





/* Write the keywords in the fitsheaderll linked list to the FITS
   file. Every keyword that is written is freed, that is why we need
   the pointer to the linked list (to correct it after we finish). */
void
updatekeys(fitsfile *fptr, struct fitsheaderll **keylist)
{
  int status=0;
  struct fitsheaderll *tmp, *ttmp;

  tmp=*keylist;
  while(tmp!=NULL)
    {
      /* Write the information: */
      if( fits_update_key(fptr, tmp->datatype, tmp->keyname, tmp->value,
			  tmp->comment, &status) )
	fitsioerror(status, NULL);
      if(tmp->unit && fits_write_key_unit(fptr, tmp->keyname,
					  tmp->unit, &status) )
	fitsioerror(status, NULL);

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
copyrightandend(fitsfile *fptr, char *spack_string)
{
  size_t i;
  int status=0;
  char version[20];
  char startblank[]="              / ";
  char *cp, *cpf, blankrec[80], titlerec[80];

  /* Set the last element of the blank array. */
  cpf=blankrec+79;
  *cpf='\0';
  titlerec[79]='\0';
  cp=blankrec; do *cp=' '; while(++cp<cpf);

  /*Print the starting information for the header.  */
  fits_write_record(fptr, blankrec, &status);
  sprintf(titlerec, "%s%s:", startblank, spack_string);
  for(i=strlen(titlerec);i<79;++i) titlerec[i]=' ';
  fits_write_record(fptr, titlerec, &status);
  fitsioerror(status, NULL);

  /* Set the version of CFITSIO as a string. */
  sprintf(version, "%-.2f", CFITSIO_VERSION);
  for(i=0;i<strlen(version);++i)
    if(version[i]==' ')
      {
	version[i]='\0';
	break;
      }

  /* Write all the information: */
  fits_write_date(fptr, &status);
  fits_update_key(fptr, TSTRING, "CFITSIO", version,
		  "Version of CFITSIO used.", &status);
  fits_write_comment(fptr, PACKAGE_STRING, &status);
  fits_write_comment(fptr, PACKAGE_URL, &status);
  fits_write_comment(fptr, SHORTCOPYRIGHT, &status);
  fits_write_comment(fptr, SHORTLICENSE, &status);
  fitsioerror(status, NULL);
}


















/*************************************************************
 ***********         FITS to array functions:      ***********
 *************************************************************/
/* Read the WCS information from the header. Unfortunately, WCS lib is
   not thread safe, so it needs a mutex. In case you are not using
   multiple threads, just pass a NULL pointer as the mutex.

   After you finish with this WCS, you should free the space with:

   status = wcsvfree(&nwcs,&wcs);

   ===================================
   WARNING: wcspih IS NOT THREAD SAFE!
   ===================================
   Don't call this function within a thread or use a mutex.
*/
void
readwcs(fitsfile *fptr, int *nwcs, struct wcsprm **wcs)
{
  /* Declaratins: */
  char *fullheader;
  int nkeys=0, status=0;
  int relax    = WCSHDR_all; /* Macro: use all informal WCS extensions. */
  int ctrl     = 0;          /* Don't report why a keyword wasn't used. */
  int nreject  = 0;          /* Number of keywords rejected for syntax. */

  /* CFITSIO function: */
  if( fits_hdr2str(fptr, 1, NULL, 0, &fullheader, &nkeys, &status) )
    fitsioerror(status, NULL);

  /* WCSlib function */
  if (wcspih(fullheader, nkeys, relax, ctrl, &nreject, nwcs, wcs) )
    error(EXIT_FAILURE, 0, "wcspih ERROR %d: %s.",
	  status, wcs_errmsg[status]);

  /* Set the internal structure: */
  if (wcsset(*wcs))
    error(EXIT_FAILURE, 0, "wcsset ERROR %d: %s.",
	  status, wcs_errmsg[status]);

  /* Initialize the wcsprm struct
  if ((status = wcsset(*wcs)))
    error(EXIT_FAILURE, 0, "wcsset ERROR %d: %s.\n", status,
	  wcs_errmsg[status]);
  */
  free(fullheader);
}





/* Read a FITS image into an array corresponding to fitstype and also
   save the size of the array. If the image has any null pixels, their
   number is returned by this function. */
int
fitsimgtoarray(char *filename, char *hdu, void *bitnul, int *bitpix,
	       void **array, size_t *s0, size_t *s1)
{
  fitsfile *fptr;
  long naxes[2], fpixel[]={1,1};
  int status=0, anynul=0, freebitnul=0;

  /* Check HDU for realistic conditions: */
  readfitshdu(filename, hdu, IMAGE_HDU, &fptr);

  /* Get the bitpix and size of the image: */
  imgbitpixsize(fptr, bitpix, naxes);
  *s0=naxes[1];
  *s1=naxes[0];

  /* Allocate space for the array. */
  *array=bitpixalloc(*s0 * *s1, *bitpix);

  /* Read the image into the allocated array: */
  if(bitnul==NULL)
    {
      bitnul=bitpixnull(*bitpix);
      freebitnul=1;
    }
  if ( fits_read_pix(fptr, bitpixtodtype(*bitpix), fpixel, *s0 * *s1,
		     bitnul, *array, &anynul, &status) )
    fitsioerror(status, NULL);
  if(freebitnul) free(bitnul);

  /* Close the FITS file: */
  fits_close_file(fptr, &status);
  fitsioerror(status, NULL);

  /* Return the number of "blank" pixels. */
  return anynul;
}




















/*************************************************************
 ******************      Array to FITS      ******************
 *************************************************************/
void
arraytofitsimg(char *filename, char *hdu, int bitpix, void *array,
	       size_t s0, size_t s1, struct wcsprm *wcs,
	       char *spack_string)
{
  int nkeyrec;
  fitsfile *fptr;
  char *wcsheader;
  int status=0, datatype;
  long fpixel=1, naxis=2, nelements, naxes[]={s1,s0};

  datatype=bitpixtodtype(bitpix);
  nelements=naxes[0]*naxes[1];

  if(access(filename,F_OK) != -1 )
    fits_open_file(&fptr,filename, READWRITE, &status);
  else
    fits_create_file(&fptr, filename, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  fits_write_img(fptr, datatype, fpixel, nelements, array, &status);

  fits_write_key(fptr, TSTRING, "EXTNAME", hdu, "", &status);
  fitsioerror(status, NULL);

  /* Delete the comments that CFITSIO puts immediately after making
     the file, we will add the CFITSIO version in the end so the user
     can be assured that everything was done with CFITSIO whcih: */
  fits_delete_key(fptr, "COMMENT", &status);
  fits_delete_key(fptr, "COMMENT", &status);
  fitsioerror(status, NULL);

  if(wcs)
    {
      /* Convert the WCS information to text. */
      status=wcshdo(WCSHDO_safe, wcs, &nkeyrec, &wcsheader);
      if(status)
	error(EXIT_FAILURE, 0, "wcshdo ERROR %d: %s.", status,
	      wcs_errmsg[status]);
      addwcstoheader(fptr, wcsheader, nkeyrec);
    }
  copyrightandend(fptr, spack_string);

  fits_close_file(fptr, &status);
  fitsioerror(status, NULL);
}





void
atofcorrectwcs(char *filename, char *hdu, int bitpix, void *array,
	       size_t s0, size_t s1, char *wcsheader, int wcsnkeyrec,
	       double *crpix, char *spack_string)
{
  fitsfile *fptr;
  int status=0, datatype;
  long fpixel=1, naxis=2, nelements, naxes[]={s1,s0};

  datatype=bitpixtodtype(bitpix);
  nelements=naxes[0]*naxes[1];

  if(access(filename,F_OK) != -1 )
    fits_open_file(&fptr,filename, READWRITE, &status);
  else
    fits_create_file(&fptr, filename, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  fits_write_img(fptr, datatype, fpixel, nelements, array, &status);

  fits_write_key(fptr, TSTRING, "EXTNAME", hdu, "", &status);
  fitsioerror(status, NULL);

  fits_delete_key(fptr, "COMMENT", &status);
  fits_delete_key(fptr, "COMMENT", &status);
  fitsioerror(status, NULL);

  if(wcsheader)
    {
      addwcstoheader(fptr, wcsheader, wcsnkeyrec);
      if(crpix)
	{
	  fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix[0], NULL, &status);
	  fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix[1], NULL, &status);
	  fitsioerror(status, NULL);
	}
    }

  copyrightandend(fptr, spack_string);

  fits_close_file(fptr, &status);
  fitsioerror(status, NULL);
}
