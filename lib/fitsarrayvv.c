/*********************************************************************
Functions to convert a FITS array to a C array and vice versa.
This is part of GNU Astronomy Utilities (gnuastro) package.

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
	error(EXIT_FAILURE, 0, "%s", message);
      else
	error(EXIT_FAILURE, 0, "%s", defmessage);
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
  if ( ( len>=4 && strcmp(&name[len-4], "fits") == 0 )
       || ( len>=7 && strcmp(&name[len-7], "fits.gz") == 0 )
       || ( len>=6 && strcmp(&name[len-6], "fits.Z") == 0 )
       || ( len>=3 && strcmp(&name[len-3], "imh") == 0 ) )
    return 1;
  else
    return 0;
}





int
nameisfitssuffix(char *name)
{
 if (strcmp(name, "fits") == 0 || strcmp(name, ".fits") == 0
      || strcmp(name, "fits.gz") == 0 || strcmp(name, ".fits.gz") == 0
      || strcmp(name, "fits.Z") == 0 || strcmp(name, ".fits.Z") == 0
      || strcmp(name, "imh") == 0 || strcmp(name, ".imh") == 0)
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
    error(EXIT_FAILURE, 0, "Currently only a 2 dimensional image array "
	  "is supported. Your array is %d dimension(s). %s", naxis,
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
  unsigned char *b;
  short *s;
  long *l;
  LONGLONG *L;
  float *f;
  double *d;

  errno=0;
  switch(bitpix)
    {
    case BYTE_IMG:
      b=malloc(sizeof(unsigned char));
      if(b==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(unsigned char));
      *b=FITSBYTENUL;
      return b;

    case SHORT_IMG:
      s=malloc(sizeof(short));
      if(s==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(short));
      *s=FITSSHORTNUL;
      return s;

    case LONG_IMG:
      l=malloc(sizeof(long));
      if(l==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(long));
      *l=FITSLONGNUL;
      return l;

    case LONGLONG_IMG:
      L=malloc(sizeof(LONGLONG));
      if(L==NULL)
	error(EXIT_FAILURE, errno, "%lu bytes", sizeof(LONGLONG));
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
      size*=sizeof(unsigned char);
      break;

    case SHORT_IMG:
      size*=sizeof(short);
      break;

    case LONG_IMG:
      size*=sizeof(long);
      break;

    case LONGLONG_IMG:
      size*=sizeof(LONGLONG);
      break;

    case FLOAT_IMG:
      if(sizeof(float)!=4)
	error(EXIT_FAILURE, 0,
              "`float` is not 32bits on this machine. The FITS standard "
              "Requires this size.");
      size*=sizeof(float);
      break;

    case DOUBLE_IMG:
      if(sizeof(double)!=8)
	error(EXIT_FAILURE, 0,
              "`double` is not 64bits on this machine. The FITS standard "
              "requires this size.");
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
  unsigned char *b, *bf, bv=*(uint8_t *) value; /* Value will only be read*/
  short *s, *sf, sv=*(int16_t *) value;   /* from one of these based on    */
  long *l, *lf, lv=*(int32_t *) value; /* bitpix. Which the caller assigned.*/
  LONGLONG *L, *Lf, Lv=*(int64_t *) value; /* If there is any problem, it  */
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
changetype(void *in, int inbitpix, size_t size, size_t numnul,
           void **out, int outbitpix)
{
  size_t i=0;
  unsigned char *b, *bf, *ib=in, *iib=in;
  short *s, *sf, *is=in, *iis=in;
  long *l, *lf, *il=in, *iil=in;
  LONGLONG *L, *Lf, *iL=in, *iiL=in;
  float *f, *ff, *iif=in, *iiif=in;
  double *d, *df, *id=in, *iid=in;

  /* Allocate space for the output and start filling it. */
  *out=bitpixalloc(size, outbitpix);
  switch(outbitpix)
    {
    case BYTE_IMG:
      switch(inbitpix)
	{
	case BYTE_IMG:
	  bf=(b=*out)+size; do *b=*ib++; while(++b<bf); return;
	case SHORT_IMG:
	  bf=(b=*out)+size; do *b=*is++; while(++b<bf);
          if(numnul)
            {b=*out; do {b[i]=(iis[i]==FITSSHORTNUL)?FITSBYTENUL:b[i];}
              while(++i!=size);}
          return;
	case LONG_IMG:
	  bf=(b=*out)+size; do *b=*il++; while(++b<bf);
          if(numnul)
            {b=*out; do {b[i]=(iil[i]==FITSLONGNUL)?FITSBYTENUL:b[i];}
              while(++i!=size);}
          return;
	case LONGLONG_IMG:
	  bf=(b=*out)+size; do *b=*iL++; while(++b<bf);
          if(numnul)
            {b=*out; do {b[i]=(iiL[i]==FITSLLONGNUL)?FITSBYTENUL:b[i];}
              while(++i!=size);}
          return;
	case FLOAT_IMG:
	  bf=(b=*out)+size; do *b=roundf(*iif++); while(++b<bf);
          if(numnul)
            {b=*out; do {b[i]=isnan(iiif[i])?FITSBYTENUL:b[i];}
              while(++i!=size);}
          return;
	case DOUBLE_IMG:
	  bf=(b=*out)+size; do *b=round(*id++); while(++b<bf);
          if(numnul)
            {b=*out; do {b[i]=isnan(iid[i])?FITSBYTENUL:b[i];}
              while(++i!=size);}
          return;
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
	  sf=(s=*out)+size; do *s=*ib++; while(++s<sf);
          if(numnul)
            {s=*out; do {s[i]=(iib[i]==FITSBYTENUL)?FITSSHORTNUL:s[i];}
              while(++i!=size);}
          return;
	case SHORT_IMG:
	  sf=(s=*out)+size; do *s=*is++; while(++s<sf); return;
	case LONG_IMG:
	  sf=(s=*out)+size; do *s=*il++; while(++s<sf);
          if(numnul)
            {s=*out; do {s[i]=(iil[i]==FITSLONGNUL)?FITSSHORTNUL:s[i];}
              while(++i!=size);}
          return;
	case LONGLONG_IMG:
	  sf=(s=*out)+size; do *s=*iL++; while(++s<sf);
          if(numnul)
            {s=*out; do {s[i]=(iiL[i]==FITSLLONGNUL)?FITSSHORTNUL:s[i];}
              while(++i!=size);}
          return;
	case FLOAT_IMG:
	  sf=(s=*out)+size; do *s=roundf(*iif++); while(++s<sf);
          if(numnul)
            {s=*out; do {s[i]=isnan(iiif[i])?FITSSHORTNUL:s[i];}
              while(++i!=size);}
          return;
	case DOUBLE_IMG:
	  sf=(s=*out)+size; do *s=round(*id++); while(++s<sf);
          if(numnul)
            {s=*out; do {s[i]=isnan(iid[i])?FITSSHORTNUL:s[i];}
              while(++i!=size);}
          return;
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
	  lf=(l=*out)+size; do *l=*ib++; while(++l<lf);
          if(numnul)
            {l=*out; do {l[i]=(iib[i]==FITSBYTENUL)?FITSLONGNUL:l[i];}
              while(++i!=size);}
          return;
	case SHORT_IMG:
	  lf=(l=*out)+size; do *l=*is++; while(++l<lf);
          if(numnul)
            {l=*out; do {l[i]=(iis[i]==FITSSHORTNUL)?FITSLONGNUL:l[i];}
              while(++i!=size);}
          return;
	case LONG_IMG:
	  lf=(l=*out)+size; do *l=*il++; while(++l<lf); return;
	case LONGLONG_IMG:
	  lf=(l=*out)+size; do *l=*iL++; while(++l<lf);
          if(numnul)
            {l=*out; do {l[i]=(iiL[i]==FITSLLONGNUL)?FITSLONGNUL:l[i];}
              while(++i!=size);}
          return;
	case FLOAT_IMG:
	  lf=(l=*out)+size; do *l=roundf(*iif++); while(++l<lf);
          if(numnul)
            {l=*out; do {l[i]=isnan(iiif[i])?FITSLONGNUL:l[i];}
              while(++i!=size);}
          return;
	case DOUBLE_IMG:
	  lf=(l=*out)+size; do *l=round(*id++); while(++l<lf);
          if(numnul)
            {l=*out; do {l[i]=isnan(iid[i])?FITSLONGNUL:l[i];}
              while(++i!=size);}
          return;
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
	  Lf=(L=*out)+size; do *L=*ib++; while(++L<Lf);
          if(numnul)
            {L=*out; do {L[i]=(iib[i]==FITSBYTENUL)?FITSLLONGNUL:L[i];}
              while(++i!=size);}
          return;
	case SHORT_IMG:
	  Lf=(L=*out)+size; do *L=*is++; while(++L<Lf);
          if(numnul)
            {L=*out; do {L[i]=(iis[i]==FITSSHORTNUL)?FITSLLONGNUL:L[i];}
              while(++i!=size);}
          return;
	case LONG_IMG:
	  Lf=(L=*out)+size; do *L=*il++; while(++L<Lf);
          if(numnul)
            {L=*out; do {L[i]=(iil[i]==FITSLONGNUL)?FITSLLONGNUL:L[i];}
              while(++i!=size);}
          return;
	case LONGLONG_IMG:
	  Lf=(L=*out)+size; do *L=*iL++; while(++L<Lf); return;
	case FLOAT_IMG:
	  Lf=(L=*out)+size; do *L=roundf(*iif++); while(++L<Lf);
          if(numnul)
            {L=*out; do {L[i]=isnan(iiif[i])?FITSLLONGNUL:L[i];}
              while(++i!=size);}
          return;
	case DOUBLE_IMG:
	  Lf=(L=*out)+size; do *L=round(*id++); while(++L<Lf);
          if(numnul)
            {L=*out; do {L[i]=isnan(iid[i])?FITSLLONGNUL:L[i];}
              while(++i!=size);}
          return;
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
	  ff=(f=*out)+size; do *f=*ib++; while(++f<ff);
          if(numnul)
            {f=*out; do {f[i]=iib[i]==FITSBYTENUL?FITSFLOATNUL:f[i];}
              while(++i!=size);}
          return;
	case SHORT_IMG:
	  ff=(f=*out)+size; do *f=*is++; while(++f<ff);
          if(numnul)
            {f=*out; do {f[i]=iis[i]==FITSSHORTNUL?FITSFLOATNUL:f[i];}
              while(++i!=size);}
          return;
	case LONG_IMG:
	  ff=(f=*out)+size; do *f=*il++; while(++f<ff);
          if(numnul)
            {f=*out; do {f[i]=iil[i]==FITSLONGNUL?FITSFLOATNUL:f[i];}
              while(++i!=size);}
          return;
	case LONGLONG_IMG:
	  ff=(f=*out)+size; do *f=*iL++; while(++f<ff);
          if(numnul)
            {f=*out; do {f[i]=iiL[i]==FITSLLONGNUL?FITSFLOATNUL:f[i];}
              while(++i!=size);}
          return;
	case FLOAT_IMG:
	  ff=(f=*out)+size; do *f=*iif++; while(++f<ff); return;
	case DOUBLE_IMG:
	  ff=(f=*out)+size; do *f=*id++; while(++f<ff); return;
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
	  df=(d=*out)+size; do *d=*ib++; while(++d<df);
          if(numnul)
            {d=*out; do {d[i]=iib[i]==FITSBYTENUL?FITSFLOATNUL:d[i];}
              while(++i!=size);}
          return;
	case SHORT_IMG:
	  df=(d=*out)+size; do *d=*is++; while(++d<df);
          if(numnul)
            {d=*out; do {d[i]=iis[i]==FITSSHORTNUL?FITSFLOATNUL:d[i];}
              while(++i!=size);}
          return;
	case LONG_IMG:
	  df=(d=*out)+size; do *d=*il++; while(++d<df);
          if(numnul)
            {d=*out; do {d[i]=iil[i]==FITSLONGNUL?FITSFLOATNUL:d[i];}
              while(++i!=size);}
          return;
	case LONGLONG_IMG:
	  df=(d=*out)+size; do *d=*iL++; while(++d<df);
          if(numnul)
            {d=*out; do {d[i]=iiL[i]==FITSLLONGNUL?FITSFLOATNUL:d[i];}
              while(++i!=size);}
          return;
	case FLOAT_IMG:
	  df=(d=*out)+size; do *d=*iif++; while(++d<df); return;
	case DOUBLE_IMG:
	  df=(d=*out)+size; do *d=*id++; while(++d<df); return;
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
  /*
  fits_write_comment(fptr, SHORTCOPYRIGHT, &status);
  fits_write_comment(fptr, SHORTLICENSE, &status);
  */
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

   If the WCS structure is not recognized, then this function will
   return a NULL pointer for the wcsprm structure and a zero for
   nwcs. It will also report the fact to the user in stderr.

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
  status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, nwcs, wcs);
  free(fullheader);
  if(status)
    {
      fprintf(stderr, "\n##################\n"
              "WCSLIB Warning: wcspih ERROR %d: %s.\n"
              "##################\n",
              status, wcs_errmsg[status]);
      *wcs=NULL; *nwcs=0;
    }

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
readfitswcs(char *filename, char *hdu, int *nwcs, struct wcsprm **wcs)
{
  int status=0;
  fitsfile *fptr;

  /* Check HDU for realistic conditions: */
  readfitshdu(filename, hdu, IMAGE_HDU, &fptr);

  /* Read the WCS information: */
  readwcs(fptr, nwcs, wcs);

  /* Close the FITS file: */
  fits_close_file(fptr, &status);
  fitsioerror(status, NULL);
}






/* Read a FITS image into an array corresponding to fitstype and also
   save the size of the array.

   If the image has any null pixels, their number is returned by this
   function. The value that is placed for those pixels is defined by
   the macros in fitsarrayvv.h and depends on the type of the data.*/
size_t
fitsimgtoarray(char *filename, char *hdu, int *bitpix, void **array,
               size_t *s0, size_t *s1)
{
  void *bitnul;
  fitsfile *fptr;
  int status=0, anynul=0;
  long naxes[2], fpixel[]={1,1};

  /* Check HDU for realistic conditions: */
  readfitshdu(filename, hdu, IMAGE_HDU, &fptr);

  /* Get the bitpix and size of the image: */
  imgbitpixsize(fptr, bitpix, naxes);
  *s0=naxes[1];
  *s1=naxes[0];

  /* Allocate space for the array. */
  bitnul=bitpixnull(*bitpix);
  *array=bitpixalloc(*s0 * *s1, *bitpix);

  /* Read the image into the allocated array: */
  fits_read_pix(fptr, bitpixtodtype(*bitpix), fpixel, *s0 * *s1,
		bitnul, *array, &anynul, &status);
  if(status) fitsioerror(status, NULL);
  free(bitnul);

  /* Close the FITS file: */
  fits_close_file(fptr, &status);
  fitsioerror(status, NULL);

  /* Return the number of nul pixels: */
  return (size_t)anynul;
}




















/*************************************************************
 ******************      Array to FITS      ******************
 *************************************************************/
void
arraytofitsimg(char *filename, char *hdu, int bitpix, void *array,
	       size_t s0, size_t s1, size_t numblank, struct wcsprm *wcs,
	       char *spack_string)
{
  int nkeyrec;
  void *blank;
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

  if(numblank)
    if(bitpix==BYTE_IMG || bitpix==SHORT_IMG
       || bitpix==LONG_IMG || bitpix==LONGLONG_IMG)
      {
        blank=bitpixnull(bitpix);
        if(fits_write_key(fptr, datatype, "BLANK", blank,
                          "Pixels with no data.", &status) )
          fitsioerror(status, "Adding the BLANK keyword.");
      }

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
