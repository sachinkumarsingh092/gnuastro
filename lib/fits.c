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
 **************      Acceptable FITS names     ***************
 *************************************************************/
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



















/*************************************************************
 **************      BITPIX Dependancies       ***************
 *************************************************************/
/* Set datatype (in CFITSIO) based on BITPIX. */
int
gal_fits_bitpix_to_datatype(int bitpix)
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
    case SBYTE_IMG:
      return TSBYTE;
    case USHORT_IMG:
      return TUSHORT;
    case ULONG_IMG:
      return TULONG;
    default:
      error(EXIT_FAILURE, 0, "bitpix value of %d not recognized",
            bitpix);
    }
  return 0;
}





/* The values to the TFORM header keyword are single letter capital
   letters, but that is useless in identifying the data type of the
   column. So this function will do the conversion based on the CFITSIO
   manual.*/
int
gal_fits_tform_to_datatype(char tform)
{
  switch(tform)
    {
    case 'X':
      return TBIT;
    case 'B':
      return TBYTE;
    case 'L':
      return TLOGICAL;
    case 'A':
      return TSTRING;
    case 'I':
      return TSHORT;
    case 'J':
      return TLONG;
    case 'K':
      return TLONGLONG;
    case 'E':
      return TFLOAT;
    case 'D':
      return TDOUBLE;
    case 'C':
      return TCOMPLEX;
    case 'M':
      return TDBLCOMPLEX;
    case 'S':
      return TSBYTE;
    case 'V':
      return TUINT;
    case 'U':
      return TUSHORT;
    default:
      error(EXIT_FAILURE, 0, "'%c' is not a recognized CFITSIO value for "
            "the TFORMn header keyword(s).", tform);
    }

  error(EXIT_FAILURE, 0, "A bug! Please contact us so we can fix this. "
        "For some reason, control has reached to the end of the "
        "gal_fits_tform_to_datatype function in fits.c.");
  return -1;
}





void
gal_fits_img_bitpix_size(fitsfile *fptr, int *bitpix, long *naxes)
{
  int status=0, maxdim=10, naxis;

  if( fits_get_img_param(fptr, maxdim, bitpix, &naxis, naxes, &status) )
    gal_fits_io_error(status, NULL);

  if(naxis!=2)
    error(EXIT_FAILURE, 0, "currently only a 2 dimensional image array "
          "is supported. Your array is %d dimension(s). %s", naxis,
          naxis ? "Please contact us to add this feature." : "This might "
          "Be due to the fact that the data in images with multiple "
          "extensions are sometimes put on the second extension. If this "
          "is the case, try changing the hdu (maybe to --hdu=1)");
}





void *
gal_fits_datatype_blank(int datatype)
{
  /* Define the pointers, note that we are ordering them based on the
     CFITSIO manual to be more easily comparable. */
  unsigned char *b;
  char *c;
  char **str;
  short *s;
  long *l;
  LONGLONG *L;
  float *f;
  double *d;
  gsl_complex_float *cx;
  gsl_complex *dcx;
  int *i;
  unsigned int *ui;
  unsigned short *us;
  unsigned long *ul;

  errno=0;
  switch(datatype)
    {
    case TBIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support TBIT "
            "datatype, please get in touch with us to implement it.");

    case TBYTE:
      b=malloc(sizeof *b);
      if(b==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TBYTE",
              sizeof *b);
      *b=GAL_FITS_BYTE_BLANK;
      return b;

      /* CFITSIO says "int for keywords, char for table columns". Here we
         are only assuming table columns. So in practice this also applies
         to TSBYTE.*/
    case TLOGICAL: case TSBYTE:
      c=malloc(sizeof *c);
      if(c==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TLOGICAL, or TSBYTE",
              sizeof *c);
      *c=GAL_FITS_LOGICAL_BLANK;
      return c;

    case TSTRING:
      str=malloc(sizeof *str);
      if(str==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TSTRING",
              sizeof *s);
      *str=GAL_FITS_STRING_BLANK;
      return str;

    case TSHORT:
      s=malloc(sizeof *s);
      if(s==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TSHORT",
              sizeof *s);
      *s=GAL_FITS_SHORT_BLANK;
      return s;

    case TLONG:
      l=malloc(sizeof *l);
      if(l==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TLONG",
              sizeof *l);
      *l=GAL_FITS_LONG_BLANK;
      return l;

    case TLONGLONG:
      L=malloc(sizeof *L);
      if(L==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TLONGLONG",
              sizeof *L);
      *L=GAL_FITS_LLONG_BLANK;
      return L;

    case TFLOAT:
      f=malloc(sizeof *f);
      if(f==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TFLOAT",
              sizeof *f);
      *f=GAL_FITS_FLOAT_BLANK;
      return f;

    case TDOUBLE:
      d=malloc(sizeof *d);
      if(d==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TDOUBLE",
              sizeof *d);
      *d=GAL_FITS_DOUBLE_BLANK;
      return d;

    case TCOMPLEX:
      cx=malloc(sizeof *cx);
      if(cx==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TCOMPLEX",
              sizeof *cx);
      GSL_SET_COMPLEX(cx,GAL_FITS_FLOAT_BLANK,GAL_FITS_FLOAT_BLANK);
      return cx;

    case TDBLCOMPLEX:
      dcx=malloc(sizeof *dcx);
      if(dcx==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TDBLCOMPLEX",
              sizeof *dcx);
      GSL_SET_COMPLEX(dcx,GAL_FITS_DOUBLE_BLANK,GAL_FITS_DOUBLE_BLANK);
      return dcx;

    case TINT:
      i=malloc(sizeof *i);
      if(i==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TINT",
              sizeof *i);
      *i=GAL_FITS_INT_BLANK;
      return i;

    case TUINT:
      ui=malloc(sizeof *ui);
      if(ui==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TUINT",
              sizeof *ui);
      *ui=GAL_FITS_UINT_BLANK;
      return ui;

    case TUSHORT:
      us=malloc(sizeof *us);
      if(us==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TUSHORT",
              sizeof *us);
      *us=GAL_FITS_USHORT_BLANK;
      return us;

    case TULONG:
      ul=malloc(sizeof *ul);
      if(ul==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for blank TULONG",
              sizeof *ul);
      *ul=GAL_FITS_ULONG_BLANK;
      return ul;


    default:
      error(EXIT_FAILURE, 0, "datatype value of %d not recognized",
            datatype);
    }

  return NULL;
}





/* Allocate an array based on the value of bitpix. Note that the argument
   `size' is the number of elements, necessary in the array, the number of
   bytes each element needs will be determined internaly by this function
   using the datatype argument, so you don't have to worry about it. */
void *
gal_fits_datatype_alloc(size_t size, int datatype)
{
  void *array;

  /* Allocate space for the array to keep the image. */
  switch(datatype)
    {
    case TBIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support TBIT "
            "datatype, please get in touch with us to implement it.");

      /* The parenthesis after sizeof is not a function, it is actually a
         type cast, so we have put a space between size of and the
         parenthesis to highlight this. In C, `sizeof' is an operator, not
         a function.*/
    case TBYTE:
      size *= sizeof (unsigned char);
      break;

    case TLOGICAL: case TSBYTE:
      size *= sizeof (char);
      break;

    case TSTRING:
      size *= sizeof (char *);
      break;

    case TSHORT:
      size *= sizeof (short);
      break;

    case TLONG:
      size *= sizeof (long);
      break;

    case TLONGLONG:
      size *= sizeof (LONGLONG);
      break;

    case TFLOAT:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0,
              "`float` is not 32bits on this machine. The FITS standard "
              "Requires this size");
      size *= sizeof (float);
      break;

    case TDOUBLE:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0,
              "`double` is not 64bits on this machine. The FITS standard "
              "requires this size");
      size *= sizeof (double);
      break;

    case TCOMPLEX:
      if( sizeof (float) != 4 )
        error(EXIT_FAILURE, 0,
              "`float` is not 32bits on this machine. The FITS standard "
              "Requires this size");
      size *= sizeof (gsl_complex_float);
      break;

    case TDBLCOMPLEX:
      if( sizeof (double) != 8 )
        error(EXIT_FAILURE, 0,
              "`double` is not 64bits on this machine. The FITS standard "
              "requires this size");
      size *= sizeof (gsl_complex);
      break;

    case TINT:
      size *= sizeof (int);
      break;

    case TUINT:
      size *= sizeof (unsigned int);
      break;

    case TUSHORT:
      size *= sizeof (unsigned short);
      break;

    case TULONG:
      size *= sizeof (unsigned long);
      break;

    default:
      error(EXIT_FAILURE, 0, "datatype value of %d not recognized in "
            "gal_fits_datatype_alloc", datatype);
    }

  errno=0;
  array=malloc(size);
  if(array==NULL)
    error(EXIT_FAILURE, errno,
          "array of %lu bytes in gal_fits_datatype_alloc", size);

  return array;
}





void
gal_fits_blank_to_value(void *array, int datatype, size_t size, void *value)
{
  /* 'value' will only be read from one of these based on the
     datatype. Which the caller assigned. If there is any problem, it is
     their responsability, not this function's.*/
  unsigned char *b, *bf,        bv = *(uint8_t *) value;
  char *c, *cf,                 cv = *(char *) value;
  char **str, **strf,        *strv = *(char **) value;
  short *s, *sf,                sv = *(int16_t *) value;
  long *l, *lf,                 lv = *(int32_t *) value;
  LONGLONG *L, *Lf,             Lv = *(int64_t *) value;
  float   *f, *ff,              fv = *(float *) value;
  double  *d, *df,              dv = *(double *) value;
  gsl_complex_float *cx, *cxf, cxv = *(gsl_complex_float *) value;
  gsl_complex *dcx, *dcxf,    dcxv = *(gsl_complex *) value;
  int *in, *inf,               inv = *(int *) value;
  unsigned int *ui, *uif,      uiv = *(unsigned int *) value;
  unsigned short *us, *usf,    usv = *(unsigned short *) value;
  unsigned long *ul, *ulf,     ulv = *(unsigned long *) value;

  switch(datatype)
    {
    case TBIT:
      error(EXIT_FAILURE, 0, "Currently Gnuastro doesn't support TBIT "
            "datatype, please get in touch with us to implement it.");

    case TBYTE:
      bf=(b=array)+size;
      do if(*b==GAL_FITS_BYTE_BLANK) *b++=bv; while(b<bf);
      break;


    case TLOGICAL: case TSBYTE:
      cf=(c=array)+size;
      do if(*c==GAL_FITS_LOGICAL_BLANK) *c++=cv; while(c<cf);
      break;


    case TSTRING:
      strf=(str=array)+size;
      do if(*str==GAL_FITS_STRING_BLANK) *str++=strv; while(str<strf);
      break;


    case TSHORT:
      sf=(s=array)+size;
      do if(*s==GAL_FITS_SHORT_BLANK) *s++=sv; while(s<sf);
      break;


    case TLONG:
      lf=(l=array)+size;
      do if(*l==GAL_FITS_LONG_BLANK) *l++=lv; while(l<lf);
      break;


    case TLONGLONG:
      Lf=(L=array)+size;
      do if(*L==GAL_FITS_LLONG_BLANK) *L++=Lv; while(L<Lf);
      break;


      /* Note that a NaN value is not equal to another NaN value, so we
         can't use the easy check for cases were the blank value is
         NaN. Also note that `isnan' is actually a macro, so it works for
         both float and double types.*/
    case TFLOAT:
      ff=(f=array)+size;
      if(isnan(GAL_FITS_FLOAT_BLANK))
        do if(isnan(*f)) *f++=fv; while(f<ff);
      else
        do if(*f==GAL_FITS_FLOAT_BLANK) *f++=fv; while(f<ff);
      break;


    case TDOUBLE:
      df=(d=array)+size;
      if(isnan(GAL_FITS_DOUBLE_BLANK))
        do if(isnan(*d)) *d++=dv; while(d<df);
      else
        do if(*d==GAL_FITS_FLOAT_BLANK) *d++=dv; while(d<df);
      break;


    case TCOMPLEX:
      cxf=(cx=array)+size;
      if(isnan(GAL_FITS_FLOAT_BLANK))
          do
            if(isnan(GSL_COMPLEX_P_REAL(cx))
               && isnan(GSL_COMPLEX_P_IMAG(cx)) )
              GSL_SET_COMPLEX(cx, GSL_COMPLEX_P_REAL(&cxv),
                              GSL_COMPLEX_P_IMAG(&cxv));
          while(++cx<cxf);
      else
        do
          if( GSL_COMPLEX_P_REAL(cx) == GAL_FITS_FLOAT_BLANK
              && GSL_COMPLEX_P_IMAG(cx) == GAL_FITS_FLOAT_BLANK)
            GSL_SET_COMPLEX(cx, GSL_COMPLEX_P_REAL(&cxv),
                            GSL_COMPLEX_P_IMAG(&cxv));
        while(++cx<cxf);
      break;


    case TDBLCOMPLEX:
      dcxf=(dcx=array)+size;
      if(isnan(GAL_FITS_DOUBLE_BLANK))
          do
            if(isnan(GSL_COMPLEX_P_REAL(dcx))
               && isnan(GSL_COMPLEX_P_IMAG(dcx)) )
              GSL_SET_COMPLEX(dcx, GSL_COMPLEX_P_REAL(&dcxv),
                              GSL_COMPLEX_P_IMAG(&dcxv));
          while(++dcx<dcxf);
      else
        do
          if( GSL_COMPLEX_P_REAL(dcx) == GAL_FITS_FLOAT_BLANK
              && GSL_COMPLEX_P_IMAG(dcx) == GAL_FITS_FLOAT_BLANK)
            GSL_SET_COMPLEX(dcx, GSL_COMPLEX_P_REAL(&dcxv),
                            GSL_COMPLEX_P_IMAG(&dcxv));
        while(++dcx<dcxf);
      break;


    case TINT:
      inf=(in=array)+size;
      do if(*in==GAL_FITS_INT_BLANK) *in++=inv; while(in<inf);
      break;


    case TUINT:
      uif=(ui=array)+size;
      do if(*ui==GAL_FITS_UINT_BLANK) *ui++=uiv; while(ui<uif);
      break;


    case TUSHORT:
      usf=(us=array)+size;
      do if(*us==GAL_FITS_USHORT_BLANK) *us++=usv; while(us<usf);
      break;


    case TULONG:
      ulf=(ul=array)+size;
      do if(*ul==GAL_FITS_ULONG_BLANK) *ul++=ulv; while(ul<ulf);
      break;

    default:
      error(EXIT_FAILURE, 0, "a bug! datatype value of %d not recognized. "
            "This should not happen here (blanktovalue in fitsarrayvv.c). "
            "Please contact us at %s to see how this happened", datatype,
            PACKAGE_BUGREPORT);
    }
}





void
gal_fits_change_type(void *in, int inbitpix, size_t size, int anyblank,
                     void **out, int outbitpix)
{
  size_t i=0;
  unsigned char *b, *bf,  *ib=in,  *iib=in;
  short         *s, *sf,  *is=in,  *iis=in;
  long          *l, *lf,  *il=in,  *iil=in;
  LONGLONG      *L, *Lf,  *iL=in,  *iiL=in;
  float         *f, *ff, *iif=in, *iiif=in;
  double        *d, *df,  *id=in,  *iid=in;

  /* Allocate space for the output and start filling it. */
  *out=gal_fits_datatype_alloc(size, gal_fits_bitpix_to_datatype(outbitpix) );
  switch(outbitpix)
    {
    case BYTE_IMG:
      switch(inbitpix)
        {
        case BYTE_IMG:
          bf=(b=*out)+size; do *b=*ib++; while(++b<bf); return;
        case SHORT_IMG:
          bf=(b=*out)+size; do *b=*is++; while(++b<bf);
          if(anyblank)
            {b=*out; do {b[i]=(iis[i]==GAL_FITS_SHORT_BLANK)
                         ?GAL_FITS_BYTE_BLANK:b[i];}
              while(++i!=size);}
          return;
        case LONG_IMG:
          bf=(b=*out)+size; do *b=*il++; while(++b<bf);
          if(anyblank)
            {b=*out; do {b[i]=(iil[i]==GAL_FITS_LONG_BLANK)
                         ?GAL_FITS_BYTE_BLANK:b[i];}
              while(++i!=size);}
          return;
        case LONGLONG_IMG:
          bf=(b=*out)+size; do *b=*iL++; while(++b<bf);
          if(anyblank)
            {b=*out; do {b[i]=(iiL[i]==GAL_FITS_LLONG_BLANK)
                         ?GAL_FITS_BYTE_BLANK:b[i];}
              while(++i!=size);}
          return;
        case FLOAT_IMG:
          bf=(b=*out)+size; do *b=roundf(*iif++); while(++b<bf);
          if(anyblank)
            {b=*out; do {b[i]=isnan(iiif[i])?GAL_FITS_BYTE_BLANK:b[i];}
              while(++i!=size);}
          return;
        case DOUBLE_IMG:
          bf=(b=*out)+size; do *b=round(*id++); while(++b<bf);
          if(anyblank)
            {b=*out; do {b[i]=isnan(iid[i])?GAL_FITS_BYTE_BLANK:b[i];}
              while(++i!=size);}
          return;
        default:
          error(EXIT_FAILURE, 0, "a bug!  In gal_fits_change_type "
                "(fitsarrayvv.c). BITPIX=%d of input not recognized. "
                "Please contact us so we can fix it", inbitpix);
        }
      break;

    case SHORT_IMG:
      switch(inbitpix)
        {
        case BYTE_IMG:
          sf=(s=*out)+size; do *s=*ib++; while(++s<sf);
          if(anyblank)
            {s=*out; do {s[i]=(iib[i]==GAL_FITS_BYTE_BLANK)
                         ?GAL_FITS_SHORT_BLANK:s[i];}
              while(++i!=size);}
          return;
        case SHORT_IMG:
          sf=(s=*out)+size; do *s=*is++; while(++s<sf); return;
        case LONG_IMG:
          sf=(s=*out)+size; do *s=*il++; while(++s<sf);
          if(anyblank)
            {s=*out; do {s[i]=(iil[i]==GAL_FITS_LONG_BLANK)
                         ?GAL_FITS_SHORT_BLANK:s[i];}
              while(++i!=size);}
          return;
        case LONGLONG_IMG:
          sf=(s=*out)+size; do *s=*iL++; while(++s<sf);
          if(anyblank)
            {s=*out; do {s[i]=(iiL[i]==GAL_FITS_LLONG_BLANK)
                         ?GAL_FITS_SHORT_BLANK:s[i];}
              while(++i!=size);}
          return;
        case FLOAT_IMG:
          sf=(s=*out)+size; do *s=roundf(*iif++); while(++s<sf);
          if(anyblank)
            {s=*out; do {s[i]=isnan(iiif[i])?GAL_FITS_SHORT_BLANK:s[i];}
              while(++i!=size);}
          return;
        case DOUBLE_IMG:
          sf=(s=*out)+size; do *s=round(*id++); while(++s<sf);
          if(anyblank)
            {s=*out; do {s[i]=isnan(iid[i])?GAL_FITS_SHORT_BLANK:s[i];}
              while(++i!=size);}
          return;
        default:
          error(EXIT_FAILURE, 0, "a bug!  In gal_fits_change_type "
                "(fitsarrayvv.c).  BITPIX=%d of input not recognized. "
                "Please contact us so we can fix it", inbitpix);
        }
      break;

    case LONG_IMG:
      switch(inbitpix)
        {
        case BYTE_IMG:
          lf=(l=*out)+size; do *l=*ib++; while(++l<lf);
          if(anyblank)
            {l=*out; do {l[i]=(iib[i]==GAL_FITS_BYTE_BLANK)
                         ?GAL_FITS_LONG_BLANK:l[i];}
              while(++i!=size);}
          return;
        case SHORT_IMG:
          lf=(l=*out)+size; do *l=*is++; while(++l<lf);
          if(anyblank)
            {l=*out; do {l[i]=(iis[i]==GAL_FITS_SHORT_BLANK)
                         ?GAL_FITS_LONG_BLANK:l[i];}
              while(++i!=size);}
          return;
        case LONG_IMG:
          lf=(l=*out)+size; do *l=*il++; while(++l<lf); return;
        case LONGLONG_IMG:
          lf=(l=*out)+size; do *l=*iL++; while(++l<lf);
          if(anyblank)
            {l=*out; do {l[i]=(iiL[i]==GAL_FITS_LLONG_BLANK)
                         ?GAL_FITS_LONG_BLANK:l[i];}
              while(++i!=size);}
          return;
        case FLOAT_IMG:
          lf=(l=*out)+size; do *l=roundf(*iif++); while(++l<lf);
          if(anyblank)
            {l=*out; do {l[i]=isnan(iiif[i])?GAL_FITS_LONG_BLANK:l[i];}
              while(++i!=size);}
          return;
        case DOUBLE_IMG:
          lf=(l=*out)+size; do *l=round(*id++); while(++l<lf);
          if(anyblank)
            {l=*out; do {l[i]=isnan(iid[i])?GAL_FITS_LONG_BLANK:l[i];}
              while(++i!=size);}
          return;
        default:
          error(EXIT_FAILURE, 0, "a bug!  In gal_fits_change_type "
                "(fitsarrayvv.c).  BITPIX=%d of input not recognized. "
                "Please contact us so we can fix it", inbitpix);
        }
      break;

    case LONGLONG_IMG:
      switch(inbitpix)
        {
        case BYTE_IMG:
          Lf=(L=*out)+size; do *L=*ib++; while(++L<Lf);
          if(anyblank)
            {L=*out; do {L[i]=(iib[i]==GAL_FITS_BYTE_BLANK)
                         ?GAL_FITS_LLONG_BLANK:L[i];}
              while(++i!=size);}
          return;
        case SHORT_IMG:
          Lf=(L=*out)+size; do *L=*is++; while(++L<Lf);
          if(anyblank)
            {L=*out; do {L[i]=(iis[i]==GAL_FITS_SHORT_BLANK)
                         ?GAL_FITS_LLONG_BLANK:L[i];}
              while(++i!=size);}
          return;
        case LONG_IMG:
          Lf=(L=*out)+size; do *L=*il++; while(++L<Lf);
          if(anyblank)
            {L=*out; do {L[i]=(iil[i]==GAL_FITS_LONG_BLANK)
                         ?GAL_FITS_LLONG_BLANK:L[i];}
              while(++i!=size);}
          return;
        case LONGLONG_IMG:
          Lf=(L=*out)+size; do *L=*iL++; while(++L<Lf); return;
        case FLOAT_IMG:
          Lf=(L=*out)+size; do *L=roundf(*iif++); while(++L<Lf);
          if(anyblank)
            {L=*out; do {L[i]=isnan(iiif[i])?GAL_FITS_LLONG_BLANK:L[i];}
              while(++i!=size);}
          return;
        case DOUBLE_IMG:
          Lf=(L=*out)+size; do *L=round(*id++); while(++L<Lf);
          if(anyblank)
            {L=*out; do {L[i]=isnan(iid[i])?GAL_FITS_LLONG_BLANK:L[i];}
              while(++i!=size);}
          return;
        default:
          error(EXIT_FAILURE, 0, "a bug!  In gal_fits_change_type "
                "(fitsarrayvv.c).  BITPIX=%d of input not recognized. "
                "Please contact us so we can fix it", inbitpix);
        }
      break;

    case FLOAT_IMG:
      switch(inbitpix)
        {
        case BYTE_IMG:
          ff=(f=*out)+size; do *f=*ib++; while(++f<ff);
          if(anyblank)
            {f=*out; do {f[i]=iib[i]==GAL_FITS_BYTE_BLANK
                         ?GAL_FITS_FLOAT_BLANK:f[i];}
              while(++i!=size);}
          return;
        case SHORT_IMG:
          ff=(f=*out)+size; do *f=*is++; while(++f<ff);
          if(anyblank)
            {f=*out; do {f[i]=iis[i]==GAL_FITS_SHORT_BLANK
                         ?GAL_FITS_FLOAT_BLANK:f[i];}
              while(++i!=size);}
          return;
        case LONG_IMG:
          ff=(f=*out)+size; do *f=*il++; while(++f<ff);
          if(anyblank)
            {f=*out; do {f[i]=iil[i]==GAL_FITS_LONG_BLANK
                         ?GAL_FITS_FLOAT_BLANK:f[i];}
              while(++i!=size);}
          return;
        case LONGLONG_IMG:
          ff=(f=*out)+size; do *f=*iL++; while(++f<ff);
          if(anyblank)
            {f=*out; do {f[i]=iiL[i]==GAL_FITS_LLONG_BLANK
                         ?GAL_FITS_FLOAT_BLANK:f[i];}
              while(++i!=size);}
          return;
        case FLOAT_IMG:
          ff=(f=*out)+size; do *f=*iif++; while(++f<ff); return;
        case DOUBLE_IMG:
          ff=(f=*out)+size; do *f=*id++; while(++f<ff); return;
        default:
          error(EXIT_FAILURE, 0, "a bug!  In gal_fits_change_type "
                "(fitsarrayvv.c).  BITPIX=%d of input not recognized. "
                "Please contact us so we can fix it", inbitpix);
        }
      break;

    case DOUBLE_IMG:
      switch(inbitpix)
        {
        case BYTE_IMG:
          df=(d=*out)+size; do *d=*ib++; while(++d<df);
          if(anyblank)
            {d=*out; do {d[i]=iib[i]==GAL_FITS_BYTE_BLANK
                         ?GAL_FITS_FLOAT_BLANK:d[i];}
              while(++i!=size);}
          return;
        case SHORT_IMG:
          df=(d=*out)+size; do *d=*is++; while(++d<df);
          if(anyblank)
            {d=*out; do {d[i]=iis[i]==GAL_FITS_SHORT_BLANK
                         ?GAL_FITS_FLOAT_BLANK:d[i];}
              while(++i!=size);}
          return;
        case LONG_IMG:
          df=(d=*out)+size; do *d=*il++; while(++d<df);
          if(anyblank)
            {d=*out; do {d[i]=iil[i]==GAL_FITS_LONG_BLANK
                         ?GAL_FITS_FLOAT_BLANK:d[i];}
              while(++i!=size);}
          return;
        case LONGLONG_IMG:
          df=(d=*out)+size; do *d=*iL++; while(++d<df);
          if(anyblank)
            {d=*out; do {d[i]=iiL[i]==GAL_FITS_LLONG_BLANK
                         ?GAL_FITS_FLOAT_BLANK:d[i];}
              while(++i!=size);}
          return;
        case FLOAT_IMG:
          df=(d=*out)+size; do *d=*iif++; while(++d<df); return;
        case DOUBLE_IMG:
          df=(d=*out)+size; do *d=*id++; while(++d<df); return;
        default:
          error(EXIT_FAILURE, 0, "a bug!  In gal_fits_change_type "
                "(fitsarrayvv.c).  BITPIX=%d of input not recognized. "
                "Please contact us so we can fix it", inbitpix);
        }
      break;


    default:
      error(EXIT_FAILURE, 0, "a bug! Output Bitpix value of %d is not "
            "recognized. This should not happen here "
            "(gal_fits_change_type in fitsarrayvv.c). Please "
            "contact us to see how this happened", outbitpix);
    }
}
















/*************************************************************
 **************      Number of extensions:     ***************
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




















/**************************************************************/
/**********       Check FITS image HDUs:      ************/
/**************************************************************/
static char *
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
      error(EXIT_FAILURE, 0, "HDU code %d in CFITSIO not recognized",
            hdutype);
    }
  return NULL;
}





/* Check the desired HDU in a FITS image and also if it has the
   desired type. */
void
gal_fits_read_hdu(char *filename, char *hdu, unsigned char img0_tab1,
                  fitsfile **outfptr)
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
    gal_fits_io_error(status, "reading this FITS file");
  fptr=*outfptr;

  /* Check the Type of the given HDU: */
  if (fits_get_hdu_type(fptr, &hdutype, &status) )
    gal_fits_io_error(status, NULL);


  /* Check if the type of the HDU is the expected type. We could have
     written these as && conditions, but this is easier to read, it makes
     no meaningful difference to the compiler. */
  if(img0_tab1)
    {
      if(hdutype==IMAGE_HDU)
        error(EXIT_FAILURE, 0, "%s: HDU %s is an image, not a table",
              filename, hdu);
    }
  else
    {
      if(hdutype!=IMAGE_HDU)
        error(EXIT_FAILURE, 0, "%s: HDU %s is %s, not an image",
              filename, hdu, hdutypestring(hdutype));
    }

  /* Clean up. */
  free(ffname);
}





/* Read keywords from a FITS file. The gal_fits_key pointer is an array of
   gal_fits_key structures, which keep the basic information for each
   keyword that is to be read and also stores the value in the appropriate
   type.

   ABOUT THE STRING VALUES:

   The space for a string value is statically allocated within the
   `gal_fits_key' structure (to be `FLEN_VALUE' characters, `FLEN_VALUE' is
   defined by CFITSIO). So if the value is necessary where `gal_fits_key'
   is no longer available, then you have to allocate space dynamically and
   copy the string there.
*/
void
gal_fits_read_keywords(char *filename, char *hdu, struct gal_fits_key *keys,
                       size_t num)
{
  int status=0;
  char *ffname;
  size_t i, len;
  fitsfile *fptr;
  void *valueptr=NULL;

  /* Add hdu to filename: */
  errno=0;
  len=strlen(filename)+strlen(hdu)+4;
  ffname=malloc(len*sizeof *ffname);
  if(ffname==NULL)
    error(EXIT_FAILURE, errno, "%lu characters", len);
  sprintf(ffname, "%s[%s#]", filename, hdu);

  /* Open the FITS file: */
  if( fits_open_file(&fptr, ffname, READONLY, &status) )
    gal_fits_io_error(status, "reading this FITS file");

  /* Get the desired keywords. */
  for(i=0;i<num;++i)
    {
      /* Initialize the status: */
      keys[i].status=0;

      /* Set the value-pointer based on the required type. */
      switch(keys[i].datatype)
        {
        case TSTRING:
          valueptr=keys[i].str;
          break;
        case TBYTE:
          valueptr=&keys[i].u;
          break;
        case TSHORT:
          valueptr=&keys[i].s;
          break;
        case TLONG:
          valueptr=&keys[i].l;
          break;
        case TLONGLONG:
          valueptr=&keys[i].L;
          break;
        case TFLOAT:
          valueptr=&keys[i].f;
          break;
        case TDOUBLE:
          valueptr=&keys[i].d;
          break;
        default:
          error(EXIT_FAILURE, 0, "the value of keys[%lu].datatype (=%d) "
                "is not recognized", i, keys[i].datatype);
        }

      /* Read the keyword and place its value in the poitner. */
      fits_read_key(fptr, keys[i].datatype, keys[i].keyname,
                    valueptr, NULL, &keys[i].status);


      /* In some cases, the caller might be fine with some kinds of errors,
         so we will only report an error here is the situation is not
         expected. For example, the caller might have alternatives for a
         keyword if it doesn't exist, or the non-existence of a keyword
         might itself be meaningful. So when the key doesn't exist, this
         function will not abort, it will just keep the status.

         The reason only non-existance is acceptable is this: if the
         keyword does exist, but CFITSIO cannot read it due to some
         technical difficulty, then the user probably wanted to give the
         value. But is not aware of the technical problem.
       */
      if(keys[i].status!=0 && keys[i].status!=KEY_NO_EXIST)
        gal_fits_io_error(keys[i].status, "reading the keyword");
    }

  /* Close the FITS file. */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);

  /* Clean up. */
  free(ffname);
}




















/*************************************************************
 ******************         Header          ******************
 *************************************************************/
/* Add on keyword to the list of header keywords that need to be added
   to a FITS file. In the end, the keywords will have to be freed, so
   it is important to know before hand if they were allocated or
   not. If not, they don't need to be freed. */
void
gal_fits_add_to_key_ll(struct gal_fits_key_ll **list, int datatype,
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
  newnode->datatype=datatype;
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
gal_fits_add_to_key_ll_end(struct gal_fits_key_ll **list, int datatype,
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
  newnode->datatype=datatype;
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
          gal_fits_add_to_key_ll_end(list, TSTRING, keyname, 1,
                                     value, 1, NULL, 0, NULL);
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
            error(EXIT_FAILURE, 0, "the filename `%sP has at least one "
                  "span of %lu characters without a `/`. It cannot be "
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
          if( fits_update_key(fptr, tmp->datatype, tmp->keyname,
                              tmp->value, tmp->comment, &status) )
            gal_fits_io_error(status, NULL);
        }
      else
        {
          if(fits_update_key_null(fptr, tmp->keyname, tmp->comment,
                                  &status))
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
#ifdef HAVE_WCSLIB_VERSION
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
#ifdef HAVE_WCSLIB_VERSION
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
  gal_fits_read_hdu(filename, hdu, 0, &fptr);

  /* Read the WCS information: */
  gal_fits_read_wcs_from_pointer(fptr, nwcs, wcs, hstartwcs, hendwcs);

  /* Close the FITS file: */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}






/* Read a FITS image into an array corresponding to fitstype and also
   save the size of the array.

   If the image has any null pixels, their number is returned by this
   function. The value that is placed for those pixels is defined by
   the macros in fitsarrayvv.h and depends on the type of the data.*/
int
gal_fits_hdu_to_array(char *filename, char *hdu, int *bitpix,
                      void **array, size_t *s0, size_t *s1)
{
  void *bitblank;
  fitsfile *fptr;
  long naxes[2], fpixel[]={1,1};
  int status=0, anyblank=0, datatype;

  /* Check HDU for realistic conditions: */
  gal_fits_read_hdu(filename, hdu, 0, &fptr);

  /* Get the bitpix and size of the image: */
  gal_fits_img_bitpix_size(fptr, bitpix, naxes);
  *s0=naxes[1];
  *s1=naxes[0];

  /* Allocate space for the array. */
  datatype=gal_fits_bitpix_to_datatype(*bitpix);
  bitblank=gal_fits_datatype_blank(datatype);
  *array=gal_fits_datatype_alloc(*s0 * *s1, datatype);

  /* Read the image into the allocated array: */
  fits_read_pix(fptr, gal_fits_bitpix_to_datatype(*bitpix), fpixel,
                *s0 * *s1, bitblank, *array, &anyblank, &status);
  if(status) gal_fits_io_error(status, NULL);
  free(bitblank);

  /* Close the FITS file: */
  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);

  /* Return the number of blank pixels: */
  return anyblank;
}




















/*************************************************************
 ******************      Array to FITS      ******************
 *************************************************************/
void
gal_fits_array_to_file(char *filename, char *extname, int bitpix, void *array,
                       size_t s0, size_t s1, int anyblank, struct wcsprm *wcs,
                       struct gal_fits_key_ll *headers, char *spack_string)
{
  int nkeyrec;
  void *blank;
  fitsfile *fptr;
  char *wcsheader;
  int status=0, datatype;
  long fpixel=1, naxis=2, nelements, naxes[]={s1,s0};

  datatype=gal_fits_bitpix_to_datatype(bitpix);
  nelements=naxes[0]*naxes[1];

  if(access(filename,F_OK) != -1 )
    fits_open_file(&fptr,filename, READWRITE, &status);
  else
    fits_create_file(&fptr, filename, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  fits_write_img(fptr, datatype, fpixel, nelements, array, &status);

  if(anyblank)
    if(bitpix==BYTE_IMG || bitpix==SHORT_IMG
       || bitpix==LONG_IMG || bitpix==LONGLONG_IMG)
      {
        blank=gal_fits_datatype_blank( gal_fits_bitpix_to_datatype(bitpix) );
        if(fits_write_key(fptr, datatype, "BLANK", blank,
                          "Pixels with no data.", &status) )
          gal_fits_io_error(status, "adding the BLANK keyword");
        free(blank);
      }

  fits_write_key(fptr, TSTRING, "EXTNAME", extname, "", &status);
  gal_fits_io_error(status, NULL);

  if(wcs)
    {
      /* Convert the WCS information to text. */
      status=wcshdo(WCSHDO_safe, wcs, &nkeyrec, &wcsheader);
      if(status)
        error(EXIT_FAILURE, 0, "wcshdo ERROR %d: %s", status,
              wcs_errmsg[status]);
      gal_fits_add_wcs_to_header(fptr, wcsheader, nkeyrec);
    }

  gal_fits_write_keys_version(fptr, headers, spack_string);

  fits_close_file(fptr, &status);
  gal_fits_io_error(status, NULL);
}





/* `atof' stands for "array to file". This is essentially the same as
   `gal_fits_array_to_file' except that the WCS structure's CRPIX values
   have changed. */
void
gal_fits_atof_correct_wcs(char *filename, char *hdu, int bitpix, void *array,
                          size_t s0, size_t s1, char *wcsheader, int wcsnkeyrec,
                          double *crpix, char *spack_string)
{
  fitsfile *fptr;
  int status=0, datatype;
  long fpixel=1, naxis=2, nelements, naxes[]={s1,s0};

  datatype=gal_fits_bitpix_to_datatype(bitpix);
  nelements=naxes[0]*naxes[1];

  if(access(filename,F_OK) != -1 )
    fits_open_file(&fptr,filename, READWRITE, &status);
  else
    fits_create_file(&fptr, filename, &status);

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  fits_write_img(fptr, datatype, fpixel, nelements, array, &status);

  fits_write_key(fptr, TSTRING, "EXTNAME", hdu, "", &status);
  gal_fits_io_error(status, NULL);

  fits_delete_key(fptr, "COMMENT", &status);
  fits_delete_key(fptr, "COMMENT", &status);
  gal_fits_io_error(status, NULL);

  if(wcsheader)
    {
      gal_fits_add_wcs_to_header(fptr, wcsheader, wcsnkeyrec);
      if(crpix)
        {
          fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix[0],
                          NULL, &status);
          fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix[1],
                          NULL, &status);
          gal_fits_io_error(status, NULL);
        }
    }

  gal_fits_write_keys_version(fptr, NULL, spack_string);

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
        return ASCII_TBL;
      else if(!strcmp(value, "BINTABLE"))
        return BINARY_TBL;
      else
        error(EXIT_FAILURE, 0, "The `XTENSION' keyword of this FITS file "
              "doesn't have a standard value (`%s')", value);
    }
  else
    {
      if(status==KEY_NO_EXIST)
        error(EXIT_FAILURE, 0, "the `gal_fits_table_type' function was "
              "called on a FITS extension which is not a table. As part "
              "of a utility, this is bug, so please contact us at %s so "
              "we can fix it.", PACKAGE_BUGREPORT);
      else
        gal_fits_io_error(status, NULL);
    }

  error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix it. "
        "for some reason, the control of `gal_fits_table_type' has reached "
        "the end of the function! This must not happen", PACKAGE_BUGREPORT);
  return -1;
}




















/**************************************************************/
/**********          Check prepare file            ************/
/**************************************************************/
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





/* The user has specified an input file and a mask file. In the
   processing, all masked pixels will be converted to NaN pixels in
   the input image so we only have to deal with one array. Also since
   all processing is done on floating point arrays, the input is
   converted to floating point, irrespective of its input type. The
   original input bitpix will be stored so if you want to, you can
   return it back to the input type if you please. */
void
gal_fits_file_to_float(char *inputname, char *maskname, char *inhdu,
                       char *mhdu, float **img, int *inbitpix,
                       int *anyblank, size_t *ins0, size_t *ins1)
{
  void *array;
  int maskbitpix;
  float *mask, *f, *ff, *fp;
  size_t maskanyblank, s0, s1;

  /* Read the input array and convert it to float. */
  *anyblank=gal_fits_hdu_to_array(inputname, inhdu, inbitpix,
                                  &array, ins0, ins1);
  if(*inbitpix==FLOAT_IMG)
    *img=array;
  else
    {
      gal_fits_change_type(array, *inbitpix, *ins0 * *ins1, *anyblank,
                                (void **)img, FLOAT_IMG);
      free(array);
    }

  /* If a mask was specified, read it as a float image, then set all
     the corresponding pixels of the input image to NaN. */
  if(maskname)
    {
      maskanyblank=gal_fits_hdu_to_array(maskname, mhdu, &maskbitpix,
                                         &array, &s0, &s1);

      if(maskbitpix==FLOAT_IMG || maskbitpix==DOUBLE_IMG)
        fprintf(stderr, "WARNING: the mask image (%s, hdu: %s) has a %s "
                "precision floating point data type (BITPIX=%d). The mask "
                "image is usually an integer type. Therefore this might "
                "be due to a mistake in the inputs and the results might "
                "not be what you intended. However, the program will not "
                "abort and continue working only with zero valued pixels "
                "in the given masked image.", maskname, mhdu,
                maskbitpix==FLOAT_IMG ? "single" : "double", maskbitpix);

      if(s0!=*ins0 || s1!=*ins1)
        error(EXIT_FAILURE, 0, "the input image %s (hdu: %s) has size: "
              "%lu x %lu. The mask image %s (hdu: %s) has size %lu x %lu. "
              "The two images have to have the same size", inputname,
              inhdu, *ins1, *ins0, maskname, mhdu, s1, s0);

      if(maskbitpix==FLOAT_IMG)
        mask=array;
      else
        {
          gal_fits_change_type(array, maskbitpix, *ins0 * *ins1,
                               maskanyblank, (void **)(&mask), FLOAT_IMG);
          free(array);
        }

      ff=mask;
      fp=(f=*img)+s0*s1;
      do if(*ff++!=0.0f) {*f=NAN; ++(*anyblank);} while(++f<fp);
      free(mask);
    }
}




/* Similar to filetofloat, but for double type */
void
gal_fits_file_to_double(char *inputname, char *maskname, char *inhdu,
                        char *mhdu, double **img, int *inbitpix,
                        int *anyblank, size_t *ins0, size_t *ins1)
{
  void *array;
  int maskbitpix;
  double *mask, *f, *ff, *fp;
  size_t maskanyblank, s0, s1;

  /* Read the input array and convert it to double. */
  *anyblank=gal_fits_hdu_to_array(inputname, inhdu, inbitpix,
                                  &array, ins0, ins1);
  if(*inbitpix==DOUBLE_IMG)
    *img=array;
  else
    {
      gal_fits_change_type(array, *inbitpix, *ins0 * *ins1, *anyblank,
                                (void **)img, DOUBLE_IMG);
      free(array);
    }

  /* If a mask was specified, read it as a double image, then set all
     the corresponding pixels of the input image to NaN. */
  if(maskname)
    {
      maskanyblank=gal_fits_hdu_to_array(maskname, mhdu, &maskbitpix,
                                         &array, &s0, &s1);

      if(maskbitpix==FLOAT_IMG || maskbitpix==DOUBLE_IMG)
        fprintf(stderr, "WARNING: the mask image (%s, hdu: %s) has a %s "
                "precision floating point data type (BITPIX=%d). The mask "
                "image is usually an integer type. Therefore this might "
                "be due to a mistake in the inputs and the results might "
                "not be what you intended. However, the program will not "
                "abort and continue working only with zero valued pixels in "
                "the given masked image.", maskname, mhdu,
                maskbitpix==FLOAT_IMG ? "single" : "double", maskbitpix);

      if(s0!=*ins0 || s1!=*ins1)
        error(EXIT_FAILURE, 0, "the input image %s (hdu: %s) has size: "
              "%lu x %lu. The mask image %s (hdu: %s) has size %lu x %lu. "
              "The two images have to have the same size", inputname,
              inhdu, *ins1, *ins0, maskname, mhdu, s1, s0);

      if(maskbitpix==DOUBLE_IMG)
        mask=array;
      else
        {
          gal_fits_change_type(array, maskbitpix, *ins0 * *ins1,
                                    maskanyblank, (void **)(&mask),
                                    DOUBLE_IMG);
          free(array);
        }

      ff=mask;
      fp=(f=*img)+s0*s1;
      do if(*ff++!=0.0f) {*f=NAN; ++(*anyblank);} while(++f<fp);
      free(mask);
    }
}





void
gal_fits_file_to_long(char *inputname, char *inhdu, long **img,
                      int *inbitpix, int *anyblank, size_t *ins0,
                      size_t *ins1)
{
  void *array;

  /* Read the input array and convert it to float. */
  *anyblank=gal_fits_hdu_to_array(inputname, inhdu, inbitpix,
                                  &array, ins0, ins1);
  if(*inbitpix==LONG_IMG)
    *img=array;
  else
    {
      gal_fits_change_type(array, *inbitpix, *ins0 * *ins1, *anyblank,
                                (void **)img, LONG_IMG);
      free(array);
    }
}





void
gal_fits_prep_float_kernel(char *inputname, char *inhdu, float **outkernel,
                           size_t *ins0, size_t *ins1)
{
  size_t i, size;
  double sum=0.0f;
  int bitpix, anyblank;
  float *f, *fp, *kernel, tmp;

  /* Read the kernel as a float array: */
  gal_fits_file_to_float(inputname, NULL, inhdu, NULL, outkernel, &bitpix,
                              &anyblank, ins0, ins1);
  size = *ins0 * *ins1;
  kernel=*outkernel;

  /* Simple sanity check: */
  if(*ins0%2==0 || *ins1%2==0)
    error(EXIT_FAILURE, 0, "the kernel image has to have an odd number "
          "of pixels on both sides (there has to be on pixel in the "
          "center). %s (hdu: %s) is %lu by %lu", inputname, inhdu,
          *ins1, *ins0);

  /* If there are any NaN pixels, set them to zero and normalize it.*/
  fp=(f=kernel)+size;
  do
    {
      if(isnan(*f)) *f=0.0f;
      else          sum+=*f;
    }
  while(++f<fp);
  f=kernel; do *f++ *= 1/sum; while(f<fp);

  /* Flip the kernel: */
  for(i=0;i<size/2;++i)
    {
      tmp=kernel[i];
      kernel[i]=kernel[size-i-1];
      kernel[size-i-1]=tmp;
    }
}
