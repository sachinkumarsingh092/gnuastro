/*********************************************************************
changetype -- Changing the array of one data type to another.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>

#include <gnuastro/data.h>

#include <changetype.h>










#define BASIC_DEFINITIONS                                          \
  unsigned char  *iuc=in->array, *iiuc=in->array;                  \
  char           *ic=in->array,  *iic=in->array;                   \
  unsigned short *ius=in->array, *iius=in->array;                  \
  short          *is=in->array,  *iis=in->array;                   \
  unsigned int   *iui=in->array, *iiui=in->array;                  \
  int            *ii=in->array,  *iii=in->array;                   \
  unsigned long  *iul=in->array, *iiul=in->array;                  \
  long           *il=in->array,  *iil=in->array;                   \
  LONGLONG       *iL=in->array,  *iiL=in->array;                   \
  float          *iif=in->array;                                   \
  double         *id=in->array;





#define IN_INTEGER                                                 \
    case GAL_DATA_TYPE_UCHAR:                                      \
      do *o=*iuc++; while(++o<of);                                 \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iiuc++ == GAL_DATA_BLANK_UCHAR ? blank : *ob;  \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_CHAR:                                       \
      do *o=*ic++; while(++o<of);                                  \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iic++ == GAL_DATA_BLANK_CHAR ? blank : *ob;    \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_USHORT:                                     \
      do *o=*ius++; while(++o<of);                                 \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iius++ == GAL_DATA_BLANK_USHORT ? blank : *ob; \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_SHORT:                                      \
      do *o=*is++; while(++o<of);                                  \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iis++ == GAL_DATA_BLANK_SHORT ? blank : *ob;   \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_UINT:                                       \
      do *o=*iui++; while(++o<of);                                 \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iiui++ == GAL_DATA_BLANK_UINT ? blank : *ob;   \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_INT:                                        \
      do *o=*ii++; while(++o<of);                                  \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iii++ == GAL_DATA_BLANK_INT ? blank : *ob;     \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_ULONG:                                      \
      do *o=*iul++; while(++o<of);                                 \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iiul++ == GAL_DATA_BLANK_ULONG ? blank : *ob;  \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_LONG:                                       \
      do *o=*il++; while(++o<of);                                  \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iil++ == GAL_DATA_BLANK_LONG ? blank : *ob;    \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_LONGLONG:                                   \
      do *o=*iL++; while(++o<of);                                  \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = *iiL++==GAL_DATA_BLANK_LONGLONG ? blank : *ob;  \
          while(++ob<of); }                                        \
      return;





#define OUT_INT_IN_FLOAT                                           \
    case GAL_DATA_TYPE_FLOAT:                                      \
      do *o=roundf(*iif++); while(++o<of);                         \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = isnan(*iiif++) ? GAL_DATA_BLANK_UCHAR : *ob;    \
          while(++ob<of); }                                        \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_DOUBLE:                                     \
      do *o=round(*id++); while(++o<of);                           \
      if(in->anyblank && in->type!=out->type)                      \
        { do *ob = isnan(*iid++) ? GAL_DATA_BLANK_UCHAR : *ob;     \
          while(++ob<of); }                                        \
      return;





#define OUT_FLOAT_IN_FLOAT                                         \
    case GAL_DATA_TYPE_FLOAT:                                      \
      do *o=*iif++; while(++o<of);                                 \
      return;                                                      \
                                                                   \
    case GAL_DATA_TYPE_DOUBLE:                                     \
      do *o=*id++; while(++o<of);                                  \
      return;





#define DEFAULT_PRINT_ERROR                                        \
    case GAL_DATA_TYPE_STRING:                                     \
      error(EXIT_FAILURE, 0, "type conversion can't be done on "   \
            "string arrays.");                                     \
                                                                   \
    default:                                                       \
      error(EXIT_FAILURE, 0, "type %d was not recognized in "      \
            "`change_type_out_is_%s' (data.c)", in->type, tforerr);





#define CHANGE_OUT_IS_INTEGER                                      \
  BASIC_DEFINITIONS;                                               \
  float *iiif=in->array;                                           \
  double *iid=in->array;                                           \
                                                                   \
  switch(in->type)                                                 \
    {                                                              \
      IN_INTEGER;                                                  \
                                                                   \
      OUT_INT_IN_FLOAT;                                            \
                                                                   \
      DEFAULT_PRINT_ERROR;                                         \
    }





#define CHANGE_OUT_IS_FLOATING                                     \
  BASIC_DEFINITIONS;                                               \
                                                                   \
  switch(in->type)                                                 \
    {                                                              \
      IN_INTEGER;                                                  \
                                                                   \
      OUT_FLOAT_IN_FLOAT;                                          \
                                                                   \
      DEFAULT_PRINT_ERROR;                                         \
    }

















void
gal_changetype_out_is_uchar(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="uchar";
  unsigned char *o = out->array, *ob = out->array;
  unsigned char blank=GAL_DATA_BLANK_UCHAR, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}





void
gal_changetype_out_is_char(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="char";
  char *o = out->array, *ob = out->array;
  char blank=GAL_DATA_BLANK_CHAR, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}




void
gal_changetype_out_is_ushort(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="ushort";
  unsigned short *o = out->array, *ob = out->array;
  unsigned short blank=GAL_DATA_BLANK_USHORT, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}





void
gal_changetype_out_is_short(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="short";
  short *o = out->array, *ob = out->array;
  short blank=GAL_DATA_BLANK_SHORT, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}




void
gal_changetype_out_is_uint(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="uint";
  unsigned int *o = out->array, *ob = out->array;
  unsigned int blank=GAL_DATA_BLANK_UINT, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}





void
gal_changetype_out_is_int(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="int";
  int *o = out->array, *ob = out->array;
  int blank=GAL_DATA_BLANK_INT, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}





void
gal_changetype_out_is_ulong(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="ulong";
  unsigned long *o = out->array, *ob = out->array;
  unsigned long blank=GAL_DATA_BLANK_LONG, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}





void
gal_changetype_out_is_long(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="long";
  long *o = out->array, *ob = out->array;
  long blank=GAL_DATA_BLANK_LONG, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}





void
gal_changetype_out_is_longlong(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="longlong";
  LONGLONG *o = out->array, *ob = out->array;
  LONGLONG blank=GAL_DATA_BLANK_LONGLONG, *of = o + out->size;

  CHANGE_OUT_IS_INTEGER;
}





void
gal_changetype_out_is_float(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="float";
  float *o = out->array, *ob = out->array;
  float blank=GAL_DATA_BLANK_FLOAT, *of = o + out->size;

  CHANGE_OUT_IS_FLOATING;
}





void
gal_changetype_out_is_double(gal_data_t *in, gal_data_t *out)
{
  char *tforerr="double";
  double *o = out->array, *ob = out->array;
  double blank=GAL_DATA_BLANK_DOUBLE, *of = o + out->size;

  CHANGE_OUT_IS_FLOATING;
}
