/*********************************************************************
txt -- functions to deal with plain text files.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/txt.h>
#include <gnuastro/list.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>

#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>









/************************************************************************/
/***************           Get table information          ***************/
/************************************************************************/

/* Format of table. Currently these constants are internal to this library,
   we don't want to crowd the name space of the user by having them in the
   header file. */
enum txt_formats_code
{
    TXT_FORMAT_INVALID,

    TXT_FORMAT_TABLE,
    TXT_FORMAT_IMAGE,
};





/* Return one of the 'txt_line_stat' constant values. */
int
gal_txt_line_stat(char *line)
{
  while(*line!='\n')
    {
      switch(*line)
        {
          /* Characters to ignore. */
        case ' ': case ',': case '\t':
          break;
        case '#':
          return GAL_TXT_LINESTAT_COMMENT;
        default:
          return GAL_TXT_LINESTAT_DATAROW;
        }
      ++line;
    }
  return GAL_TXT_LINESTAT_BLANK;
}





/* Remove the spaces around the values, and if the final/trimmed string has
   no length, return NULL. */
char *
gal_txt_trim_space(char *str)
{
  char *end;

  /* If str doesn't point to anything, just return the NULL pointer. */
  if(str==NULL) return NULL;

  /* Remove the spaces before the start of the string. */
  while(isspace(*str)) ++str;

  /* If there was nothing in the string, return NULL. */
  if(*str=='\0') return NULL;

  /* Remove the spaces at the end, and write a possibly new '\0'. */
  end = str + strlen(str) - 1;
  while(end>str && isspace(*end)) --end;
  *(end+1)='\0';

  /* Return the string. */
  return *str=='\0' ? NULL : str;
}





/* Each information comment should have a format like this (replace
   'Column' with 'Image' for 2D arrays):

      # Column N: NAME [UNITS, TYPE, BLANK] COMMENT

  TYPE has pre-defined values, and N must be an integer, but the rest can
  contain any characters (including whitespace characters). The UNITS,
  TYPE, BLANK tokens are optional, if not given, default values will be
  set. But if there are comments, then the brackets themselves are required
  to separate the name from the comments.

  Any white space characters before or after the delimiters (':', '[', ']',
  ',') is ignored, but spaces within the values are kept. For example, in
  the two following lines, NAME will be set to 'col name' (even though
  there are extra spaces in the second line, The column unit will be
  set to 'col unit'.

      # Column 2: col name
      # Column 2 :  col name     [ col unit, type ] Column comments.

  When the column type is a string, the number of characters in the string
  is also necessary, for example 'str10'. Without an integer attached, the
  line will be ignored.

  In the case of an error or mis-match, the line will be ignored.

  This function will make a linked list of information about each column
  that has information in the comments. The information on each column
  doesn't have to be in order, for example the information of column 10 can
  be before column 7.
*/
static void
txt_info_from_comment(char *in_line, gal_data_t **datall, char *comm_start,
                      int inplace)
{
  gal_data_t *tmp;
  int index, strw=0;
  char *line, *aline, *tailptr;
  size_t len=strlen(comm_start);
  int type=GAL_TYPE_FLOAT64;                     /* Default type. */
  char *number=NULL, *name=NULL, *comment=NULL;
  char *inbrackets=NULL, *unit=NULL, *typestr=NULL, *blank=NULL;

  /* Make a copy of the input line if 'inplace==0'. */
  if(inplace) line=aline=in_line;
  else
    {
      /* Because the 'line' pointer will change, we need a pointer to the
         start of the originally allocated lines. This is the purpose of
         'aline' (allocated-line). */
      gal_checkset_allocate_copy(in_line, &aline);
      line=aline;
    }


  /* Only read this comment line if it follows the convention: */
  if( !strncmp(line, comm_start, len) )
    {
      /* Set 'name', 'inbrackets', and 'comment' in the first pass through
         the line. */
      number=line+len;
      while(*line!='\0')
        {
          switch(*line)
            {
            case ':':
              if(name==NULL) { *line='\0'; name=line+1; }
              break;

            case '[':
              if(name && inbrackets==NULL) { *line='\0'; inbrackets=line+1; }
              break;

            case ']':
              if(inbrackets && comment==NULL) { *line='\0'; comment=line+1; }
              break;

            case '\n':
              *line='\0';
              break;
            }
          ++line;
        }


      /* Read the column number as an integer. If it can't be read as an
         integer, or is zero or negative then just return without adding
         anything to this line. */
      index=strtol(number, &tailptr, 0);
      if(*tailptr!='\0' || index<=0) return;


      /* If there was no name (the line is just '# Column N:'), then ignore
         the line. Relying on the column count from the first line is more
         robust and less prone to human error, for example typing a number
         larger than the total number of columns.  */
      name=gal_txt_trim_space(name);
      if(name==NULL) return;


      /* If this is a repeated index, ignore it. */
      for(tmp=*datall; tmp!=NULL; tmp=tmp->next)
        if(tmp->status==index)
          return;


      /* If there were brackets, then break it up. */
      if(inbrackets)
        {
          unit=inbrackets;
          while(*inbrackets!='\0')
            {
              if(*inbrackets==',')
                {
                  *inbrackets='\0';
                  if     (typestr==NULL)  typestr = inbrackets+1;
                  else if(blank==NULL)    blank   = inbrackets+1;
                }
              ++inbrackets;
            }
        }


      /* If 'typestr' was given, then check if this is a standard type. If
         'typestr' wasn't specified, then the default double type code will
         be used (see the variable definitions above). If the given type
         isn't a standard type then ignore the line. Just note that if we
         are dealing with the string type, we have to pull out the number
         part first. If there is no number for a string type, then ignore
         the line. */
      if(typestr && *typestr!='\0')
        {
          typestr=gal_txt_trim_space(typestr);
          if( !strncmp(typestr, "str", 3) )
            {
              type=GAL_TYPE_STRING;
              strw=strtol(typestr+3, &tailptr, 0);
              if(*tailptr!='\0' || strw<0) return;
            }
          else
            {
              type=gal_type_from_name(typestr);
              if(type==GAL_TYPE_INVALID) return;
            }
        }


      /* Add this column's information into the columns linked list. We
         will define the data structur's array to have zero dimensions (no
         array) by default. If there is a blank value its value will be put
         into the array by 'gal_table_read_blank'. To keep the name, unit,
         and comment strings, trim the white space before and after each
         before using them here.  */
      gal_list_data_add_alloc(datall, NULL, type, 0, NULL, NULL, 0, -1, 1,
                              name, gal_txt_trim_space(unit),
                              gal_txt_trim_space(comment) );


      /* Put the number of this column into the status variable of the data
         structure. If the type is string, then also copy the width into
         the structure. */
      (*datall)->status=index;
      (*datall)->disp_width = type==GAL_TYPE_STRING ? strw : 0;


      /* Write the blank value into the array. Note that this is not the
         final column, we are just collecting information now. */
      gal_tableintern_read_blank(*datall, gal_txt_trim_space(blank));
    }

  /* Clean up. */
  if(in_line!=aline) free(aline);
}





/* In the case of a table, the input might not have had information in its
   comments, or the information might not have been complete. So we need to
   go through the first row of data also. In the case of the image, this is
   necessary, because we need to find the second dimension value.

   This function will return the number of tokens in the first row of the
   given text file. If the file is a text table with string columns, the
   contents of the string column will be counted as one token.*/
static size_t
txt_info_from_first_row(char *in_line, gal_data_t **datall, int format,
                        int inplace)
{
  gal_data_t *col, *prev, *tmp;
  size_t n=0, maxcnum=0, numtokens;
  char *line, *token, *end, *aline=NULL;

  /* Make a copy of the input line if necessary. */
  if(inplace) line=in_line;
  else
    {
      gal_checkset_allocate_copy(in_line, &line);
      aline=line; /* We are going to change 'line' during this function. */
    }
  end=line+strlen(line);

  /* Remove the line termination character(s) from the end of the line. In
     Unix, the line terminator is just the new-line character, however, in
     some operating systems (like MS Windows), it is two characters:
     carriage return and new-line. To be able to deal with both, we will be
     checking the second last character first, the ASCII code for carriage
     return is 13.

     If the last column is a string, and the given length is larger than
     the available space on the line, we don't want to have the line's
     new-line character. Its better for it to actually be shorter than the
     space. */
  if( end>line+2 && *(end-2)==13 ) *(end-2)='\0';
  else if( *(end-1)=='\n' )        *(end-1)='\0';

  /* Get the maximum number of columns read from the comment
     information. */
  for(col=*datall; col!=NULL; col=col->next)
    maxcnum = maxcnum>col->status ? maxcnum : col->status;

  /* Go over the line check/fill the column information. */
  while(++n)
    {
      /* When we are dealing with a text table, check if there is
         information for this column. For a text image, only the number of
         tokens is important (as the second dimension of the image), so
         just assume there no information. */
      if(format==TXT_FORMAT_TABLE)
        for(col=*datall; col!=NULL; col=col->next) {if(col->status==n) break;}
      else
        col=NULL;


      /* If there is information for this column, then check if it is a
         string, and if so, don't use 'strtok_r' (because it might have
         delimiters). So manually go ahead in the line till you get to the
         start of the string, then increment the line until the end of the
         space set for the strings. */
      if(col)
        {
          if( col->type==GAL_TYPE_STRING )
            {
              /* Remove all delimiters before the string starts. */
              while(isspace(*line) || *line==',') ++line;

              /* Increment line to the end of the string. */
              line = (token=line) + col->disp_width;

              /* If we haven't reached the end of the line, then set a NULL
                 character where the string ends, so we can use the
                 token. VERY IMPORTANT: this should not be '<=end'. If the
                 given width is larger than line, there is no problem, the
                 '\0' of the line will also be used to end this last
                 column.*/
              if(line<end)
                {
                  *line++='\0';
                  /* printf(" col %zu: -%s-\n", i, token); */
                }
              else break;
            }
          else
            {
              token=strtok_r(n==1?line:NULL, GAL_TXT_DELIMITERS, &line);
              if(token==NULL) break;
              /* printf(" col %zu: =%s=\n", i, token); */
            }
        }
      else
        {
          /* Make sure a token exists in this undefined column. */
          token=strtok_r(n==1?line:NULL, GAL_TXT_DELIMITERS, &line);
          if(token==NULL) break;

          /* A token exists. For a table, define a new element in the
             linked list and set the column to the default double type with
             no information, then set its status value to the column
             number. So, for a table, this should be done on every
             column. But for an image, this should only be done once (when
             'datall' has not been defined yet, for example in the column
             information). */
          if( *datall==NULL || format==TXT_FORMAT_TABLE )
            {
              gal_list_data_add_alloc(datall, NULL, GAL_TYPE_FLOAT64, 0,
                                      NULL, NULL, 0, -1, 1, NULL, NULL, NULL);
              (*datall)->status=n;
            }
        }
    }


  /* When looking at a text table, 'n' is the number of columns (elements
     in the linked list). But when looking at an image, it is the size of
     the second dimension. To unify things from this step forwards, we will
     thus keep the value of 'n' until this point in another variable (that
     will be returned finally), and for an image, change 'n' to 1. This is
     necsesary in case the user has for example given two column
     information comments on an image plain text file.

     Note that 'n' counts from 1, so the total number of tokens is one less
     than 'n'.*/
  numtokens=n-1;
  if(format==TXT_FORMAT_IMAGE) n=1;

  /* If the number of columns/images given by the comments is larger than
     the actual number of lines, remove those that have larger numbers from
     the linked list before things get complicated outside of this
     function. */
  if(maxcnum>n)
    {
      prev=NULL;
      col=*datall;
      while(col!=NULL)
        {
          if(col->status > n) /* Column has no data (was only in comments) */
            {
              /* This column has to be removed/freed. But we have to make
                 some corrections before freeing it:

                  - When 'prev==NULL', then we still haven't got to the
                    first valid element yet and must free this one, but if
                    we do that, then the main pointer to the start of the
                    list will be lost (we will loose all connections with
                    the chain after leaving this loop). So we need to set
                    that to the next element.

                  - When there actually was a previous element
                    ('prev!=NULL'), then we must correct it's next
                    pointer. Otherwise we will break up the chain.*/
              if(prev) prev->next=col->next; else *datall=col->next;
              tmp=col->next;
              gal_data_free(col);
              col=tmp;
            }
          else                /* Column has data.                          */
            {
              prev=col;
              col=col->next;
            }
        }
    }

  /* Return the total number of columns/second-img-dimension. */
  if(inplace==0) free(aline);
  return numtokens;
}





/* In the steps above, we read/set the information for each column. But to
   enforce minimum standard requirements on the user, things were allowed
   to be read very loosely, for example some columns can be not defined
   (and will thus be read as a double type), or they don't necessarily have
   to be given in the same order as the table (in which case, the first
   non-commented line provides basic information like how many columns
   there are). So we just pushed each new read/set column into a linked
   list.

   With this function, we convert that badly orderd linked list into a
   clean and ordered array for much more easier random access during the
   selection/reading of the data in the columns.

   After this function, the list is freed. */
static gal_data_t *
txt_infoll_to_array(gal_data_t *datall, size_t *numdata)
{
  size_t numc=0, ind;
  gal_data_t *data, *dataarr;

  /* First find the total number of columns. */
  for(data=datall; data!=NULL; data=data->next) ++numc;

  /* Conversion to an arry is only necessary when there is more than one
     element in the list. */
  if(numc>1)
    {
      /* Now, allocate the array and put in the values. */
      dataarr=gal_data_array_calloc(numc);

      /* Put each dataset/column into its proper place in the array.  */
      while(datall!=NULL)
        {
          /* Pop the top element. */
          data=gal_list_data_pop(&datall);

          /* The 'status' value is the number of the column (counting from
             1, not 0). */
          ind=data->status-1;

          /* Put all the information from 'data' into the respective part
             of the array. About the pointers, instead of having to
             allocate them again, we will just set them to NULL so
             'gal_data_free' doesn't remove them.*/
          dataarr[ind].name       = data->name;    data->name=NULL;
          dataarr[ind].unit       = data->unit;    data->unit=NULL;
          dataarr[ind].array      = data->array;   data->array=NULL;
          dataarr[ind].dsize      = data->dsize;   data->dsize=NULL;
          dataarr[ind].comment    = data->comment; data->comment=NULL;

          dataarr[ind].type       = data->type;
          dataarr[ind].ndim       = data->ndim;
          dataarr[ind].size       = data->size;
          dataarr[ind].disp_width = data->disp_width;

          /* Clean up. */
          gal_data_free(data);
        }
    }
  else
    dataarr=datall;

  /* Return the array of all column information and put the number of
     columns into the given pointer. */
  *numdata=numc;
  return dataarr;
}





static void
txt_get_info_line(char *line, gal_data_t **datall, char *comm_start,
                  int *firstlinedone, int format, size_t *dsize, int inplace)
{
  size_t numtokens;

  switch( gal_txt_line_stat(line) )
    {
      /* Line is a comment, see if it has formatted information. */
    case GAL_TXT_LINESTAT_COMMENT:
      txt_info_from_comment(line, datall, comm_start, inplace);
      break;

      /* Line is actual data, use it to fill in the gaps.  */
    case GAL_TXT_LINESTAT_DATAROW:
      ++dsize[0];
      if(*firstlinedone==0)
        {
          *firstlinedone=1;
          numtokens=txt_info_from_first_row(line, datall, format, inplace);
          if(format==TXT_FORMAT_IMAGE) dsize[1]=numtokens;
        }
      break;

      /* We also have the case of GAL_TXT_LINESTAT_BLANK, but we don't
         need to do anything about it. */
    }
}





/* Return the information about a text file table. If there were no
   readable rows, it will return NULL.*/
static gal_data_t *
txt_get_info(char *filename, gal_list_str_t *lines, int format,
             size_t *numdata, size_t *dsize)
{
  FILE *fp;
  gal_list_str_t *tmp;
  gal_data_t *datall=NULL;
  int test, firstlinedone=0;
  char *line, *format_err="empty", *comm_start;
  size_t linelen=10; /* 'linelen' will be increased by 'getline'. */

  /* 'filename' and 'lines' cannot both be non-NULL. */
  test = (filename!=NULL) + (lines!=NULL);
  if( test!=1 )
    error(EXIT_FAILURE, 0, "%s: one of the 'filename' and 'lines' "
          "arguments must be NULL, but they are both %s", __func__,
          test==2 ? "non-NULL" : "NULL");

  /* Set the constant strings */
  switch(format)
    {
    case TXT_FORMAT_TABLE: format_err="table"; comm_start="# Column "; break;
    case TXT_FORMAT_IMAGE: format_err="image"; comm_start="# Image ";  break;
    default:
      error(EXIT_FAILURE, 0, "%s: code %d not recognized",
            __func__, format);
    }

  /* Initialize the first 'dsize' element. */
  dsize[0]=0;

  /* Parse the file or go over the lines. */
  if(filename)
    {
      /* Open the file. */
      errno=0;
      fp=fopen(filename, "r");
      if(fp==NULL)
        error(EXIT_FAILURE, errno, "%s: couldn't open to read as a plain "
              "text %s (from Gnuastro's '%s')", filename, format_err,
              __func__);


      /* Allocate the space necessary to keep each line as we parse
         it. Note that 'getline' is going to later 'realloc' this space to
         fit the line length. */
      errno=0;
      line=malloc(linelen*sizeof *line);
      if(line==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for line",
              __func__, linelen*sizeof *line);


      /* Read the comments of the line for possible information about the
         lines, but also confirm/complete the info by parsing the first
         uncommented line. */
      while( getline(&line, &linelen, fp) != -1 )
        txt_get_info_line(line, &datall, comm_start, &firstlinedone, format,
                          dsize, 1);


      /* Clean up and close the file. */
      free(line);
      errno=0;
      if(fclose(fp))
        error(EXIT_FAILURE, errno, "%s: couldn't close file after reading "
              "plain text %s information in %s", filename, format_err,
              __func__);
    }
  else
    {
      for(tmp=lines; tmp!=NULL; tmp=tmp->next)
        txt_get_info_line(tmp->v, &datall, comm_start, &firstlinedone,
                          format, dsize, 0);
    }

  /* The final dataset linked list can have any order (depending on how the
     user gave column information in tables for example). So here, we will
     convert the list into a nicely sorted array, note that this function
     frees list as part of the process. */
  return txt_infoll_to_array(datall, numdata);
}





/* Get the information of each column in a text file */
gal_data_t *
gal_txt_table_info(char *filename, gal_list_str_t *lines, size_t *numcols,
                   size_t *numrows)
{
  return txt_get_info(filename, lines, TXT_FORMAT_TABLE, numcols, numrows);
}





/* Get the information of a 2D array in a text file. */
gal_data_t *
gal_txt_image_info(char *filename, gal_list_str_t *lines, size_t *numimg,
                   size_t *dsize)
{
  return txt_get_info(filename, lines, TXT_FORMAT_IMAGE, numimg, dsize);
}





















/************************************************************************/
/***************             Read a txt table             ***************/
/************************************************************************/
static void
txt_read_token(gal_data_t *data, gal_data_t *info, char *token,
               size_t i, char *filename, size_t lineno, size_t colnum)
{
  char   *tailptr;
  char     **str = data->array, **strb;
  uint8_t    *uc = data->array,   *ucb;
  int8_t      *c = data->array,    *cb;
  uint16_t   *us = data->array,   *usb;
  int16_t     *s = data->array,    *sb;
  uint32_t   *ui = data->array,   *uib;
  int32_t    *ii = data->array,    *ib;
  uint64_t   *ul = data->array,   *ulb;
  int64_t     *l = data->array,    *lb;
  float       *f = data->array,    *fb;
  double      *d = data->array,    *db;

  /* Read the proper token into the column. */
  switch(data->type)
    {
    case GAL_TYPE_STRING:
      gal_checkset_allocate_copy(gal_txt_trim_space(token), &str[i]);
      if( (strb=info->array) && !strcmp( *strb, str[i] ) )
        {
          free(str[i]);
          gal_checkset_allocate_copy(GAL_BLANK_STRING, &str[i]);
        }
      break;

    case GAL_TYPE_UINT8:
      uc[i]=strtol(token, &tailptr, 0);
      if( (ucb=info->array) && *ucb==uc[i] )
        uc[i]=GAL_BLANK_UINT8;
      break;

    case GAL_TYPE_INT8:
      c[i]=strtol(token, &tailptr, 0);
      if( (cb=info->array) && *cb==c[i] )
        c[i]=GAL_BLANK_INT8;
      break;

    case GAL_TYPE_UINT16:
      us[i]=strtol(token, &tailptr, 0);
      if( (usb=info->array) && *usb==us[i] )
        us[i]=GAL_BLANK_UINT16;
      break;

    case GAL_TYPE_INT16:
      s[i]=strtol(token, &tailptr, 0);
      if( (sb=info->array) && *sb==s[i] )
        s[i]=GAL_BLANK_INT16;
      break;

    case GAL_TYPE_UINT32:
      ui[i]=strtol(token, &tailptr, 0);
      if( (uib=info->array) && *uib==ui[i] )
        ui[i]=GAL_BLANK_UINT32;
      break;

    case GAL_TYPE_INT32:
      ii[i]=strtol(token, &tailptr, 0);
      if( (ib=info->array) && *ib==ii[i] )
        ii[i]=GAL_BLANK_INT32;
      break;

    case GAL_TYPE_UINT64:
      ul[i]=strtoul(token, &tailptr, 0);
      if( (ulb=info->array) && *ulb==ul[i] )
        ul[i]=GAL_BLANK_UINT64;
      break;

    case GAL_TYPE_INT64:
      l[i]=strtol(token, &tailptr, 0);
      if( (lb=info->array) && *lb==l[i] )
        l[i]=GAL_BLANK_INT64;
      break;

      /* For the blank value of floating point types, we need to make
         sure it isn't a NaN, because a NaN value will fail on any
         condition check (even '=='). If it isn't NaN, then we can
         compare the values. */
    case GAL_TYPE_FLOAT32:
      f[i]=strtod(token, &tailptr);
      if( (fb=info->array)
          && ( (isnan(*fb) && isnan(f[i])) || *fb==f[i] ) )
        f[i]=GAL_BLANK_FLOAT64;
      break;

    case GAL_TYPE_FLOAT64:
      d[i]=strtod(token, &tailptr);
      if( (db=info->array)
          && ( (isnan(*db) && isnan(d[i])) || *db==d[i] ) )
        d[i]=GAL_BLANK_FLOAT64;
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, data->type);
    }

  /* If a number couldn't be read properly, then report an error. */
  if(data->type!=GAL_TYPE_STRING && *tailptr!='\0')
    error_at_line(EXIT_FAILURE, 0, filename, lineno, "column %zu "
                  "('%s') couldn't be read as a '%s' number",
                  colnum, token, gal_type_name(data->type, 1) );
}





static void
txt_fill(char *in_line, char **tokens, size_t maxcolnum, gal_data_t *info,
         gal_data_t *out, size_t rowind, char *filename, size_t lineno,
         int inplace, int format)
{
  size_t i, n=0;
  gal_data_t *data;
  int notenoughcols=0;
  char *end, *line, *aline=NULL;

  /* Make a copy of the input line if necessary. */
  if(inplace) line=in_line;
  else
    {
      gal_checkset_allocate_copy(in_line, &line);
      aline=line; /* We are going to change 'line' during this function. */
    }
  end=line+strlen(line);

  /* See explanations in 'txt_info_from_first_row'. */
  if( end>line+2 && *(end-2)==13 ) *(end-2)='\0';
  else if( *(end-1)=='\n' )        *(end-1)='\0';

  /* Start parsing the line. Note that 'n' and 'maxcolnum' start from
     one. */
  while(++n)
    {
      /* Break out of the parsing if we don't need the columns any
         more. The table might contain many more columns, but when they
         aren't needed, there is no point in tokenizing them. */
      if(n>maxcolnum) break;

      /* Set the pointer to the start of this token/column. See
         explanations in 'txt_info_from_first_row'. Note that an image has
         a single 'info' element for the whole array, while a table has one
         for each column. */
      if( info[format==TXT_FORMAT_TABLE ? n-1 : 0].type == GAL_TYPE_STRING )
        {
          /* Remove any delimiters and stop at the first non-delimiter. If
             we have reached the end of the line then its an error, because
             we were expecting a column here. */
          while(isspace(*line) || *line==',') ++line;
          if(*line=='\0') {notenoughcols=1; break;}

          /* Everything is good, set the pointer and increment the line to
             the end of the allocated space for this string. */
          line = (tokens[n]=line) + info[n-1].disp_width;
          if(line<end) *line++='\0';
        }
      else
        {
          /* If we have reached the end of the line, then 'strtok_r' will
             return a NULL pointer. */
          tokens[n]=strtok_r(n==1?line:NULL, GAL_TXT_DELIMITERS, &line);
          if(tokens[n]==NULL) {notenoughcols=1; break;}
        }
    }

  /* Report an error if there weren't enough columns. */
  if(notenoughcols)
    error_at_line(EXIT_FAILURE, 0, filename, lineno, "not enough columns in "
                  "this line. Previous (uncommented) lines in this file had "
                  "%zu columns, but this line has %zu columns", maxcolnum,
                  n-1); /* This must be 'n-1' (since n starts from 1). */

  /* For a sanity check:
  printf("row: %zu: ", rowind+1);
  for(n=1;n<=maxcolnum;++n) printf("-%s-, ", tokens[n]);
  printf("\n");
  */

  /* Read the desired tokens into the columns that need them. Note that
     when a blank value is defined for the column, the column's array
     pointer ('info[col->status-1]') is not NULL and points to the blank
     value. For strings, this will actually be a string. */
  switch(out->ndim)
    {
    case 1:
      for(data=out; data!=NULL; data=data->next)
        txt_read_token(data, &info[data->status-1], tokens[data->status],
                       rowind, filename, lineno, data->status);
      break;

    case 2:
      for(i=0;i<out->dsize[1];++i)
        txt_read_token(out, info, tokens[i+1], rowind * out->dsize[1] + i,
                       filename, lineno, i+1);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: currently only 1 and 2 dimensional "
            "datasets acceptable", __func__);
    }

  /* Clean up. */
  if(inplace==0) free(aline);
}





static gal_data_t *
txt_read(char *filename, gal_list_str_t *lines, size_t *dsize,
         gal_data_t *info, gal_list_sizet_t *indexll, size_t minmapsize,
         int quietmmap, int format)
{
  FILE *fp;
  int test;
  char *line;
  char **tokens;
  gal_list_str_t *tmp;
  gal_data_t *out=NULL;
  gal_list_sizet_t *ind;
  size_t one=1, maxcolnum=0, rowind=0, lineno=0, ndim;
  size_t linelen=10;        /* 'linelen' will be increased by 'getline'. */

  /* 'filename' and 'lines' cannot both be non-NULL. */
  test = (filename!=NULL) + (lines!=NULL);
  if( test!=1 )
    error(EXIT_FAILURE, 0, "%s: one of the 'filename' and 'lines' "
          "arguments must be NULL, but they are both %s", __func__,
          test==2 ? "non-NULL" : "NULL");

  /* Allocate the space necessary to keep a copy of each line as we parse
     it. Note that 'getline' is going to later 'realloc' this space to fit
     the line length. */
  errno=0;
  line=malloc(linelen*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'line'",
          __func__, linelen*sizeof *line);

  /* Allocate all the desired columns for output. We will be reading the
     text file line by line, and writing in the necessary values of each
     row individually. */
  switch(format)
    {

    /* This is a table. */
    case TXT_FORMAT_TABLE:
      for(ind=indexll; ind!=NULL; ind=ind->next)
        {
          /* Allocate the necessary space. We are setting a 1-element array
             to avoid any allocation errors. Then we are freeing the
             allocated spaces and correcting the sizes.*/
          ndim=1;
          maxcolnum = maxcolnum>ind->v+1 ? maxcolnum : ind->v+1;
          gal_list_data_add_alloc(&out, NULL, info[ind->v].type, ndim,
                                  dsize[0]?dsize:&one, NULL, 0, minmapsize,
                                  quietmmap, info[ind->v].name,
                                  info[ind->v].unit, info[ind->v].comment);
          out->disp_width=info[ind->v].disp_width;
          out->status=ind->v+1;

          /* If there were no actual rows (dsize[0]==0), free the allocated
             spaces and correct the size. */
          if(dsize[0]==0)
            {
              out->size=0;
              free(out->array);
              free(out->dsize);
              out->dsize=out->array=NULL;
            }
        }
      break;


    /* This is an image. */
    case TXT_FORMAT_IMAGE:
      if(info->next)
        error(EXIT_FAILURE, 0, "%s: currently reading only one image (2d "
              "array) from a text file is possible, the 'info' input has "
              "more than one element", __func__);
      ndim=2;
      maxcolnum=dsize[1];
      out=gal_data_alloc(NULL, info->type, ndim, dsize, NULL, 0, minmapsize,
                         quietmmap, info->name, info->unit, info->comment);
      break;


    /* Not recognized. */
    default:
      error(EXIT_FAILURE, 0, "%s: format code %d not recognized",
            __func__, format);
    }

  /* Allocate the space to keep the pointers to each token in the
     line. This is done here to avoid having to allocate/free this array
     for each line in 'txt_fill_columns'. Note that the column numbers are
     counted from one (unlike indexes that are counted from zero), so we
     need 'maxcolnum+1' elements in the array of tokens.*/
  errno=0;
  tokens=malloc((maxcolnum+1)*sizeof *tokens);
  if(tokens==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'tokens'",
          __func__, (maxcolnum+1)*sizeof *tokens);

  if(filename)
    {
      /* Open the file. */
      errno=0;
      fp=fopen(filename, "r");
      if(fp==NULL)
        error(EXIT_FAILURE, errno, "%s: couldn't open to read as a text "
              "table in %s", filename, __func__);

      /* Read the data columns. */
      while( getline(&line, &linelen, fp) != -1 )
        {
          ++lineno;
          if( gal_txt_line_stat(line) == GAL_TXT_LINESTAT_DATAROW )
            txt_fill(line, tokens, maxcolnum, info, out, rowind++,
                     filename, lineno, 1, format);
        }

      /* Clean up and close the file. */
      errno=0;
      if(fclose(fp))
        error(EXIT_FAILURE, errno, "%s: couldn't close file after reading "
              "ASCII table information in %s", filename, __func__);
      free(line);
    }
  else
    for(tmp=lines; tmp!=NULL; tmp=tmp->next)
      {
        ++lineno;
        if( gal_txt_line_stat(tmp->v) == GAL_TXT_LINESTAT_DATAROW )
          txt_fill(tmp->v, tokens, maxcolnum, info, out, rowind++,
                   filename, lineno, 0, format);
      }

  /* Clean up and return. */
  free(tokens);
  return out;
}





gal_data_t *
gal_txt_table_read(char *filename, gal_list_str_t *lines, size_t numrows,
                   gal_data_t *colinfo, gal_list_sizet_t *indexll,
                   size_t minmapsize, int quietmmap)
{
  return txt_read(filename, lines, &numrows, colinfo, indexll, minmapsize,
                  quietmmap, TXT_FORMAT_TABLE);
}





gal_data_t *
gal_txt_image_read(char *filename, gal_list_str_t *lines, size_t minmapsize,
                   int quietmmap)
{
  size_t numimg, dsize[2];
  gal_data_t *img, *imginfo;
  gal_list_sizet_t *indexll=NULL;

  /* Get the image information. */
  imginfo=gal_txt_image_info(filename, lines, &numimg, dsize);

  /* Read the table. */
  img=txt_read(filename, lines, dsize, imginfo, indexll, minmapsize,
               quietmmap, TXT_FORMAT_IMAGE);

  /* Clean up and return. */
  gal_data_free(imginfo);
  return img;
}




/* See if there is anything in the standard input already. This function is
   modeled on the solution provided in:

   https://stackoverflow.com/questions/3711830/set-a-timeout-for-reading-stdin */
static int
txt_stdin_has_contents(long timeout_microsec)
{
  fd_set fds;
  struct timeval tv;

  /* Set the timeout time. */
  tv.tv_sec  = 0;
  tv.tv_usec = timeout_microsec;

  /* Initialize 'fd_set'. */
  FD_ZERO(&fds);

  /* Set standard input (STDIN_FILENO is 0) as the FD that must be read. */
  FD_SET(STDIN_FILENO, &fds);

  /* 'select' takes the last file descriptor value + 1 in the fdset to
     check, the fdset for reads, writes, and errors.  We are only passing
     in reads.  the last parameter is the timeout.  select will return if
     an FD is ready or the timeout has occurred. */
  select(STDIN_FILENO+1, &fds, NULL, NULL, &tv);

  // return 0 if STDIN is not ready to be read.
  return FD_ISSET(STDIN_FILENO, &fds);
}




/* Read each line of the standard input into a linked list of strings. */
gal_list_str_t *
gal_txt_stdin_read(long timeout_microsec)
{
  char *line;
  gal_list_str_t *out=NULL;
  size_t lineno=0, linelen=10;/* 'linelen' will be increased by 'getline'. */

  /* If there is nothing  */
  if( txt_stdin_has_contents(timeout_microsec) )
    {
      /* Allocate the space necessary to keep a copy of each line as we
         parse it. Note that 'getline' is going to later 'realloc' this
         space to fit the line length. */
      errno=0;
      line=malloc(linelen*sizeof *line);
      if(line==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'line'",
              __func__, linelen*sizeof *line);

      /* Read the whole standard input. We are using getline because it can
         deal with a 'NULL' in the input, while also handing allocation
         issues while reading (allocating by line, not by a fixed buffer
         size). */
      while( getline(&line, &linelen, stdin) != -1 )
        {
          /* To help in reporting (when necessary), keep a count of how
             many lines we have. */
          ++lineno;

          /* Add the line to the output list. */
          gal_list_str_add(&out, line, 1);
        }

      /* Reverse the list (to be the same order as input). */
      gal_list_str_reverse(&out);

      /* Clean up. */
      free(line);
    }

  /* Return the result. */
  return out;
}


















/************************************************************************/
/***************              Write to txt                ***************/
/************************************************************************/
/* Make an array of 3 strings for each column (in practice a two
   dimensional array with 3 columns in a row for each input column). The
   columns are:

     Column 0: Printf format string.
     Column 1: Gnuastro type string (in plain text format).
     Column 2: Blank value string.
*/
#define FMTS_COLS 3
static char **
make_fmts_for_printf(gal_data_t *datall, int leftadjust, size_t *len)
{
  char **fmts;
  gal_data_t *data;
  size_t i=0, num=0;
  char fmt[2], lng[3];


  /* Allocate space for the output. */
  for(data=datall;data!=NULL;data=data->next) ++num;
  errno=0;
  fmts=malloc(FMTS_COLS*num*sizeof *fmts);
  if(fmts==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for fmts",
          __func__, FMTS_COLS*num*sizeof *fmts);


  /* Initialize the length to 0. */
  *len=0;


  /* Go over all the columns and make their formats. */
  for(data=datall;data!=NULL;data=data->next)
    {
      /* First allocate the necessary space to keep the string. */
      errno=0;
      fmts[ i*FMTS_COLS   ] = malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts);
      fmts[ i*FMTS_COLS+1 ] = malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts);
      if(fmts[i*FMTS_COLS]==NULL || fmts[i*FMTS_COLS+1]==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for fmts[%zu] "
              "or fmts[%zu]", __func__, GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts,
              i*FMTS_COLS, i*FMTS_COLS+1);


      /* If we have a blank value, get the blank value as a string and
         adjust the width */
      fmts[ i*FMTS_COLS+2 ] = ( gal_blank_present(data, 0)
                                ? gal_blank_as_string(data->type, 0)
                                : NULL );


      /* Fill in the printing paramters. */
      gal_tableintern_col_print_info(data, GAL_TABLE_FORMAT_TXT, fmt, lng);


      /* Adjust the width if a blank string was defined. */
      if(fmts[i*FMTS_COLS+2])
        data->disp_width = ( strlen(fmts[i*FMTS_COLS+2]) > data->disp_width
                             ? strlen(fmts[i*FMTS_COLS+2])
                             : data->disp_width );


      /* Print the result into the allocated string and add its length to
         the final length of the overall format statement. The space in the
         end of 'fmts[i*2]' is to ensure that the columns don't merge, even
         if the printed string is larger than the expected width. */
      if(data->disp_precision > 0)
        *len += 1 + sprintf(fmts[i*FMTS_COLS], "%%%s%d.%d%s%s ",
                            leftadjust ? "-" : "", data->disp_width,
                            data->disp_precision, lng, fmt);
      else
        *len += 1 + sprintf(fmts[i*FMTS_COLS], "%%%s%d%s%s ",
                            leftadjust ? "-" : "", data->disp_width,
                            lng, fmt);


      /* Set the string for the Gnuastro type. For strings, we also need to
         write the maximum number of characters.*/
      if(data->type==GAL_TYPE_STRING)
        sprintf(fmts[i*FMTS_COLS+1], "%s%d", gal_type_name(data->type, 0),
                data->disp_width);
      else
        strcpy(fmts[i*FMTS_COLS+1], gal_type_name(data->type, 0));


      /* Increment the column counter. */
      ++i;
    }

  /* Return the array. */
  return fmts;
}





static void
txt_print_value(FILE *fp, void *array, int type, size_t ind, char *fmt)
{
  switch(type)
    {
      /* Numerical types. */
    case GAL_TYPE_UINT8:   fprintf(fp, fmt, ((uint8_t *) array)[ind]); break;
    case GAL_TYPE_INT8:    fprintf(fp, fmt, ((int8_t *)  array)[ind]); break;
    case GAL_TYPE_UINT16:  fprintf(fp, fmt, ((uint16_t *)array)[ind]); break;
    case GAL_TYPE_INT16:   fprintf(fp, fmt, ((int16_t *) array)[ind]); break;
    case GAL_TYPE_UINT32:  fprintf(fp, fmt, ((uint32_t *)array)[ind]); break;
    case GAL_TYPE_INT32:   fprintf(fp, fmt, ((int32_t *) array)[ind]); break;
    case GAL_TYPE_UINT64:  fprintf(fp, fmt, ((uint64_t *)array)[ind]); break;
    case GAL_TYPE_INT64:   fprintf(fp, fmt, ((int64_t *) array)[ind]); break;
    case GAL_TYPE_FLOAT32: fprintf(fp, fmt, ((float *)   array)[ind]); break;
    case GAL_TYPE_FLOAT64: fprintf(fp, fmt, ((double *)  array)[ind]); break;

      /* Special consideration for strings. */
    case GAL_TYPE_STRING:
      if( !strcmp( ((char **)array)[ind], GAL_BLANK_STRING ) )
        fprintf(fp, fmt, GAL_BLANK_STRING);
      else
        fprintf(fp, fmt, ((char **)array)[ind]);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, type);
    }
}





static void
txt_write_metadata(FILE *fp, gal_data_t *datall, char **fmts)
{
  gal_data_t *data;
  char *tmp, *nstr;
  size_t i, j, num=0;
  int nlen, nw=0, uw=0, tw=0, bw=0;

  /* Get the maximum width for each information field. */
  for(data=datall;data!=NULL;data=data->next)
    {
      ++num;
      if( data->name && strlen(data->name)>nw ) nw=strlen(data->name);
      if( data->unit && strlen(data->unit)>uw ) uw=strlen(data->unit);
    }
  for(i=0;i<num;++i)
    {
      if( (tmp=fmts[ i*FMTS_COLS+1 ]) )            /* If it isn't NULL. */
        tw = strlen(tmp) > tw ? strlen(tmp) : tw;
      if( (tmp=fmts[ i*FMTS_COLS+2 ]) )            /* If it isn't NULL. */
        bw = strlen(tmp) > bw ? strlen(tmp) : bw;
    }


  /* When there are more than 9 columns, we don't want to have cases
     like '# Column 1 :' (note the space between '1' and ':', this
     space won't exist for the 2 digit colum numbers).

     To do this, we are first allocating and printing a string long
     enough to keep the final column's 'N:'. Then, for each column, we
     print only the number into the allocated space and put the ':' in
     manually immediately after the number. Note that the initial
     'asprintf' put a '\0' in the allocated space, so we can safely
     over-write the one that 'sprintf' puts with a ':' for the columns
     that have the same number of digits as the final column. */
  i=0;
  if( asprintf(&nstr, "%zu:", num)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  nlen=strlen(nstr);
  for(data=datall; data!=NULL; data=data->next)
    {
      /* Print the number into the number string, then add the ':'
         immediately after the number. */
      sprintf(nstr, "%zu", i+1);
      for(j=1;j<nlen;++j)
        if(!isdigit(nstr[j])) nstr[j] = isdigit(nstr[j-1]) ? ':' : ' ';

      /* Now print the full information. */
      fprintf(fp, "# %s %s %-*s [%-*s,%-*s,%-*s] %s\n",
              datall->ndim==1 ? "Column" : "Image", nstr,
              nw, data->name ? data->name    : "",
              uw, data->unit ? data->unit    : "",
              tw, fmts[i*FMTS_COLS+1] ? fmts[i*FMTS_COLS+1] : "",
              bw, fmts[i*FMTS_COLS+2] ? fmts[i*FMTS_COLS+2] : "",
              data->comment ? data->comment : "");
      ++i;
    }


  /* Clean up and return. */
  free(nstr);
}





void
gal_txt_write(gal_data_t *input, gal_list_str_t *comment, char *filename,
              uint8_t colinfoinstdout)
{
  FILE *fp;
  char **fmts;
  gal_list_str_t *strt;
  size_t i, j, num=0, fmtlen;
  gal_data_t *data, *next2d=NULL;

  /* Make sure input is valid. */
  if(input==NULL) error(EXIT_FAILURE, 0, "%s: input is NULL", __func__);


  /* Currently only 1 and 2 dimension datasets are acceptable. */
  if( input->ndim!=1 && input->ndim!=2 )
    error(EXIT_FAILURE, 0, "%s: only 1 and 2 dimensional datasets are "
          "currently supported. The input dataset has %zu dimensions",
          __func__, input->ndim);


  /* For a 2D dataset, we currently don't accept a list, we can only print
     one column. So keep the next pointer separately and restore it after
     the job of this function is finished. */
  if(input->ndim==2)
    {
      next2d=input->next;
      input->next=NULL;
    }


  /* Find the number of columns, do a small sanity check, and get the
     maximum width of the name and unit string if they are present. */
  for(data=input;data!=NULL;data=data->next)
    {
      /* Count. */
      ++num;

      /* Check if the dimensionality and size is the same for all the
         elements. */
      if( input!=data && gal_dimension_is_different(input, data) )
        error(EXIT_FAILURE, 0, "%s: the input list of datasets must have the "
              "same sizes (dimensions and length along each dimension)",
              __func__);
    }


  /* Prepare the necessary formats for each column, then allocate the space
     for the full list and concatenate all the separate inputs into it. */
  fmts=make_fmts_for_printf(input, 1, &fmtlen);


  /* Set the output FILE pointer: if it isn't NULL, its an actual file,
     otherwise, its the standard output. */
  if(filename)
    {
      /* Make sure the file doesn't already exist. */
      if( gal_checkset_check_file_return(filename) )
        error(EXIT_FAILURE, 0, "%s: %s already exists. For safety, this "
              "function will not over-write an existing file. Please delete "
              "it before calling this function", __func__, filename);

      /* Open the output file. */
      errno=0;
      fp=fopen(filename, "w");
      if(fp==NULL)
        error(EXIT_FAILURE, errno, "%s: couldn't be open to write text "
              "table by %s", filename, __func__);

      /* Write the comments if there were any. */
      for(strt=comment; strt!=NULL; strt=strt->next)
        fprintf(fp, "# %s\n", strt->v);
    }
  else
    fp=stdout;

  /* Write the meta-data if necessary. */
  if(filename ? 1 : colinfoinstdout)
    txt_write_metadata(fp, input, fmts);

  /* Print the dataset */
  switch(input->ndim)
    {
    case 1:
      /* When the dataset is bring printed on standard output and its a
         single number, don't print the column structure, because it will
         add white-space characters which can be annoying when used in an
         automatic script. */
      if(fp==stdout && input->size==1 && input->next==NULL)
        fprintf(fp, "%s\n",
                gal_type_to_string(input->array, input->type, 0));

      /* Dataset has more than one row AND more than one column, so follow
         the basic text formatting (like extra white space to keep the
         columns under each other). */
      else
        for(i=0;i<input->size;++i)                        /* Row.    */
          {
            j=0;
            for(data=input;data!=NULL;data=data->next)    /* Column. */
              txt_print_value(fp, data->array, data->type, i,
                            fmts[j++ * FMTS_COLS]);
            fprintf(fp, "\n");
          }
      break;


    case 2:
      for(i=0;i<input->dsize[0];++i)
        {
          for(j=0;j<input->dsize[1];++j)
            txt_print_value(fp, input->array, input->type,
                            i*input->dsize[1]+j, fmts[0]);
          fprintf(fp, "\n");
        }
      break;


    default:
      error(EXIT_FAILURE, 0, "%s: a bug! input->ndim=%zu is not recognized",
            __func__, input->ndim);
    }



  /* Clean up. */
  for(i=0;i<num;++i)
    {
      free(fmts[i*FMTS_COLS]);
      free(fmts[i*FMTS_COLS+1]);
      free(fmts[i*FMTS_COLS+2]);
    }
  free(fmts);


  /* Close the output file. */
  if(filename)
    {
      errno=0;
      if(fclose(fp))
        error(EXIT_FAILURE, errno, "%s: couldn't close file after writing "
              "of text table in %s", filename, __func__);
    }

  /* Restore the next pointer for a 2D dataset. */
  if(input->ndim==2) input->next=next2d;
}
