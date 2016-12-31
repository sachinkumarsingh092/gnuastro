/*********************************************************************
txt -- functions to deal with plain text files.
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

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/txt.h>
#include <gnuastro/table.h>

#include <checkset.h>










/************************************************************************/
/***************           Get table information          ***************/
/************************************************************************/
/* Status of a line: */
enum txt_line_stat
{
  TXT_LINESTAT_BLANK,
  TXT_LINESTAT_ISCOMMENT,
  TXT_LINESTAT_DATAROW,
};





/* Return one of the `txt_line_stat' constant values. */
static int
get_line_stat(char *line)
{
  while(*line!='\n')
    {
      switch(*line)
        {
          /* Characters to ignore. */
        case ' ': case ',': case '\t':
          break;
        case '#':
          return TXT_LINESTAT_ISCOMMENT;
        default:
          return TXT_LINESTAT_DATAROW;
        }
      ++line;
    }
  return TXT_LINESTAT_BLANK;
}





/* Remove the spaces around the values, and if the final/trimmed string has
   no length, return NULL. */
static char *
txt_trim_space(char *str)
{
  char *end;

  /* If str doesn't point to anything, just return the NULL pointer. */
  if(str==NULL) return NULL;

  /* Remove the spaces before the start of the string. */
  while(isspace(*str)) ++str;

  /* If there was nothing in the string, then just return the ending `\0'
     character. */
  if(*str=='\0') return NULL;

  /* Remove the spaces at the end, and write a possibly new `\0'. */
  end = str + strlen(str) - 1;
  while(end>str && isspace(*end)) --end;
  *(end+1)='\0';

  /* Return the string. */
  return *str=='\0' ? NULL : str;
}





/* Use the input `blank' string and the input column to put the blank value
   in the column's array. If no blank string is given, then free the
   column's array. */
static void
txt_read_blank(gal_data_t *col, char *blank)
{
  double d;
  long long L;
  char *tailptr;

  /* If there is nothing to use as blank, then free the array. */
  if(blank==NULL) return;

  /* Allocate space to keep the blank value. */
  col->ndim=col->size=1;
  col->array=gal_data_malloc_array(col->type, col->size);

  /* Set the dsize variable. */
  errno=0;
  col->dsize=malloc(sizeof *col->dsize);
  if(col->dsize==NULL)
    error(EXIT_FAILURE, 0, "%zu bytes for `col->dsize' in `txt_read_blank' ",
          sizeof *col->dsize);

  /* String type. Copy the string.*/
  if(col->type==GAL_DATA_TYPE_STRING)
    gal_checkset_allocate_copy(blank, col->array);

  /* Floating point: Read it as a double or long, then put it in the
     array. When the conversion can't be done (the string isn't a number
     for example), then just assume no blank value was given. */
  else if(col->type==GAL_DATA_TYPE_FLOAT || col->type==GAL_DATA_TYPE_DOUBLE)
    {
      d=strtod(blank, &tailptr);
      if(*tailptr!='\0') free(col->array);
      else
        {
          if(col->type==GAL_DATA_TYPE_FLOAT) *(float *) col->array=d;
          else                              *(double *) col->array=d;
        }
    }

  /* Integers. */
  else
    {
      L=strtoll(blank, &tailptr, 0);
      if(*tailptr!='\0') free(col->array);
      else
        switch(col->type)
          {
          case GAL_DATA_TYPE_UCHAR:   *(unsigned char *) col->array=L; break;
          case GAL_DATA_TYPE_CHAR:             *(char *) col->array=L; break;
          case GAL_DATA_TYPE_USHORT: *(unsigned short *) col->array=L; break;
          case GAL_DATA_TYPE_SHORT:           *(short *) col->array=L; break;
          case GAL_DATA_TYPE_UINT:     *(unsigned int *) col->array=L; break;
          case GAL_DATA_TYPE_INT:               *(int *) col->array=L; break;
          case GAL_DATA_TYPE_ULONG:   *(unsigned long *) col->array=L; break;
          case GAL_DATA_TYPE_LONG:             *(long *) col->array=L; break;
          case GAL_DATA_TYPE_LONGLONG:     *(LONGLONG *) col->array=L; break;
          default:
            error(EXIT_FAILURE, 0, "type code %d not recognized in "
                  "`txt_str_to_blank'", col->type);
          }
    }
}






/* Each column information comment should have a format like this:

      # Column N: NAME [UNITS, TYPE, BLANK] COMMENT

  TYPE has pre-defined values, and N must be an integer, but the rest can
  contain any characters (including whitespace characters). The UNITS,
  TYPE, BLANK tokens are optional, if not given, default values will be
  set. But if there are comments, then the brackets themselves are required
  to separate the name from the comments.

  Any white space characters before or after the delimiters (`:', `[', `]',
  `,') is ignored, but spaces within the values are kept. For example, in
  the two following lines, NAME will be set to `col name' (even though
  there are extra spaces in the second line, The column unit will be
  set to `col unit'.

      # Column 2: col name
      # Column 2 :  col name     [ col unit, type ] Column comments.

  When the column type is a string, the number of characters in the string
  is also necessary, for example `str10'. Without an integer attached, the
  line will be ignored.

  In the case of an error or mis-match, the line will be ignored.

  This function will make a linked list of information about each column
  that has information in the comments. The information on each column
  doesn't have to be in order, for example the information of column 10 can
  be before column 7.
*/
static void
txt_info_from_comment(char *line, gal_data_t **colsll)
{
  char *tailptr;
  gal_data_t *tmp;
  int index, type, strw=0;
  char *number=NULL, *name=NULL, *comment=NULL;
  char *inbrackets=NULL, *unit=NULL, *typestr=NULL, *blank=NULL;

  /* Only read this comment if it follows the convention: */
  if( !strncmp(line, "# Column ", 9) )
    {
      /* Set the name, inbrackets, and comments string in the first pass
         through the line. */
      number=line+9;
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

      /* Read the column number as an integer. If it can't be read as an
         integer, or is zero or negative then just return without adding
         anything to this line. */
      index=strtol(number, &tailptr, 0);
      if(*tailptr!='\0' || index<=0) return;

      /* See if the type is a standard type, if so, then set the type,
         otherwise, return and ignore this line. Just note that if we are
         dealing with the string type, we have to pull out the number part
         first. If there is no number, there will be an error.*/
      typestr=txt_trim_space(typestr);
      if( !strncmp(typestr, "str", 3) )
        {
          type=GAL_DATA_TYPE_STRING;
          strw=strtol(typestr+3, &tailptr, 0);
          if(*tailptr!='\0' || strw<0) return;
        }
      else
        {
          type=gal_data_string_as_type(typestr);
          if(type==-1) return;
        }

      /* If this is a repeated index, ignore it. */
      for(tmp=*colsll; tmp!=NULL; tmp=tmp->next)
        if(tmp->status==index)
          return;

      /* Add this column's information into the columns linked list. We
         will define the array to have one element to keep the blank
         value. To keep the name, unit, and comment strings, trim the white
         space before and after each before using them here.  */
      gal_data_add_to_ll(colsll, NULL, type, 0, NULL, NULL, 0, -1,
                         txt_trim_space(name), txt_trim_space(unit),
                         txt_trim_space(comment) );

      /* Put the number of this column into the status variable of the data
         structure. If the type is string, then also copy the width into
         the structure. */
      (*colsll)->status=index;
      (*colsll)->disp_width = type==GAL_DATA_TYPE_STRING ? strw : 0;

      /* Write the blank value into the array. Note that this is not the
         final column, we are just collecting information now. */
      txt_read_blank(*colsll, txt_trim_space(blank));
    }
}





/* The input ASCII table might not have had information in its comments, or
   the information might not have been complete. So we need to go through
   the first row of data also. */
void
txt_info_from_row(char *line, gal_data_t **colsll)
{
  size_t n=0;
  gal_data_t *col;
  char *token, *end=line+strlen(line);

  /* Remove the new line character from the end of the line. If the last
     column is a string, and the given length is larger than the available
     space on the line, we don't want to have the line's new-line
     character. Its better for it to actually be shorter than the space. */
  *(end-1)='\0';

  /* Go over the line check/fill the column information. */
  while(++n)
    {
      /* Check if there is information for this column. */
      for(col=*colsll; col!=NULL; col=col->next) if(col->status==n) break;

      /* If there is information for this column, then check if it is a
         string, and if so, don't use `strtok_r' (because it might have
         delimiters). So manually go ahead in the line till you get to the
         start of the string, then increment the line until the end of the
         space set for the strings. */
      if(col)
        {
          if( col->type==GAL_DATA_TYPE_STRING )
            {
              /* Remove all delimiters before the string starts. */
              while(isspace(*line) || *line==',') ++line;

              /* Increment line to the end of the string. */
              line = (token=line) + col->disp_width;

              /* If we haven't reached the end of the line, then set a NULL
                 character where the string ends, so we can use the
                 token. VERY IMPORTANT: this should not be `<=end'. If the
                 given width is larger than line, there is no problem, the
                 `\0' of the line will also be used to end this last
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
          /* printf(" col %zu: *%s*\n", i, token); */

          /* A token exists, so set this column to the default double type
             with no information, then set its status value to the column
             number. */
          gal_data_add_to_ll(colsll, NULL, GAL_DATA_TYPE_DOUBLE, 0, NULL,
                             NULL, 0, -1, NULL, NULL, NULL);
          (*colsll)->status=n;
        }
    }
}





/* In the steps above, we read/set the information for each column. But to
   enforce minimum standard requirements on the user, things were allowed
   to be read very loosely, for example some columns can be not defined
   (and will thus be read as a double type), or they don't necessarily have
   to be given in the same order as the table. So we just pushed each new
   read/set column into a linked list. Now the job is done, and we want to
   convert that linked list into an array of data structures for more
   easier random access during the selection of the columns. */
static gal_data_t *
txt_infoll_to_array(gal_data_t *colsll, size_t *numcols)
{
  size_t numc=0;
  gal_data_t *col, *allcols;

  /* First find the total number of columns. */
  for(col=colsll;col!=NULL;col=col->next)
    numc = numc > col->status ? numc : col->status;

  /* Now, allocate the array and put in the values. */
  allcols=gal_data_calloc_dataarray(numc);

  /* Put each column into its proper place in the array. After the copy,
     all the (possibly) allocated spaces in the linked list are set to
     NULL, because we didn't initialize the array of column information in
     the allocation above and we don't want to re-allocate everything
     (because freeing the linked list will free them also). */
  for(col=colsll;col!=NULL;col=col->next)
    {
      /* Note that the status value counts from 1. */
      allcols[col->status-1].name       = col->name;     col->name=NULL;
      allcols[col->status-1].unit       = col->unit;     col->unit=NULL;
      allcols[col->status-1].array      = col->array;    col->array=NULL;
      allcols[col->status-1].dsize      = col->dsize;    col->dsize=NULL;
      allcols[col->status-1].comment    = col->comment;  col->comment=NULL;

      allcols[col->status-1].type       = col->type;
      allcols[col->status-1].ndim       = col->ndim;
      allcols[col->status-1].size       = col->size;
      allcols[col->status-1].disp_width = col->disp_width;
    }

  /* Return the array of all column information and put the number of
     columns into the given pointer. */
  *numcols=numc;
  return allcols;
}





/* Return the information about a text file table. */
gal_data_t *
gal_txt_table_info(char *filename, size_t *numcols, size_t *numrows)
{
  FILE *fp;
  char *line;
  int firstlinedone=0;
  gal_data_t *colsll=NULL, *allcols;
  size_t linelen=10; /* `linelen' will be increased by `getline'. */


  /* Open the file. */
  errno=0;
  fp=fopen(filename, "r");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s: could't open to read as a text table",
          filename);


  /* Allocate the space necessary to keep each line as we parse it. Note
     that `getline' is going to later `realloc' this space to fit the line
     length. */
  errno=0;
  line=malloc(linelen*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for line in `gal_txt_table_info'",
          linelen*sizeof *line);


  /* Read the comments of the line for possible information about the
     lines, but also confirm the info by trying to read the first
     uncommented line. */
  *numrows=0;
  while( getline(&line, &linelen, fp) != -1 )
    {
      /* Line is a comment, see if it has formatted information. */
      if( get_line_stat(line) == TXT_LINESTAT_ISCOMMENT )
        txt_info_from_comment(line, &colsll);

      /* Line is actual data, use it to fill in the gaps.  */
      if( get_line_stat(line) == TXT_LINESTAT_DATAROW )
        {
          ++(*numrows);
          if(firstlinedone==0)
            {
              firstlinedone=1;
              txt_info_from_row(line, &colsll);
            }
        }
    }


  /* Write the unorganized gathered information (linked list) into an
     organized array for easy processing by later steps. */
  allcols=txt_infoll_to_array(colsll, numcols);

  /* Clean up and close the file. */
  errno=0;
  if(fclose(fp))
    error(EXIT_FAILURE, errno, "%s: couldn't close file after reading ASCII "
          "table information", filename);
  gal_data_free_ll(colsll);
  free(line);


  /* Return the array of column information. */
  return allcols;
}




















/************************************************************************/
/***************             Read a txt table             ***************/
/************************************************************************/
static void
txt_fill_columns(char *line, char **tokens, size_t maxcolnum,
                 gal_data_t *colinfo, gal_data_t *out, size_t lineind,
                 size_t lineno, char *filename)
{
  size_t n=0;
  gal_data_t *col;
  char *tailptr, *end=line+strlen(line);

  char         **str, **strb;
  unsigned char  *uc,   *ucb;
  char            *c,    *cb;
  unsigned short *us,   *usb;
  short           *s,    *sb;
  unsigned int   *ui,   *uib;
  int             *i,    *ib;
  unsigned long  *ul,   *ulb;
  long            *l,    *lb;
  LONGLONG        *L,    *Lb;
  float           *f,    *fb;
  double          *d,    *db;

  /* See explanations in `txt_info_from_row'. */
  *(end-1)='\0';

  /* Start parsing the line. Note that `n' and `maxcolnum' start from
     one. So we need column `maxcolnum'.*/
  while(++n)
    {
      /* Break out of the parsing if we don't need the columns any more. */
      if(n>maxcolnum) break;

      /* Set the pointer to the start of this token/column. See
         explanations in `txt_info_from_row'. */
      if( colinfo[n-1].type == GAL_DATA_TYPE_STRING )
        {
          while(isspace(*line) || *line==',') ++line;
          line = (tokens[n]=line) + colinfo[n-1].disp_width;
          if(line<end) *line++='\0';
        }
      else
        tokens[n]=strtok_r(n==1?line:NULL, GAL_TXT_DELIMITERS, &line);
    }

  /* For a sanity check:
  printf("row: %zu: ", lineind+1);
  for(n=1;n<=maxcolnum;++n) printf("-%s-, ", tokens[n]);
  printf("\n");
  */

  /* Read the desired tokens into the columns that need them. Note that
     when a blank value is defined for the column, the column's array
     pointer (`colinfo[col->status-1]') is not NULL and points to the blank
     value. For strings, this will actually be a string. */
  for(col=out; col!=NULL; col=col->next)
    {
      /* Read the proper token into the column. */
      switch(col->type)
        {
        case GAL_DATA_TYPE_STRING:
          str=col->array;
          gal_checkset_allocate_copy(txt_trim_space(tokens[col->status]),
                                     &str[lineind]);
          if( (strb=colinfo[col->status-1].array)
              && !strcmp( *strb, str[lineind] ) )
            {
              free(str[lineind]);
              str[lineind]=GAL_DATA_BLANK_STRING;
            }
          break;

        case GAL_DATA_TYPE_UCHAR:
          uc=col->array;
          uc[lineind]=strtol(tokens[col->status], &tailptr, 0);
          if( (ucb=colinfo[col->status-1].array) && *ucb==uc[lineind] )
            uc[lineind]=GAL_DATA_BLANK_UCHAR;
          break;

        case GAL_DATA_TYPE_CHAR:
          c=col->array;
          c[lineind]=strtol(tokens[col->status], &tailptr, 0);
          if( (cb=colinfo[col->status-1].array) && *cb==c[lineind] )
            c[lineind]=GAL_DATA_BLANK_CHAR;
          break;

        case GAL_DATA_TYPE_USHORT:
          us=col->array;
          us[lineind]=strtol(tokens[col->status], &tailptr, 0);
          if( (usb=colinfo[col->status-1].array) && *usb==us[lineind] )
            us[lineind]=GAL_DATA_BLANK_USHORT;
          break;

        case GAL_DATA_TYPE_SHORT:
          s=col->array;
          s[lineind]=strtol(tokens[col->status], &tailptr, 0);
          if( (sb=colinfo[col->status-1].array) && *sb==s[lineind] )
            s[lineind]=GAL_DATA_BLANK_SHORT;
          break;

        case GAL_DATA_TYPE_UINT:
          ui=col->array;
          ui[lineind]=strtol(tokens[col->status], &tailptr, 0);
          if( (uib=colinfo[col->status-1].array) && *uib==ui[lineind] )
            ui[lineind]=GAL_DATA_BLANK_UINT;
          break;

        case GAL_DATA_TYPE_INT:
          i=col->array;
          i[lineind]=strtol(tokens[col->status], &tailptr, 0);
          if( (ib=colinfo[col->status-1].array) && *ib==i[lineind] )
            i[lineind]=GAL_DATA_BLANK_INT;
          break;

        case GAL_DATA_TYPE_ULONG:
          ul=col->array;
          ul[lineind]=strtoul(tokens[col->status], &tailptr, 0);
          if( (ulb=colinfo[col->status-1].array) && *ulb==ul[lineind] )
            ul[lineind]=GAL_DATA_BLANK_ULONG;
          break;

        case GAL_DATA_TYPE_LONG:
          l=col->array;
          l[lineind]=strtol(tokens[col->status], &tailptr, 0);
          if( (lb=colinfo[col->status-1].array) && *lb==l[lineind] )
            l[lineind]=GAL_DATA_BLANK_LONG;
          break;

        case GAL_DATA_TYPE_LONGLONG:
          L=col->array;
          L[lineind]=strtoll(tokens[col->status], &tailptr, 0);
          if( (Lb=colinfo[col->status-1].array) && *Lb==L[lineind] )
            L[lineind]=GAL_DATA_BLANK_LONGLONG;
          break;

        /* For the blank value of floating point types, we need to make
           sure it isn't a NaN, because a NaN value will fail on any
           condition check (even `=='). If it isn't NaN, then we can
           compare the values. */
        case GAL_DATA_TYPE_FLOAT:
          f=col->array;
          f[lineind]=strtof(tokens[col->status], &tailptr);
          if( (fb=colinfo[col->status-1].array)
              && ( (isnan(*fb) && isnan(f[lineind])) || *fb==f[lineind] ) )
            f[lineind]=GAL_DATA_BLANK_FLOAT;
          break;

        case GAL_DATA_TYPE_DOUBLE:
          d=col->array;
          d[lineind]=strtod(tokens[col->status], &tailptr);
          if( (db=colinfo[col->status-1].array)
              && ( (isnan(*db) && isnan(d[lineind])) || *db==d[lineind] ) )
            d[lineind]=GAL_DATA_BLANK_DOUBLE;
          break;

        default:
          error(EXIT_FAILURE, 0, "type code %d not recognized in "
                "`txt_fill_columns'", col->type);
        }

      /* If a number couldn't be read properly, then report an error. */
      if(col->type!=GAL_DATA_TYPE_STRING && *tailptr!='\0')
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "column %d "
                      "(`%s') couldn't be read as a `%s' number",
                      col->status, tokens[col->status],
                      gal_data_type_as_string(col->type, 1) );
    }
}





gal_data_t *
gal_txt_table_read(char *filename, size_t numrows, gal_data_t *colinfo,
                   struct gal_linkedlist_sll *indexll, int minmapsize)
{
  FILE *fp;
  char *line;
  long dsize;
  char **tokens;
  gal_data_t *out=NULL;
  struct gal_linkedlist_sll *ind;
  size_t maxcolnum=0, lineind=0, lineno=0;
  size_t linelen=10; /* `linelen' will be increased by `getline'. */


  /* Open the file. */
  errno=0;
  fp=fopen(filename, "r");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s: could't open to read as a text table",
          filename);


  /* Allocate the space necessary to keep a copy of each line as we parse
     it. Note that `getline' is going to later `realloc' this space to fit
     the line length. */
  errno=0;
  line=malloc(linelen*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `line' in `gal_txt_table_read'",
          linelen*sizeof *line);


  /* Allocate all the desired columns for output. We will be reading the
     text file line by line, and writing in the necessary values of each
     row individually. */
  for(ind=indexll; ind!=NULL; ind=ind->next)
    {
      dsize=numrows;
      maxcolnum = maxcolnum>ind->v+1 ? maxcolnum : ind->v+1;
      gal_data_add_to_ll(&out, NULL, colinfo[ind->v].type, 1, &dsize, NULL,
                         0, minmapsize, colinfo[ind->v].name,
                         colinfo[ind->v].unit, colinfo[ind->v].comment);
      out->disp_width=colinfo[ind->v].disp_width;
      out->status=ind->v+1;
    }


  /* Allocate the space to keep the pointers to each token in the
     line. This is done here to avoid having to allocate/free this array
     for each line in `txt_fill_columns'. Note that the column numbers are
     counted from one (unlike indexes that are counted from zero), so we
     need `maxcolnum+1' elements in the array of tokens.*/
  errno=0;
  tokens=malloc((maxcolnum+1)*sizeof *tokens);
  if(tokens==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `tokens' in "
          "`gal_txt_table_read'", (maxcolnum+1)*sizeof *tokens);


  /* Read the data columns. */
  while( getline(&line, &linelen, fp) != -1 )
    {
      ++lineno;
      if( get_line_stat(line) == TXT_LINESTAT_DATAROW )
        txt_fill_columns(line, tokens, maxcolnum, colinfo, out, lineind++,
                         lineno, filename);
    }


  /* Clean up and close the file. */
  errno=0;
  if(fclose(fp))
    error(EXIT_FAILURE, errno, "%s: couldn't close file after reading ASCII "
          "table information", filename);
  free(tokens);
  free(line);

  /* Return the array of column information. */
  return out;
}


















/************************************************************************/
/***************          Write a txt table               ***************/
/************************************************************************/
/* Make an array of two strings for each column (in practice a two
   dimensional array with 2 columns and a row for each input column). */
#define FMTS_COLS 3
static char **
make_fmts_for_printf(gal_data_t *cols, size_t numcols, int leftadjust,
                     size_t *len)
{
  char **fmts;
  size_t i=0, j;
  gal_data_t *col;
  char *fmt=NULL, *lng, **strarr;
  char bfmt[GAL_TXT_MAX_FMT_LENGTH];
  int width=0, precision=0, maxstrlen;


  /* Allocate space for the output. */
  errno=0;
  fmts=malloc(FMTS_COLS*numcols*sizeof *fmts);
  if(fmts==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for fmts in `gal_txt_write'",
          FMTS_COLS*numcols*sizeof *fmts);


  /* Initialize the length to 0. */
  *len=0;


  /* Go over all the columns and make their formats. */
  for(col=cols;col!=NULL;col=col->next)
    {
      /* Initialize */
      lng="";


      /* First allocate the necessary space to keep the string. */
      errno=0;
      fmts[ i*FMTS_COLS   ] = malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts);
      fmts[ i*FMTS_COLS+1 ] = malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts);
      if(fmts[i*FMTS_COLS]==NULL || fmts[i*FMTS_COLS+1]==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for fmts[%zu] or fmts[%zu] "
              "in `make_fmts_for_printf' of txt.c",
              GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts,
              i*FMTS_COLS, i*FMTS_COLS+1);


      /* Write the proper format. */
      switch(col->type)
        {



        case GAL_DATA_TYPE_BIT:
          error(EXIT_FAILURE, 0, "printing of bit types is currently "
                "not supported");
          break;



        case GAL_DATA_TYPE_STRING:
          /* Set the basic information. */
          fmt="s";

          /* If `disp_width' was not set (is negative), go through all the
             strings in the column and find the maximum length to use as
             printing width when no value was given. */
          if(col->disp_width<=0)
            {
              maxstrlen=-1;
              strarr=col->array;
              for(j=0;j<col->size;++j)
                if(strarr[j])
                  maxstrlen = ( strlen(strarr[j]) > maxstrlen
                                ? strlen(strarr[j])
                                : maxstrlen );
              width = maxstrlen==-1 ? GAL_TABLE_DEF_STR_WIDTH : maxstrlen ;
            }
          else
            width = col->disp_width;
          break;



        case GAL_DATA_TYPE_UCHAR:
        case GAL_DATA_TYPE_USHORT:
        case GAL_DATA_TYPE_UINT:
        case GAL_DATA_TYPE_ULONG:

          /* Set the final printing format. */
          switch(col->disp_fmt)
            {
            case GAL_TABLE_DISPLAY_FMT_UDECIMAL: fmt="u"; break;
            case GAL_TABLE_DISPLAY_FMT_OCTAL:    fmt="o"; break;
            case GAL_TABLE_DISPLAY_FMT_HEX:      fmt="X"; break;
            default:                             fmt="u";
            }

          /* If we have a long type, then make changes. */
          if(col->type==GAL_DATA_TYPE_ULONG)
            {
              lng="l";
              width=( col->disp_width<=0 ? GAL_TABLE_DEF_LINT_WIDTH
                      : col->disp_width );
            }
          else width=( col->disp_width<=0 ? GAL_TABLE_DEF_INT_WIDTH
                       : col->disp_width );
          precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_INT_PRECISION
                      : col->disp_precision );
          break;



        case GAL_DATA_TYPE_CHAR:
        case GAL_DATA_TYPE_LOGICAL:
        case GAL_DATA_TYPE_SHORT:
        case GAL_DATA_TYPE_INT:
          fmt="d";
          width=( col->disp_width<=0 ? GAL_TABLE_DEF_INT_WIDTH
                  : col->disp_width );
          precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_INT_PRECISION
                      : col->disp_precision );
          break;



        case GAL_DATA_TYPE_LONG:
        case GAL_DATA_TYPE_LONGLONG:
          fmt="d";
          lng = col->type==GAL_DATA_TYPE_LONG ? "l" : "ll";
          width=( col->disp_width<=0 ? GAL_TABLE_DEF_LINT_WIDTH
                  : col->disp_width );
          precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_INT_PRECISION
                      : col->disp_precision );
          break;



        case GAL_DATA_TYPE_FLOAT:
        case GAL_DATA_TYPE_DOUBLE:
          switch(col->disp_fmt)
            {
            case GAL_TABLE_DISPLAY_FMT_FLOAT:    fmt="f"; break;
            case GAL_TABLE_DISPLAY_FMT_EXP:      fmt="e"; break;
            case GAL_TABLE_DISPLAY_FMT_GENERAL:  fmt="g"; break;
            default:                             fmt="f";
            }
          width=( col->disp_width<=0
                  ? ( col->type==GAL_DATA_TYPE_FLOAT
                      ? GAL_TABLE_DEF_FLT_WIDTH
                      : GAL_TABLE_DEF_DBL_WIDTH )
                  : col->disp_width );
          precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_FLT_PRECISION
                      : col->disp_precision );
          break;



        default:
          error(EXIT_FAILURE, 0, "type code %d not recognized for output "
                "column %zu (counting from 1)", col->type, i+1);
        }


      /* Print the blank value if there is any blank values in this
         column. */
      if(gal_data_has_blank(col))
        {
          sprintf(bfmt, "%%%s%s", lng, fmt);
          switch(col->type)
            {
            case GAL_DATA_TYPE_STRING:
              gal_checkset_allocate_copy(GAL_TXT_STRING_BLANK,
                                         &fmts[i*FMTS_COLS+2]);
              break;
            case GAL_DATA_TYPE_UCHAR:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (unsigned char)GAL_DATA_BLANK_UCHAR);
              break;
            case GAL_DATA_TYPE_CHAR:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (char)GAL_DATA_BLANK_CHAR);
              break;
            case GAL_DATA_TYPE_USHORT:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (unsigned short)GAL_DATA_BLANK_USHORT);
              break;
            case GAL_DATA_TYPE_SHORT:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (short)GAL_DATA_BLANK_SHORT);
              break;
            case GAL_DATA_TYPE_UINT:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (unsigned int)GAL_DATA_BLANK_UINT);
              break;
            case GAL_DATA_TYPE_INT:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (int)GAL_DATA_BLANK_INT);
              break;
            case GAL_DATA_TYPE_ULONG:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (unsigned long)GAL_DATA_BLANK_ULONG);
              break;
            case GAL_DATA_TYPE_LONG:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (long)GAL_DATA_BLANK_LONG);
              break;
            case GAL_DATA_TYPE_LONGLONG:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt,
                       (LONGLONG)GAL_DATA_BLANK_LONGLONG);
              break;
            case GAL_DATA_TYPE_FLOAT:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt, GAL_DATA_BLANK_FLOAT);
              break;
            case GAL_DATA_TYPE_DOUBLE:
              asprintf(&fmts[i*FMTS_COLS+2], bfmt, GAL_DATA_BLANK_DOUBLE);
              break;
            }

          /* Adjust the width based on the blank value. */
          width = ( strlen(fmts[i*FMTS_COLS+2]) > width
                    ? strlen(fmts[i*FMTS_COLS+2]) : width );
        }
      else fmts[i*FMTS_COLS+2]=NULL;


      /* Print the result into the allocated string and add its length to
         the final length of the overall format statement. The space in the
         end of `fmts[i*2]' is to ensure that the columns don't merge, even
         if the printed string is larger than the expected width. */
      if(precision<=0)
        *len += 1 + sprintf(fmts[i*FMTS_COLS], "%%%s%d.%d%s%s ",
                            leftadjust ? "-" : "", width, precision,
                            lng, fmt);
      else
        *len += 1 + sprintf(fmts[i*FMTS_COLS], "%%%s%d%s%s ",
                            leftadjust ? "-" : "", width, lng, fmt);


      /* Set the string for the Gnuastro type. For strings, we also need to
         write the maximum number of characters.*/
      if(col->type==GAL_DATA_TYPE_STRING)
        sprintf(fmts[i*FMTS_COLS+1], "%s%d",
                gal_data_type_as_string(col->type, 0), width);
      else
        strcpy(fmts[i*FMTS_COLS+1], gal_data_type_as_string(col->type, 0));


      /* Update the column width and precision: */
      col->disp_width = width;
      col->disp_precision = precision;


      /* Increment the column counter. */
      ++i;
    }

  /* Return the array. */
  return fmts;
}





void
gal_txt_table_write(gal_data_t *cols, char *comment, char *filename,
                    int dontdelete)
{
  FILE *fp;
  gal_data_t *col;
  char **fmts, *tmp;
  size_t i, j, numcols=0, fmtlen;
  int iw=0, nw=0, uw=0, tw=0, bw=0;


  /* Find the number of columns, do a small sanity check, and get the
     maximum width of the name and unit string if they are present. */
  for(col=cols;col!=NULL;col=col->next)
    {
      /* Count. */
      ++numcols;

      /* Make sure the columns are 1 dimensional. */
      if(cols->ndim!=1)
        error(EXIT_FAILURE, 0, "columns to print as an ASCII file must have "
              "only one dimension. column %zu of the given set has %zu "
              "dimensions", numcols, cols->ndim);

      /* Make sure sizes match. */
      if(cols->size!=col->size)
        error(EXIT_FAILURE, 0, "to print a set of columns, as an ASCII "
              "table, they must currently all have the same number of "
              "elements/rows. The inputs to `gal_txt_write' have different "
              "sizes: the first column has %zu, while column %zu as %zu "
              "elements", cols->size, numcols, col->size);
      if( col->name && strlen(col->name)>nw ) nw=strlen(col->name);
      if( col->unit && strlen(col->unit)>uw ) uw=strlen(col->unit);
    }


  /* Set the output FILE pointer: if it isn't NULL, its an actual file,
     otherwise, its the standard output. */
  if(filename)
    {
      gal_checkset_check_remove_file(filename, dontdelete);
      errno=0;
      fp=fopen(filename, "w");
      if(fp==NULL)
        error(EXIT_FAILURE, errno, "%s: couldn't be open to write text "
              "table", filename);
    }
  else
    fp=stdout;


  /* Prepare the necessary formats for each column, then allocate the space
     for the full list and concatenate all the separate inputs into it. */
  fmts=make_fmts_for_printf(cols, numcols, 1, &fmtlen);
  for(i=0;i<numcols;++i)
    {
      if( (tmp=fmts[ i*FMTS_COLS+1 ]) )
        tw = strlen(tmp) > tw ? strlen(tmp) : tw;
      if( (tmp=fmts[ i*FMTS_COLS+2 ]) )
        bw = strlen(tmp) > bw ? strlen(tmp) : bw;
    }


  /* Write the comments if there are any */
  if(comment) fprintf(fp, "%s\n", comment);


  /* Write the information for each column */
  i=0;
  iw=log10(numcols)+1;
  for(col=cols;col!=NULL;col=col->next)
    {
      fprintf(fp, "# Column %-*zu: %-*s [%-*s,%-*s,%-*s] %s\n",
              iw, i+1,
              nw, col->name    ? col->name    : "",
              uw, col->unit    ? col->unit    : "",
              tw, fmts[i*FMTS_COLS+1] ? fmts[i*FMTS_COLS+1] : "",
              bw, fmts[i*FMTS_COLS+2] ? fmts[i*FMTS_COLS+2] : "",
              col->comment ? col->comment : "");
      ++i;
    }


  /* Print the output */
  for(i=0;i<cols->size;++i)                   /* Loop over each row.    */
    {
      j=0;
      for(col=cols;col!=NULL;col=col->next)   /* Loop over each column. */
        {
          switch(col->type)
            {
            /* Integer types */
            case GAL_DATA_TYPE_UCHAR:
              fprintf(fp, fmts[j*FMTS_COLS],
                      ((unsigned char *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_CHAR:
            case GAL_DATA_TYPE_LOGICAL:
              fprintf(fp, fmts[j*FMTS_COLS], ((char *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_USHORT:
              fprintf(fp, fmts[j*FMTS_COLS],
                      ((unsigned short *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_SHORT:
              fprintf(fp, fmts[j*FMTS_COLS], ((short *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_UINT:
              fprintf(fp, fmts[j*FMTS_COLS], ((unsigned int *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_INT:
              fprintf(fp, fmts[j*FMTS_COLS], ((int *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_ULONG:
              fprintf(fp, fmts[j*FMTS_COLS],
                      ((unsigned long *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_LONG:
              fprintf(fp, fmts[j*FMTS_COLS], ((long *)col->array)[i]);
              break;
            case GAL_DATA_TYPE_LONGLONG:
              fprintf(fp, fmts[j*FMTS_COLS], ((LONGLONG *)col->array)[i]);
              break;

            /* Special consideration (string, float, double). */
            case GAL_DATA_TYPE_STRING:
              if(((char **)col->array)[i])
                fprintf(fp, fmts[j*FMTS_COLS], ((char **)col->array)[i]);
              else
                fprintf(fp, "%-*s ", col->disp_width, GAL_TXT_STRING_BLANK);
              break;

            case GAL_DATA_TYPE_FLOAT:
              fprintf(fp, fmts[j*FMTS_COLS], ((float *)col->array)[i]);
              break;

            case GAL_DATA_TYPE_DOUBLE:
              fprintf(fp, fmts[j*FMTS_COLS], ((double *)col->array)[i]);
              break;
            default:
              error(EXIT_FAILURE, 0, "type code %d not recognized for "
                    "col->type in `gal_txt_write'", col->type);
            }
          ++j;
        }
      fprintf(fp, "\n");
    }


  /* Clean up, close the input file and return. For the fmts[i*FMTS_COLS]
     elements, the reason is that fmts[i*FMTS_COLS+1] are literal strings,
     not allocated. So they don't need freeing.*/
  for(i=0;i<numcols;++i)
    {
      free(fmts[i*FMTS_COLS]);
      free(fmts[i*FMTS_COLS+1]);
      free(fmts[i*FMTS_COLS+2]);
    }
  free(fmts);
  if(filename)
    {
      errno=0;
      if(fclose(fp))
        error(EXIT_FAILURE, errno, "%s: couldn't close file after writing "
              "of text table", filename);
    }
}
