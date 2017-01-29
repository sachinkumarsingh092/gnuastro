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
/* Return one of the `txt_line_stat' constant values. */
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
  int index, strw=0;
  int type=GAL_DATA_TYPE_DOUBLE; /* Default type. */
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

      /* Read the column number as an integer. If it can't be read as an
         integer, or is zero or negative then just return without adding
         anything to this line. */
      index=strtol(number, &tailptr, 0);
      if(*tailptr!='\0' || index<=0) return;

      /* If there was no name (the line is just `# Column N:'), then ignore
         the line. Relying on the column count from the first line is more
         robust and less prone to human error, for example typing a number
         larger than the total number of columns.  */
      name=txt_trim_space(name);
      if(name==NULL) return;

      /* If this is a repeated index, ignore it. */
      for(tmp=*colsll; tmp!=NULL; tmp=tmp->next)
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

      /* If `typestr' was given, then check if this is a standard type. If
         `typestr' wasn't specified, then the default double type code will
         be used (see the variable definitions above). If the given type
         isn't a standard type then ignore the line. Just note that if we
         are dealing with the string type, we have to pull out the number
         part first. If there is no number for a string type, then ignore
         the line. */
      if(typestr)
        {
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
              if(type==GAL_DATA_TYPE_INVALID) return;
            }
        }

      /* Add this column's information into the columns linked list. We
         will define the data structur's array to have zero dimensions (no
         array) by default. If there is a blank value its value will be put
         into the array by `gal_table_read_blank'. To keep the name, unit,
         and comment strings, trim the white space before and after each
         before using them here.  */
      gal_data_add_to_ll(colsll, NULL, type, 0, NULL, NULL, 0, -1, name,
                         txt_trim_space(unit), txt_trim_space(comment) );

      /* Put the number of this column into the status variable of the data
         structure. If the type is string, then also copy the width into
         the structure. */
      (*colsll)->status=index;
      (*colsll)->disp_width = type==GAL_DATA_TYPE_STRING ? strw : 0;

      /* Write the blank value into the array. Note that this is not the
         final column, we are just collecting information now. */
      gal_table_read_blank(*colsll, txt_trim_space(blank));
    }
}





/* The input ASCII table might not have had information in its comments, or
   the information might not have been complete. So we need to go through
   the first row of data also. */
void
txt_info_from_first_row(char *line, gal_data_t **colsll)
{
  size_t n=0, maxcnum=0;
  gal_data_t *col, *prev, *tmp;
  char *token, *end=line+strlen(line);

  /* Remove the new line character from the end of the line. If the last
     column is a string, and the given length is larger than the available
     space on the line, we don't want to have the line's new-line
     character. Its better for it to actually be shorter than the space. */
  *(end-1)='\0';

  /* Get the maximum number of columns read from the comment
     information. */
  for(col=*colsll; col!=NULL; col=col->next)
    maxcnum = maxcnum>col->status ? maxcnum : col->status;

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

          /* A token exists, so set this column to the default double type
             with no information, then set its status value to the column
             number. */
          gal_data_add_to_ll(colsll, NULL, GAL_DATA_TYPE_DOUBLE, 0, NULL,
                             NULL, 0, -1, NULL, NULL, NULL);
          (*colsll)->status=n;
        }
    }

  /* If the number of columns given by the comments is larger than the
     actual number of lines, remove those that have larger numbers from the
     linked list before things get complicated outside of this function. */
  if(maxcnum>n)
    {
      prev=NULL;
      col=*colsll;
      while(col!=NULL)
        {
          if(col->status > n) /* Column has no data (was only in comments) */
            {
              /* This column has to be removed/freed. But we have to make
                 some corrections before freeing it:

                  - When `prev==NULL', then we still haven't got to the
                    first valid element yet and must free this one, but if
                    we do that, then the main pointer to the start of the
                    list will be lost (we will loose all connections with
                    the chain after leaving this loop). So we need to set
                    that to the next element.

                  - When there actually was a previous element
                    (`prev!=NULL'), then we must correct it's next
                    pointer. Otherwise we will break up the chain.*/
              if(prev) prev->next=col->next; else *colsll=col->next;
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





/* Return the information about a text file table. If there were no
   readable rows, it will return NULL.*/
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
    switch( gal_txt_line_stat(line) )
      {
      case GAL_TXT_LINESTAT_COMMENT:
        /* Line is a comment, see if it has formatted information. */
        txt_info_from_comment(line, &colsll);
        break;

        /* Line is actual data, use it to fill in the gaps.  */
      case GAL_TXT_LINESTAT_DATAROW:
        ++(*numrows);
        if(firstlinedone==0)
          {
            firstlinedone=1;
            txt_info_from_first_row(line, &colsll);
          }
        break;

        /* We also have the case of GAL_TXT_LINESTAT_BLANK */
      }


  /* If there were rows in the file, then write the unorganized gathered
     information (linked list) into an organized array for easy processing
     by later steps.  */
  allcols = *numrows ? txt_infoll_to_array(colsll, numcols) : NULL;


  /* Clean up. Note that even if there were no usable columns, there might
     have been meta-data comments, so we need to free `colsll' in any
     case. If the list is indeed empty, then `gal_data_free_ll' won't do
     anything. */
  free(line);
  gal_data_free_ll(colsll);


  /* Close the file. */
  errno=0;
  if(fclose(fp))
    error(EXIT_FAILURE, errno, "%s: couldn't close file after reading ASCII "
          "table information", filename);


  /* Return the array of column information. */
  return allcols;
}




















/************************************************************************/
/***************             Read a txt table             ***************/
/************************************************************************/
static void
txt_fill_columns(char *line, char **tokens, size_t maxcolnum,
                 gal_data_t *colinfo, gal_data_t *out, size_t rowind,
                 size_t lineno, char *filename)
{
  size_t n=0;
  int notenoughcols=0;
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
      /* Break out of the parsing if we don't need the columns any
         more. The table might contain many more columns, but when they
         aren't needed, there is no point in tokenizing them. */
      if(n>maxcolnum) break;

      /* Set the pointer to the start of this token/column. See
         explanations in `txt_info_from_row'. */
      if( colinfo[n-1].type == GAL_DATA_TYPE_STRING )
        {
          /* Remove any delimiters and stop at the first non-delimiter. If
             we have reached the end of the line then its an error, because
             we were expecting a column here. */
          while(isspace(*line) || *line==',') ++line;
          if(*line=='\0') {notenoughcols=1; break;}

          /* Everything is good, set the pointer and increment the line to
             the end of the allocated space for this string. */
          line = (tokens[n]=line) + colinfo[n-1].disp_width;
          if(line<end) *line++='\0';
        }
      else
        {
          /* If we have reached the end of the line, then `strtok_r' will
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
                  n-1); /* This must be `n-1' (since n starts from 1). */

  /* For a sanity check:
  printf("row: %zu: ", rowind+1);
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
                                     &str[rowind]);
          if( (strb=colinfo[col->status-1].array)
              && !strcmp( *strb, str[rowind] ) )
            {
              free(str[rowind]);
              gal_checkset_allocate_copy(GAL_DATA_BLANK_STRING,
                                         &str[rowind]);
            }
          break;

        case GAL_DATA_TYPE_UCHAR:
          uc=col->array;
          uc[rowind]=strtol(tokens[col->status], &tailptr, 0);
          if( (ucb=colinfo[col->status-1].array) && *ucb==uc[rowind] )
            uc[rowind]=GAL_DATA_BLANK_UCHAR;
          break;

        case GAL_DATA_TYPE_CHAR:
          c=col->array;
          c[rowind]=strtol(tokens[col->status], &tailptr, 0);
          if( (cb=colinfo[col->status-1].array) && *cb==c[rowind] )
            c[rowind]=GAL_DATA_BLANK_CHAR;
          break;

        case GAL_DATA_TYPE_USHORT:
          us=col->array;
          us[rowind]=strtol(tokens[col->status], &tailptr, 0);
          if( (usb=colinfo[col->status-1].array) && *usb==us[rowind] )
            us[rowind]=GAL_DATA_BLANK_USHORT;
          break;

        case GAL_DATA_TYPE_SHORT:
          s=col->array;
          s[rowind]=strtol(tokens[col->status], &tailptr, 0);
          if( (sb=colinfo[col->status-1].array) && *sb==s[rowind] )
            s[rowind]=GAL_DATA_BLANK_SHORT;
          break;

        case GAL_DATA_TYPE_UINT:
          ui=col->array;
          ui[rowind]=strtol(tokens[col->status], &tailptr, 0);
          if( (uib=colinfo[col->status-1].array) && *uib==ui[rowind] )
            ui[rowind]=GAL_DATA_BLANK_UINT;
          break;

        case GAL_DATA_TYPE_INT:
          i=col->array;
          i[rowind]=strtol(tokens[col->status], &tailptr, 0);
          if( (ib=colinfo[col->status-1].array) && *ib==i[rowind] )
            i[rowind]=GAL_DATA_BLANK_INT;
          break;

        case GAL_DATA_TYPE_ULONG:
          ul=col->array;
          ul[rowind]=strtoul(tokens[col->status], &tailptr, 0);
          if( (ulb=colinfo[col->status-1].array) && *ulb==ul[rowind] )
            ul[rowind]=GAL_DATA_BLANK_ULONG;
          break;

        case GAL_DATA_TYPE_LONG:
          l=col->array;
          l[rowind]=strtol(tokens[col->status], &tailptr, 0);
          if( (lb=colinfo[col->status-1].array) && *lb==l[rowind] )
            l[rowind]=GAL_DATA_BLANK_LONG;
          break;

        case GAL_DATA_TYPE_LONGLONG:
          L=col->array;
          L[rowind]=strtoll(tokens[col->status], &tailptr, 0);
          if( (Lb=colinfo[col->status-1].array) && *Lb==L[rowind] )
            L[rowind]=GAL_DATA_BLANK_LONGLONG;
          break;

        /* For the blank value of floating point types, we need to make
           sure it isn't a NaN, because a NaN value will fail on any
           condition check (even `=='). If it isn't NaN, then we can
           compare the values. */
        case GAL_DATA_TYPE_FLOAT:
          f=col->array;
          f[rowind]=strtod(tokens[col->status], &tailptr);
          if( (fb=colinfo[col->status-1].array)
              && ( (isnan(*fb) && isnan(f[rowind])) || *fb==f[rowind] ) )
            f[rowind]=GAL_DATA_BLANK_FLOAT;
          break;

        case GAL_DATA_TYPE_DOUBLE:
          d=col->array;
          d[rowind]=strtod(tokens[col->status], &tailptr);
          if( (db=colinfo[col->status-1].array)
              && ( (isnan(*db) && isnan(d[rowind])) || *db==d[rowind] ) )
            d[rowind]=GAL_DATA_BLANK_DOUBLE;
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
  size_t dsize;
  char **tokens;
  gal_data_t *out=NULL;
  struct gal_linkedlist_sll *ind;
  size_t maxcolnum=0, rowind=0, lineno=0;
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
      if( gal_txt_line_stat(line) == GAL_TXT_LINESTAT_DATAROW )
        txt_fill_columns(line, tokens, maxcolnum, colinfo, out, rowind++,
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
  size_t i=0;
  char **fmts;
  gal_data_t *col;
  char fmt[2], lng[3];


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
      /* First allocate the necessary space to keep the string. */
      errno=0;
      fmts[ i*FMTS_COLS   ] = malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts);
      fmts[ i*FMTS_COLS+1 ] = malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts);
      if(fmts[i*FMTS_COLS]==NULL || fmts[i*FMTS_COLS+1]==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for fmts[%zu] or fmts[%zu] "
              "in `make_fmts_for_printf' of txt.c",
              GAL_TXT_MAX_FMT_LENGTH*sizeof **fmts,
              i*FMTS_COLS, i*FMTS_COLS+1);


      /* If we have a blank value, get the blank value as a string and
         adjust the width */
      if(gal_data_has_blank(col)==0)
        fmts[i*FMTS_COLS+2]=NULL;
      else
        {
          /* Set the blank value. */
          if(col->type==GAL_DATA_TYPE_STRING)
            gal_checkset_allocate_copy(GAL_DATA_BLANK_STRING,
                                       &fmts[i*FMTS_COLS+2]);
          else
            fmts[i*FMTS_COLS+2]=gal_data_blank_as_string(col->type, 0);
        }


      /* Fill in the printing paramters. */
      gal_table_col_print_info(col, GAL_TABLE_FORMAT_TXT, fmt, lng);


      /* Adjust the width if a blank string was defined. */
      if(fmts[i*FMTS_COLS+2])
        col->disp_width = ( strlen(fmts[i*FMTS_COLS+2]) > col->disp_width
                            ? strlen(fmts[i*FMTS_COLS+2])
                            : col->disp_width );


      /* Print the result into the allocated string and add its length to
         the final length of the overall format statement. The space in the
         end of `fmts[i*2]' is to ensure that the columns don't merge, even
         if the printed string is larger than the expected width. */
      if(col->disp_precision > 0)
        *len += 1 + sprintf(fmts[i*FMTS_COLS], "%%%s%d.%d%s%s ",
                            leftadjust ? "-" : "", col->disp_width,
                            col->disp_precision, lng, fmt);
      else
        *len += 1 + sprintf(fmts[i*FMTS_COLS], "%%%s%d%s%s ",
                            leftadjust ? "-" : "", col->disp_width, lng, fmt);


      /* Set the string for the Gnuastro type. For strings, we also need to
         write the maximum number of characters.*/
      if(col->type==GAL_DATA_TYPE_STRING)
        sprintf(fmts[i*FMTS_COLS+1], "%s%d",
                gal_data_type_as_string(col->type, 0), col->disp_width);
      else
        strcpy(fmts[i*FMTS_COLS+1], gal_data_type_as_string(col->type, 0));


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
  char *nstr;
  gal_data_t *col;
  char **fmts, *tmp;
  size_t i, j, numcols=0, fmtlen;
  int nlen, nw=0, uw=0, tw=0, bw=0;


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


      /* Write the comments if there were any. */
      if(comment) fprintf(fp, "%s\n", comment);


      /* Write the column information if the output is a file. When the
         output is directed to standard output (the command-line), it is
         most probably intended for piping into another program (for
         example AWK for further processing, or sort, or anything) so the
         user already has the column information and is probably going to
         change them, so they are just a nuisance.

         When there are more than 9 columns, we don't want to have cases
         like `# Column 1 :' (note the space between `1' and `:', this
         space won't exist for the 2 digit colum numbers).

         To do this, we are first allocating and printing a string long
         enough to keep the final column's `N:'. Then, for each column, we
         print only the number into the allocated space and put the `:' in
         manually immediately after the number. Note that the initial
         `asprintf' put a `\0' in the allocated space, so we can safely
         over-write the one that `sprintf' puts with a `:' for the columns
         that have the same number of digits as the final column. */
      i=0;
      asprintf(&nstr, "%zu:", numcols);
      nlen=strlen(nstr);
      for(col=cols; col!=NULL; col=col->next)
        {
          /* Print the number into the number string, then add the `:'
             immediately after the number. */
          sprintf(nstr, "%zu", i+1);
          for(j=1;j<nlen;++j)
            if(!isdigit(nstr[j])) nstr[j] = isdigit(nstr[j-1]) ? ':' : ' ';

          /* Now print the full column information. */
          fprintf(fp, "# Column %s %-*s [%-*s,%-*s,%-*s] %s\n",
                  nstr,
                  nw, col->name    ? col->name    : "",
                  uw, col->unit    ? col->unit    : "",
                  tw, fmts[i*FMTS_COLS+1] ? fmts[i*FMTS_COLS+1] : "",
                  bw, fmts[i*FMTS_COLS+2] ? fmts[i*FMTS_COLS+2] : "",
                  col->comment ? col->comment : "");
          ++i;
        }


      /* Clean up */
      free(nstr);
    }
  else      /* Output wasn't a file, so set it to standard output */
    fp=stdout;


  /* Print the output */
  for(i=0;i<cols->size;++i)                   /* Loop over each row.    */
    {
      j=0;
      for(col=cols;col!=NULL;col=col->next)   /* Loop over each column. */
        {
          switch(col->type)
            {

            /* Numerical types. */
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

            case GAL_DATA_TYPE_FLOAT:
              fprintf(fp, fmts[j*FMTS_COLS], ((float *)col->array)[i]);
              break;

            case GAL_DATA_TYPE_DOUBLE:
              fprintf(fp, fmts[j*FMTS_COLS], ((double *)col->array)[i]);
              break;

            /* Special consideration for strings. */
            case GAL_DATA_TYPE_STRING:
              if( !strcmp( ((char **)col->array)[i], GAL_DATA_BLANK_STRING ) )
                fprintf(fp, "%-*s ", col->disp_width, GAL_DATA_BLANK_STRING);
              else
                fprintf(fp, fmts[j*FMTS_COLS], ((char **)col->array)[i]);
              break;

            default:
              error(EXIT_FAILURE, 0, "type code %d not recognized for "
                    "col->type in `gal_txt_write'", col->type);
            }
          ++j;
        }
      fprintf(fp, "\n");
    }


  /* Clean up. */
  for(i=0;i<numcols;++i)
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
              "of text table", filename);
    }
}
