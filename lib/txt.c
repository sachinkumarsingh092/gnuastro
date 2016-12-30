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
txt_no_space_before_after(char *str)
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
  char **strarr, *tailptr;

  /* If there is nothing to use as blank, then free the array. */
  if(blank==NULL)
    {
      free(col->array);
      return;
    }

  /* String type. Copy the string.*/
  if(col->type==GAL_DATA_TYPE_STRING)
    {
      strarr=col->array;
      gal_checkset_allocate_copy(blank, &strarr[0]);
    }

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
  long dsize=1;
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
      typestr=txt_no_space_before_after(typestr);
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
      gal_data_add_to_ll(colsll, NULL, type, 1, &dsize, NULL, 0, -1,
                         txt_no_space_before_after(name),
                         txt_no_space_before_after(unit),
                         txt_no_space_before_after(comment) );

      /* Put the number of this column into the status variable of the data
         structure. If the type is string, then also copy the width into
         the structure. */
      (*colsll)->status=index;
      if(type==GAL_DATA_TYPE_STRING) (*colsll)->disp_width=strw;

      /* Write the blank value into the array. Note that this is not the
         final column, we are just collecting information now. */
      txt_read_blank(*colsll, txt_no_space_before_after(blank));
    }
}





/* The input ASCII table might not have had information in its comments, or
   the information might not have been complete. So we need to go through
   the first row of data also. */
void
txt_info_from_row(char *line, gal_data_t **colsll)
{
  size_t i=0;
  gal_data_t *col;
  char *token, *end=line+strlen(line);

  /* Remove the new line character from the end of the line. If the last
     column is a string, and the given length is larger than the
     available space on the line, we don't want to  */
  *(end-1)=' ';

  /* Go over the line check/fill the column information. */
  while(++i)
    {
      /* Check if there is information for this column. */
      for(col=*colsll; col!=NULL; col=col->next) if(col->status==i) break;

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
              token=strtok_r(i==1?line:NULL, GAL_TXT_DELIMITERS, &line);
              if(token==NULL) break;
              /* printf(" col %zu: =%s=\n", i, token); */
            }
        }
      else
        {
          /* Make sure a token exists in this undefined column. */
          token=strtok_r(i==1?line:NULL, GAL_TXT_DELIMITERS, &line);
          if(token==NULL) break;
          /* printf(" col %zu: *%s*\n", i, token); */

          /* A token exists, so set this column to the default double type
             with no information, then set its status value to the column
             number. */
          gal_data_add_to_ll(colsll, NULL, GAL_DATA_TYPE_DOUBLE, 0, NULL,
                             NULL, 0, -1, NULL, NULL, NULL);
          (*colsll)->status=i;
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
txt_infoll_to_array(gal_data_t *colsll)
{
  size_t numcols=0;
  gal_data_t *col, *allcols;

  /* First find the total number of columns. */
  for(col=colsll;col!=NULL;col=col->next)
    numcols = numcols > col->status ? numcols : col->status;

  /* Now, allocate the array and put in the values. */
  errno=0;
  allcols=calloc(numcols, sizeof *allcols);
  if(allcols==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `allcols' in "
          "`txt_infoll_to_array'", numcols*sizeof *allcols);

  /* Put each column into its proper place in the array. We are setting all
     the allocated spaces in the linked list elements to NULL, because we
     didn't initialize the array of column information in the allocation
     above and we don't want to re-allocate everything (because freeing the
     linked list will free them also). */
  for(col=colsll;col!=NULL;col=col->next)
    {
      allcols[col->status].name=col->name;          col->name=NULL;
      allcols[col->status].unit=col->unit;          col->unit=NULL;
      allcols[col->status].array=col->array;        col->array=NULL;
      allcols[col->status].dsize=col->dsize;        col->dsize=NULL;
      allcols[col->status].comment=col->comment;    col->comment=NULL;

      allcols[col->status].type=col->type;
      allcols[col->status].ndim=col->ndim;
      allcols[col->status].size=col->size;
      allcols[col->status].disp_width=col->disp_width;
    }

  /* Return the array of all column information. */
  return allcols;
}





/* Return the information about a text file table. */
gal_data_t *
gal_txt_table_info(char *filename, size_t *numcols)
{
  FILE *fp;
  char *line;
  gal_data_t *colsll=NULL, *allcols;
  size_t linelen=10; /* `linelen' will be increased by `getline'. */


  /* Open the file. */
  errno=0;
  fp=fopen(filename, "r");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s: could't open to read as a text table",
          filename);


  /* Get the maximum line length and allocate the space necessary to keep
     copies of all lines as we parse them. Note that `getline' is going to
     put the string NULL character also, so we need one more character. */
  errno=0;
  line=malloc(linelen*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for line in `gal_txt_table_info'",
          linelen*sizeof *line);


  /* Read the comments of the line for possible information about the
     lines, but also confirm the info by trying to read the first
     uncommented line. */
  while( getline(&line, &linelen, fp) != -1 )
    {
      /* Line is a comment, see if it has formatted information. */
      if( get_line_stat(line) == TXT_LINESTAT_ISCOMMENT )
        txt_info_from_comment(line, &colsll);

      /* Line is actual data, use it to fill in the gaps.  */
      if( get_line_stat(line) == TXT_LINESTAT_DATAROW )
        {
          txt_info_from_row(line, &colsll);
          break;
        }
    }


  /* Write the unorganized gathered information (linked list) into an
     organized array for easy processing by later steps. */
  allcols=txt_infoll_to_array(colsll);


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
/***************          Write a txt table               ***************/
/************************************************************************/
/* Make an array of two strings for each column (in practice a two
   dimensional array with 2 columns and a row for each input column). */
static char **
make_fmts_for_printf(gal_data_t *cols, size_t numcols, int leftadjust,
                     size_t *len)
{
  char **fmts;
  gal_data_t *tmp;
  int width=0, precision=0;
  size_t i=0, j, maxstrlen=0;
  char *fmt=NULL, *lng, **strarr;

  /* Allocate space for the output. */
  errno=0;
  fmts=malloc(2*numcols*sizeof *fmts);
  if(fmts==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for fmts in `gal_txt_write'",
          numcols*sizeof *fmts);

  /* Initialize the length to 0. */
  *len=0;

  /* Go over all the columns and make their formats. */
  for(tmp=cols;tmp!=NULL;tmp=tmp->next)
    {
      /* Initialize */
      lng="";

      /* First allocate the necessary space to keep the string. */
      errno=0;
      fmts[i*2]=malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof *fmts[i*2]);
      fmts[i*2+1]=malloc(GAL_TXT_MAX_FMT_LENGTH*sizeof *fmts[i*2+1]);
      if(fmts[i*2]==NULL || fmts[i*2+1]==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for fmts[%zu] or fmts[%zu] "
              "in `make_fmts_for_printf' of txt.c",
              GAL_TXT_MAX_FMT_LENGTH*sizeof *fmts[i], i*2, i*2+1);

      /* Write the proper format. */
      switch(tmp->type)
        {

        case GAL_DATA_TYPE_BIT:
          error(EXIT_FAILURE, 0, "printing of bit types is currently "
                "not supported");
          break;

        case GAL_DATA_TYPE_STRING:
          /* Set the basic information. */
          fmt="s";
          width=( tmp->disp_width<0 ? GAL_TABLE_DEF_STR_WIDTH
                  : tmp->disp_width );
          precision=( tmp->disp_precision<0 ? GAL_TABLE_DEF_STR_PRECISION
                      : tmp->disp_precision );

          /* For strings, we also need the maximum length of all the
             columns, so go through all the strings in the column and find
             the maximum length. */
          strarr=tmp->array;
          for(j=0;j<tmp->size;++j)
            maxstrlen = ( strlen(strarr[j]) > maxstrlen
                          ? strlen(strarr[j])
                          : maxstrlen );
          break;


        case GAL_DATA_TYPE_UCHAR:
        case GAL_DATA_TYPE_USHORT:
        case GAL_DATA_TYPE_UINT:
        case GAL_DATA_TYPE_ULONG:

          /* If we have a long type, then make changes. */
          if(tmp->type==GAL_DATA_TYPE_ULONG)
            {
              lng="l";
              width=( tmp->disp_width<0 ? GAL_TABLE_DEF_LINT_WIDTH
                      : tmp->disp_width );
            }
          else width=( tmp->disp_width<0 ? GAL_TABLE_DEF_INT_WIDTH
                       : tmp->disp_width );
          precision=( tmp->disp_precision<0 ? GAL_TABLE_DEF_INT_PRECISION
                      : tmp->disp_precision );

          /* Set the final printing format. */
          switch(tmp->disp_fmt)
            {
            case GAL_TABLE_DISPLAY_FMT_UDECIMAL: fmt="u"; break;
            case GAL_TABLE_DISPLAY_FMT_OCTAL:    fmt="o"; break;
            case GAL_TABLE_DISPLAY_FMT_HEX:      fmt="X"; break;
            default:                             fmt="u";
            }
          break;


        case GAL_DATA_TYPE_CHAR:
        case GAL_DATA_TYPE_LOGICAL:
        case GAL_DATA_TYPE_SHORT:
        case GAL_DATA_TYPE_INT:
          fmt="d";
          width=( tmp->disp_width<0 ? GAL_TABLE_DEF_INT_WIDTH
                  : tmp->disp_width );
          precision=( tmp->disp_precision<0 ? GAL_TABLE_DEF_INT_PRECISION
                      : tmp->disp_precision );
          break;


        case GAL_DATA_TYPE_LONG:
        case GAL_DATA_TYPE_LONGLONG:
          fmt="d";
          lng = tmp->type==GAL_DATA_TYPE_LONG ? "l" : "ll";
          width=( tmp->disp_width<0 ? GAL_TABLE_DEF_LINT_WIDTH
                  : tmp->disp_width );
          precision=( tmp->disp_precision<0 ? GAL_TABLE_DEF_INT_PRECISION
                      : tmp->disp_precision );
          break;


        case GAL_DATA_TYPE_FLOAT:
        case GAL_DATA_TYPE_DOUBLE:
          switch(tmp->disp_fmt)
            {
            case GAL_TABLE_DISPLAY_FMT_FLOAT:    fmt="f"; break;
            case GAL_TABLE_DISPLAY_FMT_EXP:      fmt="e"; break;
            case GAL_TABLE_DISPLAY_FMT_GENERAL:  fmt="g"; break;
            default:                             fmt="f";
            }
          width=( tmp->disp_width<0
                  ? ( tmp->type==GAL_DATA_TYPE_FLOAT
                      ? GAL_TABLE_DEF_FLT_WIDTH
                      : GAL_TABLE_DEF_DBL_WIDTH )
                  : tmp->disp_width );
          precision=( tmp->disp_precision<0 ? GAL_TABLE_DEF_FLT_PRECISION
                      : tmp->disp_precision );
          break;

        default:
          error(EXIT_FAILURE, 0, "type code %d not recognized for output "
                "column %zu (counting from 1)", tmp->type, i+1);
        }

      /* Print the result into the allocated string and add its length to
         the final length of the overall format statement. The space in the
         end of `fmts[i*2]' is to ensure that the columns don't merge, even
         if the printed string is larger than the expected width. */
      *len += 1 + sprintf(fmts[i*2], "%%%s%d.%d%s%s ", leftadjust ? "-" : "",
                          width, precision, lng, fmt);

      /* Set the string for the Gnuastro type. For strings, we also need to
         write the maximum number of characters.*/
      if(tmp->type==GAL_DATA_TYPE_STRING)
        sprintf(fmts[i*2+1], "%s%zu", gal_data_type_as_string(tmp->type, 0),
                maxstrlen);
      else
        strcpy(fmts[i*2+1], gal_data_type_as_string(tmp->type, 0));


      /* Increment the column counter. */
      ++i;
    }

  /* Return the array. */
  return fmts;
}





void
gal_txt_write(gal_data_t *cols, char *comment, char *filename,
              int dontdelete)
{
  FILE *fp;
  char **fmts;
  gal_data_t *tmp;
  int iw=0, nw=0, uw=0, tw=0;
  size_t i, j, numcols=0, fmtlen;


  /* Find the number of columns, do a small sanity check, and get the
     maximum width of the name and unit string if they are present. */
  for(tmp=cols;tmp!=NULL;tmp=tmp->next)
    {
      /* Count. */
      ++numcols;

      /* Make sure the columns are 1 dimensional. */
      if(cols->ndim!=1)
        error(EXIT_FAILURE, 0, "columns to print as an ASCII file must have "
              "only one dimension. column %zu of the given set has %zu "
              "dimensions", numcols, cols->ndim);

      /* Make sure sizes match. */
      if(cols->size!=tmp->size)
        error(EXIT_FAILURE, 0, "to print a set of columns, as an ASCII "
              "table, they must currently all have the same number of "
              "elements/rows. The inputs to `gal_txt_write' have different "
              "sizes: the first column has %zu, while column %zu as %zu "
              "elements", cols->size, numcols, tmp->size);
      if( tmp->name && strlen(tmp->name)>nw ) nw=strlen(tmp->name);
      if( tmp->unit && strlen(tmp->unit)>uw ) uw=strlen(tmp->unit);
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
     for the full list and concatenate all the separate input into it. */
  fmts=make_fmts_for_printf(cols, numcols, 1, &fmtlen);
  for(i=0;i<numcols;++i)
    tw = strlen( fmts[ i*2+1 ] ) > tw ? strlen( fmts[ i*2+1 ] ) : tw;


  /* Write the comments if there are any */
  if(comment) fprintf(fp, "%s\n", comment);


  /* Write the given information. */
  i=0;
  iw=log10(numcols)+1;
  for(tmp=cols;tmp!=NULL;tmp=tmp->next)
    {
      fprintf(fp, "# Column %-*zu: %-*s [%-*s, %-*s] %s\n",
              iw, i+1,
              nw, tmp->name    ? tmp->name    : "",
              uw, tmp->unit    ? tmp->unit    : "",
              tw, fmts[i*2+1],
              tmp->comment ? tmp->comment : "");
      ++i;
    }


  /* Print the output */
  for(i=0;i<cols->size;++i)                   /* Loop over each row.    */
    {
      j=0;
      for(tmp=cols;tmp!=NULL;tmp=tmp->next)   /* Loop over each column. */
        {
          switch(tmp->type)
            {
            case GAL_DATA_TYPE_UCHAR:
              fprintf(fp, fmts[j*2], ((unsigned char *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_CHAR:
            case GAL_DATA_TYPE_LOGICAL:
              fprintf(fp, fmts[j*2], ((char *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_STRING:
              fprintf(fp, fmts[j*2], ((char **)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_USHORT:
              fprintf(fp, fmts[j*2], ((unsigned short *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_SHORT:
              fprintf(fp, fmts[j*2], ((short *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_UINT:
              fprintf(fp, fmts[j*2], ((unsigned int *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_INT:
              fprintf(fp, fmts[j*2], ((int *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_ULONG:
              fprintf(fp, fmts[j*2], ((unsigned long *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_LONG:
              fprintf(fp, fmts[j*2], ((long *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_LONGLONG:
              fprintf(fp, fmts[j*2], ((LONGLONG *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_FLOAT:
              fprintf(fp, fmts[j*2], ((float *)tmp->array)[i]);
              break;
            case GAL_DATA_TYPE_DOUBLE:
              fprintf(fp, fmts[j*2], ((double *)tmp->array)[i]);
              break;
            default:
              error(EXIT_FAILURE, 0, "type code %d not recognized for "
                    "tmp->type in `gal_txt_write'", tmp->type);
            }
          ++j;
        }
      fprintf(fp, "\n");
    }


  /* Clean up, close the input file and return. For the fmts[i*2] elements,
     the reason is that fmts[i*2+1] are literal strings, not allocated. So
     they don't need freeing.*/
  for(i=0;i<numcols;++i)
    {
      free(fmts[i*2]);
      free(fmts[i*2+1]);
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
