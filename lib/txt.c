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
  TXT_LINESTAT_NOTCOMMENT,
};





/* Return one of the `txt_line_stat' constant values. */
int
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
          return TXT_LINESTAT_NOTCOMMENT;
        }
      ++line;
    }
  return TXT_LINESTAT_BLANK;
}





/* Return the information about a text file table. */
gal_data_t *
gal_txt_table_info(char *filename, size_t *numcols)
{
  FILE *fp;
  char *line;
  size_t linelen=10;  /* This will be increased later by `getline'. */

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
    switch(get_line_stat(line))
      {
      case TXT_LINESTAT_BLANK:
        printf("blank\n");
        break;

      case TXT_LINESTAT_ISCOMMENT:
        printf("comment\n");
        break;

      case TXT_LINESTAT_NOTCOMMENT:
        printf("not comment\n");
        break;

      default:
        error(EXIT_FAILURE, 0, "linestatus code %d not recognized in "
              "`gal_txt_table_info'",
              linestat);
      }

  /* Clean up, close the file and return. */
  free(line);
  errno=0;
  if(fclose(fp))
    error(EXIT_FAILURE, errno, "%s: couldn't close file after reading ASCII "
          "table information", filename);

  printf("\n----end of tableinfo----\n");
  exit(0);
  return NULL;
}




















/************************************************************************/
/***************          Write a txt table               ***************/
/************************************************************************/
static char **
make_fmts_for_printf(gal_data_t *cols, size_t numcols, int leftadjust,
                     size_t *len)
{
  size_t i=0;
  char **fmts;
  gal_data_t *tmp;
  char *fmt=NULL, *lng;
  int width=0, precision=0;

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
      if(fmts[i*2]==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for fmts[%zu] in "
              "`make_fmts_for_printf' of txt.c",
              GAL_TXT_MAX_FMT_LENGTH*sizeof *fmts[i], i);

      /* Write the proper format. */
      switch(tmp->type)
        {

        case GAL_DATA_TYPE_BIT:
          error(EXIT_FAILURE, 0, "printing of bit types is currently "
                "not supported");
          break;

        case GAL_DATA_TYPE_STRING:
          fmt="s";
          width=( tmp->disp_width<0 ? GAL_TABLE_DEF_STR_WIDTH
                  : tmp->disp_width );
          precision=( tmp->disp_precision<0 ? GAL_TABLE_DEF_STR_PRECISION
                      : tmp->disp_precision );
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

      /* Print the result into the allocated string. The space in the end
         here is to ensure that if the printed string is larger than the
         expected width, it the columns will not merger and at least one
         space character will be put between them.*/
      *len += 1 + sprintf(fmts[i*2], "%%%s%d.%d%s%s ", leftadjust ? "-" : "",
                          width, precision, lng, fmt);
      fmts[i*2+1] = gal_data_type_string(tmp->type, 0);

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
  size_t i, j, numcols=0, fmtlen;
  int iw=0, fw=0, nw=0, uw=0, tw=0;


  /* Find the number of columns, do a small sanity check, and get the
     maximum width of the name and unit string if they are present. */
  for(tmp=cols;tmp!=NULL;tmp=tmp->next)
    {
      ++numcols;
      if(cols->size!=tmp->size)
        error(EXIT_FAILURE, 0, "to print a set of columns, into a file, they "
              "must all have the same number of elements/rows. The inputs to "
              "`gal_txt_write' have different sizes: the first column has "
              "%zu, while column %zu as %zu elements", cols->size, numcols,
              tmp->size);
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
    {
      fw = strlen( fmts[ i*2   ] ) > fw ? strlen( fmts[ i*2   ] ) : fw;
      tw = strlen( fmts[ i*2+1 ] ) > tw ? strlen( fmts[ i*2+1 ] ) : tw;
    }


  /* Write the comments if there as any */
  if(comment) fprintf(fp, "%s\n", comment);


  /* Write the given information. */
  i=0;
  iw=log10(numcols)+1;
  for(tmp=cols;tmp!=NULL;tmp=tmp->next)
    {
      fprintf(fp, "# Column %-*zu %-*s [ %-*s , %-*s , %-*s] %s\n",
              iw, i+1,
              nw, tmp->name    ? tmp->name    : "",
              uw, tmp->unit    ? tmp->unit    : "",
              tw, fmts[i*2+1],
              fw, fmts[i*2],
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
  for(i=0;i<numcols;++i) free(fmts[i*2]);
  free(fmts);
  if(filename)
    {
      errno=0;
      if(fclose(fp))
        error(EXIT_FAILURE, errno, "%s: couldn't close file after writing "
              "of text table", filename);
    }
}
