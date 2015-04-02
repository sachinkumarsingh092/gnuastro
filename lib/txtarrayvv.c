/*********************************************************************
txtarrayvv -- Convert a text file table to a C array.
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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "txtarrayvv.h"











/**********************************************************************/
/*****************       Read a text table        *********************/
/**********************************************************************/
void
txttablesize(char *filename, size_t *outs0, size_t *outs1)
{
  FILE *fp;
  size_t len=200, s0, s1;
  char *line, *firsttoken;

  /* Initialize the sizes any way. */
  s0=s1=0;

  /* Allocate some space for `line` with `len` elements so it can
     easily be freed later on. The value of `len` is arbitarary at
     this point, during the run, getline will change it along with the
     pointer to line. */
  errno=0;
  line=malloc(len*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes in readdefaults",
	  len * sizeof *line);

  /* Open the file: */
  errno=0;
  fp=fopen(filename, "r");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s", filename);

  while( getline(&line, &len, fp) != -1 )
    {
      /* Read the first token: */
      firsttoken=strtok(line, DELIMITERS);

      /* If there are no non-delimters in the line (can happen if the
	 line is a blank line in the end of the file). */
      if(firsttoken==NULL)
	continue;

      /* Check if it is a comment or new line character: */
      if(firsttoken[0]=='#')
	continue;

      /* We are now in a data row: */
      if(s0==0)  /* We are on the first row of data, find s1. */
	{
	  s1=1;
	  while( strtok(NULL, DELIMITERS) != NULL )
	    ++s1;
	}
      ++s0;
    }
  free(line);

  errno=0;
  if(fclose(fp)==EOF)
    error(EXIT_FAILURE, errno, "%s", filename);

  if(s0==0 && s1==0)
    error(EXIT_FAILURE, 0, "No table could be read in %s.", filename);

  *outs0=s0;
  *outs1=s1;
}





void
savetolog(FILE **log, char *filename, size_t lineno, size_t s0,
	  size_t s1, char *token)
{
  if(*log)			/* The file is already open. */
    fprintf(*log, "%-10lu%-10lu%-10lu%s\n", lineno, s0, s1, token);
  else				/* Not yet created.          */
    {
      errno=0;
      *log=fopen(TXTARRAYVVLOG, "w");
      if(*log==NULL)
	error(EXIT_FAILURE, errno, "%s", filename);
      fprintf(*log, "# Elements in %s which could not be read as a \n"
	      "# number. They are saved as nan in the array.\n"
	      "# The columns in the table below are:\n"
	      "# 0: Line number in file.\n"
	      "# 1: Row number in table (without commented or blank "
              "lines).\n"
	      "# 2: Column number in table.\n"
	      "# 3: The string that could not be converted to a number.\n"
	      "# Note that counting starts from zero.\n",
	      filename);
      fprintf(*log, "%-10lu%-10lu%-10lu%s\n", lineno-1, s0, s1, token);
    }
}





/* Macro functions: */
#define CONVERTANDSAVE {						\
  errno=0; tailptr=NULL;						\
  array[ts0*s1+ts1]=strtod(token, &tailptr);				\
  if(errno)								\
    error_at_line(EXIT_FAILURE, errno, "%s", lineno, token, filename);  \
  if(*tailptr!='\0')							\
    {									\
      savetolog(&log, filename, lineno, ts0, ts1, token);		\
      array[ts0*s1+ts1]=NAN;						\
    }									\
  ++ts1;								\
  }





void
filltable(char *filename, double *array, size_t s0, size_t s1)
{
  FILE *fp, *log=NULL;
  char *line=NULL, *token, *tailptr;
  size_t len=200, lineno=0, ts0, ts1;

  /* Open the file: */
  errno=0;
  fp=fopen(filename, "r");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s", filename);

  /* Allocate some space for `line` with `len` elements so it can
     easily be freed later on. The value of `len` is arbitarary at
     this point, during the run, getline will change it along with the
     pointer to line. */
  errno=0;
  line=malloc(len*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes in readdefaults",
	  len * sizeof *line);

  ts0=0;
  while( getline(&line, &len, fp) != -1 )
    {
      ts1=0;
      ++lineno;

      /* Read the first token: */
      token=strtok(line, DELIMITERS);

      /* If there are no non-delimters in the line (can happen if the
	 line is a blank line in the end of the file). */
      if(token==NULL)
	continue;

      /* Check if it is a comment or new line character: */
      if(token[0]=='#')
	continue;

      /* Convert the first token and put it into the array: */
      CONVERTANDSAVE;

      /* Read the rest of the tokens: */
      while( (token=strtok(NULL, DELIMITERS))!=NULL )
	{
	  if(ts1>=s1)
	    error_at_line(EXIT_FAILURE, 0, filename, lineno,
			  "Too many columns on this line. The number of "
			  "columns should be the same as the first row "
			  "of the table.");
	  CONVERTANDSAVE;
	}
      if(ts1<s1-1)		/* It should be s1-1. */
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "Not enough columns on this line. The number of "
		      "columns should be the same as the first row "
		      "of the table.");
      ++ts0;
    }
  free(line);

  errno=0;
  if(fclose(fp)==EOF)
    error(EXIT_FAILURE, errno, "%s", filename);
  if(log)
    {
      errno=0;
      if(fclose(log)==EOF)
	error(EXIT_FAILURE, errno, "%s", TXTARRAYVVLOG);
    }
}





void
txttoarray(char *filename, double **array, size_t *s0, size_t *s1)
{
  /* Find the size of the table and allocate space for it: */
  errno=0;
  txttablesize(filename, s0, s1);
  if( (*array=malloc(*s0 * *s1 * sizeof **array)) == NULL)
    error(EXIT_FAILURE, errno, "txttoarray: space for array with %lu "
	  "elements", *s0 * *s1);

  /* Fill in the table with the contents of the text file: */
  filltable(filename, *array, *s0, *s1);
}




















/**********************************************************************/
/*****************       Write a text table        ********************/
/**********************************************************************/
/* This function gets the formatting settings of the array as required
   by writeasciitable and makes an array of formatting conditions that
   is suitable for printing.  */
void
doformatting(int numcols, char **fmt, int *int_cols, int *accu_cols,
	     int *space, int *prec, char forg)
{
  int i,j, found=0;

  /* Initialize the format array: */
  for (i=0;i<numcols;++i)
    {
      /* Allocate space for the format string of each column: */
      errno=0;
      fmt[i]=malloc(FMTLENGTH * sizeof(char));
      if(fmt[i]==NULL)
	error(EXIT_FAILURE, errno, "txtarrayvv, space for format "
	      "string %d, with %d elements", i, FMTLENGTH);

      /* See if this is an int column. */
      found=0;
      for(j=0;j<numcols;++j)
        {
	  if (int_cols[j]<0) break;
	  if (i==int_cols[j])
            {
	      sprintf(fmt[i], "%%-%d.0%c", space[0], forg);
	      found=1;break;
            }
        }
      if (found==1) continue;

      /* See if this is an extra precision column. */
      found=0;
      for(j=0;j<numcols;++j)
        {
	  if (accu_cols[j]<0) break;
	  if (i==accu_cols[j])
            {
	      sprintf(fmt[i], "%%-%d.%d%c", space[2], prec[1], forg);
	      found=1;break;
            }
        }
      if (found==1) continue;

      /* It is neither of the above, so it is a normal precision column. */
      sprintf(fmt[i], "%%-%d.%d%c", space[1], prec[0], forg);
    }
}






/* Write an array to a text file. It is assumed that we have three
   types of input: Those who don't have any decimal point (ints),
   those that need a small number of decimal accuracy and finally,
   those that need to be very accurate. You specify the number of
   digits after the decimal point along with all the other settings
   with these input arrays:

   int_cols: An array of integers showing the columns that are to be
     printed with no decimal point. It should end with a negative
     value to mark the end of the column numbers.

   accu_cols: Similar to int_cols, but for columns which need extra
     accuracy.

   space: A three element array, which shows how much space should be
     given to the three different types (the value immediately after %
     in the format string).

   prec: A two element array, showing the number of decimal points to
     print the less and more accurate columns.

 */
void
arraytotxt(double *array, size_t s0, size_t s1, char *comments,
	   int *int_cols, int *accu_cols, int *space, int *prec,
	   char forg, const char *filename)
{
  int i,j;
  FILE *fp;
  char **fmt;

  /* Do a small sanity check: */
  for(i=0;int_cols[i]>0;++i)
    if(int_cols[i]>=s1)
      error(EXIT_FAILURE, 0, "arraytotxt: In int_cols[], %d is "
	    "larger than the number of columns: %lu.", int_cols[i], s1);
  for(i=0;accu_cols[i]>0;++i)
    if(accu_cols[i]>=s1)
      error(EXIT_FAILURE, 0, "arraytotxt: In accu_cols[], %d is "
	    "larger than the number of columns: %lu.", accu_cols[i], s1);
  for(i=0;i<3;++i)
    if(space[i]<=0)
      error(EXIT_FAILURE, 0, "arraytotxt: In space[], %d is "
	    "smaller or equal to zero.", space[i]);
  for(i=0;i<2;++i)
    if(prec[i]<0)
      error(EXIT_FAILURE, 0, "arraytotxt: In prec[], %d is "
	    "smaller than zero.", space[i]);

  /* Allocate the spaces: */
  errno=0;
  fmt=malloc(s1 * sizeof(char *));
  if(fmt==NULL)
    error(EXIT_FAILURE, errno, "txtarrayvv, formatting of each "
	  "column with %lu elements", s1);

  /* Prepare the formatting for each column */
  doformatting(s1, fmt, int_cols, accu_cols, space, prec, forg);

  /* Open the output file: */
  errno=0;
  fp=fopen(filename, "w");
  if (fp==NULL)
    error(EXIT_FAILURE, errno, "%s", filename);

  /* Print the headers to file: */
  fprintf(fp, "%s\n", comments);

  /* Print the data to file: */
  for(i=0;i<s0;++i)
    {
      for(j=0;j<s1;++j)
	fprintf(fp, fmt[j], array[i*s1+j]);
      fprintf(fp, "\n");
    }

  /* Close the file and free all pointers: */
  errno=0;
  if( fclose(fp) == EOF )
    error(EXIT_FAILURE, errno, "%s", filename);
  for(i=0;i<s1;++i) free(fmt[i]);
  free(fmt);
}
