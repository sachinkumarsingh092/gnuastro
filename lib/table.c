/*********************************************************************
table -- Functions for I/O on tables.
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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/git.h>
#include <gnuastro/txt.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>









/************************************************************************/
/***************         Information about a table        ***************/
/************************************************************************/
/* Store the information of each column in a table (either as a text file
   or as a FITS table) into an array of data structures with 'numcols'
   structures (one data structure for each column). The number of rows is
   stored in 'numrows'. The type of the table (e.g., ascii text file, or
   FITS binary or ASCII table) will be put in 'tableformat' (macros defined
   in 'gnuastro/table.h'.

   Note that other than the character strings (column name, units and
   comments), nothing in the data structure(s) will be allocated by this
   function for the actual data (e.g., the 'array' or 'dsize' elements). */
gal_data_t *
gal_table_info(char *filename, char *hdu, gal_list_str_t *lines,
               size_t *numcols, size_t *numrows, int *tableformat)
{
  /* Get the table format and size (number of columns and rows). */
  if(filename && gal_fits_name_is_fits(filename))
    return gal_fits_tab_info(filename, hdu, numcols, numrows, tableformat);
  else
    {
      *tableformat=GAL_TABLE_FORMAT_TXT;
      return gal_txt_table_info(filename, lines, numcols, numrows);
    }

  /* Abort with an error if we get to this point. */
  error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s so we can fix "
        "the problem. Control must not have reached the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NULL;
}





void
gal_table_print_info(gal_data_t *allcols, size_t numcols, size_t numrows)
{
  size_t i;
  int Nw=3, nw=4, uw=5, tw=4;   /* Initial width from label's width */
  char *name, *unit, *comment;

  /* Set the widths to print the column information. The width for the
     column number can easily be identified from the logarithm of the
     number of columns. */
  Nw=log10(numcols)+1;
  for(i=0;i<numcols;++i)
    {
      if(allcols[i].name && strlen(allcols[i].name)>nw)
        nw=strlen(allcols[i].name);
      if(allcols[i].unit && strlen(allcols[i].unit)>uw)
        uw=strlen(allcols[i].unit);
      if(allcols[i].type
         && strlen(gal_type_name(allcols[i].type, 1))>tw)
        tw=strlen(gal_type_name(allcols[i].type, 1));
    }

  /* We want one column space between the columns for readability, not the
     exact length, so increment all the numbers. */
  Nw+=2; nw+=2; uw+=2; tw+=2;

  /* Print these column names. */
  printf("%-*s%-*s%-*s%-*s%s\n", Nw, "---", nw, "----", uw,
         "-----", tw, "----", "-------");
  printf("%-*s%-*s%-*s%-*s%s\n", Nw, "No.", nw, "Name", uw,
         "Units", tw, "Type", "Comment");
  printf("%-*s%-*s%-*s%-*s%s\n", Nw, "---", nw, "----", uw,
         "-----", tw, "----", "-------");

  /* For each column, print the information, then free them. */
  for(i=0;i<numcols;++i)
    {
      name    = allcols[i].name;       /* Just defined for easier     */
      unit    = allcols[i].unit;       /* readability. The compiler   */
      comment = allcols[i].comment;    /* optimizer will remove them. */
      printf("%-*zu%-*s%-*s%-*s%s\n", Nw, i+1,
             nw, name ? name : GAL_BLANK_STRING ,
             uw, unit ? unit : GAL_BLANK_STRING ,
             tw, gal_type_name(allcols[i].type, 1),
             comment ? comment : GAL_BLANK_STRING);
    }

  /* Print the number of rows. */
  printf("--------\nNumber of rows: %zu\n--------\n", numrows);
}




















/************************************************************************/
/***************               Read a table               ***************/
/************************************************************************/

/* Function to print regular expression error. This is taken from the GNU C
   library manual, with small modifications to fit out style, */
static void
table_regexerrorexit(int errcode, regex_t *compiled, char *input)
{
  char *regexerrbuf;
  size_t length = regerror (errcode, compiled, NULL, 0);

  errno=0;
  regexerrbuf=malloc(length);
  if(regexerrbuf==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for regexerrbuf",
          __func__, length);
  (void) regerror(errcode, compiled, regexerrbuf, length);

  error(EXIT_FAILURE, 0, "%s: regular expression error: %s in value to "
        "'--column' ('-c'): '%s'", __func__, regexerrbuf, input);
}





/* Macro to set the string to search in */
static char *
table_set_strcheck(gal_data_t *col, int searchin)
{
  switch(searchin)
    {
    case GAL_TABLE_SEARCH_NAME:
      return col->name;

    case GAL_TABLE_SEARCH_UNIT:
      return col->unit;

    case GAL_TABLE_SEARCH_COMMENT:
      return col->comment;

    default:
      error(EXIT_FAILURE, 0, "%s: the code %d to searchin was not "
            "recognized", __func__, searchin);
    }

  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can "
        "address the problem. Control must not have reached the end of "
        "this function", __func__, PACKAGE_BUGREPORT);
  return NULL;
}





gal_list_sizet_t *
gal_table_list_of_indexs(gal_list_str_t *cols, gal_data_t *allcols,
                         size_t numcols, int searchin, int ignorecase,
                         char *filename, char *hdu, size_t *colmatch)
{
  long tlong;
  int regreturn;
  regex_t *regex;
  gal_list_str_t *tmp;
  gal_list_sizet_t *indexll=NULL;
  size_t i, nummatch, colcount=0, len;
  char *str, *strcheck, *tailptr, *errorstring;

  /* Go over the given columns.  */
  if(cols)
    for(tmp=cols; tmp!=NULL; tmp=tmp->next)
      {
        /* Counter for number of columns matched, and length of name. */
        nummatch=0;
        len=strlen(tmp->v);

        /* REGULAR EXPRESSION: the first and last characters are '/'. */
        if( tmp->v[0]=='/' && tmp->v[len-1]=='/' )
          {
            /* Remove the slashes, note that we don't want to change
               'tmp->v' (because it should be freed later). So first we set
               the last character to '\0', then define a new string from
               the first element. */
            tmp->v[len-1]='\0';
            str = tmp->v + 1;

            /* Allocate the regex_t structure: */
            errno=0;
            regex=malloc(sizeof *regex);
            if(regex==NULL)
              error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for regex",
                    __func__, sizeof *regex);

            /* First we have to "compile" the string into the regular
               expression, see the "POSIX Regular Expression Compilation"
               section of the GNU C Library.

               About the case of the string: the FITS standard says: "It is
               _strongly recommended_ that every field of the table be
               assigned a unique, case insensitive name with this
               keyword..."  So the column names can be case-sensitive.

               Here, we don't care about the details of a match, the only
               important thing is a match, so we are using the REG_NOSUB
               flag.*/
            regreturn=0;
            regreturn=regcomp(regex, str, ( ignorecase
                                            ? RE_SYNTAX_AWK | REG_ICASE
                                            : RE_SYNTAX_AWK ) );
            if(regreturn)
              table_regexerrorexit(regreturn, regex, str);


            /* With the regex structure "compile"d you can go through all
               the column names. Just note that column names are not
               mandatory in the FITS standard, so some (or all) columns
               might not have names, if so 'p->tname[i]' will be NULL. */
            for(i=0;i<numcols;++i)
              {
                strcheck=table_set_strcheck(&allcols[i], searchin);
                if(strcheck && regexec(regex, strcheck, 0, 0, 0)==0)
                  {
                    ++nummatch;
                    gal_list_sizet_add(&indexll, i);
                  }
              }

            /* Free the regex_t structure: */
            regfree(regex);

            /* Put the '/' back into the input string. This is done because
               after this function, the calling program might want to
               inform the user of their exact input string. */
            tmp->v[len-1]='/';
          }


        /* Not regular expression. */
        else
          {
            tlong=strtol(tmp->v, &tailptr, 0);

            /* INTEGER: If the string is an integer, then tailptr should
               point to the null character. If it points to anything else,
               it shows that we are not dealing with an integer (usable as
               a column number). So floating point values are also not
               acceptable. Since it is possible for the users to give zero
               for the column number, we need to read the string as a
               number first, then check it here. */
            if(*tailptr=='\0')
              {
                /* Make sure the number is larger than zero! */
                if(tlong<=0)
                  error(EXIT_FAILURE, 0, "%s: column numbers must be "
                        "positive (not zero or negative). You have asked "
                        "for column number %ld", __func__, tlong);

                /* Check if the given value is not larger than the number
                   of columns in the input catalog (note that the user is
                   counting from 1, not 0!) */
                if(tlong>numcols)
                  error(EXIT_FAILURE, 0, "%s: has %zu columns, but you "
                        "have asked for column number %ld",
                        gal_fits_name_save_as_string(filename, hdu),
                        numcols, tlong);

                /* Everything seems to be fine, put this column number in
                   the output column numbers linked list. Note that
                   internally, the column numbers start from 0, not 1.*/
                gal_list_sizet_add(&indexll, tlong-1);
                ++nummatch;
              }



            /* EXACT MATCH: */
            else
              {
                /* Go through all the desired column information and add
                   the column number when there is a match. */
                for(i=0;i<numcols;++i)
                  {
                    /* Check if this column actually has any
                       information. Then do a case-sensitive or insensitive
                       comparison of the strings. */
                    strcheck=table_set_strcheck(&allcols[i], searchin);
                    if(strcheck && ( ignorecase
                                     ? !strcasecmp(tmp->v, strcheck)
                                     : !strcmp(tmp->v, strcheck) ) )
                      {
                        ++nummatch;
                        gal_list_sizet_add(&indexll, i);
                      }
                  }
              }
          }


        /* If there was no match, then report an error. This can only happen
           for string matches, not column numbers, for numbers, the checks
           are done (and program is aborted) before this step. */
        if(nummatch==0)
          {
            if( asprintf(&errorstring, "'%s' didn't match any of the "
                         "column %ss.", tmp->v,
                         gal_tableintern_searchin_as_string(searchin))<0 )
              error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
            gal_tableintern_error_col_selection(filename, hdu, errorstring);
          }


        /* Keep the value of 'nummatch' if the user requested it. */
        if(colmatch) colmatch[colcount++]=nummatch;
      }

  /* cols==NULL */
  else
    for(i=0;i<numcols;++i)
      gal_list_sizet_add(&indexll, i);

  /* Reverse the list. */
  gal_list_sizet_reverse(&indexll);

  /* For a check.
  gal_list_sizet_print(indexll);
  exit(0);
  */

  /* Return the list. */
  return indexll;
}





/* Read the specified columns in a table (named 'filename') into a linked
   list of data structures. If the file is FITS, then 'hdu' will also be
   used, otherwise, 'hdu' is ignored. The information to search for columns
   should be specified by the 'cols' linked list as string values in each
   node of the list, the strings in each node can be a number, an exact
   match to a column name, or a regular expression (in GNU AWK format)
   enclosed in '/ /'. The 'searchin' value comes from the
   'gal_table_where_to_search' enumerator and has to be one of its given
   types. If 'cols' is NULL, then this function will read the full table.

   The output is a linked list with the same order of the cols linked
   list. Note that one column node in the 'cols' list might give multiple
   columns, in this case, the order of output columns that correspond to
   that one input, are in order of the table (which column was read first).
   So the first requested column is the first popped data structure and so
   on. */
gal_data_t *
gal_table_read(char *filename, char *hdu, gal_list_str_t *lines,
               gal_list_str_t *cols, int searchin, int ignorecase,
               size_t minmapsize, int quietmmap, size_t *colmatch)
{
  int tableformat;
  gal_list_sizet_t *indexll;
  size_t i, numcols, numrows;
  gal_data_t *allcols, *out=NULL;

  /* First get the information of all the columns. */
  allcols=gal_table_info(filename, hdu, lines, &numcols, &numrows,
                         &tableformat);

  /* If there was no actual data in the file, then return NULL. */
  if(allcols==NULL) return NULL;

  /* Get the list of indexs in the same order as the input list. */
  indexll=gal_table_list_of_indexs(cols, allcols, numcols, searchin,
                                   ignorecase, filename, hdu, colmatch);

  /* Depending on the table format, read the columns into the output
     structure. Note that the functions here pop each index, read/store the
     desired column and pop the next, so after these functions, the output
     linked list will have the opposite order of its input 'indexll'
     list. So before calling any of them, we will first reverse the
     'indexll' list, so the output data structure list will have the same
     order as the input list of desired columns. Also note that after these
     functions, the 'indexll' will be all freed (each popped element is
     actually freed).*/
  gal_list_sizet_reverse(&indexll);
  switch(tableformat)
    {
    case GAL_TABLE_FORMAT_TXT:
      out=gal_txt_table_read(filename, lines, numrows, allcols, indexll,
                             minmapsize, quietmmap);
      break;

    case GAL_TABLE_FORMAT_AFITS:
    case GAL_TABLE_FORMAT_BFITS:
      out=gal_fits_tab_read(filename, hdu, numrows, allcols, indexll,
                            minmapsize, quietmmap);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: table format code %d not recognized for "
            "'tableformat'", __func__, tableformat);
    }

  /* Clean up. */
  for(i=0;i<numcols;++i)
    gal_data_free_contents(&allcols[i]);
  free(allcols);
  gal_list_sizet_free(indexll);

  /* Return the final linked list. */
  return out;
}




















/************************************************************************/
/***************              Write a table               ***************/
/************************************************************************/
/* Write the basic information that is necessary by each program into the
   comments field. Note that the 'comments' has to be already sorted in the
   proper order. */
void
gal_table_comments_add_intro(gal_list_str_t **comments, char *program_string,
                             time_t *rawtime)
{
  char gitdescribe[100], *tmp;

  /* Get the Git description in the running folder. */
  tmp=gal_git_describe();
  if(tmp) { sprintf(gitdescribe, " from %s,", tmp); free(tmp); }
  else      gitdescribe[0]='\0';

  /* Git version and time of program's starting, this will be the second
     line. Note that ctime puts a '\n' at the end of its string, so we'll
     have to remove that. Also, note that since we are allocating 'msg', we
     are setting the allocate flag of 'gal_list_str_add' to 0. */
  if( asprintf(&tmp, "Created%s on %s", gitdescribe, ctime(rawtime))<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  tmp[ strlen(tmp)-1 ]='\0';
  gal_list_str_add(comments, tmp, 0);

  /* Program name: this will be the top of the list (first line). We will
     need to set the allocation flag for this one, because program_string
     is usually statically allocated.*/
  if(program_string)
    gal_list_str_add(comments, program_string, 1);
}





/* The input is a linked list of data structures and some comments. The
   table will then be written into 'filename' with a format that is
   specified by 'tableformat'. */
void
gal_table_write(gal_data_t *cols, gal_list_str_t *comments,
                int tableformat, char *filename, char *extname,
                uint8_t colinfoinstdout)
{
  /* If a filename was given, then the tableformat is relevant and must be
     used. When the filename is empty, a text table must be printed on the
     standard output (on the command-line). */
  if(filename)
    {
      if(gal_fits_name_is_fits(filename))
        gal_fits_tab_write(cols, comments, tableformat, filename, extname);
      else
        gal_txt_write(cols, comments, filename, colinfoinstdout);
    }
  else
    /* Write to standard output. */
    gal_txt_write(cols, comments, filename, colinfoinstdout);
}





void
gal_table_write_log(gal_data_t *logll, char *program_string,
                    time_t *rawtime, gal_list_str_t *comments,
                    char *filename, int quiet)
{
  char *msg;

  /* Write all the comments into */
  gal_table_comments_add_intro(&comments, program_string, rawtime);

  /* Write the log file to disk */
  gal_table_write(logll, comments, GAL_TABLE_FORMAT_TXT, filename, "LOG", 0);

  /* In verbose mode, print the information. */
  if(!quiet)
    {
      if( asprintf(&msg, "%s created.", filename)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(NULL, msg, 1);
      free(msg);
    }
}
