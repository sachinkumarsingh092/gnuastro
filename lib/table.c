/*********************************************************************
table -- Functions for I/O on tables.
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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/txt.h>
#include <gnuastro/table.h>










/************************************************************************/
/***************         Information about a table        ***************/
/************************************************************************/
/* Store the information of each column in a table (either as a text file
   or as a FITS table) into an array of data structures with `numcols'
   structures (one data structure for each column). The number of rows is
   stored as the `size' element of each data structure. The type of the
   table (e.g., ascii text file, or FITS binary or ASCII table) will be put
   in `tabletype' (macros defined in `gnuastro/table.h'.

   Note that other than the character strings (column name, units and
   comments), nothing in the data structure(s) will be allocated by this
   function for the actual data (e.g., the `array' or `dsize' elements). */
gal_data_t *
gal_table_info(char *filename, char *hdu, size_t *numcols, size_t *numrows,
               int *tabletype)
{
  /* Get the table type and size (number of columns and rows). */
  if(gal_fits_name_is_fits(filename))
    return gal_fits_table_info(filename, hdu, numcols, numrows, tabletype);
  else
    {
      *tabletype=GAL_TABLE_TYPE_TXT;
      return gal_txt_table_info(filename, numcols, numrows);
    }

  /* Abort with an error if we get to this point. */
  error(EXIT_FAILURE, 0, "A bug! please contact us at %s so we can fix "
        "the problem. For some reason, control has reached the end of "
        "`gal_table_info'", PACKAGE_BUGREPORT);
  return NULL;
}




















/************************************************************************/
/***************               Read a table               ***************/
/************************************************************************/
/* In programs, the `searchin' variable is much more easier to type in as a
   description than an integer (which is what `gal_table_read_cols'
   needs). This function will check the string value and give the
   corresponding integer value.*/
int
gal_table_searchin_from_str(char *searchin_str)
{
  if(strcmp(searchin_str, "name")==0)
    return GAL_TABLE_SEARCH_NAME;
  else if(strcmp(searchin_str, "unit")==0)
    return GAL_TABLE_SEARCH_UNIT;
  else if(strcmp(searchin_str, "comment")==0)
    return GAL_TABLE_SEARCH_COMMENT;
  else
    error(EXIT_FAILURE, 0, "`--searchin' only recognizes the values "
          "`name', `unit', and `comment', you have asked for `%s'",
          searchin_str);

  /* Report an error control reaches here. */
  error(EXIT_FAILURE, 0, "A bug! please contact us at %s so we can address "
        "the problem. For some reason control has reached the end of "
        "`gal_table_searchin_from_str'", PACKAGE_BUGREPORT);
  return -1;
}





/* Function to print regular expression error. This is taken from the GNU C
   library manual, with small modifications to fit out style, */
void
regexerrorexit(int errcode, regex_t *compiled, char *input)
{
  char *regexerrbuf;
  size_t length = regerror (errcode, compiled, NULL, 0);

  errno=0;
  regexerrbuf=malloc(length);
  if(regexerrbuf==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for regexerrbuf", length);
  (void) regerror(errcode, compiled, regexerrbuf, length);

  error(EXIT_FAILURE, 0, "Regular expression error: %s in value to "
        "`--column' (`-c'): `%s'", regexerrbuf, input);
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
      error(EXIT_FAILURE, 0, "the code %d to searchin was not "
            "recognized in `table_set_strcheck'", searchin);
    }

  error(EXIT_FAILURE, 0, "A bug! Please contact us at %s, For some reason "
        "control has reached the end of `table_set_strcheck'. This must "
        "not have happened", PACKAGE_BUGREPORT);
  return NULL;
}





static struct gal_linkedlist_sll *
make_list_of_indexs(struct gal_linkedlist_stll *cols, gal_data_t *allcols,
                    size_t numcols, int searchin, int ignorecase,
                    char *filename, char *hdu)
{
  long tlong;
  int regreturn;
  regex_t *regex;
  size_t i, nummatch;
  struct gal_linkedlist_stll *tmp;
  char *str, *colts, *strcheck, *tailptr;
  struct gal_linkedlist_sll *indexll=NULL;

  for(tmp=cols; tmp!=NULL; tmp=tmp->next)
    {
      /* Counter for number of columns matched. */
      nummatch=0;

      /* REGULAR EXPRESSION: When the first and last characters are `/'. */
      if( tmp->v[0]=='/' && tmp->v[strlen(tmp->v)-1]=='/' )
        {
          /* Remove the slashes, note that we don't want to change `tmp->v'
             (because it should be freed later). So first we set the last
             character to `\0', then define a new string from the first
             element. */
          tmp->v[strlen(tmp->v)-1]='\0';
          str = tmp->v + 1;

          /* Allocate the regex_t structure: */
          errno=0; regex=malloc(sizeof *regex);
          if(regex==NULL)
            error(EXIT_FAILURE, errno, "%zu bytes for regex", sizeof *regex);

          /* First we have to "compile" the string into the regular
             expression, see the "POSIX Regular Expression Compilation"
             section of the GNU C Library.

             About the case of the string: the FITS standard says: "It is
             _strongly recommended_ that every field of the table be
             assigned a unique, case insensitive name with this keyword..."
             So the column names can be case-sensitive.

             Here, we don't care about the details of a match, the only
             important thing is a match, so we are using the REG_NOSUB
             flag.*/
          regreturn=0;
          regreturn=regcomp(regex, str, ( ignorecase
                                          ? RE_SYNTAX_AWK | REG_ICASE
                                          : RE_SYNTAX_AWK ) );
          if(regreturn)
            regexerrorexit(regreturn, regex, str);


          /* With the regex structure "compile"d you can go through all the
             column names. Just note that column names are not mandatory in
             the FITS standard, so some (or all) columns might not have
             names, if so `p->tname[i]' will be NULL. */
          for(i=0;i<numcols;++i)
            {
              strcheck=table_set_strcheck(&allcols[i], searchin);
              if(strcheck && regexec(regex, strcheck, 0, 0, 0)==0)
                {
                  ++nummatch;
                  gal_linkedlist_add_to_sll(&indexll, i);
                }
            }

          /* Free the regex_t structure: */
          regfree(regex);
        }


      /* Not regular expression. */
      else
        {
          tlong=strtol(tmp->v, &tailptr, 0);

          /* INTEGER: If the string is an integer, then tailptr should
             point to the null character. If it points to anything else, it
             shows that we are not dealing with an integer (usable as a
             column number). So floating point values are also not
             acceptable. Since it is possible for the users to give zero
             for the column number, we need to read the string as a number
             first, then check it here. */
          if(*tailptr=='\0')
            {
              /* Make sure the number is larger than zero! */
              if(tlong<=0)
                error(EXIT_FAILURE, 0, "the column numbers given to the "
                      "must not be zero, or negative. You have asked for "
                      "column %ld", tlong);

              /* Check if the given value is not larger than the number of
                 columns in the input catalog (note that the user is
                 counting from 1, not 0!) */
              if(tlong>numcols)
                {
                  if(gal_fits_name_is_fits(filename))
                    error(EXIT_FAILURE, 0, "%s (hdu %s): has %zu columns, "
                          "but you have asked for column number %zu",
                          filename, hdu, numcols, tlong);
                  else
                    error(EXIT_FAILURE, 0, "%s: has %zu columns, but you "
                          "have asked for column number %zu", filename,
                          numcols, tlong);
                }

              /* Everything seems to be fine, put this column number in the
                 output column numbers linked list. Note that internally,
                 the column numbers start from 0, not 1.*/
              gal_linkedlist_add_to_sll(&indexll, tlong-1);
              ++nummatch;
            }



          /* EXACT MATCH: */
          else
            {
              /* Go through all the desired column information and add the
                 column number when there is a match. */
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
                      gal_linkedlist_add_to_sll(&indexll, i);
                    }
                }
            }
        }


      /* If there was no match, then report an error. This can only happen
         for string matches, not column numbers, for numbers, the checks
         are done before the reading.*/
      if(nummatch==0)
        {
          colts = ( searchin==GAL_TABLE_SEARCH_NAME ? "names"
                    : ( searchin==GAL_TABLE_SEARCH_UNIT ? "units"
                        : "comments") );
          error(EXIT_FAILURE, 0, "`%s' didn't match with any of the column "
                "%s in `%s'%s%s%s. You can check the available column "
                "information by running `asttable %s%s%s --information'. "
                "%s", tmp->v, colts, filename,
                gal_fits_name_is_fits(filename) ? " (hdu: " : "",
                gal_fits_name_is_fits(filename) ? hdu : "",
                gal_fits_name_is_fits(filename) ? ")" : "",
                filename,
                gal_fits_name_is_fits(filename) ? " --hdu" : "",
                gal_fits_name_is_fits(filename) ? hdu : "",
                ignorecase ? "" : "For a case-insensitive match, "
                "run again with the `--ignorecase' option. " );
        }
    }

  /* Reverse the list. */
  gal_linkedlist_reverse_sll(&indexll);
  return indexll;
}





/* Read the specified columns in a text file (named `filename') into a
   linked list of data structures. If the file is FITS, then `hdu' will
   also be used, otherwise, `hdu' is ignored. The information to search for
   columns should be specified by the `cols' linked list as string values
   in each node of the list, the strings in each node can be a number, an
   exact match to a column name, or a regular expression (in GNU AWK
   format) enclosed in `/ /'. The `tosearch' value comes from the
   `gal_table_where_to_search' enumerator and has to be one of its given
   types. If `cols' is NULL, then this function will read the full table.

   The output is a linked list with the same order of the cols linked
   list. Note that one column node in the `cols' list might give multiple
   columns, in this case, the order of output columns that correspond to
   that one input, are in order of the table (which column was read first).
   So the first requested column is the first popped data structure and so
   on. */
gal_data_t *
gal_table_read(char *filename, char *hdu, struct gal_linkedlist_stll *cols,
               int searchin, int ignorecase, int minmapsize)
{
  int tabletype;
  size_t i, numcols, numrows;
  gal_data_t *allcols, *out=NULL;
  struct gal_linkedlist_sll *indexll;

  /* First get the information of all the columns. */
  allcols=gal_table_info(filename, hdu, &numcols, &numrows, &tabletype);

  /* Get the list of indexs in the same order as the input list */
  indexll=make_list_of_indexs(cols, allcols, numcols, searchin,
                              ignorecase, filename, hdu);

  /* If no columns could be selected, just return NULL. */
  if(indexll==NULL) return NULL;

  /* Depending on the table type, read the columns into the output
     structure. Note that the functions here pop each index, read/store the
     desired column and pop the next, so after these functions, the output
     linked list will have the opposite order of its input `indexll'
     list. So before calling any of them, we will first reverse the
     `indexll' list, so the output data structure list will have the same
     order as the input list of desired columns. Also note that after these
     functions, the `indexll' will be all freed (each popped element is
     actually freed).*/
  gal_linkedlist_reverse_sll(&indexll);
  switch(tabletype)
    {
    case GAL_TABLE_TYPE_TXT:
      out=gal_txt_table_read(filename, numrows, allcols, indexll,
                             minmapsize);
      break;

    case GAL_TABLE_TYPE_AFITS:
    case GAL_TABLE_TYPE_BFITS:
      out=gal_fits_table_read(filename, hdu, numrows, allcols, indexll,
                              minmapsize);
      break;

    default:
      error(EXIT_FAILURE, 0, "table type code %d not recognized for "
            "`tabletype' in `gal_table_read_cols'", tabletype);
    }

  /* Clean up. */
  for(i=0;i<numcols;++i)
    {
      allcols[i].wcs=NULL;
      allcols[i].mmapname=NULL;
      gal_data_free(&allcols[i], 1);
    }
  free(allcols);
  gal_linkedlist_free_sll(indexll);

  /* Return the final linked list. */
  return out;
}




















/************************************************************************/
/***************              Write a table               ***************/
/************************************************************************/
void
gal_table_write(gal_data_t *cols, char *comments, int tabletype,
                char *filename, int dontdelete)
{
  /* If a filename was given, then the tabletype is relevant and must be
     used. When the filename is empty, a text table must be printed on the
     standard output (on the command-line). */
  if(filename)
    switch(tabletype)
      {
      case GAL_TABLE_TYPE_TXT:
        gal_txt_write(cols, comments, filename, dontdelete);
        break;

      case GAL_TABLE_TYPE_AFITS:
      case GAL_TABLE_TYPE_BFITS:
        gal_fits_table_write(cols, comments, tabletype, filename, dontdelete);
        break;
      }
  else
    gal_txt_write(cols, comments, filename, dontdelete);
}
