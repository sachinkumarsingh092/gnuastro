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

#include <gnuastro/git.h>
#include <gnuastro/txt.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>

#include <timing.h>
#include <checkset.h>






/************************************************************************/
/***************              Error messages              ***************/
/************************************************************************/
void
gal_table_error_col_selection(char *filename, char *hdu, char *errorstring)
{
  char *c, *name, *command;

  /* Set the proper pointers. */
  if(gal_fits_name_is_fits(filename))
    {
      asprintf(&name, "%s (hdu: %s)", filename, hdu);
      c=hdu; while(*c!='\0') if(isspace(*c++)) break;
      asprintf(&command, *c=='\0' ? "%s --hdu=%s" : "%s --hdu=\"%s\"",
               filename, hdu);
    }
  else command=name=filename;

  /* Abort with with the proper error. */
  error(EXIT_FAILURE, 0, "%s: %s\n\nFor more information on selecting "
        "columns in Gnuastro, please run the following command (press "
        "`SPACE' to go down and `q' to return to the command-line):\n\n"
        "    $ info gnuastro \"Selecting table columns\"\n\n"
        "To define a better column selection criteria, you can see "
        "the list of column meta-data in this table, with the following "
        "command:\n\n"
        "    $ asttable %s --info\n", name, errorstring, command);
}




















/************************************************************************/
/***************                 Formats                  ***************/
/************************************************************************/
/* Return the type of desired table based on a standard string. */
uint8_t
gal_table_string_to_format(char *string)
{
  if(string)                    /* Its not NULL. */
    {
      if(!strcmp(string, "txt"))              return GAL_TABLE_FORMAT_TXT;
      else if(!strcmp(string,"fits-ascii"))   return GAL_TABLE_FORMAT_AFITS;
      else if(!strcmp(string, "fits-binary")) return GAL_TABLE_FORMAT_BFITS;
      else                                    return GAL_TABLE_FORMAT_INVALID;
    }
  else                                        return GAL_TABLE_FORMAT_INVALID;
}





char *
gal_table_format_as_string(uint8_t format)
{
  switch(format)
    {
    case GAL_TABLE_FORMAT_TXT:    return "txt";
    case GAL_TABLE_FORMAT_AFITS:  return "fits-ascii";
    case GAL_TABLE_FORMAT_BFITS:  return "fits-binary";
    default:
      error(EXIT_FAILURE, 0, "code %d not recognized as a valid format "
            "in `gal_table_format_as_string'", format);
      return NULL;
    }
}






/* In programs, the `searchin' variable is much more easier to format in as
   a description than an integer (which is what `gal_table_read_cols'
   needs). This function will check the string value and give the
   corresponding integer value.*/
uint8_t
gal_table_string_to_searchin(char *string)
{
  if(string)                    /* Its not NULL. */
    {
      if(!strcmp(string, "name"))          return GAL_TABLE_SEARCH_NAME;
      else if(!strcmp(string, "unit"))     return GAL_TABLE_SEARCH_UNIT;
      else if(!strcmp(string, "comment"))  return GAL_TABLE_SEARCH_COMMENT;
      else                                 return GAL_TABLE_SEARCH_INVALID;
    }
  else                                     return GAL_TABLE_SEARCH_INVALID;
}





char *
gal_table_searchin_as_string(uint8_t searchin)
{
  switch(searchin)
    {
    case GAL_TABLE_SEARCH_NAME:    return "name";
    case GAL_TABLE_SEARCH_UNIT:    return "unit";
    case GAL_TABLE_SEARCH_COMMENT: return "comment";
    default:
      error(EXIT_FAILURE, 0, "code %d not recognized as a valid search "
            "field in `gal_table_searchin_as_string'", searchin);
      return NULL;
    }
}





/* For programs that output tables, the `--tableformat' option will be used
   to specify what format the output table should be in. When the output
   file is a FITS file, there are multiple formats, so to simplify the
   coding in each program, this function will do a sanity check on the
   value given to the `--tableformat' parameter. */
void
gal_table_check_fits_format(char *filename, int tableformat)
{
  if( filename && gal_fits_name_is_fits(filename) )
    {
      /* When `--tableformat' was not given. */
      if(tableformat==GAL_TABLE_FORMAT_INVALID)
        error(EXIT_FAILURE, 0, "`%s' (output file) is a FITS file but the "
              "desired format of the FITS table has not been specified with "
              "the `--tableformat' option. For FITS tables, this option can "
              "take two values: `fits-ascii', or `fits-binary'", filename);

      /* When `--tableformat' didn't have the correct value. */
      if( tableformat != GAL_TABLE_FORMAT_AFITS
          && tableformat != GAL_TABLE_FORMAT_BFITS )
        error(EXIT_FAILURE, 0, "`%s' (output file) is a FITS file but "
              "is not a recognized FITS table format. For FITS tables, "
              "`--tableformat' can take two values: `fits-ascii', or "
              "`fits-binary'", filename);
    }
}




















/************************************************************************/
/***************          Printing information            ***************/
/************************************************************************/
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
         && strlen(gal_type_to_string(allcols[i].type, 1))>tw)
        tw=strlen(gal_type_to_string(allcols[i].type, 1));
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
             tw, gal_type_to_string(allcols[i].type, 1),
             comment ? comment : GAL_BLANK_STRING);
    }

  /* Print the number of rows. */
  printf("--------\nNumber of rows: %zu\n--------\n", numrows);
}





/* Fill in/adjust the basic information necessary to print a column. This
   information can be used for printing a plain text file or for FITS ASCII
   tables. The `fmt' and `lng' should point to pre-allocated arrays. The
   best way is: `char fmt[2], lng[3];' in the same function calling this.

   The width and precision, which are also necessary for printing, are
   updated in the data structure's `disp_width' and `disp_precision'
   elements.
*/
void
gal_table_col_print_info(gal_data_t *col, int tableformat, char *fmt,
                         char *lng)
{
  size_t j;
  char **strarr;
  int maxstrlen, width=0, precision=0;


  /* First do a sanity check, so we can safly stop checking in the steps
     below. */
  switch(tableformat)
    {
    case GAL_TABLE_FORMAT_TXT:
    case GAL_TABLE_FORMAT_AFITS:
      break;
    default:
      error(EXIT_FAILURE, 0, "tableformat code %d not recognized in "
            "`gal_table_col_print_info'", tableformat);
    }



  /* Set the formats and widths based on the type of the column. Initialize
     the characters and blank pointer. The long prefix is not necessary for
     most types, so just initialize it once up here.*/
  fmt[0]=fmt[1]=lng[0]=lng[1]=lng[2]='\0';
  switch(col->type)
    {
    case GAL_TYPE_BIT:
      error(EXIT_FAILURE, 0, "printing of bit types is currently "
            "not supported");
      break;




    case GAL_TYPE_STRING:

      /* Set the basic information. */
      fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 's' : 'A';

      /* Go through all the strings in the column and find the maximum
         length to use as printing. If the user asked for a larger width
         (through the data structure's disp_width element), then set
         that. */
      maxstrlen=0;
      strarr=col->array;
      for(j=0;j<col->size;++j)
        maxstrlen = ( (int)strlen(strarr[j]) > maxstrlen
                      ? (int)strlen(strarr[j]) : maxstrlen );
      width = col->disp_width>maxstrlen ? col->disp_width : maxstrlen;
      break;




    case GAL_TYPE_UINT8:
    case GAL_TYPE_UINT16:
    case GAL_TYPE_UINT32:
    case GAL_TYPE_UINT64:

      /* For the FITS ASCII table, there is only one format for all
         integers.  */
      if(tableformat==GAL_TABLE_FORMAT_AFITS)
        fmt[0]='I';
      else
        switch(col->disp_fmt)
          {
          case GAL_TABLE_DISPLAY_FMT_UDECIMAL: fmt[0]='u'; break;
          case GAL_TABLE_DISPLAY_FMT_OCTAL:    fmt[0]='o'; break;
          case GAL_TABLE_DISPLAY_FMT_HEX:      fmt[0]='X'; break;
          default:                             fmt[0]='u';
          }

      /* If we have a long type, then make changes. */
      if(col->type==GAL_TYPE_UINT64)
        {
          lng[0]='l';
          width=( col->disp_width<=0 ? GAL_TABLE_DEF_LINT_WIDTH
                  : col->disp_width );
        }
      else width=( col->disp_width<=0 ? GAL_TABLE_DEF_INT_WIDTH
                    : col->disp_width );
      precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_INT_PRECISION
                  : col->disp_precision );
      break;




    case GAL_TYPE_INT8:
    case GAL_TYPE_INT16:
    case GAL_TYPE_INT32:
      fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'd' : 'I';
      width = ( col->disp_width<=0 ? GAL_TABLE_DEF_INT_WIDTH
                : col->disp_width );
      precision = ( col->disp_precision<=0 ? GAL_TABLE_DEF_INT_PRECISION
                    : col->disp_precision );
      break;




    case GAL_TYPE_INT64:
      lng[0] = 'l';
      fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'd' : 'I';
      width=( col->disp_width<=0 ? GAL_TABLE_DEF_LINT_WIDTH
              : col->disp_width );
      precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_INT_PRECISION
                  : col->disp_precision );
      break;



    /* We need a default value (because in most cases, it won't be set. */
    case GAL_TYPE_FLOAT32:
    case GAL_TYPE_FLOAT64:
      switch(col->disp_fmt)
        {
        case GAL_TABLE_DISPLAY_FMT_FLOAT:
          fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'f' : 'F'; break;
        case GAL_TABLE_DISPLAY_FMT_EXP:
          fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'e' : 'E'; break;
        case GAL_TABLE_DISPLAY_FMT_GENERAL:
          fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'g' : 'E'; break;
        default:
          fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'g' : 'E'; break;
        }
      width = ( col->disp_width<=0
                ? ( col->type==GAL_TYPE_FLOAT32
                    ? GAL_TABLE_DEF_FLT_WIDTH
                    : GAL_TABLE_DEF_DBL_WIDTH )
                : col->disp_width );
      precision = ( col->disp_precision<=0 ? GAL_TABLE_DEF_FLT_PRECISION
                    : col->disp_precision );
      break;



    default:
      error(EXIT_FAILURE, 0, "type code %d not recognized in "
            "`gal_table_col_print_info'", col->type);
    }

  /* Write the final width and precision into the column's data structure. */
  col->disp_width=width;
  col->disp_precision=precision;
}





/* Use the input `blank' string and the input column to put the blank value
   in the column's array. This function should later be generalized into a
   function to read a string into a given data type (see
   `gal_data_string_to_array_elem'). It is only here temporarily. */
void
gal_table_read_blank(gal_data_t *col, char *blank)
{
  void *colarr=col->array;

  /* If there is nothing to use as blank, then don't continue, note that
     the column data structure was initialized to mean that there is no
     blank value. */
  if(blank==NULL) return;

  /* Just for a sanity check, the ndim and array elements should be zero. */
  if(col->ndim || col->array)
    error(EXIT_FAILURE, 0, "A bug! The number of dimensions, and the pointer "
          "to the array in the data structure passed "
          "`to gal_table_read_blank' must be zero");

  /* Read the blank value as the given type. If successful, then
     `gal_data_string_to_type' will return 0. In that case, we need to
     initialize the necessary paramters to read this data structure
     correctly. */
  if( !gal_data_string_to_type(&colarr, blank, col->type) )
    {
      errno=0;
      col->dsize=malloc(sizeof *col->dsize);
      if(col->dsize==NULL)
        error(EXIT_FAILURE, 0, "%zu bytes for `col->dsize' in "
              "`txt_read_blank' ", sizeof *col->dsize);
      col->dsize[0]=col->ndim=col->size=1;
    }
}




















/************************************************************************/
/***************         Information about a table        ***************/
/************************************************************************/
/* Store the information of each column in a table (either as a text file
   or as a FITS table) into an array of data structures with `numcols'
   structures (one data structure for each column). The number of rows is
   stored as the `size' element of each data structure. The type of the
   table (e.g., ascii text file, or FITS binary or ASCII table) will be put
   in `tableformat' (macros defined in `gnuastro/table.h'.

   Note that other than the character strings (column name, units and
   comments), nothing in the data structure(s) will be allocated by this
   function for the actual data (e.g., the `array' or `dsize' elements). */
gal_data_t *
gal_table_info(char *filename, char *hdu, size_t *numcols, size_t *numrows,
               int *tableformat)
{
  /* Get the table format and size (number of columns and rows). */
  if(gal_fits_name_is_fits(filename))
    return gal_fits_tab_info(filename, hdu, numcols, numrows, tableformat);
  else
    {
      *tableformat=GAL_TABLE_FORMAT_TXT;
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

/* Function to print regular expression error. This is taken from the GNU C
   library manual, with small modifications to fit out style, */
static void
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
  size_t i, nummatch, len;
  struct gal_linkedlist_stll *tmp;
  struct gal_linkedlist_sll *indexll=NULL;
  char *str, *strcheck, *tailptr, *errorstring;

  for(tmp=cols; tmp!=NULL; tmp=tmp->next)
    {
      /* Counter for number of columns matched, and length of name. */
      nummatch=0;
      len=strlen(tmp->v);

      /* REGULAR EXPRESSION: When the first and last characters are `/'. */
      if( tmp->v[0]=='/' && tmp->v[len-1]=='/' )
        {
          /* Remove the slashes, note that we don't want to change `tmp->v'
             (because it should be freed later). So first we set the last
             character to `\0', then define a new string from the first
             element. */
          tmp->v[len-1]='\0';
          str = tmp->v + 1;

          /* Allocate the regex_t structure: */
          errno=0;
          regex=malloc(sizeof *regex);
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

          /* Put the `/' back into the input string. This is done because
             after this function, the calling program might want to inform
             the user of their exact input string. */
          tmp->v[len-1]='/';
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
                error(EXIT_FAILURE, 0, "column numbers must be positive "
                      "(not zero or negative). You have asked for column "
                      "number %ld", tlong);

              /* Check if the given value is not larger than the number of
                 columns in the input catalog (note that the user is
                 counting from 1, not 0!) */
              if(tlong>numcols)
                error(EXIT_FAILURE, 0, "%s: has %zu columns, but you "
                      "have asked for column number %zu",
                      gal_fits_name_save_as_string(filename, hdu),
                      numcols, tlong);

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
         are done (and program is aborted) before this step. */
      if(nummatch==0)
        {
          asprintf(&errorstring, "`%s' didn't match any of the column %ss.",
                   tmp->v, gal_table_searchin_as_string(searchin));
          gal_table_error_col_selection(filename, hdu, errorstring);
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
   format) enclosed in `/ /'. The `searchin' value comes from the
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
  int tableformat;
  size_t i, numcols, numrows;
  gal_data_t *allcols, *out=NULL;
  struct gal_linkedlist_sll *indexll;

  /* If the column string linked list is empty, no need to continue. */
  if(cols==NULL) return NULL;

  /* First get the information of all the columns. */
  allcols=gal_table_info(filename, hdu, &numcols, &numrows, &tableformat);

  /* If there was no actual data in the file, then return NULL. */
  if(allcols==NULL) return NULL;

  /* Get the list of indexs in the same order as the input list */
  indexll=make_list_of_indexs(cols, allcols, numcols, searchin,
                              ignorecase, filename, hdu);

  /* Depending on the table format, read the columns into the output
     structure. Note that the functions here pop each index, read/store the
     desired column and pop the next, so after these functions, the output
     linked list will have the opposite order of its input `indexll'
     list. So before calling any of them, we will first reverse the
     `indexll' list, so the output data structure list will have the same
     order as the input list of desired columns. Also note that after these
     functions, the `indexll' will be all freed (each popped element is
     actually freed).*/
  gal_linkedlist_reverse_sll(&indexll);
  switch(tableformat)
    {
    case GAL_TABLE_FORMAT_TXT:
      out=gal_txt_table_read(filename, numrows, allcols, indexll,
                             minmapsize);
      break;

    case GAL_TABLE_FORMAT_AFITS:
    case GAL_TABLE_FORMAT_BFITS:
      out=gal_fits_tab_read(filename, hdu, numrows, allcols, indexll,
                            minmapsize);
      break;

    default:
      error(EXIT_FAILURE, 0, "table format code %d not recognized for "
            "`tableformat' in `gal_table_read_cols'", tableformat);
    }

  /* Clean up. */
  for(i=0;i<numcols;++i)
    gal_data_free_contents(&allcols[i]);
  free(allcols);
  gal_linkedlist_free_sll(indexll);

  /* Return the final linked list. */
  return out;
}




















/************************************************************************/
/***************              Write a table               ***************/
/************************************************************************/
/* Write the basic information that is necessary by each program into the
   comments field. Note that the `comments' has to be already sorted in the
   proper order. */
void
gal_table_comments_add_intro(struct gal_linkedlist_stll **comments,
                             char *program_string, time_t *rawtime)
{
  char gitdescribe[100], *tmp;

  /* Get the Git description in the running folder. */
  tmp=gal_git_describe();
  if(tmp) sprintf(gitdescribe, " from %s,", tmp);
  else    gitdescribe[0]='\0';
  free(tmp);


  /* Git version and time of program's starting, this will be the second
     line. Note that ctime puts a `\n' at the end of its string, so we'll
     have to remove that. Also, note that since we are allocating `msg', we
     are setting the allocate flag of `gal_linkedlist_add_to_stll' to 0. */
  asprintf(&tmp, "Created%s on %s", gitdescribe, ctime(rawtime));
  tmp[ strlen(tmp)-1 ]='\0';
  gal_linkedlist_add_to_stll(comments, tmp, 0);


  /* Program name: this will be the top of the list (first line). We will
     need to set the allocation flag for this one, because program_string
     is usually statically allocated.*/
  gal_linkedlist_add_to_stll(comments, program_string, 1);
}





/* The input is a linked list of data structures and some comments. The
   table will then be written into `filename' with a format that is
   specified by `tableformat'. If it already exists, and `dontdelete' has a
   value of 1, then it won't be deleted and an error will be printed. */
void
gal_table_write(gal_data_t *cols, struct gal_linkedlist_stll *comments,
                int tableformat, char *filename, int dontdelete)
{
  /* If a filename was given, then the tableformat is relevant and must be
     used. When the filename is empty, a text table must be printed on the
     standard output (on the command-line). */
  if(filename)
    {
      if(gal_fits_name_is_fits(filename))
        gal_fits_tab_write(cols, comments, tableformat, filename,
                           dontdelete);
      else
        gal_txt_write(cols, comments, filename, dontdelete);
    }
  else
    gal_txt_write(cols, comments, filename, dontdelete);
}





void
gal_table_write_log(gal_data_t *logll, char *program_string,
                    time_t *rawtime, struct gal_linkedlist_stll *comments,
                    char *filename, int dontdelete, int quiet)
{
  char *msg;

  /* Write all the comments into */
  gal_table_comments_add_intro(&comments, program_string, rawtime);

  /* Write the log file to disk */
  gal_table_write(logll, comments, GAL_TABLE_FORMAT_TXT, filename,
                  dontdelete);

  /* In verbose mode, print the information. */
  if(!quiet)
    {
      asprintf(&msg, "%s created.", filename);
      gal_timing_report(NULL, msg, 1);
      free(msg);
    }
}
