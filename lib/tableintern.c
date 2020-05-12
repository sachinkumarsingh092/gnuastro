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
#include <gnuastro/pointer.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>






/************************************************************************/
/***************              Error messages              ***************/
/************************************************************************/
void
gal_tableintern_error_col_selection(char *filename, char *hdu,
                                    char *errorstring)
{
  char *c, *name, *command;

  /* Set the proper pointers. */
  if(gal_fits_name_is_fits(filename))
    {
      if( asprintf(&name, "%s (hdu: %s)", filename, hdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      c=hdu; while(*c!='\0') if(isspace(*c++)) break;
      if( asprintf(&command, *c=='\0' ? "%s --hdu=%s" : "%s --hdu=\"%s\"",
                   filename, hdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }
  else command=name=filename?filename:"stdin";

  /* Abort with with the proper error. */
  error(EXIT_FAILURE, 0, "%s: %s\n\nFor more information on selecting "
        "columns in Gnuastro, please run the following command (press "
        "'SPACE' to go down and 'q' to return to the command-line):\n\n"
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
gal_tableintern_string_to_format(char *string)
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
gal_tableintern_format_as_string(uint8_t tableformat)
{
  switch(tableformat)
    {
    case GAL_TABLE_FORMAT_TXT:    return "txt";
    case GAL_TABLE_FORMAT_AFITS:  return "fits-ascii";
    case GAL_TABLE_FORMAT_BFITS:  return "fits-binary";
    default:
      error(EXIT_FAILURE, 0, "%s: code %d not recognized", __func__,
            tableformat);
      return NULL;
    }
}






/* In programs, the 'searchin' variable is much more easier to format in as
   a description than an integer (which is what 'gal_table_read_cols'
   needs). This function will check the string value and give the
   corresponding integer value.*/
uint8_t
gal_tableintern_string_to_searchin(char *string)
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
gal_tableintern_searchin_as_string(uint8_t searchin)
{
  switch(searchin)
    {
    case GAL_TABLE_SEARCH_NAME:    return "name";
    case GAL_TABLE_SEARCH_UNIT:    return "unit";
    case GAL_TABLE_SEARCH_COMMENT: return "comment";
    default:
      error(EXIT_FAILURE, 0, "%s: code %d not recognized as a valid search "
            "field", __func__, searchin);
      return NULL;
    }
}





/* For programs that output tables, the '--tableformat' option will be used
   to specify what format the output table should be in. When the output
   file is a FITS file, there are multiple formats, so to simplify the
   coding in each program, this function will do a sanity check on the
   value given to the '--tableformat' parameter. */
void
gal_tableintern_check_fits_format(char *filename, int tableformat)
{
  if( filename && gal_fits_name_is_fits(filename) )
    {
      /* When '--tableformat' was not given. */
      if(tableformat==GAL_TABLE_FORMAT_INVALID)
        error(EXIT_FAILURE, 0, "'%s' (output file) is a FITS file but the "
              "desired format of the FITS table has not been specified with "
              "the '--tableformat' option. For FITS tables, this option can "
              "take two values: 'fits-ascii', or 'fits-binary'", filename);

      /* When '--tableformat' didn't have the correct value. */
      if( tableformat != GAL_TABLE_FORMAT_AFITS
          && tableformat != GAL_TABLE_FORMAT_BFITS )
        error(EXIT_FAILURE, 0, "'%s' (output file) is a FITS file but "
              "is not a recognized FITS table format. For FITS tables, "
              "'--tableformat' can take two values: 'fits-ascii', or "
              "'fits-binary'", filename);
    }
}




















/************************************************************************/
/***************          Printing information            ***************/
/************************************************************************/
/* Fill in/adjust the basic information necessary to print a column. This
   information can be used for printing a plain text file or for FITS ASCII
   tables. The 'fmt' and 'lng' should point to pre-allocated arrays. The
   best way is: 'char fmt[2], lng[3];' in the same function calling this.

   The width and precision, which are also necessary for printing, are
   updated in the data structure's 'disp_width' and 'disp_precision'
   elements. */
void
gal_tableintern_col_print_info(gal_data_t *col, int tableformat,
                               char *fmt, char *lng)
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
      error(EXIT_FAILURE, 0, "%s: is only for plain text or FITS ASCII "
            "tables. The input 'tableformat' code %d not recognized",
            __func__, tableformat);
    }



  /* Set the formats and widths based on the type of the column. Initialize
     the characters and blank pointer. The long prefix is not necessary for
     most types, so just initialize it once up here.*/
  fmt[0]=fmt[1]=lng[0]=lng[1]=lng[2]='\0';
  switch(col->type)
    {
    case GAL_TYPE_BIT:
      error(EXIT_FAILURE, 0, "%s: printing of bit types is currently "
            "not supported", __func__);
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
          width=( col->disp_width<=0 ? GAL_TABLE_DEF_WIDTH_LINT
                  : col->disp_width );
        }
      else width=( col->disp_width<=0 ? GAL_TABLE_DEF_WIDTH_INT
                    : col->disp_width );
      precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_PRECISION_INT
                  : col->disp_precision );
      break;




    case GAL_TYPE_INT8:
    case GAL_TYPE_INT16:
    case GAL_TYPE_INT32:
      fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'd' : 'I';
      width = ( col->disp_width<=0 ? GAL_TABLE_DEF_WIDTH_INT
                : col->disp_width );
      precision = ( col->disp_precision<=0 ? GAL_TABLE_DEF_PRECISION_INT
                    : col->disp_precision );
      break;




    case GAL_TYPE_INT64:
      lng[0] = 'l';
      fmt[0] = tableformat==GAL_TABLE_FORMAT_TXT ? 'd' : 'I';
      width=( col->disp_width<=0 ? GAL_TABLE_DEF_WIDTH_LINT
              : col->disp_width );
      precision=( col->disp_precision<=0 ? GAL_TABLE_DEF_PRECISION_INT
                  : col->disp_precision );
      break;



    /* We need a default value (because in most cases, it won't be set. */
    case GAL_TYPE_FLOAT32:
    case GAL_TYPE_FLOAT64:
      /* Set the format. */
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

      /* Set the width and precision */
      switch(col->type)
        {
        case GAL_TYPE_FLOAT32:
          width     = ( col->disp_width<=0
                        ? GAL_TABLE_DEF_WIDTH_FLT : col->disp_width );
          precision = ( col->disp_precision<=0
                        ? GAL_TABLE_DEF_PRECISION_FLT : col->disp_precision );
          break;
        case GAL_TYPE_FLOAT64:
          width     = ( col->disp_width<=0
                        ? GAL_TABLE_DEF_WIDTH_DBL : col->disp_width );

          /* CFITSIO doesn't recognize the double precision defined here
             for ASCII FITS tables. */
          precision = ( col->disp_precision<=0
                        ? ( tableformat==GAL_TABLE_FORMAT_TXT
                            ? GAL_TABLE_DEF_PRECISION_DBL
                            : GAL_TABLE_DEF_PRECISION_FLT )
                        : col->disp_precision );
          break;
        }
      break;



    default:
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",
            __func__, col->type);
    }

  /* Write the final width and precision into the column's data structure. */
  col->disp_width=width;
  col->disp_precision=precision;
}





/* Use the input 'blank' string and the input column to put the blank value
   in the column's array. This function should later be generalized into a
   function to read a string into a given data type (see
   'gal_data_string_to_array_elem'). It is only here temporarily. */
void
gal_tableintern_read_blank(gal_data_t *col, char *blank)
{
  /* If there is nothing to use as blank, then don't continue, note that
     the column data structure was initialized to mean that there is no
     blank value. */
  if(blank==NULL) return;

  /* Just for a sanity check, the ndim and array elements should be zero. */
  if(col->ndim || col->array)
    error(EXIT_FAILURE, 0, "%s: the number of dimensions, and the "
          "'array' element of 'col' must be zero", __func__);

  /* Read the blank value as the given type. If successful, then
     'gal_data_string_to_type' will return 0. In that case, we need to
     initialize the necessary paramters to read this data structure
     correctly. */
  if( !gal_type_from_string((void **)(&col->array), blank, col->type) )
    {
      col->dsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, 1, 0, __func__,
                                      "col->dsize");
      col->dsize[0]=col->ndim=col->size=1;
    }
}
