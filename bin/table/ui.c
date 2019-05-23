/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2019, Free Software Foundation, Inc.

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

#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "authors-cite.h"





/**************************************************************/
/*********      Argp necessary global entities     ************/
/**************************************************************/
/* Definition parameters for the Argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "ASTRdata";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" can be used to view the "
  "information, select columns, or convert tables. The inputs and outputs "
  "can be plain text (with white-space or comma as delimiters), FITS ascii, "
  "or FITS binary tables. The output columns can either be selected by "
  "number (counting from 1), name or using regular expressions. For regular "
  "expressions, enclose the value to the `--column' (`-c') option in "
  "slashes (`\\', as in `-c\\^mag\\'). To print the selected columns on the "
  "command-line, don't specify an output file.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct tableparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;

  /* Program-specific initialization. */
  p->head                = GAL_BLANK_SIZE_T;
  p->tail                = GAL_BLANK_SIZE_T;
  p->wcstoimg            = GAL_BLANK_SIZE_T;
  p->imgtowcs            = GAL_BLANK_SIZE_T;

  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        /* Mandatory options. */
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
        case GAL_OPTIONS_KEY_TABLEFORMAT:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;

        /* Options to ignore. */
        case GAL_OPTIONS_KEY_TYPE:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }

      /* Select by group. */
      switch(cp->coptions[i].group)
        {
        case GAL_OPTIONS_GROUP_TESSELLATION:
          cp->coptions[i].doc=NULL; /* Necessary to remove title. */
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }
    }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct tableparams *p = state->input;

  /* Pass `gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then `arg`
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign (`=`). For short "
               "options, `=` should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {

    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      if(p->filename)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->filename=arg;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct tableparams *p)
{
  double *darr;
  gal_data_t *tmp;

  /* Check if the format of the output table is valid, given the type of
     the output. */
  gal_tableintern_check_fits_format(p->cp.output, p->cp.tableformat);

  /* Some checks on `--range'. */
  if(p->range)
    for(tmp=p->range;tmp!=NULL;tmp=tmp->next)
      {
        /* Range needs two input numbers. */
        if(tmp->size!=2)
          error(EXIT_FAILURE, 0, "two values (separated by comma) necessary "
                "for `--range' in this format: `--range=COLUMN,min,max'");

        /* The first must be smaller than the second. */
        darr=tmp->array;
        if( darr[0] > darr[1] )
          error(EXIT_FAILURE, 0, "first value (%g) given to `--range' must "
                "be smaller than the second (%g)", darr[0], darr[1]);
      }

  /* Make sure `--head' and `--tail' aren't given together. */
  if(p->head!=GAL_BLANK_SIZE_T && p->tail!=GAL_BLANK_SIZE_T)
    error(EXIT_FAILURE, 0, "`--head' and `--tail' options cannot be "
          "called together");
}





static void
ui_check_options_and_arguments(struct tableparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->filename)
    {
      if( gal_fits_name_is_fits(p->filename) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified. When the input is a FITS "
              "file, a HDU must also be specified, you can use the `--hdu' "
              "(`-h') option and give it the HDU number (starting from "
              "zero), extension name, or anything acceptable by CFITSIO");

    }
}


















/**************************************************************/
/***************   List of range datasets   *******************/
/**************************************************************/
static void
ui_list_range_add(struct list_range **list, gal_data_t *dataset)
{
  struct list_range *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  newnode->v=dataset;
  newnode->next=*list;
  *list=newnode;
}





static gal_data_t *
ui_list_range_pop(struct list_range **list)
{
  gal_data_t *out=NULL;
  struct list_range *tmp;
  if(*list)
    {
      tmp=*list;
      out=tmp->v;
      *list=tmp->next;
      free(tmp);
    }
  return out;
}





static void
ui_list_range_reverse(struct list_range **list)
{
  gal_data_t *thisdata;
  struct list_range *correctorder=NULL;

  /* Only do the reversal if there is more than one element. */
  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisdata=ui_list_range_pop(list);
          ui_list_range_add(&correctorder, thisdata);
        }
      *list=correctorder;
    }
}





void
ui_list_range_free(struct list_range *list, int freevalue)
{
  struct list_range *tmp;
  while(list!=NULL)
    {
      tmp=list->next;
      if(freevalue)
        gal_data_free(list->v);
      free(list);
      list=tmp;
    }
}


















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_print_info_exit(struct tableparams *p)
{
  char *tmp;
  int tableformat;
  gal_data_t *allcols;
  gal_list_str_t *lines;
  size_t numcols, numrows;

  /* Read the table information for the number of columns and rows. */
  lines=gal_options_check_stdin(p->filename, p->cp.stdintimeout, "input");
  allcols=gal_table_info(p->filename, p->cp.hdu, lines, &numcols,
                         &numrows, &tableformat);
  if(p->filename==NULL) p->filename="Standard-input";
  gal_list_str_free(lines, 1);

  /* If there was no actual data in the file, then inform the user */
  if(allcols==NULL)
    error(EXIT_FAILURE, 0, "%s: no usable data rows", p->filename);


  /* Print the file information. */
  printf("--------\n");
  tmp=gal_fits_name_save_as_string(p->filename, p->cp.hdu);
  printf("%s\n", tmp);
  free(tmp);


  /* Print each column's information. */
  gal_table_print_info(allcols, numcols, numrows);


  /* Free the information from all the columns. */
  gal_data_array_free(allcols, numcols, 0);


  /* Free the allocated spaces and exit. Otherwise, add the number of
     columns to the list if the user wanted to print the columns
     (didn't just want their information. */
  ui_free_report(p);
  exit(EXIT_SUCCESS);
}





/* WCS <-> Image conversion */
static void
ui_wcs_conversion(struct tableparams *p, char *instr,
                  gal_list_str_t **wcstoimg_ptr,
                  gal_list_str_t **imgtowcs_ptr,
                  gal_list_str_t **new)
{
  const char delimiter[]="/";
  char *c1, *c2, *c3, *optname, *saveptr;

  /* Set the basic option properties. */
  instr[8]='\0';
  optname = instr;

  /* If a WCS hasn't been read yet, read it.*/
  if(p->wcs==NULL)
    {
      /* A small sanity check. */
      if(p->wcsfile==NULL || p->wcshdu==NULL)
        error(EXIT_FAILURE, 0, "`--wcsfile' and `--wcshdu' are necessary for "
              "the `%s' conversion", optname);

      /* Read the WCS. */
      p->wcs=gal_wcs_read(p->wcsfile, p->wcshdu, 0, 0, &p->nwcs);
      if(p->wcs==NULL)
        error(EXIT_FAILURE, 0, "%s (hdu: %s): no WCS could be read by "
              "WCSLIB", p->wcsfile, p->wcshdu);
    }

  /* Make sure this conversion is only requested once. */
  if( instr[0]=='w' ? *wcstoimg_ptr : *imgtowcs_ptr )
    error(EXIT_FAILURE, 0, "`%s' can only be called once.", optname);

  /* Read the first token. */
  c1=strtok_r(instr+9, delimiter, &saveptr);
  if(c1==NULL)
    error(EXIT_FAILURE, 0, "`%s' must end with `)'", optname);
  if(*c1==')')
    error(EXIT_FAILURE, 0, "`%s' has no defining column. Please specify "
          "two column names or numbers with a slash between them for "
          "example `%s(RA/DEC)'", optname, optname);

  /* Read the second token and make sure it ends with `)'. Then set the `)'
     as the string-NULL character. */
  c2=strtok_r(NULL,        delimiter, &saveptr);
  if(c2==NULL || c2[strlen(c2)-1]!=')')
    error(EXIT_FAILURE, 0, "`%s' needs two input columns. Please "
          "specify two columns with a slash between them for "
          "example `%s(RA/DEC)'", optname, optname);
  c2[strlen(c2)-1]='\0';

  /* Make sure there aren't any more elements. */
  c3=strtok_r(NULL, delimiter, &saveptr);
  if(c3!=NULL)
    error(EXIT_FAILURE, 0, "only two values can be given to `%s'", optname);

  /* Add the column name/number to the list of columns to read, then keep
     the pointer to the node (so we can later identify it). */
  gal_list_str_add(new, c1, 1);
  gal_list_str_add(new, c2, 1);

  /* Keep this node's pointer. */
  if(instr[0]=='w') *wcstoimg_ptr=*new;
  else              *imgtowcs_ptr=*new;

  /* Clean up (the `c1' and `c2') are actually within the allocated space
     of `instr'. That is why we are setting the last argument of
     `gal_list_str_add' to `1' (to allocate space for them). */
  free(instr);
}





/* The columns can be given as comma-separated values to one option or
   multiple calls to the column option. Here, we'll check if the input list
   has comma-separated values. If they do then the list will be updated to
   be fully separate. */
static void
ui_columns_prepare(struct tableparams *p, size_t *wcstoimg, size_t *imgtowcs)
{
  size_t i;
  char **strarr;
  gal_data_t *strs;
  gal_list_str_t *tmp, *new=NULL;
  gal_list_str_t *wcstoimg_ptr=NULL, *imgtowcs_ptr=NULL;

  /* Go over the whole original list (where each node may have more than
     one value separated by a comma. */
  for(tmp=p->columns;tmp!=NULL;tmp=tmp->next)
    {
      /* Read the different comma-separated strings into an array (within a
         `gal_data_t'). */
      strs=gal_options_parse_csv_strings_raw(tmp->v, NULL, 0);
      strarr=strs->array;

      /* Go over all the given colum names/numbers. */
      for(i=0;i<strs->size;++i)
        {
          if(!strncmp(strarr[i],"wcstoimg(",9)
             || !strncmp(strarr[i],"imgtowcs(",9))
            ui_wcs_conversion(p, strarr[i], &wcstoimg_ptr, &imgtowcs_ptr, &new);
          else
            gal_list_str_add(&new, strarr[i], 0);

          /* The pointer allocated string is either being used (and later
             freed) else, or has already been freed. So its necessary to
             set it to NULL. */
          strarr[i]=NULL;
        }

      /* Clean up. */
      gal_data_free(strs);
    }

  /* Delete the old list. */
  gal_list_str_free(p->columns, 1);

  /* Reverse the new list, then put it into `p->columns'. */
  gal_list_str_reverse(&new);
  p->columns=new;

  /* If conversion is necessary, set the final column number. */
  i=0;
  if(wcstoimg_ptr || imgtowcs_ptr)
    for(tmp=p->columns;tmp!=NULL;tmp=tmp->next)
      {
        if(wcstoimg_ptr==tmp) *wcstoimg=i-1;
        if(imgtowcs_ptr==tmp) *imgtowcs=i-1;
        ++i;
      }
}





/* The users give column numbers counting from 1. But we need an `index'
   (starting from 0). So if we can read it as a number, we'll subtract one
   from it. */
static size_t
ui_check_range_sort_read_col_ind(char *string)
{
  size_t out;
  void *ptr=&out;

  if( gal_type_from_string(&ptr, string, GAL_TYPE_SIZE_T) )
    out=GAL_BLANK_SIZE_T;
  else out-=1;

  return out;
}





/* See if the `--range' and `--sort' columns should also be added. */
static void
ui_check_range_sort_before(struct tableparams *p, gal_list_str_t *lines,
                           size_t *nrange, size_t *origoutncols,
                           size_t *sortindout, size_t **rangeindout_out)
{
  size_t *rangeind=NULL;
  size_t *rangeindout=NULL;
  gal_data_t *dtmp, *allcols;
  size_t sortind=GAL_BLANK_SIZE_T;
  int tableformat, rangehasname=0;
  gal_list_sizet_t *tmp, *indexll;
  gal_list_str_t *stmp, *add=NULL;
  size_t i, j, *s, *sf, allncols, numcols, numrows;


  /* Allocate necessary spaces. */
  if(p->range)
    {
      *nrange=gal_list_data_number(p->range);
      rangeind=gal_pointer_allocate(GAL_TYPE_SIZE_T, *nrange, 0,
                                    __func__, "rangeind");
      rangeindout=gal_pointer_allocate(GAL_TYPE_SIZE_T, *nrange, 0,
                                        __func__, "rangeindout");
      sf=(s=rangeindout)+*nrange; do *s++=GAL_BLANK_SIZE_T; while(s<sf);
      *rangeindout_out=rangeindout;
    }


  /* See if the given columns are numbers or names. */
  i=0;
  if(p->sort)  sortind  = ui_check_range_sort_read_col_ind(p->sort);
  if(p->range)
    for(dtmp=p->range;dtmp!=NULL;dtmp=dtmp->next)
      {
        rangeind[i] = ui_check_range_sort_read_col_ind(dtmp->name);
        ++i;
      }


  /* Get all the column information. */
  allcols=gal_table_info(p->filename, p->cp.hdu, lines, &numcols,
                         &numrows, &tableformat);


  /* If the values are numbers, we'll check if the given value is less than
     the total number of columns. Just note that the indexs count from
     zero. */
  if(p->sort && sortind!=GAL_BLANK_SIZE_T && sortind>=numcols)
    error(EXIT_FAILURE, 0, "%s has %zu columns, less than the column "
          "number given to  `--sort' (%s)",
          gal_fits_name_save_as_string(p->filename, p->cp.hdu), numcols,
          p->sort);
  if(p->range)
    for(i=0;i<*nrange;++i)
      if(rangeind[i]!=GAL_BLANK_SIZE_T && rangeind[i]>=numcols)
        error(EXIT_FAILURE, 0, "%s has %zu columns, less than the column "
              "number given to  `--range' (%zu)",
              gal_fits_name_save_as_string(p->filename, p->cp.hdu), numcols,
              rangeind[i]);


  /* If any of the columns isn't specified by an index, go over the table
     information and set the index based on the names. */
  if(p->range)
    for(i=0;i<*nrange;++i)
      if(rangeind[i]==GAL_BLANK_SIZE_T) { rangehasname=1; break; }
  if( (p->sort && sortind==GAL_BLANK_SIZE_T) || rangehasname )
    {
      /* For `--sort', go over all the columns if an index hasn't been set
         yet. If the input columns have a name, see if their names matches
         the name given to `sort'. */
      if(p->sort && sortind==GAL_BLANK_SIZE_T)
        for(i=0;i<numcols;++i)
          if( allcols[i].name && !strcasecmp(allcols[i].name, p->sort) )
            { sortind=i; break; }

      /* Same for `--range'. Just note that here we may have multiple calls
         to `--range'. It is thus important to loop over the values given
         to range first, then loop over the column names from the start for
         each new `--ran */
      i=0;
      if(p->range)
        for(dtmp=p->range;dtmp!=NULL;dtmp=dtmp->next)
          {
           if(rangeind[i]==GAL_BLANK_SIZE_T)
             for(j=0;j<numcols;++j)
               if( allcols[j].name
                   && !strcasecmp(allcols[j].name, dtmp->name) )
                 { rangeind[i]=j; break; }
           ++i;
          }
    }


  /* Both columns must be good indexs now, if they aren't the user didn't
     specify the name properly and the program must abort. */
  if( p->sort && sortind==GAL_BLANK_SIZE_T )
    error(EXIT_FAILURE, 0, "%s: no column named `%s' (value to `--sort') "
          "you can either specify a name or number",
          gal_fits_name_save_as_string(p->filename, p->cp.hdu), p->sort);
  if(p->range)
    {
      i=0;
      for(dtmp=p->range;dtmp!=NULL;dtmp=dtmp->next)
        {
          if(rangeind[i]==GAL_BLANK_SIZE_T)
            error(EXIT_FAILURE, 0, "%s: no column named `%s' (value to "
                  "`--range') you can either specify a name or number",
                  gal_fits_name_save_as_string(p->filename, p->cp.hdu),
                  dtmp->name);
          ++i;
        }
    }


  /* See which columns the user has asked for. */
  indexll=gal_table_list_of_indexs(p->columns, allcols, numcols,
                                   p->cp.searchin, p->cp.ignorecase,
                                   p->filename, p->cp.hdu, NULL);
  allncols=*origoutncols=gal_list_sizet_number(indexll);


  /* See if the requested columns are already on the to-read list. If so,
     keep the counter. */
  i=0;
  for(tmp=indexll; tmp!=NULL; tmp=tmp->next)
    {
      if(p->sort  && *sortindout==GAL_BLANK_SIZE_T  && tmp->v == sortind)
        *sortindout=i;
      if(p->range)
        for(j=0;j<*nrange;++j)
          if(rangeindout[j]==GAL_BLANK_SIZE_T && tmp->v==rangeind[j])
            rangeindout[j]=i;
      ++i;
    }


  /* See if any of the necessary columns (for `--sort' and `--range')
     aren't requested as an output by the user. If there is any, such
     columns, keep them here. */
  if( p->sort && *sortindout==GAL_BLANK_SIZE_T )
    { *sortindout=allncols++;  gal_list_str_add(&add, p->sort, 0); }


  /* Note that the sorting and range may be requested on the same
     column. In this case, we don't want to read the same column twice. */
  if(p->range)
    {
      i=0;
      for(dtmp=p->range;dtmp!=NULL;dtmp=dtmp->next)
        {
          if(*sortindout!=GAL_BLANK_SIZE_T
             && rangeindout[i]==*sortindout)
            rangeindout[i]=*sortindout;
          else
            {
              if( rangeindout[i]==GAL_BLANK_SIZE_T )
                {
                  rangeindout[i]=allncols++;
                  gal_list_str_add(&add, dtmp->name, 0);
                }
            }
          ++i;
        }
    }


  /* Add the possibly new set of columns to read. */
  if(add)
    {
      gal_list_str_reverse(&add);
      for(stmp=p->columns; stmp!=NULL; stmp=stmp->next)
        if(stmp->next==NULL) { stmp->next=add; break; }
    }


  /* Clean up. */
  if(rangeind) free(rangeind);
  gal_list_sizet_free(indexll);
  gal_data_array_free(allcols, numcols, 0);
}





static void
ui_check_range_sort_after(struct tableparams *p, size_t nrange,
                          size_t origoutncols, size_t sortindout,
                          size_t *rangeindout)
{
  struct list_range *rtmp;
  size_t i, j, *rangein=NULL;
  gal_data_t *tmp, *last=NULL;

  /* Allocate the necessary arrays. */
  if(p->range)
    {
      rangein=gal_pointer_allocate(GAL_TYPE_UINT8, nrange, 0,
                                   __func__, "rangein");
      p->freerange=gal_pointer_allocate(GAL_TYPE_UINT8, nrange, 1,
                                        __func__, "p->freerange");
    }


  /* Set the proper pointers. For `rangecol' we'll need to do it separately
     (because the orders can get confused).*/
  i=0;
  for(tmp=p->table; tmp!=NULL; tmp=tmp->next)
    {
      if(i==origoutncols-1)           last=tmp;
      if(p->sort && i==sortindout) p->sortcol=tmp;
      ++i;
    }


  /* Find the range columns. */
  for(i=0;i<nrange;++i)
    {
      j=0;
      for(tmp=p->table; tmp!=NULL; tmp=tmp->next)
        {
          if(j==rangeindout[i])
            {
              ui_list_range_add(&p->rangecol, tmp);
              break;
            }
          ++j;
        }
    }
  ui_list_range_reverse(&p->rangecol);


  /* Terminate the actual table where it should be terminated (by setting
     `last->next' to NULL. */
  last->next=NULL;


  /*  Also, remove any possibly existing `next' pointer for `sortcol' and
     `rangecol'. */
  if(p->sort && sortindout>=origoutncols)
    { p->sortcol->next=NULL;  p->freesort=1; }
  else p->sortin=1;
  if(p->range)
    {
      i=0;
      for(rtmp=p->rangecol;rtmp!=NULL;rtmp=rtmp->next)
        {
          if(rangeindout[i]>=origoutncols)
            {
              rtmp->v->next=NULL;
              p->freerange[i] = (rtmp->v==p->sortcol) ? 0 : 1;
            }
          else rangein[i]=1;
          ++i;
        }
    }


  /* Clean up. */
  if(rangein) free(rangein);
}






static void
ui_preparations(struct tableparams *p)
{
  gal_list_str_t *lines;
  size_t i, ncolnames, *colmatch;
  size_t nrange=0, origoutncols=0;
  struct gal_options_common_params *cp=&p->cp;
  size_t sortindout=GAL_BLANK_SIZE_T, *rangeindout=NULL;
  size_t wcstoimg=GAL_BLANK_SIZE_T, imgtowcs=GAL_BLANK_SIZE_T;

  /* If there were no columns specified or the user has asked for
     information on the columns, we want the full set of columns. */
  if(p->information)
    ui_print_info_exit(p);


  /* Prepare the column names. */
  ui_columns_prepare(p, &wcstoimg, &imgtowcs);


  /* If the input is from stdin, save it as `lines'. */
  lines=gal_options_check_stdin(p->filename, p->cp.stdintimeout, "input");


  /* If sort or range are given, see if we should read them also. */
  if(p->range || p->sort)
    ui_check_range_sort_before(p, lines, &nrange, &origoutncols, &sortindout,
                               &rangeindout);


  /* If any conversions must be done, we need to know how many matches
     there were for each column. */
  ncolnames=gal_list_str_number(p->columns);
  colmatch = ( (wcstoimg!=GAL_BLANK_SIZE_T || imgtowcs!=GAL_BLANK_SIZE_T)
               ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ncolnames, 1,
                                      __func__, "colmatch")
               : NULL );


  /* Read the necessary columns. */
  p->table=gal_table_read(p->filename, cp->hdu, lines, p->columns,
                          cp->searchin, cp->ignorecase, cp->minmapsize,
                          colmatch);
  if(p->filename==NULL) p->filename="stdin";
  gal_list_str_free(lines, 1);


  /* If we have a unit conversion, find the proper column to use. */
  if(wcstoimg!=GAL_BLANK_SIZE_T)
    { p->wcstoimg=0; for(i=0;i<wcstoimg;++i) p->wcstoimg += colmatch[i]; }
  if(imgtowcs!=GAL_BLANK_SIZE_T)
    { p->imgtowcs=0; for(i=0;i<imgtowcs;++i) p->imgtowcs += colmatch[i]; }


  /* If the range and sort options are requested, keep them as separate
     datasets. */
  if(p->range || p->sort)
    ui_check_range_sort_after(p, nrange, origoutncols, sortindout,
                              rangeindout);


  /* If there was no actual data in the file, then inform the user and
     abort. */
  if(p->table==NULL)
    error(EXIT_FAILURE, 0, "%s: no usable data rows (non-commented and "
          "non-blank lines)", p->filename);


  /* Now that the data columns are ready, we can free the string linked
     list. */
  gal_list_str_free(p->columns, 1);
  p->columns=NULL;


  /* Make sure the (possible) output name is writable. */
  gal_checkset_writable_remove(p->cp.output, 0, p->cp.dontdelete);


  /* If the head or tail values are given and are larger than the number of
     rows, just set them to the number of rows (print the all the final
     rows). This is how the `head' and `tail' programs of GNU Coreutils
     operate. */
  p->head = ( ((p->head!=GAL_BLANK_SIZE_T) && (p->head > p->table->size))
              ? p->table->size
              : p->head );
  p->tail = ( ((p->tail!=GAL_BLANK_SIZE_T) && (p->tail > p->table->size))
              ? p->table->size
              : p->tail );


  /* Clean up. */
  free(colmatch);
  if(rangeindout) free(rangeindout);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct tableparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program (`args.h')
     and for the common options to all Gnuastro (`commonopts.h'). We want
     to directly put the pointers to the fields in `p' and `cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <gnuastro-internal/commonopts.h>
#include "args.h"


  /* Initialize the options and necessary information.  */
  ui_initialize_options(p, program_options, gal_commonopts_options);


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Read the configuration files and set the common values. */
  gal_options_read_config_set(&p->cp);


  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
  ui_read_check_only_options(p);


  /* Print the option values if asked. Note that this needs to be done
     after the option checks so un-sane values are not printed in the
     output state. */
  gal_options_print_state(&p->cp);


  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);


  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct tableparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
  gal_list_data_free(p->table);
}
