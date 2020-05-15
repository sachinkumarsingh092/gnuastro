/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include "arithmetic.h"
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
  "expressions, enclose the value to the '--column' ('-c') option in "
  "slashes ('\\', as in '-c\\^mag\\'). To print the selected columns on the "
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
  cp->program_struct     = p;
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;

  /* Program-specific initialization. */
  p->head                = GAL_BLANK_SIZE_T;
  p->tail                = GAL_BLANK_SIZE_T;

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

  /* Pass 'gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then 'arg'
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign ('='). For short "
               "options, '=' should not be used and for long options, "
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
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct tableparams *p)
{
  double *darr;
  gal_data_t *tmp;

  /* Check if the format of the output table is valid, given the type of
     the output. */
  gal_tableintern_check_fits_format(p->cp.output, p->cp.tableformat);

  /* Some checks on '--range'. */
  if(p->range)
    for(tmp=p->range;tmp!=NULL;tmp=tmp->next)
      {
        /* Range needs two input numbers. */
        if(tmp->size!=2)
          error(EXIT_FAILURE, 0, "two values (separated by ':' or ',') are "
                "necessary for '--range' in this format: "
                "'--range=COLUMN,min:max'");

        /* The first must be smaller than the second. */
        darr=tmp->array;
        if( darr[0] > darr[1] )
          error(EXIT_FAILURE, 0, "first value (%g) given to '--range' must "
                "be smaller than the second (%g)", darr[0], darr[1]);
      }

  /* Basic checks for '--inpolygon' or '--outpolygon'. */
  if(p->inpolygon || p->outpolygon)
    {
      if(p->inpolygon && p->outpolygon)
        error(EXIT_FAILURE, 0, "'--inpolygon' and '--outpolygon' options "
              "cannot be called together");

      if(p->inpolygon && p->inpolygon->size!=2)
        error(EXIT_FAILURE, 0, "two columns (for coordinates) can "
              "be given to '--inpolygon'");

      if(p->outpolygon && p->outpolygon->size!=2)
        error(EXIT_FAILURE, 0, "two columns (for coordinates) can "
              "be given to '--outpolygon'");

      if(p->polygon==NULL)
        error(EXIT_FAILURE, 0, "no polygon specified for '--inpolygon' or "
              "'--outpolygon'! Please provide the vertices of the desired "
              "polygon with the '--polygon' option in this format: "
              "v1x,v1y:v2x,v2y:v3x,v3y:...");
    }


  /* Make sure '--head' and '--tail' aren't given together. */
  if(p->head!=GAL_BLANK_SIZE_T && p->tail!=GAL_BLANK_SIZE_T)
    error(EXIT_FAILURE, 0, "'--head' and '--tail' options cannot be "
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
              "file, a HDU must also be specified, you can use the '--hdu' "
              "('-h') option and give it the HDU number (starting from "
              "zero), extension name, or anything acceptable by CFITSIO");

    }
}


















/**************************************************************/
/************   List of row-selection requests   **************/
/**************************************************************/
static void
ui_list_select_add(struct list_select **list, gal_data_t *col, int type)
{
  struct list_select *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  newnode->col=col;
  newnode->type=type;
  newnode->next=*list;
  *list=newnode;
}





static gal_data_t *
ui_list_select_pop(struct list_select **list, int *type)
{
  gal_data_t *out=NULL;
  struct list_select *tmp;
  if(*list)
    {
      /* Extract all the necessary components of the node. */
      tmp=*list;
      out=tmp->col;
      *type=tmp->type;
      *list=tmp->next;

      /* Delete the node. */
      free(tmp);
    }
  return out;
}





static void
ui_list_select_reverse(struct list_select **list)
{
  int thistype;
  gal_data_t *thisdata;
  struct list_select *correctorder=NULL;

  /* Only do the reversal if there is more than one element. */
  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisdata=ui_list_select_pop(list, &thistype);
          ui_list_select_add(&correctorder, thisdata, thistype);
        }
      *list=correctorder;
    }
}





void
ui_list_select_free(struct list_select *list, int freevalue)
{
  struct list_select *tmp;
  while(list!=NULL)
    {
      tmp=list->next;
      if(freevalue)
        gal_data_free(list->col);
      free(list);
      list=tmp;
    }
}




















/**************************************************************/
/***************      Packaged columns      *******************/
/**************************************************************/
/* Return the last outcols element. */
static struct column_pack *
ui_outcols_last(struct column_pack *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}





/* Allocate a clean 'out_columns' structure and put it at the top of the
   list. */
static struct column_pack *
ui_outcols_add_new_to_end(struct column_pack **list)
{
  struct column_pack *last, *node;

  /* Allocate a new node. */
  errno=0;
  node=malloc(sizeof *node);
  if(node==NULL)
    error(EXIT_FAILURE, errno, "%s: couldn't allocate new node (%zu bytes)",
          __func__, sizeof *node);

  /* Initialize its elements. */
  node->next=NULL;
  node->numsimple=0;
  node->tokens=NULL;
  node->start=GAL_BLANK_SIZE_T;

  /* If the list already has elements, go to the last node in the list and
     add this node. */
  if(*list)
    {
      last=ui_outcols_last(*list);
      last->next=node;
    }
  else
    *list=node;

  /* Return a pointer to this node (to use temporarily). */
  return node;
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





/* The columns can be given as comma-separated values to one option or
   multiple calls to the column option. Here, we'll check if the input list
   has comma-separated values. If they do then the list will be updated to
   be fully separate. */
static void
ui_columns_prepare(struct tableparams *p)
{
  char *c, **strarr;
  gal_data_t *strs;
  size_t i, totcalled=0;
  struct column_pack *node, *last;
  gal_list_str_t *tmp, *toread=NULL;

  /* Go over the whole original list (where each node may have more than
     one value separated by a comma. */
  for(tmp=p->columns;tmp!=NULL;tmp=tmp->next)
    {
      /* Remove any possibly commented new-line where we have a backslash
         followed by a new-line character (replace the two characters with
         two single space characters). This can happen with the 'arith'
         argument in a script, when it gets long (bug #58371). But to be
         general in other cases too, we'll just correct it here. */
      for(c=tmp->v;*c!='\0';++c)
        if(*c=='\\' && *(c+1)=='\n') { *c=' '; *(++c)=' '; }

      /* Read the different comma-separated strings into an array (within a
         'gal_data_t'). */
      strs=gal_options_parse_csv_strings_raw(tmp->v, NULL, 0);
      strarr=strs->array;

      /* Go over all the given colum names/numbers. */
      for(i=0;i<strs->size;++i)
        {
          /* See if this is an arithmetic column to be processed, or the
             contents should just be printed. */
          if(!strncmp(strarr[i], ARITHMETIC_CALL, ARITHMETIC_CALL_LENGTH))
            {
              /* If this is the first arithmetic operation and the user has
                 already asked for some columns, we'll need to put all
                 previously requested simply-printed columns into an
                 'outcols' structure, then add this arithmetic operation's
                 'outcols'. */
              if(p->outcols==NULL && toread)
                {
                  /* Allocate an empty structure and set the necessary
                     pointers. */
                  node=ui_outcols_add_new_to_end(&p->outcols);
                  node->start=0;
                  node->numsimple=gal_list_str_number(toread);
                  totcalled=node->numsimple;
                }

              /* Add a new column pack, then read all the tokens (while
                 specifying which columns it needs). */
              node=ui_outcols_add_new_to_end(&p->outcols);
              arithmetic_init(p, &node->tokens, &toread, &totcalled,
                              strarr[i]+ARITHMETIC_CALL_LENGTH);
              free(strarr[i]);
            }
          /* This is a simple column (no change in values). */
          else
            {
              /* Add this column to the list of columns to read. */
              gal_list_str_add(&toread, strarr[i], 0);

              /* See if we have packaged the output columns. */
              if(p->outcols)
                {
                  /* If the previous column package was an arithmetic
                     operation, allocate a new node. */
                  last=ui_outcols_last(p->outcols);
                  if(last->tokens)
                    {
                      node=ui_outcols_add_new_to_end(&p->outcols);
                      node->start=totcalled;
                      node->numsimple=1;
                    }

                  /* The previous package of columns are simple (we don't
                     need to change their value), so we can just increment
                     the number of columns there and don't need to allocate
                     a new one. */
                  else
                    last->numsimple+=1;
                }

              /* Increment the total number of called columns. */
              totcalled+=1;
            }

          /* The pointer allocated string is either being used (and later
             freed) else, or has already been freed. So its necessary to
             set it to NULL. */
          strarr[i]=NULL;
        }

      /* Clean up. */
      gal_data_free(strs);
    }

  /* For a check
  if(p->outcols)
    {
      struct column_pack *tmp;
      struct arithmetic_token *atmp;
      for(tmp=p->outcols;tmp!=NULL;tmp=tmp->next)
        {
          if(tmp->tokens)
            for(atmp=tmp->tokens;atmp!=NULL;atmp=atmp->next)
              {
                printf("Arithmetic: ");
                if(atmp->constant) printf("Constant number\n");
                else if(atmp->index) printf("Called column: %zu\n",
                                            atmp->index);
                else if(atmp->operator!=ARITHMETIC_TABLE_OP_INVALID)
                  printf("Operator: %d\n", atmp->operator);
                else
                  error(EXIT_FAILURE, 0, "%s: UNKNOWN SITUATION!",
                        __func__);
              }
          else
            printf("Simple: start: %zu, num: %zu\n", tmp->start,
                   tmp->numsimple);
        }
    }
  */


  /* Delete the old list, then reverse the 'toread' list, and put it into
     'p->columns'. */
  gal_list_str_free(p->columns, 1);
  gal_list_str_reverse(&toread);
  p->columns=toread;
}





/* The users give column numbers counting from 1. But we need an 'index'
   (starting from 0). So if we can read it as a number, we'll subtract one
   from it. */
static size_t
ui_check_select_sort_read_col_ind(char *string)
{
  size_t out;
  void *ptr=&out;

  if( gal_type_from_string(&ptr, string, GAL_TYPE_SIZE_T) )
    out=GAL_BLANK_SIZE_T;
  else out-=1;

  return out;
}





/* See if row selection or sorting needs any extra columns to be read. */
static void
ui_check_select_sort_before(struct tableparams *p, gal_list_str_t *lines,
                            size_t *nselect, size_t *origoutncols,
                            size_t *sortindout, size_t **selectindout_out,
                            size_t **selecttypeout_out)
{
  char **strarr;
  gal_list_sizet_t *tmp, *indexll;
  gal_list_str_t *stmp, *add=NULL;
  int tableformat, selecthasname=0;
  size_t one=1, sortind=GAL_BLANK_SIZE_T;
  size_t *selectind=NULL, *selecttype=NULL;
  size_t *selectindout=NULL, *selecttypeout=NULL;
  size_t i, j, k, *s, *sf, allncols, numcols, numrows;
  gal_data_t *dtmp, *allcols, *inpolytmp=NULL, *outpolytmp=NULL;

  /* Important note: these have to be in the same order as the 'enum
     select_types' in 'main.h'. We'll fill the two polygon elements
     later. */
  gal_data_t *select[SELECT_TYPE_NUMBER]={p->range, NULL, NULL,
                                          p->equal, p->notequal};


  /* The inpolygon dataset is currently a single dataset with two elements
     (strings). But we need to have it as two linked datasets with a set
     name. */
  if(p->inpolygon)
    {
      strarr=p->inpolygon->array;
      inpolytmp=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &one, NULL, 1, -1, 1,
                               strarr[0], NULL, NULL);
      inpolytmp->next=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &one, NULL,
                                     1, -1, 1, strarr[1], NULL, NULL);
      select[SELECT_TYPE_INPOLYGON]=inpolytmp;
    }
  if(p->outpolygon)
    {
      strarr=p->outpolygon->array;
      outpolytmp=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &one, NULL, 1, -1, 1,
                               strarr[0], NULL, NULL);
      outpolytmp->next=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &one, NULL,
                                     1, -1, 1, strarr[1], NULL, NULL);
      select[SELECT_TYPE_OUTPOLYGON]=outpolytmp;
    }


  /* Allocate necessary spaces. */
  if(p->selection)
    {
      *nselect = ( gal_list_data_number(p->range)
                   + gal_list_data_number(inpolytmp)
                   + gal_list_data_number(outpolytmp)
                   + gal_list_data_number(p->equal)
                   + gal_list_data_number(p->notequal) );
      selectind=gal_pointer_allocate(GAL_TYPE_SIZE_T, *nselect, 0,
                                     __func__, "selectind");
      selecttype=gal_pointer_allocate(GAL_TYPE_SIZE_T, *nselect, 0,
                                      __func__, "selecttype");
      selectindout=gal_pointer_allocate(GAL_TYPE_SIZE_T, *nselect, 0,
                                        __func__, "selectindout");
      selecttypeout=gal_pointer_allocate(GAL_TYPE_SIZE_T, *nselect, 0,
                                         __func__, "selecttypeout");
      sf=(s=selectindout)+*nselect; do *s++=GAL_BLANK_SIZE_T; while(s<sf);
      *selectindout_out=selectindout;
      *selecttypeout_out=selecttypeout;
    }


  /* See if the given columns are numbers or names. */
  i=0;
  if(p->sort) sortind=ui_check_select_sort_read_col_ind(p->sort);
  if(p->selection)
    for(k=0;k<SELECT_TYPE_NUMBER;++k)
      for(dtmp=select[k];dtmp!=NULL;dtmp=dtmp->next)
        {
          selecttype[i] = k;
          selectind[i] = ui_check_select_sort_read_col_ind(dtmp->name);
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
          "number given to  '--sort' (%s)",
          gal_fits_name_save_as_string(p->filename, p->cp.hdu), numcols,
          p->sort);
  if(p->selection)
    for(i=0;i<*nselect;++i)
      if(selectind[i]!=GAL_BLANK_SIZE_T && selectind[i]>=numcols)
        error(EXIT_FAILURE, 0, "%s has %zu columns, less than the column "
              "number given to  '--range', '--inpolygon', '--outpolygon', "
              "'--equal', or '--sort' (%zu)",
              gal_fits_name_save_as_string(p->filename, p->cp.hdu), numcols,
              selectind[i]);


  /* If any of the columns isn't specified by an index, go over the table
     information and set the index based on the names. */
  if(p->selection)
    for(i=0;i<*nselect;++i)
      if(selectind[i]==GAL_BLANK_SIZE_T) { selecthasname=1; break; }
  if( (p->sort && sortind==GAL_BLANK_SIZE_T) || selecthasname )
    {
      /* For '--sort', go over all the columns if an index hasn't been set
         yet. If the input columns have a name, see if their names matches
         the name given to 'sort'. */
      if(p->sort && sortind==GAL_BLANK_SIZE_T)
        for(i=0;i<numcols;++i)
          if( allcols[i].name && !strcasecmp(allcols[i].name, p->sort) )
            { sortind=i; break; }

      /* Same for the selection. Just note that here we may have multiple
         calls. It is thus important to loop over the values given to range
         first, then loop over the column names from the start for each new
         '--ran */
      i=0;
      for(k=0;k<SELECT_TYPE_NUMBER;++k)
        for(dtmp=select[k];dtmp!=NULL;dtmp=dtmp->next)
          {
            if(selectind[i]==GAL_BLANK_SIZE_T)
              for(j=0;j<numcols;++j)
                if( allcols[j].name
                    && !strcasecmp(allcols[j].name, dtmp->name) )
                  { selecttype[i]=k; selectind[i]=j; break; }
            ++i;
          }
    }


  /* The columns must be good indexs now, if they don't the user didn't
     specify the name properly and the program must abort. */
  if( p->sort && sortind==GAL_BLANK_SIZE_T )
    error(EXIT_FAILURE, 0, "%s: no column named '%s' (value to '--sort') "
          "you can either specify a name or number",
          gal_fits_name_save_as_string(p->filename, p->cp.hdu), p->sort);
  if(p->selection)
    {
      i=0;
      for(k=0;k<SELECT_TYPE_NUMBER;++k)
        for(dtmp=select[k];dtmp!=NULL;dtmp=dtmp->next)
          {
            if(selectind[i]==GAL_BLANK_SIZE_T)
              error(EXIT_FAILURE, 0, "%s: no column named '%s' (value to "
                    "'--%s') you can either specify a name or number",
                    gal_fits_name_save_as_string(p->filename, p->cp.hdu),
                    dtmp->name,
                    ( k==0?"range":( k==1?"equal":"notequal") ));
            ++i;
          }
    }


  /* See which columns the user has asked to output. */
  indexll=gal_table_list_of_indexs(p->columns, allcols, numcols,
                                   p->cp.searchin, p->cp.ignorecase,
                                   p->filename, p->cp.hdu, NULL);
  allncols=*origoutncols=gal_list_sizet_number(indexll);


  /* See if the requested columns are already on the to-read list. If so,
     keep the counter. */
  i=0;
  for(tmp=indexll; tmp!=NULL; tmp=tmp->next)
    {
      if(p->sort && *sortindout==GAL_BLANK_SIZE_T  && tmp->v == sortind)
        *sortindout=i;
      if(p->selection)
        for(j=0;j<*nselect;++j)
          if(selectindout[j]==GAL_BLANK_SIZE_T && tmp->v==selectind[j])
            {
              selectindout[j]=i;
              selecttypeout[j]=selecttype[j];
            }
      ++i;
    }


  /* See if any of the sorting or selection columns aren't requested as an
     output by the user. If there is, keep their new label.

     Note that the sorting and range may be requested on the same
     column. In this case, we don't want to read the same column twice. */
  if( p->sort && *sortindout==GAL_BLANK_SIZE_T )
    { *sortindout=allncols++;  gal_list_str_add(&add, p->sort, 0); }
  if(p->selection)
    {
      i=0;
      for(k=0;k<SELECT_TYPE_NUMBER;++k)
        for(dtmp=select[k];dtmp!=NULL;dtmp=dtmp->next)
          {
            if(*sortindout!=GAL_BLANK_SIZE_T && selectindout[i]==*sortindout)
              {
                selecttypeout[i]=k;
                selectindout[i]=*sortindout;
              }
            else
              {
                if( selectindout[i]==GAL_BLANK_SIZE_T )
                  {
                    selecttypeout[i]=k;
                    selectindout[i]=allncols++;
                    gal_list_str_add(&add, dtmp->name, 1);
                  }
              }
            ++i;
          }
    }


  /* If any new (not requested by the user to output) columns must be read,
     add them to the list of columns to read from the input file. */
  if(add)
    {
      gal_list_str_reverse(&add);
      for(stmp=p->columns; stmp!=NULL; stmp=stmp->next)
        if(stmp->next==NULL) { stmp->next=add; break; }
    }


  /* Clean up. */
  gal_list_sizet_free(indexll);
  if(selectind) free(selectind);
  if(selecttype) free(selecttype);
  gal_data_array_free(allcols, numcols, 0);
  if(inpolytmp) gal_list_data_free(inpolytmp);
  if(outpolytmp) gal_list_data_free(outpolytmp);
}





static void
ui_check_select_sort_after(struct tableparams *p, size_t nselect,
                           size_t origoutncols, size_t sortindout,
                           size_t *selectindout, size_t *selecttypeout)
{
  size_t i, j;
  struct list_select *rtmp;
  gal_data_t *tmp, *origlast=NULL;

  /* Allocate the necessary arrays. */
  if(p->selection)
    p->freeselect=gal_pointer_allocate(GAL_TYPE_UINT8, nselect, 1,
                                       __func__, "p->freeselect");


  /* Set some necessary pointers (last pointer of actual output table and
     pointer to the sort column). */
  i=0;
  for(tmp=p->table; tmp!=NULL; tmp=tmp->next)
    {
      if(i==origoutncols-1)        origlast=tmp;
      if(p->sort && i==sortindout) p->sortcol=tmp;
      ++i;
    }


  /* Since we can have several selection columns, we'll treat them
     differently. */
  for(i=0;i<nselect;++i)
    {
      j=0;
      for(tmp=p->table; tmp!=NULL; tmp=tmp->next)
        {
          if(j==selectindout[i])
            {
              ui_list_select_add(&p->selectcol, tmp, selecttypeout[i]);
              break;
            }
          ++j;
        }
    }
  ui_list_select_reverse(&p->selectcol);


  /* Terminate the desired output table where it should be terminated (by
     setting 'origlast->next' to NULL. */
  origlast->next=NULL;


  /*  Also, remove any possibly existing 'next' pointer for 'sortcol' and
     'selectcol'. */
  if(p->sort && sortindout>=origoutncols)
    { p->sortcol->next=NULL;  p->freesort=1; }
  else p->sortin=1;
  if(p->selection)
    {
      i=0;
      for(rtmp=p->selectcol;rtmp!=NULL;rtmp=rtmp->next)
        {
          if(selectindout[i]>=origoutncols)
            {
              rtmp->col->next=NULL;
              p->freeselect[i] = (rtmp->col==p->sortcol) ? 0 : 1;
            }
          ++i;
        }
    }
}






static void
ui_preparations(struct tableparams *p)
{
  size_t *colmatch;
  gal_list_str_t *lines;
  size_t nselect=0, origoutncols=0;
  size_t sortindout=GAL_BLANK_SIZE_T;
  struct gal_options_common_params *cp=&p->cp;
  size_t *selectindout=NULL, *selecttypeout=NULL;

  /* If there were no columns specified or the user has asked for
     information on the columns, we want the full set of columns. */
  if(p->information)
    ui_print_info_exit(p);


  /* Prepare the column names. */
  ui_columns_prepare(p);


  /* If the input is from stdin, save it as 'lines'. */
  lines=gal_options_check_stdin(p->filename, p->cp.stdintimeout, "input");


  /* If any kind of row-selection is requested set 'p->selection' to 1. */
  p->selection = ( p->range || p->inpolygon || p->outpolygon || p->equal
                   || p->notequal );


  /* If row sorting or selection are requested, see if we should read any
     extra columns. */
  if(p->selection || p->sort)
    ui_check_select_sort_before(p, lines, &nselect, &origoutncols,
                                &sortindout, &selectindout, &selecttypeout);


  /* If we have any arithmetic operations, we need to make sure how many
     columns match every given column name. */
  colmatch = ( p->outcols
               ? gal_pointer_allocate(GAL_TYPE_SIZE_T,
                                      gal_list_str_number(p->columns), 1,
                                      __func__, "colmatch")
               : NULL);


  /* Read the necessary columns. */
  p->table=gal_table_read(p->filename, cp->hdu, lines, p->columns,
                          cp->searchin, cp->ignorecase, cp->minmapsize,
                          p->cp.quietmmap, colmatch);
  if(p->filename==NULL) p->filename="stdin";
  gal_list_str_free(lines, 1);


  /* If row sorting or selection are requested, keep them as separate
     datasets.*/
  if(p->selection || p->sort)
    ui_check_select_sort_after(p, nselect, origoutncols, sortindout,
                               selectindout, selecttypeout);


  /* If there was no actual data in the file, then inform the user and
     abort. */
  if(p->table==NULL)
    error(EXIT_FAILURE, 0, "%s: no usable data rows (non-commented and "
          "non-blank lines)", p->filename);


  /* Set the final indexs. */
  if(p->outcols)
    arithmetic_indexs_final(p, colmatch);


  /* Now that the data columns are ready, we can free the string linked
     list. */
  gal_list_str_free(p->columns, 1);
  p->columns=NULL;


  /* Make sure the (possible) output name is writable. */
  gal_checkset_writable_remove(p->cp.output, 0, p->cp.dontdelete);


  /* If the head or tail values are given and are larger than the number of
     rows, just set them to the number of rows (print the all the final
     rows). This is how the 'head' and 'tail' programs of GNU Coreutils
     operate. */
  p->head = ( ((p->head!=GAL_BLANK_SIZE_T) && (p->head > p->table->size))
              ? p->table->size
              : p->head );
  p->tail = ( ((p->tail!=GAL_BLANK_SIZE_T) && (p->tail > p->table->size))
              ? p->table->size
              : p->tail );


  /* Clean up. */
  free(colmatch);
  if(selectindout) free(selectindout);
  if(selecttypeout) free(selecttypeout);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct tableparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program ('args.h')
     and for the common options to all Gnuastro ('commonopts.h'). We want
     to directly put the pointers to the fields in 'p' and 'cp', so we are
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
  if(p->colarray) free(p->colarray);
}
