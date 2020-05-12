/*********************************************************************
Match - A program to match catalogs and WCS warps
Match is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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
#include <stdlib.h>
#include <string.h>

#include <gnuastro/match.h>
#include <gnuastro/table.h>
#include <gnuastro/pointer.h>
#include <gnuastro/permutation.h>

#include <gnuastro-internal/checkset.h>

#include <main.h>



/* Number of columns in a file. */
static gal_list_str_t *
match_add_all_cols(char *filename, char *extname, gal_list_str_t *stdinlines,
                   gal_list_str_t *incols, size_t *num)
{
  char *tstr;
  int tableformat;
  gal_data_t *colinfo=NULL;
  gal_list_str_t *tmp, *finalcols=NULL;
  size_t i, numrows, numcols=GAL_BLANK_SIZE_T;

  /* Go over all the given input columns. */
  for(tmp=incols; tmp!=NULL; tmp=tmp->next)
    {
      if(!strcmp(tmp->v,"_all"))
        {
          /* Read all the column information (if it hasn't been read until
             now). */
          if( numcols == GAL_BLANK_SIZE_T )
            {
              colinfo=gal_table_info(filename, extname,
                                     filename ? NULL : stdinlines, &numcols,
                                     &numrows, &tableformat);
              gal_data_array_free(colinfo, numcols, 1);
            }

          /* Add each column number to the list of columns. */
          for(i=0;i<numcols;++i)
            {
              errno=0;
              if( asprintf(&tstr, "%zu", i+1)<0 )
                error(EXIT_FAILURE, errno, "asprintf allocation");
              gal_list_str_add(&finalcols, tstr, 0);
            }
        }
      else
        gal_list_str_add(&finalcols, tmp->v, 1);
    }

  /* If a new list of columns is ready, re-order tham and write
     them in. Note that there may be multiple '_all' terms, so we
     need to do this after parsing all the requested columns. */
  gal_list_str_reverse(&finalcols);

  /* For a check.
  gal_list_str_print(finalcols);
  exit(1);
  */

  /* Clean up and return. */
  *num=numcols;
  return finalcols;
}





static gal_data_t *
match_cat_from_coord(struct matchparams *p, gal_list_str_t *cols,
                     size_t *numcolmatch)
{
  void *rptr;
  gal_list_str_t *col;
  uint8_t read, readtype;
  size_t colcounter, counter;
  gal_data_t *tmp, *ttmp, *out=NULL;

  /* Go over the desired columns and only return the good ones. */
  colcounter=0;
  for(col=cols;col!=NULL;col=col->next)
    {
      /* In 'ui_preparations_out_cols', we have done the necessary sanity
         checks, so we can safely use the values. */
      rptr=gal_type_string_to_number(col->v, &readtype);
      if(readtype!=GAL_TYPE_UINT8)
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
              "the problem. The given string didn't have a 'uint8' type",
              __func__, PACKAGE_BUGREPORT);
      read=*((uint8_t *)rptr);

      /* Find the proper column in the second input's columns. Just note
         that column counting starts from 1.*/
      counter=1;
      for(tmp=p->cols2;tmp!=NULL;tmp=tmp->next)
        if(counter++ == read)
          {
            ttmp=gal_data_copy(tmp);
            ttmp->next=NULL;
            gal_list_data_add(&out, ttmp);
            ++numcolmatch[colcounter];
            break;
          }

      /* Increment the column counter. */
      ++colcounter;
    }

  /* Reverse the list. */
  gal_list_data_reverse(&out);

  /* Return the output columns. */
  return out;
}





/* Read the catalog in the given file and use the given permutation to keep
   the proper columns. */
static gal_data_t *
match_catalog_read_write_all(struct matchparams *p, size_t *permutation,
                             size_t nummatched, int f1s2,
                             size_t **numcolmatch)
{
  int hasall=0;
  size_t origsize;
  gal_data_t *tmp, *cat;
  gal_list_str_t *cols, *tcol;
  gal_list_void_t *arrays=NULL;

  char *hdu              = (f1s2==1) ? p->cp.hdu     : p->hdu2;
  gal_list_str_t *incols = (f1s2==1) ? p->acols      : p->bcols;
  size_t *numcols        = (f1s2==1) ? &p->anum      : &p->bnum;
  char *extname          = (f1s2==1) ? "INPUT_1"     : "INPUT_2";
  char *outname          = (f1s2==1) ? p->out1name   : p->out2name;
  char *filename         = (f1s2==1) ? p->input1name : p->input2name;

  /* If special columns are requested. */
  if(p->outcols)
    {
      /* As a special situation, the user can ask to incude all of the
         columns from one of the inputs with the special '_all' name. So,
         we'll check if that is the case and write in all the columns where
         they are requested.*/
      for(tcol=incols; tcol!=NULL; tcol=tcol->next)
        if(!strcmp(tcol->v,"_all")) { hasall=1; break; }

      /* If atleast one instance of '_all' is present, then reset the list
         of columns to include in output. */
      if(hasall)
        {
          cols=match_add_all_cols(filename, hdu, p->stdinlines, incols,
                                  numcols);
          if(f1s2==1) { gal_list_str_free(p->acols, 0); p->acols=cols; }
          else        { gal_list_str_free(p->bcols, 0); p->bcols=cols; }
        }
      else
        cols=incols;


      /* When the output contains columns from both inputs, we need to keep
         the number of columns matched against each column identifier. */
      *numcolmatch=gal_pointer_allocate(GAL_TYPE_SIZE_T,
                                        gal_list_str_number(cols), 1,
                                        __func__, "numcolmatch");
    }
  else cols=incols;


  /* Read the full table. NOTE that with '--coord', for the second input,
     both 'filename' and 'p->stdinlines' will be NULL. */
  if(filename || p->stdinlines)
    cat=gal_table_read(filename, hdu, filename ? NULL : p->stdinlines, cols,
                       p->cp.searchin, p->cp.ignorecase, p->cp.minmapsize,
                       p->cp.quietmmap, *numcolmatch);
  else
    cat=match_cat_from_coord(p, cols, *numcolmatch);
  origsize = cat ? cat->size : 0;


  /* Go over each column and permute its contents. */
  if(permutation)
    for(tmp=cat; tmp!=NULL; tmp=tmp->next)
      {
        /* Do the permutation. */
        gal_permutation_apply(tmp, permutation);

        /* Correct the size of the array so only the matching columns are
           saved as output. This is only Gnuastro's convention, it has no
           effect on later freeing of the array in the memory. */
        if(p->notmatched)
          {
            /* Add the original array pointer to a list (we need to reset it
               later). */
            gal_list_void_add(&arrays, tmp->array);

            /* Reset the data structure's array element to start where the
               non-matched elements start. */
            tmp->array=gal_pointer_increment(tmp->array, nummatched,
                                             tmp->type);

            /* Correct the size of the tile. */
            tmp->size = tmp->dsize[0] = tmp->size - nummatched;
          }
        else
          tmp->size=tmp->dsize[0]=nummatched;
      }
  /* If no match was found ('permutation==NULL'), and the matched columns
     are requested, empty all the columns that are to be written (only
     keeping the meta-data). */
  else
    if(p->notmatched==0)
      {
        for(tmp=cat; tmp!=NULL; tmp=tmp->next)
          {
            tmp->size=0;
            free(tmp->dsize); tmp->dsize=NULL;
            free(tmp->array); tmp->array=NULL;
          }
      }


  /* Write the catalog to the output. */
  if(p->outcols)
    return cat;
  else if(cat)
    {
      /* Write the catalog to a file. */
      gal_table_write(cat, NULL, p->cp.tableformat, outname, extname, 0);

      /* Correct arrays and sizes (when 'notmatched' was called). The
         'array' element has to be corrected for later freeing.

         IMPORTANT: '--notmatched' cannot be called with '--outcols'. So
         you don't have to worry about the checks here being done later. */
      if(p->notmatched)
        {
          /* Reverse the list of array pointers to write them back in. */
          gal_list_void_reverse(&arrays);

          /* Correct the array and size pointers. */
          for(tmp=cat; tmp!=NULL; tmp=tmp->next)
            {
              tmp->array=gal_list_void_pop(&arrays);
              tmp->size=tmp->dsize[0]=origsize;
              tmp->block=NULL;
            }
        }

      /* Clean up. */
      gal_list_data_free(cat);
    }

  return NULL;
}





/* When specific columns from both inputs are requested, this function
   will write them out into a single table. */
static void
match_catalog_write_one(struct matchparams *p, gal_data_t *a, gal_data_t *b,
                        size_t *acolmatch, size_t *bcolmatch)
{
  gal_data_t *cat=NULL;
  size_t i, j, k, ac=0, bc=0, npop;
  char **strarr=p->outcols->array;

  /* Go over the initial list of strings. */
  for(i=0; i<p->outcols->size; ++i)
    switch(strarr[i][0])
      {
      case 'a':
        for(j=0;j<acolmatch[ac];++j)
          {
            npop = strcmp(strarr[i]+1,"_all") ? 1 : p->anum;
            for(k=0;k<npop;++k)
              gal_list_data_add(&cat, gal_list_data_pop(&a));
          }
        ac++;
        break;

      case 'b':
        for(j=0;j<bcolmatch[bc];++j)
          {
            npop = strcmp(strarr[i]+1,"_all") ? 1 : p->bnum;
            for(k=0;k<npop;++k)
              gal_list_data_add(&cat, gal_list_data_pop(&b));
          }
        bc++;
        break;

      default:
        error(EXIT_FAILURE, 0, "a bug! Please contact us at %s to fix the "
              "problem. the value to strarr[%zu][0] (%c) is not recognized",
              PACKAGE_BUGREPORT, i, strarr[i][0]);
      }

  /* A small sanity check. */
  if(a || b)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us to fix the problem. "
          "The two 'a' and 'b' arrays must be NULL by this point: "
          "'a' %s NULL, 'b' %s NULL", __func__, a?"is not":"is",
          b?"is not":"is");

  /* Reverse the table and write it out. */
  gal_list_data_reverse(&cat);
  gal_table_write(cat, NULL, p->cp.tableformat, p->out1name, "MATCHED", 0);
}





static void
match_catalog(struct matchparams *p)
{
  uint32_t *u, *uf;
  gal_data_t *tmp, *mcols;
  gal_data_t *a=NULL, *b=NULL;
  size_t nummatched, *acolmatch=NULL, *bcolmatch=NULL;

  /* Find the matching coordinates. We are doing the processing in
     place, */
  mcols=gal_match_coordinates(p->cols1, p->cols2, p->aperture->array, 0, 1,
                              p->cp.minmapsize, p->cp.quietmmap, &nummatched);

  /* If the output is to be taken from the input columns (it isn't just the
     log), then do the job. */
  if(p->logasoutput==0)
    {
      /* Read (and possibly write) the outputs. Note that we only need to
         read the table when it is necessary for the output (the user might
         have asked for '--outcols', only with columns of one of the two
         inputs). */
      if(p->outcols==NULL || p->acols)
        a=match_catalog_read_write_all(p, mcols?mcols->array:NULL,
                                       nummatched, 1, &acolmatch);
      if(p->outcols==NULL || p->bcols)
        b=match_catalog_read_write_all(p, mcols?mcols->next->array:NULL,
                                       nummatched, 2, &bcolmatch);

      /* If one catalog (with specific columns from either of the two
         inputs) was requested, then write it out. */
      if(p->outcols)
        {
          /* Arrange the columns and write the output. */
          match_catalog_write_one(p, a, b, acolmatch, bcolmatch);

          /* Clean up. */
          if(acolmatch) free(acolmatch);
          if(bcolmatch) free(bcolmatch);
        }
    }

  /* Write the raw information in a log file if necessary.  */
  if(p->logname && mcols)
    {
      /* Note that unsigned 64-bit integers are not recognized in FITS
         tables. So if the log file is a FITS table, covert the two
         index columns to uint32. */
      tmp=gal_data_copy_to_new_type(mcols, GAL_TYPE_UINT32);
      tmp->next=mcols->next;
      tmp->size=nummatched;
      gal_data_free(mcols);
      mcols=tmp;

      /* We also want everything to be incremented by one. In a C
         program, counting starts with zero, so 'gal_match_coordinates'
         will return indexs starting from zero. But outside a C
         program, on the command-line people expect counting to start
         from 1 (for example with AWK). */
      uf = (u=mcols->array) + tmp->size; do (*u)++; while(++u<uf);

      /* Same for the second set of indexs. */
      tmp=gal_data_copy_to_new_type(mcols->next, GAL_TYPE_UINT32);
      uf = (u=tmp->array) + tmp->size; do (*u)++; while(++u<uf);
      tmp->next=mcols->next->next;
      gal_data_free(mcols->next);
      tmp->size=nummatched;
      mcols->next=tmp;

      /* Correct the comments. */
      free(mcols->comment);
      mcols->comment="Row index in first catalog (counting from 1).";
      free(mcols->next->comment);
      mcols->next->comment="Row index in second catalog (counting "
        "from 1).";

      /* Write them into the table. */
      gal_table_write(mcols, NULL, p->cp.tableformat, p->logname,
                      "LOG_INFO", 0);

      /* Set the comment pointer to NULL: they weren't allocated. */
      mcols->comment=NULL;
      mcols->next->comment=NULL;
    }

  /* Clean up. */
  gal_list_data_free(mcols);

  /* Print the number of matches if not in quiet mode. */
  if(!p->cp.quiet)
    {
      fprintf(stdout, "Number of matching rows in both catalogs: %zu\n",
              nummatched);
      if(p->out2name && strcmp(p->out1name, p->out2name))
        fprintf(stdout, "Output:\n %s\n %s\n", p->out1name, p->out2name);
      else
        fprintf(stdout, "Output: %s\n", p->out1name);
    }
}




















/*******************************************************************/
/*************            Top level function           *************/
/*******************************************************************/
void
match(struct matchparams *p)
{
  /* Do the correct type of matching. */
  switch(p->mode)
    {
    case MATCH_MODE_CATALOG: match_catalog(p); break;
    case MATCH_MODE_WCS:
      error(EXIT_FAILURE, 0, "matching by WCS is not yet supported");
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to fix "
            "the problem: %d is not a recognized mode",
            __func__, PACKAGE_BUGREPORT, p->mode);
    }

  /* Write Match's configuration as keywords into the first extension of
     the output. */
  if(gal_fits_name_is_fits(p->out1name))
    {
      gal_fits_key_write_filename("input1", ( p->input1name
                                              ? p->input1name
                                              : "Standard input" ),
                                  &p->cp.okeys, 1);
      gal_fits_key_write_filename("input2",
                                  p->input2name?p->input2name:"--coord",
                                  &p->cp.okeys, 1);
      gal_fits_key_write_config(&p->cp.okeys, "Match configuration",
                                "MATCH-CONFIG", p->out1name, "0");
    }
}
