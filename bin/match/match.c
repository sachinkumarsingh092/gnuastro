/*********************************************************************
Match - A program to match catalogs and WCS warps
Match is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2019, Free Software Foundation, Inc.

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





/* Read the catalog in the given file and use the given permutation to keep
   the proper columns. */
static gal_data_t *
match_catalog_read_write_all(struct matchparams *p, size_t *permutation,
                             size_t nummatched, int f1s2,
                             size_t **numcolmatch)
{
  size_t origsize;
  gal_data_t *tmp, *cat;
  gal_list_void_t *arrays=NULL;

  char *hdu            = (f1s2==1) ? p->cp.hdu     : p->hdu2;
  gal_list_str_t *cols = (f1s2==1) ? p->acols      : p->bcols;
  char *extname        = (f1s2==1) ? "INPUT_1"     : "INPUT_2";
  char *outname        = (f1s2==1) ? p->out1name   : p->out2name;
  char *filename       = (f1s2==1) ? p->input1name : p->input2name;

  /* When the output contains columns from both inputs, we need to keep the
     number of columns matched against each column identifier. */
  if(p->outcols)
    *numcolmatch=gal_pointer_allocate(GAL_TYPE_SIZE_T,
                                      gal_list_str_number(cols), 0,
                                      __func__, "numcolmatch");

  /* Read the full table. */
  cat=gal_table_read(filename, hdu, filename ? NULL : p->stdinlines, cols,
                     p->cp.searchin, p->cp.ignorecase, p->cp.minmapsize,
                     *numcolmatch);
  origsize=cat->size;


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

  /* If no match was found (`permutation==NULL'), and the matched columns
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
  else
    {
      /* Write the catalog to a file. */
      gal_table_write(cat, NULL, p->cp.tableformat, outname, extname, 0);

      /* Correct arrays and sizes (when `notmatched' was called). The
         `array' element has to be corrected for later freeing.

         IMPORTANT: `--notmatched' cannot be called with `--outcols'. So
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
      return NULL;
    }
}





/* When specific columns from both inputs are requested, this function
   will write them out into a single table. */
static void
match_catalog_write_one(struct matchparams *p, gal_data_t *a, gal_data_t *b,
                        size_t *acolmatch, size_t *bcolmatch)
{
  gal_data_t *cat=NULL;
  size_t i, j, ac=0, bc=0;
  char **strarr=p->outcols->array;

  /* Go over the initial list of strings. */
  for(i=0; i<p->outcols->size; ++i)
    switch(strarr[i][0])
      {
      case 'a':
        for(j=0;j<acolmatch[ac];++j)
          gal_list_data_add(&cat, gal_list_data_pop(&a));
        ac++;
        break;

      case 'b':
        for(j=0;j<bcolmatch[bc];++j)
          gal_list_data_add(&cat, gal_list_data_pop(&b));
        bc++;
        break;

      default:
        error(EXIT_FAILURE, 0, "a bug! Please contact us at %s to fix the "
              "problem. the value to strarr[%zu][0] (%c) is not recognized",
              PACKAGE_BUGREPORT, i, strarr[i][0]);
      }

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
  size_t nummatched, *acolmatch, *bcolmatch;

  /* Find the matching coordinates. We are doing the processing in
     place, */
  mcols=gal_match_coordinates(p->cols1, p->cols2, p->aperture->array, 0, 1,
                              p->cp.minmapsize, &nummatched);

  /* If the output is to be taken from the input columns (it isn't just the
     log), then do the job. */
  if(p->logasoutput==0)
    {
      /* Read (and possibly write) the outputs. Note that we only need to
         read the table when it is necessary for the output (the user might
         have asked for `--outcols', only with columns of one of the two
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
         program, counting starts with zero, so `gal_match_coordinates'
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
    fprintf(stdout, "%zu\n", nummatched);
}




















/*******************************************************************/
/*************            Top level function           *************/
/*******************************************************************/
void
match(struct matchparams *p)
{
  /* Do the correct type of matching. */
  if(p->mode==MATCH_MODE_CATALOG)
    match_catalog(p);

  /* Write Match's configuration as keywords into the first extension of
     the output. */
  if(gal_fits_name_is_fits(p->out1name))
    {
      gal_fits_key_write_filename("input1", ( p->input1name
                                              ? p->input1name
                                              : "Standard input" ),
                                  &p->cp.okeys, 1);
      gal_fits_key_write_filename("input2", p->input2name, &p->cp.okeys, 1);
      gal_fits_key_write_config(&p->cp.okeys, "Match configuration",
                                "MATCH-CONFIG", p->out1name, "0");
    }
}
