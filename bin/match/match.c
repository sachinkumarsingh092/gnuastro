/*********************************************************************
Match - A program to match catalogs and WCS warps
Match is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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
#include <gnuastro/permutation.h>

#include <gnuastro-internal/checkset.h>

#include <main.h>





/* Read the catalog in the given file and use the given permutation to keep
   the proper columns. */
static void
match_catalog_write(struct matchparams *p, char *filename, char *hdu,
                    size_t *permutation, size_t nummatched, char *outname,
                    char *extname)
{
  gal_data_t *tmp, *cat;

  /* Read the full table. */
  cat=gal_table_read(filename, hdu, NULL,p->cp.searchin, p->cp.ignorecase,
                     p->cp.minmapsize);

  /* Go over each column and permute its contents. */
  for(tmp=cat; tmp!=NULL; tmp=tmp->next)
    {
      /* Do the permutation. */
      gal_permutation_apply(tmp, permutation);

      /* Correct the size of the array so only the matching columns are
         saved as output. This is only Gnuastro's convention, it has no
         effect on later freeing of the array in the memory. */
      tmp->size=tmp->dsize[0]=nummatched;
    }

  /* Write the catalog to the output. */
  gal_table_write(cat, NULL, p->cp.tableformat, outname, extname);

  /* Clean up. */
  gal_list_data_free(cat);
}





static void
match_catalog(struct matchparams *p)
{
  uint32_t *u, *uf;
  size_t nummatched;
  gal_data_t *tmp, *mcols;

  /* Find the matching coordinates. We are doing the processing in
     place, */
  mcols=gal_match_coordinates(p->cols1, p->cols2, p->aperture->array, 0, 1,
                              p->cp.minmapsize, &nummatched);

  /* If a match was found, then make the output files. */
  if(mcols)
    {
      /* Read all the first catalog columns. */
      if(p->logasoutput==0)
        {
          match_catalog_write(p, p->input1name, p->cp.hdu, mcols->array,
                              nummatched, p->out1name, "INPUT_1");
          match_catalog_write(p, p->input2name, p->hdu2, mcols->next->array,
                              nummatched, p->out2name, "INPUT_2");
        }

      /* Write the raw information in a log file if necessary.  */
      if(p->logname)
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
                          "LOG_INFO");

          /* Set the comment pointer to NULL: they weren't allocated. */
          mcols->comment=NULL;
          mcols->next->comment=NULL;
        }
      gal_list_data_free(mcols);
    }

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
  if(p->mode==MATCH_MODE_CATALOG)
    match_catalog(p);
}
