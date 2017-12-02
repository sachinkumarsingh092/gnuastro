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
#include <stdlib.h>
#include <string.h>

#include <gnuastro/match.h>
#include <gnuastro/table.h>

#include <main.h>





/* Read the catalog in the given file and use the given permutation to keep
   the proper columns. */
static void
match_catalog_write(struct matchparams *p, char *filename, char *hdu,
                    gal_data_t *permutation, char *outname, char *extname)
{
  size_t i, *perm=permutation->array;
  gal_data_t *tmp, *cat, *newcol, *newcat=NULL;

  /* Read the full table. */
  cat=gal_table_read(filename, hdu, NULL,p->cp.searchin, p->cp.ignorecase,
                     p->cp.minmapsize);

  /* Go over each column and make a new copy. */
  for(tmp=cat; tmp!=NULL; tmp=tmp->next)
    {
      /* Allocate space for the new column. */
      newcol=gal_data_alloc(NULL, tmp->type, 1, &permutation->size,
                            NULL, 0, p->cp.minmapsize, tmp->name,
                            tmp->unit, tmp->comment);

      /* Copy the elements from the old column to the new one. We can't use
         the permute functions because the number of elements in
         `permutation', might (will probably) be different from the number
         of elements in the columns.*/
      for(i=0;i<permutation->size;++i)
        memcpy( gal_data_ptr_increment(newcol->array, i,       tmp->type),
                gal_data_ptr_increment(tmp->array,    perm[i], tmp->type),
                gal_type_sizeof(tmp->type) );

      /* Add the new column to the new catalog. */
      gal_list_data_add(&newcat, newcol);
    }

  /* Reverse the columns to match the input, free the input catalog and
     return.*/
  gal_list_data_reverse(&newcat);

  /* Write the catalog to the output. */
  gal_table_write(newcat, NULL, p->cp.tableformat, outname, extname);

  /* Clean up. */
  gal_list_data_free(cat);
  gal_list_data_free(newcat);
}





static void
match_catalog(struct matchparams *p)
{
  uint32_t *u, *uf;
  gal_data_t *tmp, *mcols;

  /* Find the matching coordinates. */
  mcols=gal_match_coordinates(p->cols1, p->cols2, p->aperture, 0, 1,
                              p->cp.minmapsize);

  /* Read all the first catalog columns. */
  if(p->logasoutput==0)
    {
      match_catalog_write(p, p->input1name, p->cp.hdu, mcols,
                          p->out1name, "INPUT_1");
      match_catalog_write(p, p->input2name, p->cp.hdu, mcols->next,
                          p->out2name, "INPUT_2");
    }

  /* Write the raw information in a log file if necessary.  */
  if(p->logname)
    {
      /* Note that unsigned 64-bit integers are not recognized in FITS
         tables. So if the log file is a FITS table, covert the two
         index columns to uint32. */
      tmp=gal_data_copy_to_new_type(mcols, GAL_TYPE_UINT32);
      tmp->next=mcols->next;
      gal_data_free(mcols);
      mcols=tmp;

      /* We also want everything to be incremented by one. In a C program,
         counting starts with zero, so `gal_match_coordinates' will return
         indexs starting from zero. But outside a C program, on the
         command-line people expect counting to start from 1 (for example
         with AWK). */
      uf = (u=mcols->array) + tmp->size; do (*u)++; while(++u<uf);

      /* Same for the second set of indexs. */
      tmp=gal_data_copy_to_new_type(mcols->next, GAL_TYPE_UINT32);
      uf = (u=tmp->array) + tmp->size; do (*u)++; while(++u<uf);
      tmp->next=mcols->next->next;
      gal_data_free(mcols->next);
      mcols->next=tmp;

      /* Correct the comments. */
      free(mcols->comment);
      mcols->comment="Row index in first catalog (counting from 1).";
      free(mcols->next->comment);
      mcols->next->comment="Row index in second catalog (counting from 1).";

      /* Write them into the table. */
      gal_table_write(mcols, NULL, p->cp.tableformat, p->logname, "LOG_INFO");

      /* Set the comment pointer to NULL: they weren't allocated. */
      mcols->comment=NULL;
      mcols->next->comment=NULL;
    }
  gal_list_data_free(mcols);
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
