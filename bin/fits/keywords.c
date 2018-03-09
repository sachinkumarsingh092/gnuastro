/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2018, Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>

#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "fits.h"








/***********************************************************************/
/******************           Preparations          ********************/
/***********************************************************************/
static void
keywords_open(struct fitsparams *p, fitsfile **fptr, int iomode)
{
  if(*fptr==NULL)
    *fptr=gal_fits_hdu_open(p->filename, p->cp.hdu, iomode);
}




















/***********************************************************************/
/******************        File manipulation        ********************/
/***********************************************************************/
static void
keywords_rename_keys(struct fitsparams *p, fitsfile **fptr, int *r)
{
  int status=0;
  char *copy, *str, *from, *to;

  /* Set the FITS file pointer. */
  keywords_open(p, fptr, READWRITE);

  /* Tokenize the */
  while(p->rename!=NULL)
    {
      /* Pop out the top element. */
      str=gal_list_str_pop(&p->rename);

      /* Take a copy of the input string for error reporting, because
         `strtok' will write into the array. */
      gal_checkset_allocate_copy(str, &copy);

      /* Tokenize the input. */
      from = strtok(str,  ", ");
      to   = strtok(NULL, ", ");

      /* Make sure both elements were read. */
      if(from==NULL || to==NULL)
        error(EXIT_FAILURE, 0, "`%s' could not be tokenized in order to "
              "complete rename. There should be a space character "
              "or a comma (,) between the two keyword names. If you have "
              "used the space character, be sure to enclose the value to "
              "the `--rename' option in double quotation marks", copy);

      /* Rename the keyword */
      fits_modify_name(*fptr, from, to, &status);
      if(status) *r=fits_has_error(p, FITS_ACTION_RENAME, from, status);
      status=0;

      /* Clean up the user's input string. Note that `strtok' just changes
         characters within the allocated string, no extra allocation is
         done. */
      free(str);
      free(copy);
    }
}





static void
keywords_write_update(struct fitsparams *p, fitsfile **fptr,
                      gal_fits_list_key_t *keyll, int u1w2)
{
  int status=0;
  gal_fits_list_key_t *tmp;

  /* Open the FITS file if it hasn't been opened yet. */
  keywords_open(p, fptr, READWRITE);

  /* Go through each key and write it in the FITS file. */
  while(keyll!=NULL)
    {
      /* Write the information: */
      if(u1w2==1)
        {
          if(keyll->value)
            {
              if( fits_update_key(*fptr,
                                  gal_fits_type_to_datatype(keyll->type),
                                  keyll->keyname, keyll->value,
                                  keyll->comment, &status) )
                gal_fits_io_error(status, NULL);
            }
          else
            {
              if(fits_write_key_null(*fptr, keyll->keyname, keyll->comment,
                                     &status))
                gal_fits_io_error(status, NULL);
            }
        }
      else if (u1w2==2)
        {
          if(keyll->value)
            {
              if( fits_write_key(*fptr,
                                 gal_fits_type_to_datatype(keyll->type),
                                 keyll->keyname, keyll->value, keyll->comment,
                                 &status) )
                gal_fits_io_error(status, NULL);
            }
          else
            {
              if(fits_write_key_null(*fptr, keyll->keyname, keyll->comment,
                                     &status))
                gal_fits_io_error(status, NULL);
            }
          if(keyll->unit
             && fits_write_key_unit(*fptr, keyll->keyname, keyll->unit,
                                    &status) )
            gal_fits_io_error(status, NULL);
        }
      else
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at `%s' so we "
              "can fix this problem. The value %d is not valid for `u1w2'",
              __func__, PACKAGE_BUGREPORT, u1w2);

      /* Add the unit (if one was given). */
      if(keyll->unit
         && fits_write_key_unit(*fptr, keyll->keyname, keyll->unit, &status) )
        gal_fits_io_error(status, NULL);

      /* Free the allocated spaces if necessary: */
      if(keyll->vfree) free(keyll->value);
      if(keyll->kfree) free(keyll->keyname);
      if(keyll->cfree) free(keyll->comment);

      /* Keep the pointer to the next keyword and free the allocated
         space for this keyword.*/
      tmp=keyll->next;
      free(keyll);
      keyll=tmp;
    }
}





static void
keywords_print_all_keys(struct fitsparams *p, fitsfile **fptr)
{
  size_t i=0;
  int nkeys, status=0;
  char *fullheader, *c, *cf;

  /* Conver the header into a contiguous string. */
  if( fits_hdr2str(*fptr, 0, NULL, 0, &fullheader, &nkeys, &status) )
    gal_fits_io_error(status, NULL);

  /* FLEN_CARD supposes that the NULL string character is in the
     end of each keyword header card. In fits_hdr2str, the NULL
     characters are removed and so the maximum length is one
     less. */
  cf=(c=fullheader)+nkeys*(FLEN_CARD-1);
  do
    {
      if(i && i%(FLEN_CARD-1)==0)
        putc('\n', stdout);
      putc(*c++, stdout);
      ++i;
    }
  while(c<cf);
  printf("\n");

  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in header.c for freeing "
                      "the memory used to keep all the headers");
}




















/***********************************************************************/
/******************           Main function         ********************/
/***********************************************************************/
/* NOTE ON CALLING keywords_open FOR EACH OPERATION:

   `keywords_open' is being called individually for each separate operation
   because the necessary permissions differ: when the user only wants to
   read keywords, they don't necessarily need write permissions. So if they
   haven't asked for any writing/editing operation, we shouldn't open in
   write-mode. Because the user might not have the permissions to write and
   they might not want to write. `keywords_open' will only open the file
   once (if the pointer is already allocated, it won't do anything). */
int
keywords(struct fitsparams *p)
{
  int status=0;
  int r=EXIT_SUCCESS;
  fitsfile *fptr=NULL;
  gal_list_str_t *tstll;


  /* Delete the requested keywords. */
  if(p->delete)
    {
      /* Open the FITS file. */
      keywords_open(p, &fptr, READWRITE);

      /* Go over all the keywords to delete. */
      for(tstll=p->delete; tstll!=NULL; tstll=tstll->next)
        {
          fits_delete_key(fptr, tstll->v, &status);
          if(status)
            r=fits_has_error(p, FITS_ACTION_DELETE, tstll->v, status);
          status=0;
        }
    }


  /* Rename the requested keywords. */
  if(p->rename)
    {
      keywords_open(p, &fptr, READWRITE);
      keywords_rename_keys(p, &fptr, &r);
    }


  /* Update the requested keywords. */
  if(p->update)
    {
      keywords_open(p, &fptr, READWRITE);
      keywords_write_update(p, &fptr, p->update_keys, 1);
    }


  /* Write the requested keywords. */
  if(p->write)
    {
      keywords_open(p, &fptr, READWRITE);
      keywords_write_update(p, &fptr, p->write_keys, 2);
    }


  /* Put in any full line of keywords as-is. */
  if(p->asis)
    {
      keywords_open(p, &fptr, READWRITE);
      for(tstll=p->asis; tstll!=NULL; tstll=tstll->next)
        {
          fits_write_record(fptr, tstll->v, &status);
          if(status) r=fits_has_error(p, FITS_ACTION_WRITE, tstll->v, status);
          status=0;
        }
    }


  /* Add the history keyword(s). */
  if(p->history)
    {
      keywords_open(p, &fptr, READWRITE);
      for(tstll=p->history; tstll!=NULL; tstll=tstll->next)
        {
          fits_write_history(fptr, tstll->v, &status);
          if(status)
            r=fits_has_error(p, FITS_ACTION_WRITE, "HISTORY", status);
          status=0;
        }
    }


  /* Add comment(s). */
  if(p->comment)
    {
      keywords_open(p, &fptr, READWRITE);
      for(tstll=p->comment; tstll!=NULL; tstll=tstll->next)
        {
          fits_write_comment(fptr, tstll->v, &status);
          if(status)
            r=fits_has_error(p, FITS_ACTION_WRITE, "COMMENT", status);
          status=0;
        }
    }


  /* Update/add the date. */
  if(p->date)
    {
      keywords_open(p, &fptr, READWRITE);
      fits_write_date(fptr, &status);
      if(status) r=fits_has_error(p, FITS_ACTION_WRITE, "DATE", status);
      status=0;
    }


  /* Print all the keywords in the extension. */
  if(p->printallkeys)
    {
      keywords_open(p, &fptr, READONLY);
      keywords_print_all_keys(p, &fptr);
    }


  /* Close the FITS file */
  if(fits_close_file(fptr, &status))
    gal_fits_io_error(status, NULL);


  /* Return. */
  return r;
}
