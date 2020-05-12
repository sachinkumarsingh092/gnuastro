/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <gnuastro-internal/timing.h>

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
         'strtok' will write into the array. */
      gal_checkset_allocate_copy(str, &copy);

      /* Tokenize the input. */
      from = strtok(str,  ", ");
      to   = strtok(NULL, ", ");

      /* Make sure both elements were read. */
      if(from==NULL || to==NULL)
        error(EXIT_FAILURE, 0, "'%s' could not be tokenized in order to "
              "complete rename. There should be a space character "
              "or a comma (,) between the two keyword names. If you have "
              "used the space character, be sure to enclose the value to "
              "the '--rename' option in double quotation marks", copy);

      /* Rename the keyword */
      fits_modify_name(*fptr, from, to, &status);
      if(status) *r=fits_has_error(p, FITS_ACTION_RENAME, from, status);
      status=0;

      /* Clean up the user's input string. Note that 'strtok' just changes
         characters within the allocated string, no extra allocation is
         done. */
      free(str);
      free(copy);
    }
}





/* Special write options don't have any value and the value has to be found
   within the script. */
static int
keywords_write_set_value(struct fitsparams *p, fitsfile **fptr,
                         gal_fits_list_key_t *keyll)
{
  int status=0;

  if( !strcasecmp(keyll->keyname,"checksum")
      || !strcasecmp(keyll->keyname,"datasum") )
    {
      /* If a value is given, then just write what the user gave. */
      if( keyll->value )
        return 1;
      else
        {
          /* Calculate and write the 'CHECKSUM' and 'DATASUM' keywords. */
          if( fits_write_chksum(*fptr, &status) )
            gal_fits_io_error(status, NULL);

          /* If the user just wanted datasum, remove the checksum
             keyword. */
          if( !strcasecmp(keyll->keyname,"datasum") )
            if( fits_delete_key(*fptr, "CHECKSUM", &status) )
              gal_fits_io_error(status, NULL);

          /* Inform the caller that everything is done. */
          return 0;
        }
    }
  else if( keyll->keyname[0]=='/' )
    {
      gal_fits_key_write_title_in_ptr(keyll->value, *fptr);
      return 0;
    }
  else
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
          "fix the problem. The 'keyname' value '%s' is not "
          "recognized as one with no value", __func__,
          PACKAGE_BUGREPORT, keyll->keyname);

  /* Function should not reach here. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix this "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return 0;
}





static void
keywords_write_update(struct fitsparams *p, fitsfile **fptr,
                      gal_fits_list_key_t *keyll, int u1w2)
{
  int status=0, continuewriting=0;
  gal_fits_list_key_t *tmp;

  /* Open the FITS file if it hasn't been opened yet. */
  keywords_open(p, fptr, READWRITE);

  /* Go through each key and write it in the FITS file. */
  while(keyll!=NULL)
    {
      /* Deal with special keywords. */
      continuewriting=1;
      if( keyll->value==NULL || keyll->keyname[0]=='/' )
        continuewriting=keywords_write_set_value(p, fptr, keyll);

      /* Write the information: */
      if(continuewriting)
        {
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
                  if(fits_write_key_null(*fptr, keyll->keyname,
                                         keyll->comment, &status))
                    gal_fits_io_error(status, NULL);
                }
            }
          else if (u1w2==2)
            {
              if(keyll->value)
                {
                  if( fits_write_key(*fptr,
                                     gal_fits_type_to_datatype(keyll->type),
                                     keyll->keyname, keyll->value,
                                     keyll->comment, &status) )
                    gal_fits_io_error(status, NULL);
                }
              else
                {
                  if(fits_write_key_null(*fptr, keyll->keyname,
                                         keyll->comment, &status))
                    gal_fits_io_error(status, NULL);
                }
              if(keyll->unit
                 && fits_write_key_unit(*fptr, keyll->keyname, keyll->unit,
                                        &status) )
                gal_fits_io_error(status, NULL);
            }
          else
            error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' so "
                  "we can fix this problem. The value %d is not valid for "
                  "'u1w2'", __func__, PACKAGE_BUGREPORT, u1w2);

          /* Add the unit (if one was given). */
          if(keyll->unit
             && fits_write_key_unit(*fptr, keyll->keyname, keyll->unit,
                                    &status) )
            gal_fits_io_error(status, NULL);
        }

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





static int
keywords_verify(struct fitsparams *p, fitsfile **fptr)
{
  int dataok, hduok, status=0;

  /* Ask CFITSIO to verify the two keywords. */
  if( fits_verify_chksum(*fptr, &dataok, &hduok, &status) )
    gal_fits_io_error(status, NULL);

  /* Print some introduction: */
  if(!p->cp.quiet)
    printf("%s\n"
           "Checking integrity of %s (hdu %s)\n"
           "%s"
           "--------\n"
           "Basic info (remove all extra info with '--quiet'):\n"
           "    - DATASUM: verifies only the data (not keywords).\n"
           "    - CHECKSUM: verifies data and keywords.\n"
           "They can be added-to/updated-in an extension/HDU with:\n"
           "    $ astfits %s -h%s --write=checksum\n"
           "--------\n", PROGRAM_STRING, p->filename, p->cp.hdu,
           ctime(&p->rawtime), p->filename, p->cp.hdu);

  /* Print the verification result. */
  printf("DATASUM:  %s\n",
         dataok==1 ? "Verified" : (dataok==0 ? "NOT-PRESENT" : "INCORRECT"));
  printf("CHECKSUM: %s\n",
         hduok==1  ? "Verified" : (hduok==0  ? "NOT-PRESENT" : "INCORRECT"));

  /* Return failure if any of the keywords are not verified. */
  return (dataok==-1 || hduok==-1) ? EXIT_FAILURE : EXIT_SUCCESS;
}






static void
keywords_copykeys(struct fitsparams *p, char *inkeys, size_t numinkeys)
{
  size_t i;
  int status=0;
  long initial;
  fitsfile *fptr;

  /* Initial sanity check. Since 'numinkeys' includes 'END' (counting from
     1, as we do here), the first keyword must not be larger OR EQUAL to
     'numinkeys'. */
  if(p->copykeysrange[0]>=numinkeys)
    error(EXIT_FAILURE, 0, "%s (hdu %s): first keyword number give to "
          "'--copykeys' (%ld) is larger than the number of keywords in this "
          "header (%zu, including the 'END' keyword)", p->filename, p->cp.hdu,
          p->copykeysrange[0], numinkeys);

  /* If the user wanted to count from the end (by giving a negative value),
     then do that. */
  if(p->copykeysrange[1]<0)
    {
      /* Set the last keyword requested. */
      initial=p->copykeysrange[1];
      p->copykeysrange[1] += numinkeys;

      /* Sanity check. */
      if(p->copykeysrange[0]>=p->copykeysrange[1])
        error(EXIT_FAILURE, 0, "%s (hdu %s): the last keyword given to "
              "'--copykeys' (%ld, or %ld after counting from the bottom) "
              "is earlier than the first (%ld)", p->filename, p->cp.hdu,
              initial, p->copykeysrange[1], p->copykeysrange[0]);
    }

  /* Final sanity check (on range limit). */
  if(p->copykeysrange[1]>=numinkeys)
    error(EXIT_FAILURE, 0, "%s (hdu %s): second keyword number give to "
          "'--copykeys' (%ld) is larger than the number of keywords in this "
          "header (%zu, including the 'END' keyword)", p->filename, p->cp.hdu,
          p->copykeysrange[1], numinkeys);


  /* Open the output HDU. */
  fptr=gal_fits_hdu_open(p->cp.output, p->outhdu, READWRITE);


  /* Copy the requested headers into the output. */
  for(i=p->copykeysrange[0]-1; i<=p->copykeysrange[1]-1; ++i)
    if( fits_write_record(fptr, &inkeys[i*80], &status ) )
      gal_fits_io_error(status, NULL);

  /* Close the output FITS file. */
  status=0;
  if(fits_close_file(fptr, &status))
    gal_fits_io_error(status, NULL);
}





static void
keywords_date_to_seconds(struct fitsparams *p, fitsfile *fptr)
{
  int status=0;
  double subsec;
  size_t seconds;
  char *subsecstr;
  char fitsdate[FLEN_KEYWORD];

  /* Read the requested FITS keyword. */
  if( fits_read_key(fptr, TSTRING, p->datetosec, &fitsdate, NULL, &status) )
    gal_fits_io_error(status, NULL);

  /* Return the number of seconds (and subseconds) that it corresponds
     to. */
  seconds=gal_fits_key_date_to_seconds(fitsdate, &subsecstr, &subsec);

  /* Print the result (for the sub-seconds, print everything after the */
  if( !p->cp.quiet )
    {
      printf("%s (hdu %s), key '%s': %s\n", p->filename, p->cp.hdu,
             p->datetosec, fitsdate);
      printf("Seconds since 1970/01/01 (00:00:00): %zu%s\n\n", seconds,
             subsecstr?subsecstr:"");
      printf("(To suppress verbose output, run with '-q')\n");
    }
  else
    printf("%zu%s\n", seconds, subsecstr?subsecstr:"");

  /* Clean up. */
  if(subsecstr) free(subsecstr);
}




















/***********************************************************************/
/******************           Main function         ********************/
/***********************************************************************/
/* NOTE ON CALLING keywords_open FOR EACH OPERATION:

   'keywords_open' is being called individually for each separate operation
   because the necessary permissions differ: when the user only wants to
   read keywords, they don't necessarily need write permissions. So if they
   haven't asked for any writing/editing operation, we shouldn't open in
   write-mode. Because the user might not have the permissions to write and
   they might not want to write. 'keywords_open' will only open the file
   once (if the pointer is already allocated, it won't do anything). */
int
keywords(struct fitsparams *p)
{
  char *inkeys=NULL;
  int r=EXIT_SUCCESS;
  fitsfile *fptr=NULL;
  gal_list_str_t *tstll;
  int status=0, numinkeys;


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


  /* Verify the CHECKSUM and DATASUM keys. */
  if(p->verify)
    {
      keywords_open(p, &fptr, READONLY);
      r=keywords_verify(p, &fptr);
    }


  /* If a range of keywords must be copied, get all the keywords as a
     single string. */
  if(p->copykeys)
    {
      keywords_open(p, &fptr, READONLY);
      if( fits_convert_hdr2str(fptr, 0, NULL, 0, &inkeys, &numinkeys,
                               &status) )
        gal_fits_io_error(status, NULL);
      status=0;
    }


  /* Convert the FITS date string into seconds. */
  if(p->datetosec)
    {
      keywords_open(p, &fptr, READONLY);
      keywords_date_to_seconds(p, fptr);
    }


  /* Close the FITS file */
  if(fits_close_file(fptr, &status))
    gal_fits_io_error(status, NULL);


  /* Write desired keywords into output. */
  if(p->copykeys)
    {
      keywords_copykeys(p, inkeys, numinkeys);
      free(inkeys);
    }

  /* Return. */
  return r;
}
