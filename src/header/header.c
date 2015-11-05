/*********************************************************************
Header - View and manipulate a data file header
Header is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

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

#include "linkedlist.h"
#include "fitsarrayvv.h"

#include "main.h"



int
haserror(struct headerparams *p, int actionid, char *string, int status)
{
  char action[20];
  int r=EXIT_SUCCESS;

  switch(actionid)
    {
    case 1:
      strcpy(action, "deleted");
      break;
    case 2:
      strcpy(action, "renamed");
      break;
    case 3:
      strcpy(action, "updated");
      break;
    case 4:
      strcpy(action, "written");
      break;
    default:
      error(EXIT_FAILURE, 0, "A bug! Please contact us at `%s' so we can fix "
            "this problem. In `header.c'. The value of actionid in "
            "`haserror' must not be %d.", PACKAGE_BUGREPORT, actionid);
    }

  if(p->quitonerror)
    {
      fits_report_error(stderr, status);
      error(EXIT_FAILURE, 0, "Not deleted: %s\n", string);
    }
  else
    {
      fprintf(stderr, "Not deleted: %s\n", string);
      r=EXIT_FAILURE;
    }
  return r;
}





void
writeupdatekeys(fitsfile *fptr, struct fitsheaderll **keylist, int u1w2)
{
  int status=0;
  struct fitsheaderll *tmp, *ttmp;

  tmp=*keylist;
  while(tmp!=NULL)
    {

      /* Write the information: */
      if(u1w2==1)
        {
          if(tmp->value)
            {
              if( fits_update_key(fptr, tmp->datatype, tmp->keyname,
                                  tmp->value, tmp->comment, &status) )
                gal_fitsarray_io_error(status, NULL);
            }
          else
            {
              if(fits_write_key_null(fptr, tmp->keyname, tmp->comment,
                                     &status))
                gal_fitsarray_io_error(status, NULL);
            }
        }
      else if (u1w2==2)
        {
          if(tmp->value)
            {
              if( fits_write_key(fptr, tmp->datatype, tmp->keyname,
                                 tmp->value, tmp->comment, &status) )
                gal_fitsarray_io_error(status, NULL);
            }
          else
            {
              if(fits_write_key_null(fptr, tmp->keyname, tmp->comment,
                                     &status))
                gal_fitsarray_io_error(status, NULL);
            }
          if(tmp->unit
             && fits_write_key_unit(fptr, tmp->keyname, tmp->unit, &status) )
            gal_fitsarray_io_error(status, NULL);
        }
      else
        error(EXIT_FAILURE, 0, "A bug! Please contact us at `%s' so we can "
              "fix this problem. In `header.c'. The value of u1w2 in "
              "writeupdatekeys must not be %d.\n", PACKAGE_BUGREPORT, u1w2);

      /* Add the unit: */
      if(tmp->unit
         && fits_write_key_unit(fptr, tmp->keyname, tmp->unit, &status) )
        gal_fitsarray_io_error(status, NULL);

      /* Free the value pointer if desired: */
      if(tmp->kfree) free(tmp->keyname);
      if(tmp->vfree) free(tmp->value);
      if(tmp->cfree) free(tmp->comment);

      /* Keep the pointer to the next keyword and free the allocated
	 space for this keyword.*/
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
  *keylist=NULL;
}





int
header(struct headerparams *p)
{
  size_t i=0;
  int r=EXIT_SUCCESS;
  int nkeys, status=0;
  char *fullheader, *c, *cf;
  struct stll *tstll, *ttstll;

  if(p->onlyview)
    {
      if( fits_hdr2str(p->fptr, 0, NULL, 0, &fullheader, &nkeys, &status) )
        gal_fitsarray_io_error(status, NULL);

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
        gal_fitsarray_io_error(status, "Problem in header.c for freeing "
                               "the memory used to keep all the headers.");
    }
  else
    {
      if(p->delete)
        {
          for(tstll=p->delete; tstll!=NULL; tstll=tstll->next)
            {
              fits_delete_key(p->fptr, tstll->v, &status);
              if(status) r=haserror(p, 1, tstll->v, status);
            }
        }
      if(p->up.rename)
        {
          ttstll=p->renameto;
          for(tstll=p->renamefrom; tstll!=NULL; tstll=tstll->next)
            {
              fits_modify_name(p->fptr, tstll->v, ttstll->v, &status);
              if(status) r=haserror(p, 2, tstll->v, status);
              ttstll=ttstll->next;
            }
        }
      if(p->update)
        writeupdatekeys(p->fptr, &p->update, 1);
      if(p->write)
        writeupdatekeys(p->fptr, &p->write, 2);
      if(p->asis)
        for(tstll=p->asis; tstll!=NULL; tstll=tstll->next)
          {
            fits_write_record(p->fptr, tstll->v, &status);
            if(status) r=haserror(p, 1, tstll->v, status);
          }
      if(p->history)
        {
          fits_write_history(p->fptr, p->history, &status);
          if(status) r=haserror(p, 4, "HISTORY", status);
        }
      if(p->comment)
        {
          fits_write_comment(p->fptr, p->comment, &status);
          if(status) r=haserror(p, 4, "COMMENT", status);
        }
      if(p->date)
        {
          fits_write_date(p->fptr, &status);
          if(status) r=haserror(p, 4, "DATE", status);
        }
    }

  return r;
}
