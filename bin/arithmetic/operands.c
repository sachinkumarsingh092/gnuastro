/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "operands.h"






size_t
operands_num(struct arithmeticparams *p)
{
  size_t counter=0;
  struct operand *tmp=NULL;
  for(tmp=p->operands;tmp!=NULL;tmp=tmp->next)
    ++counter;
  return counter;
}





void
operands_add(struct arithmeticparams *p, char *filename, gal_data_t *data)
{
  struct operand *newnode;

  /* Some operators might not actually return any dataset (data=NULL), in
     such cases filename will also be NULL (since the operand was not added
     from the command-line). So, we shouldn't add anything to the stack. */
  if(data || filename)
    {
      /* Allocate space for the new operand. */
      errno=0;
      newnode=malloc(sizeof *newnode);
      if(newnode==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `newnode'",
              __func__, sizeof *newnode);

      /* Fill in the values. */
      newnode->data=data;
      newnode->filename=filename;

      if(filename != NULL && gal_fits_name_is_fits(filename))
        {
          /* Set the HDU for this filename. */
          if(p->globalhdu)
            gal_checkset_allocate_copy(p->globalhdu, &newnode->hdu);
          else
            newnode->hdu=gal_list_str_pop(&p->hdus);

          /* Increment the FITS counter. */
          ++p->addcounter;
        }

      /* Make the link to the previous list. */
      newnode->next=p->operands;
      p->operands=newnode;
    }
}





gal_data_t *
operands_pop(struct arithmeticparams *p, char *operator)
{
  size_t i;
  gal_data_t *data;
  char *filename, *hdu;
  struct operand *operands=p->operands;

  /* If the operand linked list has finished, then give an error and
     exit. */
  if(operands==NULL)
    error(EXIT_FAILURE, 0, "not enough operands for the \"%s\" operator",
          operator);


  /* Set the dataset. If filename is present then read the file
     and fill in the array, if not then just set the array. */
  if(operands->filename)
    {
      hdu=operands->hdu;
      filename=operands->filename;

      /* Read the dataset. */
      data=gal_fits_img_read(filename, hdu, p->cp.minmapsize, 0, 0);

      /* In case this is the first image that is read, then keep the WCS
         information in the `refdata' structure. Otherwise, the WCS is not
         necessary and we can safely free it. In any case, `data' must not
         have a WCS structure. */
      if(p->popcounter==0)
        {
          p->refdata.wcs=data->wcs;
          p->refdata.nwcs=data->nwcs;
        }
      else
        wcsfree(data->wcs);
      data->wcs=NULL;
      data->nwcs=0;

      /* When the reference data structure's dimensionality is non-zero, it
         means that this is not the first image read. So, write its basic
         information into the reference data structure for future
         checks. */
      if(p->refdata.ndim)
        {
          if(gal_data_dsize_is_different(&p->refdata, data))
            error(EXIT_FAILURE, 0, "%s (hdu=%s): has a different size "
                  "compared to previous images. All the images must be "
                  "the same size in order for Arithmetic to work",
                  filename, hdu);
        }
      else
        {
          /* Set the dimensionality. */
          p->refdata.ndim=(data)->ndim;

          /* Allocate the dsize array. */
          errno=0;
          p->refdata.dsize=malloc(p->refdata.ndim * sizeof *p->refdata.dsize);
          if(p->refdata.dsize==NULL)
            error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for "
                  "p->refdata.dsize", __func__,
                  p->refdata.ndim * sizeof *p->refdata.dsize);

          /* Write the values into it. */
          for(i=0;i<p->refdata.ndim;++i)
            p->refdata.dsize[i]=data->dsize[i];
        }

      /* Report the read image if desired: */
      if(!p->cp.quiet) printf(" - %s (hdu %s) is read.\n", filename, hdu);

      /* Free the HDU string: */
      free(hdu);

      /* Add to the number of popped FITS images: */
      ++p->popcounter;
    }
  else
    data=operands->data;

  /* Remove this node from the queue, return the data structure. */
  p->operands=operands->next;
  free(operands);
  return data;
}
