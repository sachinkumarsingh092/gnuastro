/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/tiff.h>
#include <gnuastro/array.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "operands.h"





/**********************************************************************/
/************            General info on operands       ***************/
/**********************************************************************/
size_t
operands_num(struct arithmeticparams *p)
{
  size_t counter=0;
  struct operand *tmp=NULL;
  for(tmp=p->operands;tmp!=NULL;tmp=tmp->next)
    ++counter;
  return counter;
}




















/**********************************************************************/
/************                Named operands             ***************/
/**********************************************************************/
static int
operands_name_is_used_later(struct arithmeticparams *p, char *name)
{
  size_t counter=0;
  gal_list_str_t *token;

  /* If the name indeed exists afterwards, then just return 1. */
  for(token=p->tokens;token!=NULL;token=token->next)
    if( counter++ > p->tokencounter && !strcmp(token->v, name) )
      return 1;

  /* If we get to this point, it means that the name doesn't exist. */
  return 0;
}





/* Remove a name from the list of names and return the dataset it points
   to. */
static gal_data_t *
operands_remove_name(struct arithmeticparams *p, char *name)
{
  gal_data_t *tmp, *removed=NULL, *prev=NULL;

  /* Go over all the given names. */
  for(tmp=p->named;tmp!=NULL;tmp=tmp->next)
    {
      if( !strcmp(tmp->name, name) )
        {
          removed=tmp;
          if(prev) prev->next = tmp->next;
          else     p->named   = tmp->next;
        }

      /* Set this node as the 'prev' pointer. */
      prev=tmp;
    }

  /* A small sanity check. */
  if(removed==NULL)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
          "problem. 'removed' must not be NULL at this point", __func__,
          PACKAGE_BUGREPORT);

  /* Nothing in the list points to it now. So we can safely modify and
     return it. */
  free(removed->name);
  removed->next=NULL;
  removed->name=NULL;
  return removed;
}





/* Pop a dataset and keep it in the 'named' list for later use. */
void
operands_set_name(struct arithmeticparams *p, char *token)
{
  gal_data_t *tmp, *tofree;
  char *varname=&token[ OPERATOR_PREFIX_LENGTH_SET ];

  /* If a dataset with this name already exists, it will be removed/deleted
     so we can use the name for the newly designated dataset. */
  for(tmp=p->named; tmp!=NULL; tmp=tmp->next)
    if( !strcmp(varname, tmp->name) )
      {
        tofree=operands_remove_name(p, varname);
        gal_data_free(tofree);

        /* IMPORTANT: we MUST break here! 'tmp' does't point to the right
           place any more. We can define a 'prev' node and modify it on
           every attempt, but since there is only one dataset with a given
           name, that is redundant and will just make the program slow. */
        break;
      }

  /* Pop the top operand, then add it to the list of named datasets, but
     only if it is used in later tokens. If it isn't, free the popped
     dataset. The latter case (to define a name, but not use it), is
     obviously a redundant operation, but that is upto the user, we
     shouldn't worry about it here. We should just have everything in
     place, so no crashes occur or no extra memory is consumed. */
  if( operands_name_is_used_later(p, varname) )
    {
      /* Add the top popped operand to the list of names. */
      gal_list_data_add(&p->named, operands_pop(p, "set"));

      /* Write the requested name into this dataset. But note that 'name'
         MUST be already empty. So to be safe, we'll do a sanity check. */
      if(p->named->name)
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
              "the problem. The 'name' element should be NULL at this "
              "point, but it isn't", __func__, PACKAGE_BUGREPORT);
      gal_checkset_allocate_copy(varname, &p->named->name);
    }
  else
    {
      /* Pop the top operand, then free it. */
      tmp=operands_pop(p, "set");
      gal_data_free(tmp);
    }
}





/* See if a given token is the name of a variable. */
int
operands_is_name(struct arithmeticparams *p, char *token)
{
  gal_data_t *tmp;

  /* Make sure the variable name hasn't been set before. */
  for(tmp=p->named; tmp!=NULL; tmp=tmp->next)
    if( !strcmp(token, tmp->name) )
      return 1;

  /* If control reaches here, then there was no match*/
  return 0;
}





/* Return a copy of the named dataset. */
static gal_data_t *
operands_copy_named(struct arithmeticparams *p, char *name)
{
  gal_data_t *out=NULL, *tmp;

  /* Find the proper named element to use. */
  for(tmp=p->named;tmp!=NULL;tmp=tmp->next)
    if( !strcmp(tmp->name, name) )
      {
        /* If the named operand is used later, then copy it into the
           output. */
        if( operands_name_is_used_later(p, name) )
          {
            out=gal_data_copy(tmp);
            free(out->name);
            out->name=NULL;
            out->next=NULL;
          }
        /* The named operand is not used any more. Remove it from the list
           of named datasets and continue. */
        else out=operands_remove_name(p, name);
      }

  /* A small sanity check. */
  if(out==NULL)
    error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to fix the "
          "problem. The requested name '%s' couldn't be found in the list",
          __func__, PACKAGE_BUGREPORT, name);

  /* Return. */
  return out;
}




















/**********************************************************************/
/************      Adding to and popping from stack     ***************/
/**********************************************************************/
void
operands_add(struct arithmeticparams *p, char *filename, gal_data_t *data)
{
  int readwcs;
  size_t ndim, *dsize;
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
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'newnode'",
              __func__, sizeof *newnode);

      /* If the 'filename' is the name of a dataset, then use a copy of it.
         otherwise, do the basic analysis. */
      if( filename && operands_is_name(p, filename) )
        {
          newnode->filename=NULL;
          newnode->data=operands_copy_named(p, filename);
        }
      else
        {
          /* Set the basic parameters. */
          newnode->data=data;
          newnode->filename=filename;

          /* See if a HDU must be read or not. */
          if(filename != NULL
             && ( gal_fits_name_is_fits(filename)
                  || gal_tiff_name_is_tiff(filename) ) )
            {
              /* Set the HDU for this filename. */
              if(p->globalhdu)
                gal_checkset_allocate_copy(p->globalhdu, &newnode->hdu);
              else
                newnode->hdu=gal_list_str_pop(&p->hdus);

              /* If no WCS is set yet, use the WCS of this image and remove
                 possibly extra dimensions if necessary. */
              readwcs = (p->wcsfile && !strcmp(p->wcsfile,"none")) ? 0 : 1;
              if(readwcs && p->refdata.wcs==NULL)
                {
                  dsize=gal_fits_img_info_dim(filename, newnode->hdu, &ndim);
                  p->refdata.wcs=gal_wcs_read(filename, newnode->hdu, 0, 0,
                                              &p->refdata.nwcs);
                  ndim=gal_dimension_remove_extra(ndim, dsize, p->refdata.wcs);
                  if(p->refdata.wcs && !p->cp.quiet)
                    printf(" - WCS: %s (hdu %s).\n", filename, newnode->hdu);
                  free(dsize);
                }
            }
          else newnode->hdu=NULL;
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
    error(EXIT_FAILURE, 0, "not enough operands for the '%s' operator",
          operator);

  /* Set the dataset. If filename is present then read the file
     and fill in the array, if not then just set the array. */
  if(operands->filename)
    {
      /* Set the HDU and filename */
      hdu=operands->hdu;
      filename=operands->filename;

      /* Read the dataset and remove possibly extra dimensions. */
      data=gal_array_read_one_ch(filename, hdu, NULL, p->cp.minmapsize,
                                 p->cp.quietmmap);
      data->ndim=gal_dimension_remove_extra(data->ndim, data->dsize, NULL);

      /* Arithmetic changes the contents of a dataset, so the existing name
         (in the FITS 'EXTNAME' keyword) should not be passed on beyond
         this point. Also, in Arithmetic, the 'name' element is used to
         identify variables. */
      if(data->name) { free(data->name); data->name=NULL; }

      /* When the reference data structure's dimensionality is non-zero, it
         means that this is not the first image read. So, write its basic
         information into the reference data structure for future
         checks. */
      if(p->refdata.ndim==0)
        {
          /* Set the dimensionality. */
          p->refdata.ndim=(data)->ndim;

          /* Allocate the dsize array. */
          errno=0;
          p->refdata.dsize=malloc(p->refdata.ndim
                                  * sizeof *p->refdata.dsize);
          if(p->refdata.dsize==NULL)
            error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for "
                  "p->refdata.dsize", __func__,
                  p->refdata.ndim * sizeof *p->refdata.dsize);

          /* Write the values into it. */
          for(i=0;i<p->refdata.ndim;++i)
            p->refdata.dsize[i]=data->dsize[i];
        }

      /* Report the read image if desired: */
      if(!p->cp.quiet) printf(" - Read: %s (hdu %s).\n", filename, hdu);

      /* Free the HDU string: */
      if(hdu) free(hdu);

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
