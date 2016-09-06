/*********************************************************************
ImageArithmetic - Do arithmetic operations on images.
ImageArithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro/fits.h>
#include <gnuastro/checkset.h>
#include <gnuastro/arraymanip.h>
#include <gnuastro/statistics.h>

#include "main.h"
#include "arithmetic.h"            /* needs main.h.                  */




/***************************************************************/
/*************    Operand linked list functions    *************/
/***************************************************************/
size_t
num_operands(struct imgarithparams *p)
{
  size_t counter=0;
  struct operand *tmp=NULL;
  for(tmp=p->operands;tmp!=NULL;tmp=tmp->next)
    ++counter;
  return counter;
}





void
add_operand(struct imgarithparams *p, char *filename, double number,
            double *array)
{
  struct operand *newnode;

  /* Allocate space for the new operand. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "imgarith.c: Making new element in operand");

  /* Fill in the values. */
  newnode->array=array;
  newnode->number=number;
  newnode->filename=filename;

  if(filename != NULL && gal_fits_name_is_fits(filename))
    {
      /* Set the HDU for this filename. */
      gal_linkedlist_pop_from_stll(&p->hdus, &newnode->hdu);

      /* Increment the FITS counter. */
      ++p->addcounter;
    }

  /* Make the link to the previous list. */
  newnode->next=p->operands;
  p->operands=newnode;
}





void
pop_operand(struct imgarithparams *p, double *number, double **array,
            char *operator)
{
  int bitpix;
  size_t s0, s1;
  struct uiparams *up=&p->up;
  struct operand *operands=p->operands;
  char *maskname, *mhdu, *filename, *hdu;

  /* If the operand linked list has finished, then give an error and
     exit. */
  if(operands==NULL)
    error(EXIT_FAILURE, 0, "not enough operands for the \"%s\" operator",
          operator);


  /* Set the array output. If filename is present then read the file
     and fill in the array, if not then just set the array. */
  if(strlen(operands->filename))
    {
      hdu=operands->hdu;
      filename=operands->filename;

      /* In case this is the first image that is read, then read the
         WCS information and set the mask name so masked pixels can be
         set to NaN. For the other images, the mask can be completely
         ignored. */
      if(p->popcounter)         /* This is not the first FITS file. */
        {
          maskname=NULL;
          mhdu=NULL;
        }
      else
        {
          mhdu=up->mhdu;
          maskname=up->maskname;
          gal_fits_read_wcs(filename, hdu, 0, 0, &p->nwcs, &p->wcs);
        }
      gal_fits_file_to_double(filename, maskname, hdu, mhdu,
                                   array, &bitpix, &p->anyblank, &s0, &s1);

      /* If the output size was not set yet, then set it. Otherwise,
         make sure the size of this image is the same as the previous
         images. */
      if(p->s0==0 && p->s1==0)
        {
          p->s0=s0;
          p->s1=s1;
        }
      else
        {
          if(p->s0!=s0 || p->s1!=s1)
            error(EXIT_FAILURE, 0, "%s (hdu=%s): has size of %lu x %lu. "
                  "However, previous images had a size of %lu x %lu. All "
                  "the images must be the same size in order for "
                  "ImageArithmetic to work", filename, hdu, s0, s1,
                  p->s0, p->s1);
        }

      /* Free the HDU string: */
      free(hdu);

      /* Add to the number of popped FITS images: */
      ++p->popcounter;
    }
  else
    *array=operands->array;

  /* Set the number: */
  *number=operands->number;

  /* Remove this node from the queue. */
  p->operands=operands->next;
  free(operands);
}




















/***************************************************************/
/*************              Operators              *************/
/***************************************************************/
void
sum(struct imgarithparams *p)
{
  size_t size;
  char *operator="+";
  double fnum, snum;            /* First or second number.    */
  double *farr, *sarr;          /* First or second array.     */

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);
  pop_operand(p, &snum, &sarr, operator);

  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;

  /* Do the operation: */
  if(farr && sarr)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in the first
         input. */
      gal_arraymanip_dsum_arrays(farr, sarr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);

      /* Clean up. */
      free(sarr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dsum_const(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dsum_const(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, fnum+snum, NOOPTARRAY);
}





void
subtract(struct imgarithparams *p)
{
  size_t size;
  char *operator="-";
  double fnum, snum;            /* First or second number.    */
  double *farr, *sarr;          /* First or second array.     */

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);
  pop_operand(p, &snum, &sarr, operator);

  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;

  /* Do the operation: */
  if(farr && sarr)              /* Both are arrays. */
    {

      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_arraymanip_dsubtract_arrays(sarr, farr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);

      /* Clean up. */
      free(farr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dconst_subtract(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dsubtract_const(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, snum-fnum, NOOPTARRAY);
}





void
multiply(struct imgarithparams *p)
{
  size_t size;
  char *operator="*";
  double fnum, snum;            /* First or second number.    */
  double *farr, *sarr;          /* First or second array.     */

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);
  pop_operand(p, &snum, &sarr, operator);

  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;

  /* Do the operation: */
  if(farr && sarr)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in farr. */
      gal_arraymanip_dmultip_arrays(farr, sarr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);

      /* Clean up. */
      free(sarr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dmultip_const(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dmultip_const(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, fnum*snum, NOOPTARRAY);
}





void
divide(struct imgarithparams *p)
{
  size_t size;
  char *operator="/";
  double fnum, snum;            /* First or second number.    */
  double *farr, *sarr;          /* First or second array.     */

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);
  pop_operand(p, &snum, &sarr, operator);

  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;

  /* Do the operation: */
  if(farr && sarr)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_arraymanip_ddivide_arrays(sarr, farr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);

      /* Clean up. */
      free(farr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dconst_divide(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      gal_arraymanip_ddivide_const(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, snum/fnum, NOOPTARRAY);
}





void
topower(struct imgarithparams *p, char *op)
{
  size_t size;
  char *operator = op?op:"pow";
  double fnum, snum;            /* First or second number.    */
  double *farr, *sarr;          /* First or second array.     */

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);
  pop_operand(p, &snum, &sarr, operator);

  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;

  /* Do the operation: */
  if(farr && sarr)              /* Both are arrays. */
    {

      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_arraymanip_dpower_arrays(sarr, farr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);


      /* Clean up. */
      free(farr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dconst_power(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      gal_arraymanip_dpower_const(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, pow(snum, fnum), NOOPTARRAY);
}





void
alloppixs(struct imgarithparams *p, char *operator)
{
  int i, j;                  /* Integer to allow negative checks.     */
  double num;                /* Temporary number holder.              */
  double *arr;               /* Temporary array holder.               */
  int firstisarray=-1;       /* ==-1: unset. ==1: array. ==0: number. */
  double *allpixels=NULL;    /* Array for all values in one pixel.    */
  double **allarrays=NULL;   /* Array for pointers to input arrays.   */
  size_t size, numop=num_operands(p);
  double (*thisfunction)(double *, size_t)=NULL;


  /* First set the appropriate function to call. */
  if(!strcmp(operator, "min"))
    thisfunction = &gal_statistics_double_min_return;
  else if(!strcmp(operator, "max"))
    thisfunction = &gal_statistics_double_max_return;
  else if(!strcmp(operator, "median"))
    thisfunction = &gal_statistics_median_double_in_place;
  else if(!strcmp(operator, "average"))
    thisfunction = &gal_statistics_double_average;
  else
    error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we "
          "can address the problem. The value of `operator' in "
          "alloppixs (%s) is not recognized",
          PACKAGE_BUGREPORT, operator);


  /* Allocate the array of pointers to all input arrays and also the
     array to temporarily keep all values for each pixel */
  errno=0;
  allarrays=malloc(numop*sizeof *allarrays);
  if(allarrays==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for allarrays in alloppixs",
          numop*sizeof *allarrays);
  errno=0;
  allpixels=malloc(numop*sizeof *allpixels);
  if(allpixels==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for allpixels in alloppixs",
          numop*sizeof *allpixels);


  /* Prepare all the inputs, note that since it is a linked list, the
     operands pop from the last to first. Here order is not important,
     but in other cases that it might be (and also for debugging),
     here the inputs are put in order so the indexs start from the
     last to first. */
  for(i=numop-1;i>=0;--i)
    {
      /* Pop out the operand. */
      pop_operand(p, &num, &arr, operator);

      /* Do the appropriate action if it is an array or a number. */
      if(arr)
        switch(firstisarray)
          {
          case -1:
            firstisarray=1;
            allarrays[i]=arr;
            break;
          case 1:
            allarrays[i]=arr;
            break;
          case 0:
            error(EXIT_FAILURE, 0, "for the %s operator, all operands "
                  "must be either an array or number", operator);
            break;
          default:
            error(EXIT_FAILURE, 0, "a Bug! Please contact us at %s so we "
                  "can address the problem. The value of firstisarray (%d) "
                  "in the alloppixs function is not recognized",
                  PACKAGE_BUGREPORT, firstisarray);
          }
      else
        switch(firstisarray)
          {
          case -1:
            firstisarray=0;
            allpixels[i]=num;
            break;
          case 0:
            allpixels[i]=num;
            break;
          case 1:
            error(EXIT_FAILURE, 0, "for the %s operator, all operands "
                  "must be either an array or number", operator);
            break;
          default:
            error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we "
                  "can address the problem. The value of firstisarray (%d) "
                  "in the alloppixs function is not recognized",
                  PACKAGE_BUGREPORT, firstisarray);
          }
    }


  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;


  /* Do the operation and report the result: */
  if(arr)
    {
      /* Find the value and replace it with the first operand. */
      for(i=0;i<size;++i)
        {
          /* Go over all the inputs and put the appropriate pixel
             values in the allpixels array. */
          for(j=0;j<numop;++j)
            allpixels[j]=allarrays[j][i];

          /* Do the appropriate action. */
          allarrays[0][i]=(*thisfunction)(allpixels, numop);
        }

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, allarrays[0]);

      /* Free all the extra operands */
      for(i=1;i<numop;++i) free(allarrays[i]);
    }
  else
    add_operand(p, NOOPTFILENAME, (*thisfunction)(allpixels, numop),
                NOOPTARRAY);

  /* Clean up: */
  free(allarrays);
  free(allpixels);
}





void
takesqrt(struct imgarithparams *p)
{
  char *operator="sqrt";

  /* Add a 0.5 number to the operand stack */
  add_operand(p, NOOPTFILENAME, 0.5f, NOOPTARRAY);

  /* Call the power operator. */
  topower(p, operator);
}





void
takelog(struct imgarithparams *p)
{
  char *operator="log";
  double fnum, *farr;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);

  /* Do the operation: */
  if(farr)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_arraymanip_dlog_array(farr, p->s0*p->s1);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, log(fnum), NOOPTARRAY);
}





void
takelog10(struct imgarithparams *p)
{
  char *operator="log10";
  double fnum, *farr;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);

  /* Do the operation: */
  if(farr)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_arraymanip_dlog10_array(farr, p->s0*p->s1);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, log10(fnum), NOOPTARRAY);
}





void
takeabs(struct imgarithparams *p)
{
  char *operator="abs";
  double fnum, *farr;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);

  /* Do the operation: */
  if(farr)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_arraymanip_dabs_array(farr, p->s0*p->s1);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, fabs(fnum), NOOPTARRAY);
}





void
findmin(struct imgarithparams *p)
{
  char *operator="min";
  double min, fnum, *farr;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);

  /* Do the operation: */
  if(farr)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_statistics_double_min(farr, p->s0*p->s1, &min);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, min, NOOPTARRAY);

      /* Clean up. */
      free(farr);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, fnum, NOOPTARRAY);
}





void
findmax(struct imgarithparams *p)
{
  char *operator="max";
  double max, fnum, *farr;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);

  /* Do the operation: */
  if(farr)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_statistics_double_max(farr, p->s0*p->s1, &max);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, max, NOOPTARRAY);

      /* Clean up. */
      free(farr);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, fnum, NOOPTARRAY);
}





int
lessthan(double left, double right)
{ return left<right; }

int
lessequal(double left, double right)
{ return left<=right; }

int
greaterthan(double left, double right)
{ return left>right; }

int
greaterequal(double left, double right)
{ return left>=right; }

int
equal(double left, double right)
{ return left==right; }

int
notequal(double left, double right)
{ return left!=right; }





void
conditionals(struct imgarithparams *p, char *operator)
{
  size_t size;
  double fnum, snum;            /* First or second number.    */
  double *farr, *sarr;          /* First or second array.     */
  double *f, *s, *ff, *ss;
  int (*thisfunction)(double, double)=NULL;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);
  pop_operand(p, &snum, &sarr, operator);

  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;

  if(!strcmp(operator, "lt"))       thisfunction = &lessthan;
  else if(!strcmp(operator, "le"))  thisfunction = &lessequal;
  else if(!strcmp(operator, "gt"))  thisfunction = &greaterthan;
  else if(!strcmp(operator, "ge"))  thisfunction = &greaterequal;
  else if(!strcmp(operator, "eq"))  thisfunction = &equal;
  else if(!strcmp(operator, "neq")) thisfunction = &notequal;
  else
    error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we "
          "can address the problem. The value of `operator' in "
          "conditionals (%s) is not recognized",
          PACKAGE_BUGREPORT, operator);

  /* Do the operation: */
  if(farr && sarr)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      f=farr;
      ss=(s=sarr)+size;
      do *s = thisfunction(*s, *f++); while(++s<ss);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);

      /* Clean up. */
      free(farr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      ff=(f=farr)+size;
      do *f = thisfunction(snum, *f); while(++f<ff);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      ss=(s=sarr)+size;
      do *s = thisfunction(*s, fnum); while(++s<ss);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, thisfunction(snum, fnum), NOOPTARRAY);
}





/* In order to not conflict with the internal C `is...' functions, and in
   particular the `isblank' function, we are calling this function
   opisblank for operator-isblank. */
void
opisblank(struct imgarithparams *p)
{
  size_t size=p->s0*p->s1;
  char *operator="isblank";
  double *f, *ff, fnum, *farr;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, operator);

  /* Do the operation: */
  if(farr)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      ff=(f=farr)+size;
      do *f = isnan(*f); while(++f<ff);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, isnan(fnum), NOOPTARRAY);
}





/* Replace the pixels in the second popped element with the first. While
   choosing the pixels that are selected from the third The third popped
   element. The third is considered to be an array that can only be filled
   with 0 or 1. */
void
where(struct imgarithparams *p)
{
  size_t size;
  double *f, *s, *t, *ss;
  double fnum, snum, tnum;      /* First, second, or third number.    */
  double *farr, *sarr, *tarr;   /* First, second, or third array.     */

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &farr, "where");
  pop_operand(p, &snum, &sarr, "where");
  pop_operand(p, &tnum, &tarr, "where");

  /* Set the total number of pixels, note that we can't do this in the
     definition of the variable because p->s0 and p->s1 will be set in
     pop_operand for the first image. */
  size=p->s0*p->s1;

  /* Do the operation: */
  if(sarr && tarr)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      t=tarr;
      ss=(s=sarr)+size;
      if(farr)
        {
          printf("\nhere\n");
          f=farr;
          /* Note that we need to increment "f" after the check and
             replace. If the increment is inside the conditional replace,
             then when t==0, the increment won't work.*/
          do {*s = *t++ ? *f : *s; ++f;} while(++s<ss);
        }
      else
        do *s = *t++ ? fnum : *s; while(++s<ss);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);

      /* Clean up. */
      free(farr);
      free(tarr);
    }
  else
    error(EXIT_FAILURE, 0, "the first and second arguments (second and "
          "third popped elements) to `where' have to be arrays.");
}


















/***************************************************************/
/*************      Reverse Polish algorithm       *************/
/***************************************************************/
/* This function implements the reverse polish algorithm as explained
   in the Wikipedia page.

   NOTE that in ui.c, the input linked list of tokens was ordered to
   have the same order as what the user provided. */
void
reversepolish(struct imgarithparams *p)
{
  float *farray;
  double number;
  char *tokeepvalue;
  struct gal_linkedlist_stll *token;

  /* Prepare the processing: */
  p->s0=p->s1=0;
  p->operands=NULL;
  p->addcounter=p->popcounter=0;

  /* Go over each input token and do the work. */
  for(token=p->tokens;token!=NULL;token=token->next)
    {
      /* If we have a name or number, then add it to the operands
         linked list. Otherwise, pull out two members and do the
         specified operation on them. */
      if(gal_fits_name_is_fits(token->v))
        add_operand(p, token->v, NOOPTNUMBER, NOOPTARRAY);
      else if(strisdouble(token->v, &number))
        add_operand(p, NOOPTFILENAME, number, NOOPTARRAY);
      else
        {
          if     (!strcmp(token->v, "+"))       sum(p);
          else if(!strcmp(token->v, "-"))       subtract(p);
          else if(!strcmp(token->v, "*"))       multiply(p);
          else if(!strcmp(token->v, "/"))       divide(p);
          else if(!strcmp(token->v, "abs"))     takeabs(p);
          else if(!strcmp(token->v, "pow"))     topower(p, NULL);
          else if(!strcmp(token->v, "sqrt"))    takesqrt(p);
          else if(!strcmp(token->v, "log"))     takelog(p);
          else if(!strcmp(token->v, "log10"))   takelog10(p);
          else if(!strcmp(token->v, "minvalue"))findmin(p);
          else if(!strcmp(token->v, "maxvalue"))findmax(p);
          else if(!strcmp(token->v, "min")
                  || !strcmp(token->v, "max")
                  || !strcmp(token->v, "average")
                  || !strcmp(token->v, "median")) alloppixs(p, token->v);
          else if(!strcmp(token->v, "lt")
                  || !strcmp(token->v, "le")
                  || !strcmp(token->v, "gt")
                  || !strcmp(token->v, "ge")
                  || !strcmp(token->v, "eq")
                  || !strcmp(token->v, "neq")) conditionals(p, token->v);
          else if(!strcmp(token->v, "isblank")) opisblank(p);
          else if(!strcmp(token->v, "where")) where(p);
          else
            error(EXIT_FAILURE, 0, "the argument \"%s\" could not be "
                  "interpretted as an operator", token->v);
        }
    }

  /* If there is more than one node in the operands stack, then the
     user has given too many operands and there is an error. */
  if(p->operands->next!=NULL)
    error(EXIT_FAILURE, 0, "there are too many operands for the operators "
          "in the given expression");


  /* If the remaining operand is an array then save the array as a
     FITS image, if not, simply print the floating point number. */
  if(p->operands->array)
    {
      /* Internally, all arrays were double type. But the user could set
         the the output type using the type option. So if the user has
         asked for anything other than a double, we will have to convert
         the arrays.*/
      if(p->type==DOUBLE_IMG)
        gal_fits_array_to_file(p->cp.output, "astimgarith",
                               DOUBLE_IMG, p->operands->array,
                               p->s0, p->s1, p->anyblank,
                               p->wcs, NULL, SPACK_STRING);
      else
        {
          gal_fits_change_type(p->operands->array, DOUBLE_IMG,
                               p->s0*p->s1, p->anyblank,
                               (void **)(&farray), p->type);
          gal_fits_array_to_file(p->cp.output, "astimgarith",
                                 p->type, farray, p->s0, p->s1,
                                 p->anyblank, p->wcs, NULL,
                                 SPACK_STRING);
        }

    }
  else
    printf("%g\n", p->operands->number);


  /* If there are any remaining HDUs in the hdus linked list, then
     free them. */
  while(p->hdus!=NULL)
    {
      gal_linkedlist_pop_from_stll(&p->hdus, &tokeepvalue);
      free(tokeepvalue);
    }
}



















/***************************************************************/
/*************             Top function            *************/
/***************************************************************/
void
imgarith(struct imgarithparams *p)
{
  reversepolish(p);
}
