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

#include "checkset.h"
#include "arraymanip.h"
#include "fitsarrayvv.h"

#include "main.h"
#include "imgarith.h"            /* needs main.h.                  */




/***************************************************************/
/*************    Operand linked list functions    *************/
/***************************************************************/
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

  if(strlen(filename) && nameisfits(filename))
    {
      /* Set the HDU for this filename. */
      newnode->hdu=p->up.hdus[p->addcounter];

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
    error(EXIT_FAILURE, 0, "Not enough operands for the \"%s\" operator.",
          operator);


  /* Do a sanity check. The basic idea behind this is that all the
     conditionals below will evaluate to 1 or 0. So if more than one
     of them are true, then the sum will be larger than 1 and if none
     of them are true then the sum will be 0. So if the sum (check) is
     not equal to 1, then there is a bug and the user should be
     warned. The parenthesis will help in avoiding compiler
     warnings.*/
  if( (strlen(operands->filename)>0) + !(isnan(operands->number))
      + (operands->array!=NOOPTARRAY) != 1)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix the "
          "problem. For some reason, one node in the operands linked list "
          "has more than one value.", PACKAGE_BUGREPORT);


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
          readfitswcs(filename, hdu, 0, 0, &p->nwcs, &p->wcs);
        }
      filetodouble(filename, maskname, hdu, mhdu, array, &bitpix,
                   &p->anyblank, &s0, &s1);

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
                  "ImageArithmetic to work.", filename, hdu, s0, s1,
                  p->s0, p->s1);
        }

      /* Set the bitpix of the output. */
      if(bitpix==DOUBLE_IMG) p->obitpix=DOUBLE_IMG;

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
      dsumarrays(farr, sarr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);

      /* Clean up. */
      free(sarr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      dsumconst(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      dsumconst(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, fnum+snum, NOOPTARRAY);
}





void
subtract(struct imgarithparams *p)
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
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      dsubtractarrays(sarr, farr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);

      /* Clean up. */
      free(farr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      dconstsubtract(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      dsubtractconst(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, snum-fnum, NOOPTARRAY);
}





void
multiply(struct imgarithparams *p)
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
      /* Do the operation, note that the output is stored in farr. */
      dmultiparrays(farr, sarr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);

      /* Clean up. */
      free(sarr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      dmultipconst(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      dmultipconst(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, fnum*snum, NOOPTARRAY);
}





void
divide(struct imgarithparams *p)
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
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      ddividearrays(sarr, farr, size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);

      /* Clean up. */
      free(farr);
    }
  else if(farr)                 /* Only the first is an array. */
    {
      dconstdivide(farr, size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, farr);
    }
  else if(sarr)                 /* Only the first is an array. */
    {
      ddivideconst(sarr, size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, sarr);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, snum/fnum, NOOPTARRAY);
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
  struct stll *token;

  /* Prepare the processing: */
  p->s0=p->s1=0;
  p->operands=NULL;
  p->obitpix=FLOAT_IMG;
  p->addcounter=p->popcounter=0;

  /* Go over each input token and do the work. */
  for(token=p->tokens;token!=NULL;token=token->next)
    {
      /* If we have a name or number, then add it to the operands
         linked list. Otherwise, pull out two members and do the
         specified operation on them. */
      if(nameisfits(token->v))
        add_operand(p, token->v, NOOPTNUMBER, NOOPTARRAY);
      else if(strisdouble(token->v, &number))
        add_operand(p, NOOPTFILENAME, number, NOOPTARRAY);
      else
        {
          /* Currently all operators are a single character length, we
             don't have bitwise operators or functions yet, so to make
             things easy, we just do a sanity check here */
          if(strlen(token->v)!=1)
            error(EXIT_FAILURE, 0, "The argument \"%s\" could not be read "
                  "as a FITS file, number or an operator.", token->v);

          switch(token->v[0])
            {
            case '+':
              sum(p);
              break;
            case '-':
              subtract(p);
              break;
            case '*':
              multiply(p);
              break;
            case '/':
              divide(p);
              break;
            default:
              error(EXIT_FAILURE, 0, "The argument \"%s\" could not be "
                    "interpretted as an operator.", token->v);
            }
        }
    }


  /* If there is more than one node in the operands stack, then the
     user has given too many operands and there is an error. */
  if(p->operands->next!=NULL)
    error(EXIT_FAILURE, 0, "There are too many operands for the operators "
          "in the given expression.");


  /* There is only one operand left and that is the output. So write
     it to a file and finish. If none of the inputs had a double type,
     then convert the output array into a float and then save it. Note
     that the last operand must be an array. */
  if(p->obitpix==FLOAT_IMG)
    {
      changetype(p->operands->array, DOUBLE_IMG, p->s0*p->s1, p->anyblank,
                 (void **)(&farray), FLOAT_IMG);
      arraytofitsimg(p->cp.output, "astimgarith", FLOAT_IMG, farray,
                     p->s0, p->s1, p->anyblank, p->wcs, NULL, SPACK_STRING);
    }
  else
    arraytofitsimg(p->cp.output, "astimgarith", DOUBLE_IMG,
                   p->operands->array, p->s0, p->s1, p->anyblank,
                   p->wcs, NULL, SPACK_STRING);
}



















/***************************************************************/
/*************             Top function            *************/
/***************************************************************/
void
imgarith(struct imgarithparams *p)
{
  reversepolish(p);
}
