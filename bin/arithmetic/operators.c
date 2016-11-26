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

#include <gnuastro/data.h>
#include <gnuastro/array.h>
#include <gnuastro/statistics.h>

#include "main.h"

#include "operands.h"
#include "operators.h"










void
sum(struct imgarithparams *p)
{
  double fnum, snum;            /* First or second number.    */
  gal_data_t *f, *s;            /* First or second dataset.   */
  char *operator="+";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);
  pop_operand(p, &snum, &s, operator);

  /* Do the operation: */
  if(f && s)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in the first
         input. */
      gal_array_dsum_arrays(f->array, s->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);

      /* Clean up. */
      gal_data_free(s);
    }
  else if(f)                 /* Only the first is an array. */
    {
      gal_array_dsum_const(f->array, f->size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else if(s)                 /* Only the second is an array. */
    {
      gal_array_dsum_const(s->array, s->size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, fnum+snum, NOOPTDATA);
}





void
subtract(struct imgarithparams *p)
{
  double fnum, snum;            /* First or second number.    */
  gal_data_t *f, *s;            /* First or second dataset.   */
  char *operator="-";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);
  pop_operand(p, &snum, &s, operator);

  /* Do the operation: */
  if(f && s)              /* Both are arrays. */
    {

      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_array_dsubtract_arrays(s->array, f->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s->array);

      /* Clean up. */
      gal_data_free(f);
    }
  else if(f)                 /* Only the first is an array. */
    {
      gal_array_dconst_subtract(f->array, f->size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f->array);
    }
  else if(s)                 /* Only the second is an array. */
    {
      gal_array_dsubtract_const(s->array, s->size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);
    }
  else                       /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, snum-fnum, NOOPTDATA);
}





void
multiply(struct imgarithparams *p)
{
  double fnum, snum;            /* First or second number.    */
  gal_data_t *f, *s;            /* First or second dataset.   */
  char *operator="*";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);
  pop_operand(p, &snum, &s, operator);

  /* Do the operation: */
  if(f && s)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in farr. */
      gal_array_dmultip_arrays(f->array, s->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);

      /* Clean up. */
      gal_data_free(s);
    }
  else if(f)                 /* Only the first is an array. */
    {
      gal_array_dmultip_const(f->array, f->size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else if(s)                 /* Only the first is an array. */
    {
      gal_array_dmultip_const(s->array, s->size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, fnum*snum, NOOPTDATA);
}





void
divide(struct imgarithparams *p)
{
  double fnum, snum;            /* First or second number.    */
  gal_data_t *f, *s;            /* First or second dataset.   */
  char *operator="/";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);
  pop_operand(p, &snum, &s, operator);

  /* Do the operation: */
  if(f && s)              /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_array_ddivide_arrays(s->array, f->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);

      /* Clean up. */
      gal_data_free(f);
    }
  else if(f)                 /* Only the first is an array. */
    {
      gal_array_dconst_divide(f->array, f->size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else if(s)                 /* Only the first is an array. */
    {
      gal_array_ddivide_const(s->array, s->size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, snum/fnum, NOOPTDATA);
}





void
topower(struct imgarithparams *p, char *op)
{
  double fnum, snum;            /* First or second number.    */
  gal_data_t *f, *s;            /* First or second dataset.   */
  char *operator = op?op:"pow";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);
  pop_operand(p, &snum, &s, operator);

  /* Do the operation: */
  if(f && s)              /* Both are arrays. */
    {

      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_array_dpower_arrays(s->array, f->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);

      /* Clean up. */
      gal_data_free(f);
    }
  else if(f)                 /* Only the first is an array. */
    {
      gal_array_dconst_power(f->array, f->size, snum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else if(s)                 /* Only the first is an array. */
    {
      gal_array_dpower_const(s->array, s->size, fnum);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, pow(snum, fnum), NOOPTDATA);
}





void
alloppixs(struct imgarithparams *p, char *operator)
{
  int i, j;                  /* Integer to allow negative checks.     */
  double num;                /* Temporary number holder.              */
  int firstisdata=-1;        /* ==-1: unset. ==1: array. ==0: number. */
  double *allpixels=NULL;    /* Array for all values in one pixel.    */
  double **allarrays=NULL;   /* Pointer to all the arrays.            */
  gal_data_t **alldata=NULL; /* Array for pointers to input datasets. */
  size_t numop=num_operands(p);
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
  alldata=malloc(numop*sizeof *alldata);
  if(alldata==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for alldata in alloppixs",
          numop*sizeof *alldata);
  errno=0;
  allarrays=malloc(numop*sizeof *allarrays);
  if(allarrays==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for allarrays in alloppixs",
          numop*sizeof *allarrays);
  errno=0;
  allpixels=malloc(numop*sizeof *allpixels);
  if(allpixels==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for allpixels in alloppixs",
          numop*sizeof *allpixels);


  /* Prepare all the inputs, note that since it is a linked list, the
     operands pop from the last to first. Here order is not important,
     but in other cases that it might be (and also for debugging),
     here the inputs are put in order so the indexs start from the
     last to first. */
  for(i=numop-1;i>=0;--i)
    {
      /* Pop out the operand. */
      pop_operand(p, &num, &alldata[i], operator);

      /* Do the appropriate action if it is an array or a number. */
      if(alldata[i])
        switch(firstisdata)
          {
          case -1:
            firstisdata=1;
            allarrays[i]=alldata[i]->array;
            break;
          case 1:
            allarrays[i]=alldata[i]->array;
            break;
          case 0:
            error(EXIT_FAILURE, 0, "for the %s operator, all operands "
                  "must be either an array or number", operator);
            break;
          default:
            error(EXIT_FAILURE, 0, "a Bug! Please contact us at %s so we "
                  "can address the problem. The value of firstisdata (%d) "
                  "in the alloppixs function is not recognized",
                  PACKAGE_BUGREPORT, firstisdata);
          }
      else
        switch(firstisdata)
          {
          case -1:
            firstisdata=0;
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
                  "can address the problem. The value of firstisdata (%d) "
                  "in the alloppixs function is not recognized",
                  PACKAGE_BUGREPORT, firstisdata);
          }
    }


  /* Do the operation and report the result: */
  if(firstisdata==1)
    {
      /* For each pixel, find the value and replace it with the first
         operand. */
      for(i=0;i<alldata[0]->size;++i)
        {
          /* Go over all the inputs and put the appropriate pixel
             values in the allpixels array. */
          for(j=0;j<numop;++j)
            allpixels[j]=allarrays[j][i];

          /* Do the appropriate action. */
          allarrays[0][i]=(*thisfunction)(allpixels, numop);
        }

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, alldata[0]);

      /* Free all the extra operands */
      for(i=1;i<numop;++i) gal_data_free(alldata[i]);
    }
  else
    add_operand(p, NOOPTFILENAME, (*thisfunction)(allpixels, numop),
                NOOPTDATA);


  /* Clean up. Reminder: the allocated arrays were freed above as part of
     `gal_data_free'. */
  free(alldata);
  free(allarrays);
  free(allpixels);
}





void
takesqrt(struct imgarithparams *p)
{
  char *operator="sqrt";

  /* Add a 0.5 number to the operand stack */
  add_operand(p, NOOPTFILENAME, 0.5f, NOOPTDATA);

  /* Call the power operator. */
  topower(p, operator);
}





void
takelog(struct imgarithparams *p)
{
  double fnum;
  gal_data_t *f;
  char *operator="log";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);

  /* Do the operation: */
  if(f)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_array_dlog_array(f->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, log(fnum), NOOPTDATA);
}





void
takelog10(struct imgarithparams *p)
{
  double fnum;
  gal_data_t *f;
  char *operator="log10";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);

  /* Do the operation: */
  if(f)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_array_dlog10_array(f->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, log10(fnum), NOOPTDATA);
}





void
takeabs(struct imgarithparams *p)
{
  double fnum;
  gal_data_t *f;
  char *operator="abs";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);

  /* Do the operation: */
  if(f)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_array_dabs_array(f->array, f->size);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, fabs(fnum), NOOPTDATA);
}





void
findmin(struct imgarithparams *p)
{
  gal_data_t *f;
  double min, fnum;
  char *operator="min";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);

  /* Do the operation: */
  if(f)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_statistics_double_min(f->array, f->size, &min);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, min, NOOPTDATA);

      /* Clean up. */
      gal_data_free(f);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, fnum, NOOPTDATA);
}





void
findmax(struct imgarithparams *p)
{
  gal_data_t *f;
  double max, fnum;
  char *operator="max";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);

  /* Do the operation: */
  if(f)                       /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      gal_statistics_double_max(f->array, f->size, &max);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, max, NOOPTDATA);

      /* Clean up. */
      gal_data_free(f);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, fnum, NOOPTDATA);
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
  double fnum, snum;            /* First or second number.    */
  gal_data_t *f, *s;            /* First or second dataset.   */
  double *fp, *sp, *ff, *ss;
  int (*thisfunction)(double, double)=NULL;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);
  pop_operand(p, &snum, &s, operator);

  /* Set the function to use. */
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
  if(f && s)                /* Both are arrays. */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      fp=f->array;
      ss=(sp=s->array)+s->size;
      do *sp = thisfunction(*sp, *fp++); while(++sp<ss);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);

      /* Clean up. */
      gal_data_free(f);
    }
  else if(f)                 /* Only the first is an array. */
    {
      ff=(fp=f->array)+f->size;
      do *fp = thisfunction(snum, *fp); while(++fp<ff);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else if(s)                 /* Only the first is an array. */
    {
      ss=(sp=s->array)+s->size;
      do *sp = thisfunction(*sp, fnum); while(++sp<ss);
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, s);
    }
  else                          /* Both are numbers.           */
    add_operand(p, NOOPTFILENAME, thisfunction(snum, fnum), NOOPTDATA);
}





void
andor(struct imgarithparams *p, char *operator)
{
  double fnum, snum;
  gal_data_t *f, *s;
  double *fp, *sp, *ff;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);
  pop_operand(p, &snum, &s, operator);

  /* Do a small sanity check: */
  if(strcmp(operator, "and") && strcmp(operator, "or") )
    error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we "
          "can address the problem. The value of `operator' in "
          "`andor' (%s) is not recognized", PACKAGE_BUGREPORT, operator);

  /* Do the operation: */
  if(f && s)
    {
      /* Fill the first array with the result. IMPORTANT: It is important
         that the second array pointer is the first checked array, since it
         is also incremented. In the `or' operation, if `*f' is successful,
         `*s++' is never called and so not incremented. */
      sp = s->array;
      ff = (fp=f->array) + f->size;
      if(!strcmp(operator, "and"))
        do *fp = *sp++ && *fp; while(++fp<ff);
      else
        do *fp = *sp++ || *fp; while(++fp<ff);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);

      /* Clean up */
      gal_data_free(s);
    }
  else if(f==NULL || s==NULL)
    error(EXIT_FAILURE, 0, "The `and' and `or' operators need two operators "
          "of the same type: either both images or both numbers.");
  else
    add_operand(p, NOOPTFILENAME,
                !strcmp(operator, "and") ? fnum && snum : fnum || snum,
                NOOPTDATA);
}




void
notfunc(struct imgarithparams *p)
{
  double fnum;
  gal_data_t *f;
  double *fp, *ff;
  char *operator="not";

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);

  /* Do the operation: */
  if(f)
    {
      /* Fill the array with the output values. */
      ff = (fp=f->array) + f->size;
      do *fp = !(*fp); while(++fp<ff);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else
    add_operand(p, NOOPTFILENAME, !fnum, NOOPTDATA);
}





/* In order to not conflict with the internal C `is...' functions, and in
   particular the `isblank' function, we are calling this function
   opisblank for operator-isblank. */
void
opisblank(struct imgarithparams *p)
{
  gal_data_t *f;
  char *operator="isblank";
  double *fp, *ff, fnum;

  /* Pop out the number of operands needed. */
  pop_operand(p, &fnum, &f, operator);

  /* Do the operation: */
  if(f)                         /* Operand is array.        */
    {
      /* Do the operation, note that the output is stored in the first
         input. Also note that since the linked list is
         first-in-first-out, the second operand should be put first
         here. */
      ff=(fp=f->array)+f->size;
      do *fp = isnan(*fp); while(++fp<ff);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, f);
    }
  else                          /* Operand is a number.      */
    add_operand(p, NOOPTFILENAME, isnan(fnum), NOOPTDATA);
}





/* Replace the pixels in the second popped element with the first. While
   choosing the pixels that are selected from the third The third popped
   element. The third is considered to be an array that can only be filled
   with 0 or 1. */
void
where(struct imgarithparams *p)
{
  gal_data_t *n, *c, *i;        /* New, condition, input datasets.  */
  double nnum, cnum, inum;      /* First, second, or third number.  */
  double *np, *cp, *ip, *ii;

  /* ORDER IS VERY IMPORTANT HERE. Pop out the number of operands needed. */
  pop_operand(p, &nnum, &n, "where");              /* New value. */
  pop_operand(p, &cnum, &c, "where");              /* Condition. */
  pop_operand(p, &inum, &i, "where");              /* Input.     */

  /* Do the operation: */
  if(i && c)              /* Both are arrays. */
    {
      cp=c->array;
      ii=(ip=i->array)+i->size;
      if(n)                  /* `new' is an array, not number. */
        {
          np=n->array;
          /* Note that we need to increment "fp" after the check and
             replace. If the increment is inside the conditional replace,
             then when t==0, the increment won't work.*/
          do {*ip = *cp++ ? *np : *ip; ++np;} while(++ip<ii);
        }
      else                      /* `new' is a number, not array. */
        do *ip = *cp++ ? nnum : *ip; while(++ip<ii);

      /* Push the output onto the stack. */
      add_operand(p, NOOPTFILENAME, NOOPTNUMBER, i);

      /* Clean up. */
      gal_data_free(c);
      gal_data_free(n);
    }
  else if ( i==NULL && c==NULL && n==NULL )
    add_operand(p, NOOPTFILENAME, cnum ? nnum : inum, NOOPTDATA);
  else
    error(EXIT_FAILURE, 0, "the first and second arguments (second and "
          "third popped elements) to `where' have to be arrays, or all have "
          "to be numbers.");
}
