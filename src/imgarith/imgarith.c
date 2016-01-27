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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include "main.h"
#include "fitsarrayvv.h"

#include "imgarith.h"            /* needs main.h.                  */




/***************************************************************/
/*************    Operand linked list functions    *************/
/***************************************************************/
void
add_operand(struct operand **list, int isname, char *filename, double number)
{
  struct operand *newnode;

  /* Allocate space for the new operand. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "imgarith.c: Making new element in operand");

  /* Fill in the values. */
  newnode->isname=isname;
  newnode->number=number;
  newnode->filename=filename;

  /* Make the link to the previous list. */
  newnode->next=*list;
  *list=newnode;
}





void
pop_operand(struct operand **list, int *isname, char **filename,
            double *number)
{
  struct operand *tmp=*list;

  /* Return the values. */
  *isname=tmp->isname;
  *number=tmp->number;
  *filename=tmp->filename;

  /* Remove this node from the queue. */
  *list=tmp->next;
  free(tmp);
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
  struct stll *token;

  for(token=p->tokens;token!=NULL;token=token->next)
    {
      if(nameisfits(token->v))
        printf("\n%s\n", token->v);
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
