/*********************************************************************
Functions for linked lists.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <assert.h>

#include <gnuastro/linkedlist.h>










/****************************************************************
 *****************        Two doubles        ********************
 ****************************************************************/
void
gal_linkedlist_add_to_tdll(struct gal_linkedlist_tdll **list,
                           double a, double b)
{
  struct gal_linkedlist_tdll *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: New element in "
          "gal_linkedlist_tdll");

  newnode->a=a;
  newnode->b=b;
  newnode->next=*list;
  *list=newnode;
}





void
gal_linkedlist_pop_from_tdll(struct gal_linkedlist_tdll **list,
                             double *a, double *b)
{
  struct gal_linkedlist_tdll *tmp=*list;

  *a=tmp->a;
  *b=tmp->b;
  *list=tmp->next;
  free(tmp);
}





size_t
gal_linkedlist_num_int_dll(struct gal_linkedlist_tdll *list)
{
  size_t num=0;
  struct gal_linkedlist_tdll *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





void
gal_linkedlist_tdll_to_array_inv(struct gal_linkedlist_tdll *list,
                                 double **d, size_t *num)
{
  size_t i;
  double *td;
  struct gal_linkedlist_tdll *tmp;

  /* Find the number of elements: */
  if(*num==0)
    *num=gal_linkedlist_num_int_dll(list);

  /* Allocate the space (every element of the list has two
     elements.) */
  errno=0;
  td=*d=malloc(2 * *num * sizeof(double));
  if(*d==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: array of gal_linkedlist_tdll "
          "with %lu elements", *num);

  /* Fill in the array in reverse order */
  i = 2 * *num - 2;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    {
      td[i]=tmp->a;
      td[i+1]=tmp->b;
      i-=2;
    }
}





void
gal_linkedlist_free_tdll(struct gal_linkedlist_tdll *list)
{
  struct gal_linkedlist_tdll *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}




















/****************************************************************
 *****************            Float          ********************
 ****************************************************************/
void
gal_linkedlist_print_fll_array(struct gal_linkedlist_fll **afll, size_t num)
{
  size_t i;
  struct gal_linkedlist_fll *tmp;
  for(i=0;i<num;++i)
    {
      printf(" %lu:\n", i);
      for(tmp=afll[i];tmp!=NULL;tmp=tmp->next)
        printf("%f, ", tmp->v);
      printf("\n");
    }
}





void
gal_linkedlist_add_to_fll(struct gal_linkedlist_fll **list, float value)
{
  struct gal_linkedlist_fll *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: New element in "
          "gal_linkedlist_fll");

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





void
gal_linkedlist_pop_from_fll(struct gal_linkedlist_fll **list, float *value)
{
  struct gal_linkedlist_fll *tmp;
  tmp=*list;
  *value=tmp->v;
  *list=tmp->next;
  free(tmp);
}





size_t
gal_linkedlist_num_in_fll(struct gal_linkedlist_fll *list)
{
  size_t num=0;
  struct gal_linkedlist_fll *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





void
gal_linkedlist_fll_to_array(struct gal_linkedlist_fll *list,
                            float **f, size_t *num)
{
  float *tf;
  size_t i=0;
  struct gal_linkedlist_fll *tmp;

  /* Find the number of elements: */
  if(*num==0)
    *num=gal_linkedlist_num_in_fll(list);

  /* Allocate the space: */
  errno=0;
  *f=malloc(*num*sizeof(float));
  if(*f==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: array of gal_linkedlist_fll "
          "with %lu elements", *num);
  tf=*f;

  /* Fill in the array: */
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    tf[i++]=tmp->v;
}





void
gal_linkedlist_free_fll(struct gal_linkedlist_fll *list)
{
  struct gal_linkedlist_fll *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}





void
gal_linkedlist_free_fll_array(struct gal_linkedlist_fll **afll, size_t num)
{
  size_t i;
  for(i=0;i<num;++i)
    gal_linkedlist_free_fll(afll[i]);
  free(afll);
}




















/****************************************************************
 *****************           string          ********************
 ****************************************************************/
void
gal_linkedlist_add_to_stll(struct gal_linkedlist_stll **list, char *value)
{
  struct gal_linkedlist_stll *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno,
          "linkedlist: New element in gal_linkedlist_stll");

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





void
gal_linkedlist_pop_from_stll(struct gal_linkedlist_stll **list, char **value)
{
  struct gal_linkedlist_stll *tmp;
  tmp=*list;
  *value=tmp->v;
  *list=tmp->next;
  free(tmp);
}





void
gal_linkedlist_reverse_stll(struct gal_linkedlist_stll **list)
{
  char *thisstring;
  struct gal_linkedlist_stll *correctorder=NULL;

  while(*list!=NULL)
    {
      gal_linkedlist_pop_from_stll(list, &thisstring);
      gal_linkedlist_add_to_stll(&correctorder, thisstring);
    }
  *list=correctorder;
}




void
gal_linkedlist_print_stll(struct gal_linkedlist_stll *list)
{
  struct gal_linkedlist_stll *tmp;
  for(tmp=list; tmp!=NULL; tmp=tmp->next)
    printf("%s\n", tmp->v);
}





size_t
gal_linkedlist_num_in_stll(struct gal_linkedlist_stll *list)
{
  size_t num=0;
  struct gal_linkedlist_stll *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next) ++num;
  return num;
}

















/****************************************************************
 *****************           size_t          ********************
 ****************************************************************/
void
add_to_sll(struct gal_linkedlist_sll **list, size_t value)
{
  struct gal_linkedlist_sll *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: New element in gal_linkedlist_sll");

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





void
pop_from_sll(struct gal_linkedlist_sll **list, size_t *value)
{
  struct gal_linkedlist_sll *tmp;
  tmp=*list;
  *value=tmp->v;
  *list=tmp->next;
  free(tmp);
}





size_t
numinsll(struct gal_linkedlist_sll *list)
{
  size_t num=0;
  struct gal_linkedlist_sll *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





void
gal_linkedlist_sll_to_array(struct gal_linkedlist_sll *list,
                            size_t **f, size_t *num)
{
  size_t i=0, *tf;
  struct gal_linkedlist_sll *tmp;

  *num=numinsll(list);

  errno=0;
  *f=malloc(*num*sizeof(size_t));
  if(*f==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: array of gal_linkedlist_sll "
          "with %lu elements", *num);
  tf=*f;

  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    tf[i++]=tmp->v;
}





void
gal_linkedlist_print_sll(struct gal_linkedlist_sll *list)
{
  struct gal_linkedlist_sll *tmp;
  printf("\n\n");
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    printf("%lu, ", tmp->v);
  printf("\b\b.\n\n");
  return;
}





void
gal_linkedlist_free_sll(struct gal_linkedlist_sll *list)
{
  struct gal_linkedlist_sll *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}




















/****************************************************************
 ******************        Two way SLL    ***********************
 *****************           size_t          ********************
 ****************************************************************/

void
gal_linkedlist_add_to_tsll_end(struct gal_linkedlist_tsll **last, size_t value)
{
  struct gal_linkedlist_tsll *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: New element in "
          "gal_linkedlist_tsll");

  newnode->v=value;
  newnode->next=*last;
  newnode->prev=NULL;
  if(*last)                        /* If *list is not NULL */
    (*last)->prev=newnode;
  *last=newnode;
}





/* Note that start has to be initialized. */
void
gal_linkedlist_pop_from_tsll_start(struct gal_linkedlist_tsll **first,
                                   size_t *value)
{
  struct gal_linkedlist_tsll *tmp;
  tmp=*first;
  *value=tmp->v;
  *first=tmp->prev;
  free(tmp);
  if(*first)
    (*first)->next=NULL;
}




















/****************************************************************
 ******************        Ordered SLL       ********************
 *****************           size_t          ********************
 ****************************************************************/
/* We want to put the nodes in order based on the 'tosort' value of
each node. The top element should always have the smallest radius. */
void
gal_linkedlist_add_to_osll(struct gal_linkedlist_osll **list,
                           size_t value, float tosort)
{
  struct gal_linkedlist_osll *newnode, *tmp=*list, *prev=NULL;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: New element "
          "in gal_linkedlist_osll");

  newnode->v=value;
  newnode->s=tosort;

  /* *list points to the smallest value in the queue!*/
  while(tmp!=NULL)
    {
      if(tosort<tmp->s) break;
      /* No need for else, it will only come here if the condition
         above is not satisfied. */
      prev=tmp;
      tmp=tmp->next;
    }

  if(tmp==NULL)      /* This is the largest value so far. */
    {                /* '*list' only changes if it is NULL. */
      newnode->next=NULL;
      if(prev) prev->next=newnode;   /* 'prev' is not NULL! */
      else     *list=newnode;        /* Only for initial node. */
    }
  else
    {
      if(prev) prev->next=newnode;
      else     *list=newnode;        /* 'tosort' is smaller than all. */
      newnode->next=tmp;
    }
}





/* Note that the popped element is the smallest! */
void
gal_linkedlist_pop_from_osll(struct gal_linkedlist_osll **list,
                             size_t *value, float *sortvalue)
{
  struct gal_linkedlist_osll *tmp;
  tmp=*list;
  *value=tmp->v;
  *sortvalue=tmp->s;
  *list=tmp->next;
  free(tmp);
}





/* Add the elements of an gal_linkedlist_osll to a gal_linkedlist_sll. */
void
gal_linkedlist_osll_into_sll(struct gal_linkedlist_osll *in,
                             struct gal_linkedlist_sll **out)
{
  struct gal_linkedlist_osll *tmp;
  while(in!=NULL)
    {
      tmp=in->next;
      add_to_sll(out, in->v);
      free(in);
      in=tmp;
    }
}




















/****************************************************************
 ******************   Two way, Ordered SLL   ********************
 *****************           size_t          ********************
 ****************************************************************/

/* The two way ordered SLL looks something like this:

            largest pointer
            |
   NULL <-- (v0,s0) <--> (v1,s1) <--> ... (vn,sn) --> NULL
                                          |
                           smallest pointer

   Where s(n)>s(n+1) for all n.
*/

void
gal_linkedlist_print_tosll(struct gal_linkedlist_tosll *l,
                           struct gal_linkedlist_tosll *s)
{
  size_t counter=1;   /* We are not counting array elements :-D ! */
  while(l!=NULL)
    {
      printf("\t%-5lu (%lu, %.4f) \n", counter++,
             l->v, l->s);
      l=l->next;
      printf("\t\t\t\t(%lu, %.4f)\n", s->v, s->s);
      s=s->prev;
    }
  printf("\n");
}





/* Very similar to Ordered SLL, but now it is two way. */
void
gal_linkedlist_add_to_tosll_end(struct gal_linkedlist_tosll **largest,
                                struct gal_linkedlist_tosll **smallest,
                                size_t value, float tosort)
{
  struct gal_linkedlist_tosll *newnode, *tmp=*largest;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "linkedlist: New element "
          "in gal_linkedlist_tosll");

  newnode->v=value;
  newnode->s=tosort;
  newnode->prev=NULL;

  while(tmp!=NULL)
    {
      if(tosort >= tmp->s) break;
      /* No need for else, it will only come here if the condition
         above is not satisfied. */
      newnode->prev=tmp;
      tmp=tmp->next;
    }

  if(tmp==NULL)      /* This is the smallest value so far.     */
    {                /* '*largest' only changes if it is NULL. */
      newnode->next=NULL;
      *smallest=newnode;
      if(newnode->prev)         /* 'prev' is not NULL! */
        newnode->prev->next=newnode;
      else                      /* 'prev is NULL, Only first. */
        *largest=newnode;
    }
  else
    {
      if(newnode->prev)
        {
          newnode->prev->next->prev=newnode;
          newnode->prev->next=newnode;
        }
      else
        {
          (*largest)->prev=newnode;
          *largest=newnode;       /* 'tosort' is larger than all. */
        }
      newnode->next=tmp;
    }
}





/* Note that start has to be initialized. */
void
gal_linkedlist_pop_from_tosll_start(struct gal_linkedlist_tosll **largest,
                                    struct gal_linkedlist_tosll **smallest,
                                    size_t *value, float *tosort)
{
  struct gal_linkedlist_tosll *tmp=*smallest;

  *value=tmp->v;
  *tosort=tmp->s;

  *smallest=tmp->prev;
  free(tmp);
  if(*smallest)
    (*smallest)->next=NULL;
  else
    *largest=NULL;

  /*printf("Popped v: %lu, s: %f\n", *value, *tosort);*/
}





void
gal_linkedlist_smallest_tosll(struct gal_linkedlist_tosll *largest,
                              struct gal_linkedlist_tosll **smallest)
{
  struct gal_linkedlist_tosll *tmp=largest;

  while(tmp!=NULL)
    {
      if(tmp->next==NULL)
        {
          *smallest=tmp;
          break;
        }
      tmp=tmp->next;
    }

  /* If *largest wasn't NULL initially, tmp should not be NULL because
     the loop terminated before it becomes null. But if it was
     initiall NULL, it will never enter the loop, and so smallest
     should also be NULL. */
  if(tmp==NULL) *smallest=NULL;
}




void
gal_linkedlist_tosll_into_sll(struct gal_linkedlist_tosll *in,
                              struct gal_linkedlist_sll **out)
{
  struct gal_linkedlist_tosll *tmp;
  while(in!=NULL)
    {
      tmp=in->next;
      add_to_sll(out, in->v);
      free(in);
      in=tmp;
    }
}





void
gal_linkedlist_tosll_free(struct gal_linkedlist_tosll *largest)
{
  struct gal_linkedlist_tosll *tmp;
  while(largest!=NULL)
    {
      tmp=largest->next;
      free(largest);
      largest=tmp;
    }
}
