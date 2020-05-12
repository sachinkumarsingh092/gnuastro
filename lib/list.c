/*********************************************************************
Functions for linked lists.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

#include <gnuastro/list.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>









/****************************************************************
 *****************           String          ********************
 ****************************************************************/
void
gal_list_str_add(gal_list_str_t **list, char *value,
                 int allocate)
{
  gal_list_str_t *newnode;

  /* If the value is a NULL pointer, don't add to the list. */
  if(value==NULL) return;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  if(allocate)
    gal_checkset_allocate_copy(value, &newnode->v);
  else
    newnode->v=value;

  newnode->next=*list;
  *list=newnode;
}





char *
gal_list_str_pop(gal_list_str_t **list)
{
  char *out=NULL;
  gal_list_str_t *tmp;
  if(*list)
    {
      tmp=*list;
      out=tmp->v;
      *list=tmp->next;
      free(tmp);
    }
  return out;
}





size_t
gal_list_str_number(gal_list_str_t *list)
{
  size_t num=0;
  gal_list_str_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next) ++num;
  return num;
}





gal_list_str_t *
gal_list_str_last(gal_list_str_t *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}





void
gal_list_str_print(gal_list_str_t *list)
{
  gal_list_str_t *tmp;
  for(tmp=list; tmp!=NULL; tmp=tmp->next)
    printf("%s\n", tmp->v);
}





void
gal_list_str_reverse(gal_list_str_t **list)
{
  char *thisstring;
  gal_list_str_t *correctorder=NULL;

  /* Only do the reversal if there is more than one element. */
  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisstring=gal_list_str_pop(list);
          gal_list_str_add(&correctorder, thisstring, 0);
        }
      *list=correctorder;
    }
}





void
gal_list_str_free(gal_list_str_t *list, int freevalue)
{
  gal_list_str_t *tmp;
  while(list!=NULL)
    {
      tmp=list->next;
      if(freevalue)
        free(list->v);
      free(list);
      list=tmp;
    }
}




















/****************************************************************
 *****************            int            ********************
 ****************************************************************/
void
gal_list_i32_add(gal_list_i32_t **list, int32_t value)
{
  gal_list_i32_t *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





int32_t
gal_list_i32_pop(gal_list_i32_t **list)
{
  gal_list_i32_t *tmp;
  int out=GAL_BLANK_INT32;

  if(*list)
    {
      tmp=*list;
      out=tmp->v;
      *list=tmp->next;
      free(tmp);
    }
  return out;
}





size_t
gal_list_i32_number(gal_list_i32_t *list)
{
  size_t num=0;
  gal_list_i32_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





gal_list_i32_t *
gal_list_i32_last(gal_list_i32_t *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}





void
gal_list_i32_print(gal_list_i32_t *list)
{
  gal_list_i32_t *tmp;
  for(tmp=list; tmp!=NULL; tmp=tmp->next)
    printf("%"PRId32"\n", tmp->v);
}





void
gal_list_i32_reverse(gal_list_i32_t **list)
{
  int thisnum;
  gal_list_i32_t *correctorder=NULL;

  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisnum=gal_list_i32_pop(list);
          gal_list_i32_add(&correctorder, thisnum);
        }
      *list=correctorder;
    }
}





int32_t *
gal_list_i32_to_array(gal_list_i32_t *list, int reverse, size_t *num)
{
  size_t i;
  int32_t *out=NULL;
  gal_list_i32_t *tmp;

  *num=gal_list_i32_number(list);

  if(*num)
    {
      out=gal_pointer_allocate(GAL_TYPE_SIZE_T, *num, 0, __func__, "out");

      i = reverse ? *num-1: 0;
      if(reverse)
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i--]=tmp->v;
      else
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i++]=tmp->v;
    }

  return out;
}





void
gal_list_i32_free(gal_list_i32_t *list)
{
  gal_list_i32_t *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}




















/****************************************************************
 *****************           size_t          ********************
 ****************************************************************/
void
gal_list_sizet_add(gal_list_sizet_t **list, size_t value)
{
  gal_list_sizet_t *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





size_t
gal_list_sizet_pop(gal_list_sizet_t **list)
{
  gal_list_sizet_t *tmp;
  size_t out=GAL_BLANK_SIZE_T;

  if(list)
    {
      tmp=*list;
      out=tmp->v;
      *list=tmp->next;
      free(tmp);
    }
  return out;
}





size_t
gal_list_sizet_number(gal_list_sizet_t *list)
{
  size_t num=0;
  gal_list_sizet_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





gal_list_sizet_t *
gal_list_sizet_last(gal_list_sizet_t *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}





void
gal_list_sizet_print(gal_list_sizet_t *list)
{
  gal_list_sizet_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    printf("%zu\n", tmp->v);
  return;
}





void
gal_list_sizet_reverse(gal_list_sizet_t **list)
{
  size_t thisnum;
  gal_list_sizet_t *correctorder=NULL;

  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisnum=gal_list_sizet_pop(list);
          gal_list_sizet_add(&correctorder, thisnum);
        }
      *list=correctorder;
    }
}





size_t *
gal_list_sizet_to_array(gal_list_sizet_t *list, int reverse, size_t *num)
{
  size_t i, *out=NULL;
  gal_list_sizet_t *tmp;

  *num=gal_list_sizet_number(list);

  if(*num)
    {
      out=gal_pointer_allocate(GAL_TYPE_SIZE_T, *num, 0, __func__, "out");

      i = reverse ? *num-1: 0;
      if(reverse)
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i--]=tmp->v;
      else
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i++]=tmp->v;
    }

  return out;
}





void
gal_list_sizet_free(gal_list_sizet_t *list)
{
  gal_list_sizet_t *tmp, *ttmp;
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
gal_list_f32_add(gal_list_f32_t **list, float value)
{
  struct gal_list_f32_t *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





float
gal_list_f32_pop(gal_list_f32_t **list)
{
  float out=NAN;
  gal_list_f32_t *tmp;

  if(*list)
    {
      tmp=*list;
      out=tmp->v;
      *list=tmp->next;
      free(tmp);
    }
  return out;
}





size_t
gal_list_f32_number(gal_list_f32_t *list)
{
  size_t num=0;
  gal_list_f32_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





gal_list_f32_t *
gal_list_f32_last(gal_list_f32_t *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}





void
gal_list_f32_reverse(gal_list_f32_t **list)
{
  float thisnum;
  gal_list_f32_t *correctorder=NULL;

  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisnum=gal_list_f32_pop(list);
          gal_list_f32_add(&correctorder, thisnum);
        }
      *list=correctorder;
    }
}





void
gal_list_f32_print(gal_list_f32_t *list)
{
  gal_list_f32_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    printf("%f\n", tmp->v);
  return;
}





float *
gal_list_f32_to_array(gal_list_f32_t *list, int reverse, size_t *num)
{
  size_t i;
  float *out=NULL;
  gal_list_f32_t *tmp;

  /* Find the number of elements: */
  *num=gal_list_f32_number(list);

  /* If there is actually anything in the list, then allocate the array and
     fill it in. */
  if(*num)
    {
      /* Allocate the space: */
      out=gal_pointer_allocate(GAL_TYPE_FLOAT32, *num, 0, __func__, "out");

      /* Fill in the array. */
      i = reverse ? *num-1: 0;
      if(reverse)
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i--]=tmp->v;
      else
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i++]=tmp->v;
    }

  /* Return the created array. */
  return out;
}





void
gal_list_f32_free(gal_list_f32_t *list)
{
  gal_list_f32_t *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}




















/****************************************************************
 *****************          Double           ********************
 ****************************************************************/
void
gal_list_f64_add(gal_list_f64_t **list, double value)
{
  gal_list_f64_t *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





double
gal_list_f64_pop(gal_list_f64_t **list)
{
  double out=NAN;
  gal_list_f64_t *tmp;

  if(*list)
    {
      tmp=*list;
      out=tmp->v;
      *list=tmp->next;
      free(tmp);
    }
  return out;
}





size_t
gal_list_f64_number(gal_list_f64_t *list)
{
  size_t num=0;
  gal_list_f64_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





gal_list_f64_t *
gal_list_f64_last(gal_list_f64_t *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}





void
gal_list_f64_print(gal_list_f64_t *list)
{
  gal_list_f64_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    printf("%f\n", tmp->v);
  return;
}





void
gal_list_f64_reverse(gal_list_f64_t **list)
{
  double thisvalue;
  gal_list_f64_t *correctorder=NULL;

  /* Only do the reversal if there is more than one element. */
  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisvalue=gal_list_f64_pop(list);
          gal_list_f64_add(&correctorder, thisvalue);
        }
      *list=correctorder;
    }
}




double *
gal_list_f64_to_array(gal_list_f64_t *list, int reverse, size_t *num)
{
  size_t i;
  double *out=NULL;
  gal_list_f64_t *tmp;

  /* Find the number of elements: */
  *num=gal_list_f64_number(list);

  /* If there is actually anything in the list, then allocate the array and
     fill it in. */
  if(*num)
    {
      /* Allocate the space: */
      out=gal_pointer_allocate(GAL_TYPE_FLOAT64, *num, 0, __func__, "out");

      /* Fill in the array. */
      i = reverse ? *num-1: 0;
      if(reverse)
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i--]=tmp->v;
      else
        for(tmp=list;tmp!=NULL;tmp=tmp->next)
          out[i++]=tmp->v;
    }

  /* Return the created array. */
  return out;
}





void
gal_list_f64_free(gal_list_f64_t *list)
{
  gal_list_f64_t *tmp, *ttmp;
  tmp=list;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}



















/****************************************************************
 *****************          void *           ********************
 ****************************************************************/
void
gal_list_void_add(gal_list_void_t **list, void *value)
{
  gal_list_void_t *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

  newnode->v=value;
  newnode->next=*list;
  *list=newnode;
}





void *
gal_list_void_pop(gal_list_void_t **list)
{
  void *out=NULL;
  gal_list_void_t *tmp;

  if(*list)
    {
      tmp=*list;
      out=tmp->v;
      *list=tmp->next;
      free(tmp);
    }
  return out;
}





size_t
gal_list_void_number(gal_list_void_t *list)
{
  size_t num=0;
  gal_list_void_t *tmp;
  for(tmp=list;tmp!=NULL;tmp=tmp->next)
    ++num;
  return num;
}





gal_list_void_t *
gal_list_void_last(gal_list_void_t *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}





void
gal_list_void_reverse(gal_list_void_t **list)
{
  void *thisptr;
  gal_list_void_t *correctorder=NULL;

  if( *list && (*list)->next )
    {
      while(*list!=NULL)
        {
          thisptr=gal_list_void_pop(list);
          gal_list_void_add(&correctorder, thisptr);
        }
      *list=correctorder;
    }
}





void
gal_list_void_free(gal_list_void_t *list, int freevalue)
{
  gal_list_void_t *tmp=list, *ttmp;
  while(tmp!=NULL)
    {
      if(freevalue) free(tmp->v);
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}




















/****************************************************************
 ****************       Ordered size_t       ********************
 ****************************************************************/
/* We want to put the nodes in order based on the 'tosort' value of
each node. The top element should always have the smallest radius. */
void
gal_list_osizet_add(gal_list_osizet_t **list,
                    size_t value, float tosort)
{
  gal_list_osizet_t *newnode, *tmp=*list, *prev=NULL;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

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
size_t
gal_list_osizet_pop(gal_list_osizet_t **list, float *sortvalue)
{
  size_t value;
  gal_list_osizet_t *tmp=*list;

  if(*list)
    {
      value=tmp->v;
      *sortvalue=tmp->s;
      *list=tmp->next;
      free(tmp);
    }
  else
    {
      value=GAL_BLANK_SIZE_T;
      *sortvalue=NAN;
    }

  return value;
}





/* Add the elements of an gal_list_osll to a gal_list_sll. */
void
gal_list_osizet_to_sizet_free(gal_list_osizet_t *in,
                              gal_list_sizet_t **out)
{
  gal_list_osizet_t *tmp;
  while(in!=NULL)
    {
      tmp=in->next;
      gal_list_sizet_add(out, in->v);
      free(in);
      in=tmp;
    }
}




















/****************************************************************
 ******************   Two way, Ordered SLL   ********************
 *****************           size_t          ********************
 ****************************************************************/
/* Doubly-linked ordered size_t list can be visualized like this:

            largest pointer
            |
   NULL <-- (v0,s0) <--> (v1,s1) <--> ... (vn,sn) --> NULL
                                          |
                           smallest pointer

   Where s(n)>s(n+1) for all n.
*/
/* Very similar to Ordered SLL, but now it is two way. */
void
gal_list_dosizet_add(gal_list_dosizet_t **largest,
                     gal_list_dosizet_t **smallest, size_t value, float tosort)
{
  gal_list_dosizet_t *newnode, *tmp=*largest;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating new node", __func__);

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
size_t
gal_list_dosizet_pop_smallest(gal_list_dosizet_t **largest,
                              gal_list_dosizet_t **smallest, float *tosort)
{
  size_t value;
  gal_list_dosizet_t *tmp=*smallest;

  if(*smallest)
    {
      value=tmp->v;
      *tosort=tmp->s;

      *smallest=tmp->prev;
      free(tmp);
      if(*smallest)
        (*smallest)->next=NULL;
      else
        *largest=NULL;
    }
  else
    {
      /* If 'smallest' is NULL, 'largest' should also be NULL. */
      if(*largest)
        error(EXIT_FAILURE, 0, "%s: 'largest' and 'smallest' pointers must "
              "both be non-NULL or both be NULL. However, in this call, "
              "'smallest' was NULL while 'largest' isn't NULL", __func__);
      value=GAL_BLANK_SIZE_T;
      *tosort=NAN;
    }

  /*printf("Popped v: %zu, s: %f\n", *value, *tosort);*/

  return value;
}





void
gal_list_dosizet_print(gal_list_dosizet_t *largest,
                       gal_list_dosizet_t *smallest)
{
  size_t counter=1;   /* We are not counting array elements :-D ! */
  while(largest!=NULL)
    {
      printf("\t%-5zu (%zu, %.4f) \n", counter++,
             largest->v, largest->s);
      largest=largest->next;
      printf("\t\t\t\t(%zu, %.4f)\n", smallest->v, smallest->s);
      smallest=smallest->prev;
    }
  printf("\n");
}





void
gal_list_dosizet_to_sizet(gal_list_dosizet_t *in, gal_list_sizet_t **out)
{
  gal_list_dosizet_t *tmp;
  while(in!=NULL)
    {
      tmp=in->next;
      gal_list_sizet_add(out, in->v);
      free(in);
      in=tmp;
    }
}





void
gal_list_dosizet_free(gal_list_dosizet_t *largest)
{
  gal_list_dosizet_t *tmp;
  while(largest!=NULL)
    {
      tmp=largest->next;
      free(largest);
      largest=tmp;
    }
}




















/*********************************************************************/
/*************    Data structure as a linked list   ******************/
/*********************************************************************/
/* Add a new data structure to the top of an existing linked list of data
   structures. Note that if the new node is its self a list, all its nodes
   will be added to the list. */
void
gal_list_data_add(gal_data_t **list, gal_data_t *newnode)
{
  gal_data_t *tmp=newnode, *toadd;

  /* Check if newnode is itself a list or not. */
  if(newnode->next)
    {
      /* Go onto the last node in newnode's existing list. */
      while(tmp->next) tmp=tmp->next;

      /* Set the last node as the node to add to the list. */
      toadd=tmp;
    }
  else
    /* Its not a list, so just set it to 'toadd'. */
    toadd=newnode;


  /* Set the next element of toadd and update what list points to.*/
  toadd->next=*list;
  *list=newnode;
}





void
gal_list_data_add_alloc(gal_data_t **list, void *array, uint8_t type,
                        size_t ndim, size_t *dsize, struct wcsprm *wcs,
                        int clear, size_t minmapsize, int quietmmap,
                        char *name, char *unit, char *comment)
{
  gal_data_t *newnode;

  /* Put all the input information into a new data structure node. */
  newnode=gal_data_alloc(array, type, ndim, dsize, wcs, clear,
                         minmapsize, quietmmap, name, unit, comment);

  /* Add the new node to the list. */
  gal_list_data_add(list, newnode);
}





gal_data_t *
gal_list_data_pop(gal_data_t **list)
{
  gal_data_t *out=NULL;

  /* If list is not empty. */
  if(*list)
    {
      /* Keep the top pointer. */
      out=*list;

      /* Move the list pointer to the next node. */
      *list=out->next;

      /* Set the next poitner of the out pointer to NULL so it isn't
         interpretted as a list any more. */
      out->next=NULL;
    }
  return out;
}





void
gal_list_data_reverse(gal_data_t **list)
{
  gal_data_t *popped, *in=*list, *reversed=NULL;

  /* Only do the job if the list is not NULL and has more than one node. */
  if( in && in->next )
    {
      while(in!=NULL)
        {
          popped=gal_list_data_pop(&in);
          gal_list_data_add(&reversed, popped);
        }
      *list=reversed;
    }
}





gal_data_t **
gal_list_data_to_array_ptr(gal_data_t *list, size_t *num)
{
  size_t i, n;
  gal_data_t *tmp, **out;

  /* Count how many columns are necessary. */
  n=*num=gal_list_data_number(list);

  /* Allocate space for the array. */
  errno=0;
  out=malloc(n * sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, 0, "%s: couldn't allocate %zu bytes", __func__,
          n * sizeof *out);

  /* Fill up the array with the pointers and return. */
  i=0;
  for(tmp=list;tmp!=NULL;tmp=tmp->next) out[i++]=tmp;
  return out;
}





size_t
gal_list_data_number(gal_data_t *list)
{
  size_t num=0;
  while(list!=NULL)
    {
      ++num;
      list=list->next;
    }
  return num;
}





gal_data_t *
gal_list_data_last(gal_data_t *list)
{
  if(list)
    {
      while(list->next!=NULL) list=list->next;
      return list;
    }
  else return NULL;
}




void
gal_list_data_free(gal_data_t *list)
{
  struct gal_data_t *tmp;
  while(list!=NULL)
    {
      tmp=list->next;
      gal_data_free(list);
      list=tmp;
    }
}
