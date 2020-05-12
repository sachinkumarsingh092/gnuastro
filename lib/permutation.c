/*********************************************************************
Permutation -- Work on permutations (arrays of indexs).
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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
#include <string.h>
#include <stdlib.h>

#include <gnuastro/pointer.h>
#include <gnuastro/permutation.h>



/*********************************************************************/
/***************          Permutation info         *******************/
/*********************************************************************/
void
gal_permutation_check(size_t *permutation, size_t size)
{
  size_t i;
  for(i=0; i<size; ++i)
    printf("after[ %-5zu ]    =   before [ %-5zu ]\n", i, permutation[i]);
}




















/*********************************************************************/
/***************          Apply permutation        *******************/
/*********************************************************************/
/* Re-order the input dataset based on the given permutation. If
   'permutation' is NULL, then the input won't be touched (no re-ordering).

   This is a re-implementation of GSL's 'gsl_permute' function (from its
   'permutation/permute_source.c'). The reason we didn't use that function
   was that it uses system-specific types (like 'long' and 'int') which are
   not easily convertable to Gnuastro's width-based types. There is also a
   separate function for each type, heavily using macros to allow a "base"
   function to work on all the types. Thus it is hard to
   read/understand. Since we use fixed-width types, we can easily use
   'memcpy' and have a type-agnostic implementation (only needing the width
   of the type).

   As described in GSL's source code and manual, this implementation comes
   from Knuth's "Art of computer programmin" book, the "Sorting and
   Searching" chapter of Volume 3 (3rd ed). Section 5.2 Exercise 10
   (answers), p 617. See there fore more explanations. The algorithm is a
   little too abstract, but memory and CPU efficient.

   Definition of permutations:

      permute:    OUT[ i       ]   =   IN[ perm[i] ]     i = 0 .. N-1
      inverse:    OUT[ perm[i] ]   =   IN[ i       ]     i = 0 .. N-1
*/
void
gal_permutation_apply(gal_data_t *input, size_t *permutation)
{
  void *tmp;
  size_t i, k, pk, width;
  uint8_t *array=input->array;

  /* If permutation is NULL, then it is assumed that the data doesn't need
     to be re-ordered. */
  if(permutation)
    {
      /* Necessary initializations. */
      width=gal_type_sizeof(input->type);
      tmp=gal_pointer_allocate(input->type, 1, 0, __func__, "tmp");

      /* Do the permutation. */
      for(i=0;i<input->size;++i)
        {
          k=permutation[i];

          while(k>i) k=permutation[k];

          if(k>=i)
            {
              pk = permutation[k];
              if( pk != i )
                {
                  memcpy(tmp, &array[i*width], width);

                  while(pk!=i)
                    {
                      memcpy(&array[k*width], &array[pk*width], width);
                      k  = pk;
                      pk = permutation[k];
                    }

                  memcpy(&array[k*width], tmp, width);
                }
            }
        }

      /* Clean up. */
      free(tmp);
    }
}





/* Apply the inverse of given permutation on the input dataset, see
   'gal_permutation_apply_inverse'. */
void
gal_permutation_apply_inverse(gal_data_t *input, size_t *permutation)
{
  void *tmp, *ttmp;
  size_t i, k, pk, width;
  uint8_t *array=input->array;

  if(permutation)
    {
      /* Initializations */
      width=gal_type_sizeof(input->type);
      tmp=gal_pointer_allocate(input->type, 1, 0, __func__, "tmp");
      ttmp=gal_pointer_allocate(input->type, 1, 0, __func__, "ttmp");

      /* Re-order the values. */
      for(i=0;i<input->size;++i)
        {
          k=permutation[i];

          while(k>i) k=permutation[k];

          if(k>=i)
            {
              pk = permutation[k];

              if(pk!=i)
                {
                  memcpy(tmp, &array[k*width], width);

                  while(pk!=i)
                    {
                      memcpy(ttmp, &array[pk*width], width);
                      memcpy(&array[pk*width], tmp, width);
                      memcpy(tmp, ttmp, width);
                      k  = pk;
                      pk = permutation[k];
                    }

                  memcpy(&array[pk*width], tmp, width);
                }
            }
        }

      /* Clean up. */
      free(tmp);
      free(ttmp);
    }
}
