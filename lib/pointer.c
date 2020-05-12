/*********************************************************************
pointer -- facilitate working with pointers and allocation.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

#include <gnuastro/type.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>


/* Increment a give pointer depending on the given type.

   When working with the 'array' elements of 'gal_data_t', we are actually
   dealing with 'void *' pointers. Pointer arithmetic doesn't apply to
   'void *', because the system doesn't know how much space each element
   has to increment the pointer respectively.

   So, here, we will use the type information to find the increment. This
   is mainly useful when dealing with the 'block' pointer of a tile over a
   larger image. This function reads the address as a 'char *' type (note
   that 'char' is guaranteed to have a size of 1 (byte)). It then
   increments the 'char *' by 'increment*sizeof(type)' */
void *
gal_pointer_increment(void *pointer, size_t increment, uint8_t type)
{
  char *p=(char *)pointer;
  return p + increment * gal_type_sizeof(type);
}





/* Find the number of values between two void pointers with a given
   type. See the explanations before 'gal_data_ptr_increment'. */
size_t
gal_pointer_num_between(void *earlier, void *later, uint8_t type)
{
  char *e=(char *)earlier, *l=(char *)later;
  return (l-e)/gal_type_sizeof(type);
}





/* Allocate an array based on the value of type. Note that the argument
   'size' is the number of elements, necessary in the array, the number of
   bytes each element needs will be determined internaly by this function
   using the datatype argument, so you don't have to worry about it. */
void *
gal_pointer_allocate(uint8_t type, size_t size, int clear,
                     const char *funcname, const char *varname)
{
  void *array;

  errno=0;
  array = ( clear
            ? calloc( size,  gal_type_sizeof(type) )
            : malloc( size * gal_type_sizeof(type) ) );
  if(array==NULL)
    {
      if(varname)
        error(EXIT_FAILURE, errno, "%s: %zu bytes couldn't be allocated "
              "for variable '%s'", funcname ? funcname : __func__,
              size * gal_type_sizeof(type), varname);
      else
        error(EXIT_FAILURE, errno, "%s: %zu bytes couldn't be allocated",
              funcname ? funcname : __func__, size * gal_type_sizeof(type));
    }

  return array;
}





void *
gal_pointer_allocate_mmap(uint8_t type, size_t size, int clear,
                          char **filename, int quietmmap)
{
  void *out;
  int filedes;
  uint8_t uc=0;
  char *dirname=NULL;
  size_t bsize=size*gal_type_sizeof(type);


  /* Check if the '.gnuastro_mmap' folder exists, write the file there. If
     it doesn't exist, then make it. If it can't be built, we'll make a
     randomly named file in the current directory. */
  gal_checkset_allocate_copy("./.gnuastro_mmap/", &dirname);
  if( gal_checkset_mkdir(dirname) )
    {
      /* The directory couldn't be built. Free the old name. */
      free(dirname);

      /* Set 'dirname' to NULL so it knows not to write in a directory. */
      dirname=NULL;
    }


  /* Set the filename. If 'dirname' couldn't be allocated, directly make
     the memory map file in the current directory (just as a hidden
     file). */
  if( asprintf(filename, "%sXXXXXX", dirname?dirname:"./.gnuastro_mmap_")<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  if(dirname) free(dirname);


  /* Create a zero-sized file and keep its descriptor.  */
  errno=0;
  /*filedes=open(filename, O_RDWR | O_CREAT | O_EXCL | O_TRUNC );*/
  filedes=mkstemp(*filename);
  if(filedes==-1)
    error(EXIT_FAILURE, errno, "%s: %s couldn't be created", __func__,
          *filename);


  /* Make the necessary space on the file. */
  errno=0;
  if( lseek(filedes, bsize, SEEK_SET) == -1 )
    error(EXIT_FAILURE, errno, "%s: %s: unable to change file position by "
          "%zu bytes", __func__, *filename, bsize);


  /* Inform the user. */
  if(!quietmmap)
    error(EXIT_SUCCESS, 0, "%s: temporary %zu byte file (consider "
          "'--minmapsize')", *filename, bsize);


  /* Write to the newly set file position so the space is allocated. To do
     this, we are simply writing 'uc' (a byte with value 0) into the space
     we identified by 'lseek' (above). This will ensure that this space is
     set a side for this array and prepare us to use 'mmap'. */
  if( write(filedes, &uc, 1) == -1)
    error(EXIT_FAILURE, errno, "%s: %s: unable to write one byte at the "
          "%zu-th position", __func__, *filename, bsize);


  /* Map the memory. */
  errno=0;
  out=mmap(NULL, bsize, PROT_READ | PROT_WRITE, MAP_SHARED, filedes, 0);
  if(out==MAP_FAILED)
    {
      fprintf(stderr, "\n%s: WARNING: the following error may be due to "
              "many mmap allocations. Recall that the kernel only allows "
              "finite number of mmap allocations. It is recommended to use "
              "ordinary RAM allocation for smaller arrays and keep mmap'd "
              "allocation only for the large volumes.\n\n", __func__);
      error(EXIT_FAILURE, errno, "couldn't map %zu bytes into the file '%s'",
            bsize, *filename);
    }


  /* Close the file. */
  if( close(filedes) == -1 )
    error(EXIT_FAILURE, errno, "%s: %s couldn't be closed",
          __func__, *filename);


  /* If it was supposed to be cleared, then clear the memory. */
  if(clear) memset(out, 0, bsize);


  /* Return the mmap'd pointer and save the file name. */
  return out;
}
