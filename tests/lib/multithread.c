/*********************************************************************
A test program to multithreaded building using Gnuastro's helpers.

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
#include <stdio.h>
#include <stdlib.h>

#include "gnuastro/fits.h"
#include "gnuastro/threads.h"


/* This structure can keep all information you want to pass onto the worker
   function on each thread. */
struct params
{
  gal_data_t *image;            /* Dataset to print values of. */
};




/* This is the main worker function which will be called by the different
   threads. 'gal_threads_params' is defined in 'gnuastro/threads.h' and
   contains the pointer to the paramter we want. Note that its input and
   output must have 'void *' types. */
void *
worker_on_thread(void *in_prm)
{
  /* Low-level definitions to be done first. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct params *p=(struct params *)tprm->params;


  /* Subsequent definitions. */
  float *array=p->image->array;
  size_t i, index, *dsize=p->image->dsize;


  /* Go over all the pixels that were assigned to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. */
      index = tprm->indexs[i];


      /* Print the information. */
      printf("(%zu, %zu) on thread %zu: %g\n", index%dsize[1]+1,
             index/dsize[1]+1, tprm->id, array[index]);
    }


  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}




/* A simple program to open a FITS image, distributes its pixels between
   different threads and print the value of each pixel and the thread it
   was assigned to, this will test both the opening of a FITS file and also
   the multi-threaded functions. After running 'make check' you can see the
   outputs in 'tests/multithread.log'.

   Please run the following command for an explanation on easily linking
   and compiling C programs that use Gnuastro's libraries (without having
   to worry about the libraries to link to) anywhere on your system:

      $ info gnuastro "Automatic linking script"
*/
int
main(void)
{
  struct params p;
  char *filename="psf.fits", *hdu="1";
  size_t numthreads=gal_threads_number();


  /* Read the image into memory as a float32 data type. */
  p.image=gal_fits_img_read_to_type(filename, hdu, GAL_TYPE_FLOAT32, -1, 1);


  /* Print some basic information before the actual contents: */
  printf("Pixel values of %s (HDU: %s) on %zu threads.\n", filename, hdu,
         numthreads);
  printf("Used to check the compiled library's capability in opening a "
         "FITS file, and also spinning-off threads.\n");


  /* A small sanity check: this is only intended for 2D arrays. */
  if(p.image->ndim!=2)
    {
      fprintf(stderr, "only 2D images are supported.");
      exit(EXIT_FAILURE);
    }


  /* Spin-off the threads and do the processing on each thread. */
  gal_threads_spin_off(worker_on_thread, &p, p.image->size, numthreads);


  /* Clean up and return. */
  gal_data_free(p.image);
  return EXIT_SUCCESS;
}
