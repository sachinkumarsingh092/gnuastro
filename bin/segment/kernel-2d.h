/*********************************************************************
The default 2D kernel to be used in Segment.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
#ifndef __GAL_KERNEL2D_H__
#define __GAL_KERNEL2D_H__


/* How to produce this kernel
   ==========================

   Below, the steps necessary to easily create the C contents of this file
   are described. The first step can be modified to change the default
   kernel properties and put the new contents into this file.

   Make the kernel
   ---------------

   We'll first make the kernel using MakeProfiles with the following
   commands. IMPORTANT NOTE: because the kernel is so sharp, random
   sampling is going to be used for all the pixels (the center won't be
   used). So it is important to have a large number of random points to
   make the very slight differences between symmetric parts of the profile
   even less significant.

     export GSL_RNG_SEED=1
     export GSL_RNG_TYPE=ranlxs2
     astmkprof --kernel=gaussian,1.5,5 --oversample=1 --envseed --numrandom=100000

   Convert it to C code
   --------------------

   Put the following C program into a file called 'kernel.c'.

     #include <stdio.h>
     #include <stdlib.h>
     #include <gnuastro/fits.h>

     int
     main(void)
     {
       size_t i;
       float *arr;
       gal_data_t *img=gal_fits_img_read_to_type("kernel.fits", "1",
                                                 GAL_TYPE_FLOAT32, -1);

       arr=img->array;

       printf("size_t kernel_2d_dsize[2]={%zu, %zu};\n",
              img->dsize[0], img->dsize[1]);
       printf("float kernel_2d[%zu]={", img->size);
       for(i=0;i<img->size;++i)
         {
           if(i>0)
             {
               if(i % img->dsize[1] == 0 ) printf("\n\n");
             }

           // We cannot use '\b' here, since we are writing directly
           // to the command-line, so we'll first write the number,
           // then decide if any subsequent character (a comma)
           // should be written.
           printf("%.7g", arr[i]);

           // The last element doesn't need a comma. In each line,
           // the last character must not be a space, but for easy
           // readability, the elements in between need a space.
           if( i!=(img->size-1) )
             printf("%s", ((i+1)%img->dsize[1]) ? ", " : ",");
         }
       printf("};\n");

       gal_data_free(img);
       return EXIT_SUCCESS;
     }

   Run the C program
   -----------------

   We can now compile and run that C program and put the outputs in
   'kernel.c'. Once its created, copy the contents of 'kernel-2d.h' after
   these comments.

     $ astbuildprog -q kernel.c > kernel-2d.h
 */

size_t kernel_2d_dsize[2]={7, 7};
float kernel_2d[49]={0, 3.992438e-07, 8.88367e-06, 2.470061e-05, 8.96143e-06, 3.961747e-07, 0,

3.961645e-07, 8.509836e-05, 0.001905851, 0.005246491, 0.001900595, 8.399635e-05, 3.977891e-07,

8.959198e-06, 0.00190299, 0.04301567, 0.1174493, 0.0428412, 0.001911332, 8.923742e-06,

2.455387e-05, 0.005209642, 0.1172349, 0.3221542, 0.1174603, 0.005248448, 2.447141e-05,

9.018465e-06, 0.001908686, 0.04294781, 0.1173853, 0.04282322, 0.001887719, 8.985901e-06,

3.969509e-07, 8.505241e-05, 0.001909065, 0.005238522, 0.001906396, 8.491996e-05, 3.998521e-07,

0, 3.998288e-07, 9.012383e-06, 2.466673e-05, 9.072039e-06, 4.024199e-07, 0};


#endif
