/*********************************************************************
The default 2D kernel to be used in NoiseChisel.
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
     astmkprof --kernel=gaussian,2,5 --oversample=1 --envseed --numrandom=100000

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

size_t kernel_2d_dsize[2]={11, 11};
float kernel_2d[121]={0, 0, 0, 0, 0, 2.599797e-08, 0, 0, 0, 0, 0,

0, 0, 3.008479e-08, 6.938075e-07, 4.493532e-06, 8.276223e-06, 4.515019e-06, 6.947793e-07, 3.04628e-08, 0, 0,

0, 3.009687e-08, 2.556034e-06, 5.936867e-05, 0.0003808578, 0.0007126221, 0.0003827095, 5.902729e-05, 2.553342e-06, 2.978137e-08, 0,

0, 7.021852e-07, 5.912285e-05, 0.00137637, 0.008863639, 0.01648383, 0.008855942, 0.001365171, 5.925718e-05, 7.021184e-07, 0,

0, 4.490787e-06, 0.0003826718, 0.008857355, 0.05742518, 0.1062628, 0.05727194, 0.008880079, 0.0003826067, 4.478989e-06, 0,

2.595735e-08, 8.31301e-06, 0.0007113572, 0.01640853, 0.1061298, 0.1971036, 0.1062611, 0.01647962, 0.000708363, 8.379878e-06, 2.593496e-08,

0, 4.516684e-06, 0.0003846966, 0.008860709, 0.05739478, 0.1062216, 0.05725683, 0.00881713, 0.000383981, 4.473017e-06, 0,

0, 6.950547e-07, 5.920586e-05, 0.00137483, 0.00887785, 0.0164709, 0.008855232, 0.001372743, 5.939038e-05, 7.016624e-07, 0,

0, 3.006322e-08, 2.587011e-06, 5.92911e-05, 0.0003843824, 0.0007118155, 0.000386519, 5.974654e-05, 2.585581e-06, 3.048036e-08, 0,

0, 0, 3.041056e-08, 7.05225e-07, 4.497418e-06, 8.388542e-06, 4.478833e-06, 7.018358e-07, 2.995504e-08, 0, 0,

0, 0, 0, 0, 0, 2.567377e-08, 0, 0, 0, 0, 0};


#endif
