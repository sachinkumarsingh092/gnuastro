/*********************************************************************
The default 2D kernel to be used in NoiseChisel and Segment.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018, Free Software Foundation, Inc.

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

   We'll first make the kernel using MakeProfiles, then crop the outer
   layers that are zero. Put these lines into a script and run it.

     set -o errexit           # Stop if a program returns false.
     export GSL_RNG_TYPE=ranlxs2
     export GSL_RNG_SEED=1
     astmkprof --kernel=gaussian,2,5 --oversample=1 --envseed -ok2.fits
     astcrop k2.fits --section=2:*-1,2:*-1 --zeroisnotblank    \
             --mode=img --output=fwhm2.fits
     rm k2.fits

   Convert it to C code
   --------------------

   We'll use this tiny C program to convert the FITS file into the two
   important C variables:

     #include <stdio.h>
     #include <stdlib.h>
     #include <gnuastro/fits.h>

     int
     main(void)
     {
       size_t i;
       float *arr;
       gal_data_t *img=gal_fits_img_read_to_type("fwhm2.fits", "1",
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

           // We cannot use `\b' here, since we are writing directly
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
   `ready.c'. Once its created, copy the contents of `ready.c' after these
   comments.

     $ astbuildprog -q prepare.c > ready.c
 */

size_t kernel_2d_dsize[2]={11, 11};
float kernel_2d[121]={0, 0, 0, 0, 0, 2.570758e-08, 0, 0, 0, 0, 0,

0, 0, 2.981546e-08, 7.249833e-07, 4.468747e-06, 8.409227e-06, 4.554846e-06, 7.034199e-07, 3.002102e-08, 0, 0,

0, 3.054498e-08, 2.614154e-06, 5.891601e-05, 0.0003810036, 0.000708165, 0.0003842406, 5.963722e-05, 2.618934e-06, 2.990584e-08, 0,

0, 7.199899e-07, 5.801019e-05, 0.001365485, 0.009023659, 0.01638159, 0.008892864, 0.001345278, 5.920425e-05, 6.984741e-07, 0,

0, 4.584869e-06, 0.0003830431, 0.008917402, 0.05743425, 0.1061489, 0.05746412, 0.008902563, 0.0003849257, 4.448404e-06, 0,

2.572769e-08, 8.414041e-06, 0.0007008284, 0.0164456, 0.1055995, 0.19753, 0.1061855, 0.01653461, 0.0007141303, 8.41643e-06, 2.550312e-08,

0, 4.582525e-06, 0.0003775396, 0.00898499, 0.05741534, 0.1062144, 0.05700329, 0.008838926, 0.0003822096, 4.543726e-06, 0,

0, 6.883925e-07, 6.09077e-05, 0.001339333, 0.008817007, 0.01636454, 0.008995386, 0.001407854, 6.004799e-05, 7.203602e-07, 0,

0, 3.095966e-08, 2.575403e-06, 5.89859e-05, 0.0003804447, 0.0007091904, 0.0003810006, 5.903253e-05, 2.575202e-06, 2.934356e-08, 0,

0, 0, 3.040937e-08, 7.018197e-07, 4.543086e-06, 8.296753e-06, 4.434901e-06, 6.659026e-07, 3.066215e-08, 0, 0,

0, 0, 0, 0, 0, 2.603901e-08, 0, 0, 0, 0, 0};

#endif
