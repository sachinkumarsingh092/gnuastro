/*********************************************************************
ImageWarp - Warp images using projective mapping.
ImageWarp is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <string.h>
#include <fitsio.h>

#include <nproc.h>               /* From Gnulib.                   */

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/txtarray.h>

#include <timing.h>              /* Includes time.h and sys/time.h */
#include <checkset.h>
#include <commonargs.h>
#include <configfiles.h>

#include "main.h"

#include "ui.h"                  /* Needs main.h                   */
#include "args.h"                /* Needs main.h, includes argp.h. */


/* Set the file names of the places where the default parameters are
   put. */
#define CONFIG_FILE SPACK CONF_POSTFIX
#define SYSCONFIG_FILE SYSCONFIG_DIR "/" CONFIG_FILE
#define USERCONFIG_FILEEND USERCONFIG_DIR CONFIG_FILE
#define CURDIRCONFIG_FILE CURDIRCONFIG_DIR CONFIG_FILE










/**************************************************************/
/**************       Options and parameters    ***************/
/**************************************************************/
void
readconfig(char *filename, struct imgwarpparams *p)
{
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;
  char key='a';        /* Not used, just a place holder. */

  /* When the file doesn't exist or can't be opened, it is ignored. It
     might be intentional, so there is no error. If a parameter is
     missing, it will be reported after all defaults are read. */
  fp=fopen(filename, "r");
  if (fp==NULL) return;


  /* Allocate some space for `line` with `len` elements so it can
     easily be freed later on. The value of `len` is arbitarary at
     this point, during the run, getline will change it along with the
     pointer to line. */
  errno=0;
  line=malloc(len*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes in readdefaults",
          len * sizeof *line);

  /* Read the tokens in the file:  */
  while(getline(&line, &len, fp) != -1)
    {
      /* Prepare the "name" and "value" strings, also set lineno. */
      GAL_CONFIGFILES_START_READING_LINE;




      /* Inputs: */
      if(strcmp(name, "hdu")==0)
        gal_checkset_allocate_copy_set(value, &cp->hdu, &cp->hduset);
      else if(strcmp(name, "hstartwcs")==0)
        {
          if(up->hstartwcsset) continue;
          gal_checkset_sizet_el_zero(value, &p->hstartwcs, name, key,
                                     SPACK, filename, lineno);
          up->hstartwcsset=1;
        }
      else if(strcmp(name, "hendwcs")==0)
        {
          if(up->hendwcsset) continue;
          gal_checkset_sizet_el_zero(value, &p->hendwcs, name, key, SPACK,
                                     filename, lineno);
          up->hendwcsset=1;
        }





      /* Outputs */
      else if(strcmp(name, "matrix")==0)
        gal_checkset_allocate_copy_set(value, &up->matrixstring,
                                       &up->matrixstringset);

      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output, &cp->outputset);

      else if(strcmp(name, "maxblankfrac")==0)
        {
          if(up->maxblankfracset) continue;
          gal_checkset_float_l_0_s_1(value, &p->maxblankfrac, name, key,
                                     SPACK, filename, lineno);
          up->maxblankfracset=1;
        }
      else if(strcmp(name, "align")==0)
        {
          if(up->alignset) continue;
          gal_checkset_int_zero_or_one(value, &up->align, name, key,
                                       SPACK, filename, lineno);
          up->alignset=1;
        }
      else if(strcmp(name, "rotate")==0)
        {
          if(up->rotateset) continue;
          gal_checkset_any_float(value, &up->rotate, name, key, SPACK,
                                 filename, lineno);
          up->rotateset=1;
        }
      else if(strcmp(name, "scale")==0)
        {
          if(up->scaleset) continue;
          gal_checkset_any_float(value, &up->scale, name, key, SPACK,
                                 filename, lineno);
          up->scaleset=1;
        }




      /* Operating modes: */
      /* Read options common to all programs */
      GAL_CONFIGFILES_READ_COMMONOPTIONS_FROM_CONF


      else
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct imgwarpparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);
  if(up->hstartwcsset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "hstartwcs", p->hstartwcs);
  if(up->hendwcsset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "hendwcs", p->hendwcs);

  fprintf(fp, "\n# Output parameters:\n");
  if(up->matrixstringset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("matrix", up->matrixstring);
  if(up->alignset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "align", up->align);
  if(up->rotateset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "rotate", up->rotate);
  if(up->scaleset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "scale", up->scale);

  if(cp->outputset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("output", cp->output);

  if(up->maxblankfracset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "maxblankfrac", p->maxblankfrac);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct imgwarpparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");
  if(up->maxblankfracset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("maxblankfrac");

  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Prepare Matrix       *******************/
/**************************************************************/
void
readmatrixoption(struct imgwarpparams *p)
{
  size_t counter=0;
  char *t, *tailptr;

  /* Allocate the necessary space. */
  errno=0;
  p->matrix=malloc(9*sizeof *p->matrix);
  if(p->matrix==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for temporary array to keep "
          "the matrix option", 9*sizeof *p->matrix);

  /* Go over the string and set the values. */
  t=p->up.matrixstring;
  while(*t!='\0')
    {
      switch(*t)
        {
        case ' ': case '\t': case ',':
          ++t;
          break;
        default:
          errno=0;
          p->matrix[counter++]=strtod(t, &tailptr);
          if(errno) error(EXIT_FAILURE, errno, "In reading `%s`", t);
          if(tailptr==t)
            error(EXIT_FAILURE, 0, "the provided string `%s' for matrix "
                  "could not be read as a number", t);
          t=tailptr;
          if(counter>9)       /* Note that it was ++'d! */
            error(EXIT_FAILURE, 0, "there are %lu elements in `%s', there "
                  "should be 4 or 9", counter, p->up.matrixstring);
          /*printf("%f, %s\n", p->matrix[counter-1], t);*/
        }
    }

  /* Add the other necessary information: */
  if(counter==4)
    p->ms1=p->ms0=2;
  else if (counter==9)
    p->ms1=p->ms0=3;
  else
    error(EXIT_FAILURE, 0, "there are %lu numbers in the string `%s'! "
          "It should contain 4 or 9 numbers (for a 2 by 2 or 3 by 3 "
          "matrix)", counter, p->up.matrixstring);
}





/* Set the matrix so the image is aligned with the axises. Note that
   WCSLIB automatically fills the CRPI */
void
makealignmatrix(struct imgwarpparams *p)
{
  double A, dx, dy;
  double w[4]={0,0,0,0};

  /* Check if there is only two WCS axises: */
  if(p->wcs->naxis!=2)
    error(EXIT_FAILURE, 0, "the WCS structure of %s (hdu: %s) has %d "
          "axises. For the `--align' option to operate it must be 2",
          p->up.inputname, p->cp.hdu, p->wcs->naxis);

  /* Depending on the type of data, make the input matrix. Note that
     wcs->altlin is actually bit flags, not integers, so we have to compare
     with powers of two. */
  if(p->wcs->altlin==1)
    {
      w[0]=p->wcs->cdelt[0]*p->wcs->pc[0];
      w[1]=p->wcs->cdelt[0]*p->wcs->pc[1];
      w[2]=p->wcs->cdelt[1]*p->wcs->pc[2];
      w[3]=p->wcs->cdelt[1]*p->wcs->pc[3];
    }
  else
    error(EXIT_FAILURE, 0, "currently the `--align' option only recognizes "
          "PCi_j keywords, not any others");

  /* Find the pixel scale along the two dimensions. Note that we will be
     using the scale along the image X axis for both values. */
  gal_wcs_pixel_scale_deg(p->wcs, &dx, &dy);

  /* Allocate space for the matrix: */
  errno=0;
  p->matrix=malloc(4*sizeof *p->matrix);
  if(p->matrix==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->matrix in makealignmatrix",
          4*sizeof *p->matrix);
  p->ms0=p->ms1=2;

  /* Lets call the given WCS orientation `W', the rotation matrix we want
     to find as `X' and the final (aligned matrix) to have just one useful
     value: `a' (which is the pixel scale):

        x0  x1       w0  w1      -a  0
        x2  x3   *   w2  w3   =   0  a

     Let's open up the matrix multiplication, so we can find the `X'
     elements as function of the `W' elements and `a'.

        x0*w0 + x1*w2 = -a                                         (1)
        x0*w1 + x1*w3 =  0                                         (2)
        x2*w0 + x3*w2 =  0                                         (3)
        x2*w1 + x3*w3 =  a                                         (4)

     Let's bring the X with the smaller index in each equation to the left
     side:

        x0 = (-w2/w0)*x1 - a/w0                                    (5)
        x0 = (-w3/w1)*x1                                           (6)
        x2 = (-w2/w0)*x3                                           (7)
        x2 = (-w3/w1)*x3 + a/w1                                    (8)

    Using (5) and (6) we can find x0 and x1, by first eliminating x0:

       (-w2/w0)*x1 - a/w0 = (-w3/w1)*x1 -> (w3/w1 - w2/w0) * x1 = a/w0

    For easy reading/writing, let's define: A = (w3/w1 - w2/w0)

       --> x1 = a / w0 / A
       --> x0 = -1 * x1 * w3 / w1

    Similar to the above, we can find x2 and x3 from (7) and (8):

       (-w2/w0)*x3 = (-w3/w1)*x3 + a/w1 -> (w3/w1 - w2/w0) * x3 = a/w1

       --> x3 = a / w1 / A
       --> x2 = -1 * x3 * w2 / w0

   */
  A = (w[3]/w[1]) - (w[2]/w[0]);
  p->matrix[1] = dx / w[0] / A;
  p->matrix[3] = dx / w[1] / A;
  p->matrix[0] = -1 * p->matrix[1] * w[3] / w[1];
  p->matrix[2] = -1 * p->matrix[3] * w[2] / w[0];

  /* For a check:
  printf("dx: %e\n", dx);
  printf("w:\n");
  printf("  %.8e    %.8e\n", w[0], w[1]);
  printf("  %.8e    %.8e\n", w[2], w[3]);
  printf("x:\n");
  printf("  %.8e    %.8e\n", p->matrix[0], p->matrix[1]);
  printf("  %.8e    %.8e\n", p->matrix[2], p->matrix[3]);
  */
}





/* Create rotation matrix */
void
makebasicmatrix(struct imgwarpparams *p)
{
  struct uiparams *up=&p->up;

  /* Allocate space for the matrix: */
  errno=0;
  p->matrix=malloc(4*sizeof *p->matrix);
  if(p->matrix==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->matrix in makerotationmatrix",
          4*sizeof *p->matrix);
  p->ms0=p->ms1=2;

  /* Put in the elements */
  if(up->rotateset)
    {
      p->matrix[1] = sin( up->rotate*M_PI/180 );
      p->matrix[2] = p->matrix[1] * -1.0f;
      p->matrix[3] = p->matrix[0] = cos( up->rotate*M_PI/180 );
    }
  else if (up->scaleset)
    {
      p->matrix[1] = p->matrix[2] = 0.0f;
      p->matrix[0] = p->matrix[3] = up->scale;
    }
}





/* Fill in the warping matrix elements based on the options/arguments */
void
preparematrix(struct imgwarpparams *p)
{
  double *tmp;
  int mcheck=0;
  struct uiparams *up=&p->up;

  /* Make sure that none of the matrix creation options/arguments are given
     together. Note that a matrix string is optional (will only be used if
     there is no matrix file). */
  mcheck = ( up->alignset + up->rotateset + up->scaleset +
             (up->matrixname!=NULL) );
  if( mcheck > 1 )
    error(EXIT_FAILURE, 0, "More than one method to define the warping "
          "matrix has been given. Please only specify one.");


  /* Read the input image WCS structure. We are doing this here because
     some of the matrix operations might need it. */
  gal_fits_read_wcs(up->inputname, p->cp.hdu, p->hstartwcs,
                    p->hendwcs, &p->nwcs, &p->wcs);


  /* Depending on the given situation make the matrix. */
  if(up->alignset)
    makealignmatrix(p);
  else if( up->rotateset || up->scaleset)
    makebasicmatrix(p);
  else
    {
      if(up->matrixname)
        gal_txtarray_txt_to_array(up->matrixname, &p->matrix,
                                  &p->ms0, &p->ms1);
      else
        {
          /* Check if a matrix string is actually present. */
          if(up->matrixstringset==0)
            error(EXIT_FAILURE, 0, "no warping matrix string has been and "
                  "no other means of making the warping matrix (a file or "
                  "other options have been specified");
          readmatrixoption(p);
        }
  }


  /* Convert a 2 by 2 matrix into a 3 by 3 matrix and also correct it for
     the FITS definition. */
  if(p->ms0==2 && p->ms1==2)
    {
      errno=0;
      tmp=malloc(9*sizeof *tmp);
      if(tmp==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for 3 by 3 matrix",
              9*sizeof *tmp);

      /* Put the four identical elements in place. */
      tmp[0]=p->matrix[0];
      tmp[1]=p->matrix[1];
      tmp[3]=p->matrix[2];
      tmp[4]=p->matrix[3];

      /* Add the other elements. Note that we need to correct for the FITS
         standard that defines the coordinates of the center of the first
         pixel in the image to be 1. We want 0 to be the coordinates of the
         bottom corner of the image.

         1  0  0.5     a  b  0     a  b  0.5
         0  1  0.5  *  c  d  0  =  c  d  0.5
         0  0   1      0  0  1     0  0   1

         and

         a  b  0.5     1  0  -0.5     a  b  (a*-0.5)+(b*-0.5)+0.5
         c  d  0.5  *  0  1  -0.5  =  c  d  (c*-0.5)+(d*-0.5)+0.5
         0  0   1      0  0   1       0  0           1
      */
      tmp[8] = 1.0f;
      tmp[6] = tmp[7] = 0.0f;
      tmp[2] = ((p->matrix[0] + p->matrix[1]) * -0.5f) + 0.5f;
      tmp[5] = ((p->matrix[2] + p->matrix[3]) * -0.5f) + 0.5f;

      /* Free the previously allocated 2D matrix and put the put the newly
         allocated array with correct values in it. */
      free(p->matrix);
      p->matrix=tmp;

      /* Set the new sizes */
      p->ms0=p->ms1==3;
    }
  else if (p->ms0!=3 || p->ms1!=3)
    error(EXIT_FAILURE, 0, "a bug! please contact us at %s so we can "
          "address the problem. For some reason p->ms0=%lu and p->ms1=%lu! "
          "They should both have a value of 3.", PACKAGE_BUGREPORT,
          p->ms0, p->ms1);
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct imgwarpparams *p)
{
  double *d, *df, *m=p->matrix;

  /* Make sure the input file exists. */
  gal_checkset_check_file(p->up.inputname);

  /* Set the output name: */
  if(p->cp.output)
    gal_checkset_check_remove_file(p->cp.output, p->cp.dontdelete);
  else
    gal_checkset_automatic_output(p->up.inputname, "_warped.fits",
                                  p->cp.removedirinfo, p->cp.dontdelete,
                                  &p->cp.output);


  /* Check the size of the input matrix, note that it might only have
     the wrong numbers when it is read from a file. */
  if(p->up.matrixname)
     if( (p->ms0 != 2 && p->ms0 != 3) || p->ms0 != p->ms1 )
       error(EXIT_FAILURE, 0, "the given matrix in %s has %lu rows and "
             "%lu columns. Its size must be either 2x2 or 3x3",
             p->up.matrixname, p->ms0, p->ms1);

  /* Check if there are any non-normal numbers in the matrix: */
  df=(d=m)+p->ms0*p->ms1;
  do
    if(!isfinite(*d++))
      error(EXIT_FAILURE, 0, "%f is not a `normal' number", *(d-1));
  while(d<df);

  /* Check if the determinant is not zero: */
  if( ( p->ms0==2 && (m[0]*m[3] - m[1]*m[2] == 0) )
      ||
      ( p->ms0==3 &&
        (m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7]
         - m[2]*m[4]*m[6] - m[1]*m[3]*m[8] - m[0]*m[5]*m[7] == 0) ) )
      error(EXIT_FAILURE, 0, "the determinant of the given matrix "
            "is zero");

  /* Check if the transformation is spatially invariant */
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* It is important that the image names are stored in an array (for
   WCS mode in particular). We do that here. */
void
preparearrays(struct imgwarpparams *p)
{
  void *array;
  size_t numnul;
  double *inv, *m=p->matrix;

  /* Read in the input image: */
  numnul=gal_fits_hdu_to_array(p->up.inputname, p->cp.hdu,
                                     &p->inputbitpix, &array, &p->is0,
                                     &p->is1);
  if(p->inputbitpix==DOUBLE_IMG)
    p->input=array;
  else
    {
      gal_fits_change_type(array, p->inputbitpix, p->is0*p->is1, numnul,
                                (void **)&p->input, DOUBLE_IMG);
      free(array);
    }

  /* Make the inverse matrix: */
  errno=0;
  p->inverse=inv=malloc(9*sizeof *inv);
  if(inv==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for the inverse array",
          9*sizeof *inv);
  inv[0]=m[4]*m[8]-m[5]*m[7];
  inv[1]=m[2]*m[7]-m[1]*m[8];
  inv[2]=m[1]*m[5]-m[2]*m[4];
  inv[3]=m[5]*m[6]-m[3]*m[8];
  inv[4]=m[0]*m[8]-m[2]*m[6];
  inv[5]=m[2]*m[3]-m[0]*m[5];
  inv[6]=m[3]*m[7]-m[4]*m[6];
  inv[7]=m[1]*m[6]-m[0]*m[7];
  inv[8]=m[0]*m[4]-m[1]*m[3];
  /* Just for a test:
  {
    size_t i;
    printf("\nInput matrix:");
    for(i=0;i<9;++i) { if(i%3==0) printf("\n"); printf("%-10.5f", m[i]); }
    printf("\n-----------\n");
    printf("Inverse matrix:");
    for(i=0;i<9;++i) { if(i%3==0) printf("\n"); printf("%-10.5f", inv[i]); }
    printf("\n\n");
  }
  */
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct imgwarpparams *p)
{
  struct gal_commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = num_processors(NPROC_CURRENT);
  cp->removedirinfo = 1;
  p->correctwcs     = 1;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    GAL_CONFIGFILES_REPORT_PARAMETERS_SET;

  /* Read the input matrix */
  preparematrix(p);

  /* Do a sanity check. */
  sanitycheck(p);
  gal_checkset_check_remove_file(GAL_TXTARRAY_LOG, 0);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf(" Using %lu CPU thread%s\n", p->cp.numthreads,
             p->cp.numthreads==1 ? "." : "s.");
      printf(" Input image: %s\n", p->up.inputname);
      printf(" matrix:"
             "\n\t%.4f   %.4f   %.4f"
             "\n\t%.4f   %.4f   %.4f"
             "\n\t%.4f   %.4f   %.4f\n",
             p->matrix[0], p->matrix[1], p->matrix[2],
             p->matrix[3], p->matrix[4], p->matrix[5],
             p->matrix[6], p->matrix[7], p->matrix[8]);
    }

  /* Make the array of input images. */
  preparearrays(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct imgwarpparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->input);
  free(p->matrix);
  free(p->cp.hdu);
  free(p->inverse);
  free(p->cp.output);

  if(p->wcs)
    wcsvfree(&p->nwcs, &p->wcs);

  /* Print the final message. */
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in: ", 0);
}
