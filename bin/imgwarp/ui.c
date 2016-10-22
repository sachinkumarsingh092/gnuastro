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
    error(EXIT_FAILURE, errno, "ui.c: %zu bytes in readdefaults",
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
      else if(strcmp(name, "nofitscorrect")==0)
        {
          if(up->nofitscorrectset) continue;
          gal_checkset_int_zero_or_one(value, &up->nofitscorrect, name,
                                       key, SPACK, filename, lineno);
          up->nofitscorrectset=1;
        }





      /* Modular warpings */
      else if(strcmp(name, "align")==0)
        add_to_optionwapsll(&p->up.owll, ALIGN_WARP, NULL);

      else if(strcmp(name, "rotate")==0)
        add_to_optionwapsll(&p->up.owll, ROTATE_WARP, value);

      else if(strcmp(name, "scale")==0)
        add_to_optionwapsll(&p->up.owll, SCALE_WARP, value);

      else if(strcmp(name, "flip")==0)
        add_to_optionwapsll(&p->up.owll, FLIP_WARP, value);

      else if(strcmp(name, "shear")==0)
        add_to_optionwapsll(&p->up.owll, SHEAR_WARP, value);

      else if(strcmp(name, "translate")==0)
        add_to_optionwapsll(&p->up.owll, TRANSLATE_WARP, value);

      else if(strcmp(name, "project")==0)
        add_to_optionwapsll(&p->up.owll, PROJECT_WARP, value);



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
    fprintf(fp, CONF_SHOWFMT"%zu\n", "hstartwcs", p->hstartwcs);
  if(up->hendwcsset)
    fprintf(fp, CONF_SHOWFMT"%zu\n", "hendwcs", p->hendwcs);

  fprintf(fp, "\n# Output parameters:\n");
  if(up->matrixstringset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("matrix", up->matrixstring);

  if(cp->outputset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("output", cp->output);

  if(up->maxblankfracset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "maxblankfrac", p->maxblankfrac);



  fprintf(fp, "\n# Modular transformations:\n");
  if(up->nofitscorrectset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "nofitscorrect", up->nofitscorrect);



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
/**********      Modular matrix linked list       *************/
/**************************************************************/
void
add_to_optionwapsll(struct optionwarpsll **list, int type, char *value)
{
  double v1=NAN, v2=NAN;
  char *tailptr, *secondstr;
  struct optionwarpsll *newnode;

  /* Allocate the necessary space. */
  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for newnode in "
          "add_to_optionwarpsll", sizeof *newnode);

  /* Read the numbers when necessary. */
  if(value)
    {
      /* Parse the first number */
      v1=strtod(value, &tailptr);
      if(tailptr==value)
        error(EXIT_FAILURE, 0, "The start of the string `%s' could not be "
              "read as a number", value);

      /* If there is any white space characters, ignore them and make sure
         that the first character is a coma (`,'). */
      secondstr=tailptr;
      while(isspace(*secondstr)) ++secondstr;
      if(*secondstr==',')
        {
          /* If the type is rotate, then print an error since rotate only
             needs one input, not two. */
          if(type==ROTATE_WARP)
            error(EXIT_FAILURE, 0, "The `--rotate' (`-r') option only needs "
                  "one input number, not more. It was given `%s'", value);

          /* Ignore the coma. */
          ++secondstr;

          /* Read the second number: */
          v2=strtod(secondstr, &tailptr);
          if(tailptr==secondstr)
            error(EXIT_FAILURE, 0, "The second part (after the coma) of "
                  "`%s' (`%s') could not be read as a number", value,
                  secondstr);
        }

      /* If there was only one number given, secondstr will be '\0' when
         control reaches here. */
      else if(*secondstr!='\0')
        error(EXIT_FAILURE, 0, "the character between the two numbers (`%s') "
              "must be a coma (`,')\n", value);
    }

  /* Put in the values. Note that both v1 and v2 were initialized to NaN,
     so if v2 is not given, it will be NaN and the later function can
     decide what it wants to replace it with.*/
  newnode->v1=v1;
  newnode->v2=v2;
  newnode->type=type;
  newnode->next=*list;

  /* Set list to point to the new node. */
  *list=newnode;
}





/* Allocate space for a new node: */
struct optionwarpsll *
alloc_owll_node(void)
{
  struct optionwarpsll *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for newnode in "
          "add_to_optionwarpsll", sizeof *newnode);

  return newnode;
}





/* The input list of warpings are recorded in a last-in-first-out order. So
   we reverse the order and also add the transformations necessary for the
   FITS definition (where the center of the pixel has a value of 1, not its
   corner. */
void
prepare_optionwapsll(struct imgwarpparams *p)
{
  struct optionwarpsll *tmp, *next, *newnode, *prepared=NULL;

  /* Add the FITS correction for the first warp (before everything else).

     IMPORTANT: This is the last transform that will be done, so we have to
     translate the image by -0.5.*/
  if(!p->up.nofitscorrect)
    {
      newnode = alloc_owll_node();
      newnode->v1 = newnode->v2 = -0.5f;
      newnode->type = TRANSLATE_WARP;
      newnode->next = prepared;
      prepared = newnode;
    }

  /* Put in the rest of the warpings */
  tmp=p->up.owll;
  while(tmp!=NULL)
    {
      /* Allocate space for the new element, and put the values in. */
      newnode = alloc_owll_node();
      newnode->v1 = tmp->v1;
      newnode->v2 = tmp->v2;
      newnode->type = tmp->type;
      newnode->next = prepared;

      /* Now that previous nodes have been linked to next, set the
         reversed to the new node. */
      prepared = newnode;

      /* Now keep the next element and free the old allocated space. */
      next = tmp->next;
      free(tmp);
      tmp = next;
    }

  /* Add the FITS correction for the last warp (after everything else).

     IMPORTANT: This is the first transform that will be done, so we have
     to translate the image by +0.5.*/
  if(!p->up.nofitscorrect)
    {
      newnode = alloc_owll_node();
      newnode->v1 = newnode->v2 = 0.5f;
      newnode->type = TRANSLATE_WARP;
      newnode->next = prepared;
      prepared = newnode;
    }

  /* Put the pointer in the output */
  p->up.owll=prepared;
}




















/**************************************************************/
/*************      Fill temporary matrix     *****************/
/**************************************************************/
void
read_matrix(struct imgwarpparams *p)
{
  char *t, *tailptr;
  size_t i, counter=0, m0, m1;
  double *matrix=p->matrix, rmatrix[9], *fmatrix;

  /* Read the matrix either as a file or from the command-line. */
  if(p->up.matrixname)
    {
      gal_txtarray_txt_to_array(p->up.matrixname, &fmatrix, &m0, &m1);
      counter=m0*m1;
      for(i=0;i<(counter<9 ? counter : 9);++i)
        rmatrix[i]=fmatrix[i];
      free(fmatrix);
    }
  else
    {
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
              rmatrix[counter++]=strtod(t, &tailptr);
              if(errno) error(EXIT_FAILURE, errno, "reading `%s`", t);
              if(tailptr==t)
                error(EXIT_FAILURE, 0, "the provided string `%s' for "
                      "matrix could not be read as a number", t);
              t=tailptr;
              if(counter>9)       /* Note that it was incremented! */
                error(EXIT_FAILURE, 0, "there are %zu elements in `%s', "
                      "there should be 4 or 9", counter, p->up.matrixstring);
              /*printf("%f, %s\n", p->matrix[counter-1], t);*/
            }
        }
    }

  /* If there was 4 elements (a 2 by 2 matrix), put them into a 3 by 3
     matrix. */
  if(counter==4)
    {
      /* Fill in the easy 3 by 3 matrix: */
      matrix[0]=rmatrix[0];   matrix[1]=rmatrix[1];
      matrix[3]=rmatrix[2];   matrix[4]=rmatrix[3];
      matrix[6]=0.0f;         matrix[7]=0.0f;         matrix[8]=1.0f;

      /* If we need to correct for the FITS standard, then correc the last
         two elements. Recall that the coordinates of the center of the
         first pixel in the FITS standard are 1. We want 0 to be the
         coordinates of the bottom corner of the image.

         1  0  0.5      a  b  0      a  b  0.5
         0  1  0.5   *  c  d  0   =  c  d  0.5
         0  0   1       0  0  1      0  0   1

         and

         a  b  0.5     1  0  -0.5     a  b  (a*-0.5)+(b*-0.5)+0.5
         c  d  0.5  *  0  1  -0.5  =  c  d  (c*-0.5)+(d*-0.5)+0.5
         0  0   1      0  0   1       0  0           1
      */
      if(p->up.nofitscorrect)
        matrix[2] = matrix[5] = 0.0f;
      else
        {
          matrix[2] = ((rmatrix[0] + rmatrix[1]) * -0.5f) + 0.5f;
          matrix[5] = ((rmatrix[2] + rmatrix[3]) * -0.5f) + 0.5f;
        }
    }
  else if (counter==9)
    {
      matrix[0]=rmatrix[0];   matrix[1]=rmatrix[1];   matrix[2]=rmatrix[2];
      matrix[3]=rmatrix[3];   matrix[4]=rmatrix[4];   matrix[5]=rmatrix[5];
      matrix[6]=rmatrix[6];   matrix[7]=rmatrix[7];   matrix[8]=rmatrix[8];
    }
  else
    error(EXIT_FAILURE, 0, "there are %zu numbers in the string `%s'! "
          "It should contain 4 or 9 numbers (for a 2 by 2 or 3 by 3 "
          "matrix)", counter, p->up.matrixstring);
}





/* Set the matrix so the image is aligned with the axises. Note that
   WCSLIB automatically fills the CRPI */
void
makealignmatrix(struct imgwarpparams *p, double *tmatrix)
{
  double A, dx, dy;
  double amatrix[4];
  double w[4]={0,0,0,0};

  /* Check if there is only two WCS axises: */
  if(p->wcs->naxis!=2)
    error(EXIT_FAILURE, 0, "the WCS structure of %s (hdu: %s) has %d "
          "axises. For the `--align' option to operate it must be 2",
          p->up.inputname, p->cp.hdu, p->wcs->naxis);

  /* Depending on the type of data, make the input matrix. Note that
     wcs->altlin is actually bit flags, not integers, so we have to compare
     with powers of two. */
  if(p->wcs->altlin |= 1)
    {
      w[0]=p->wcs->cdelt[0]*p->wcs->pc[0];
      w[1]=p->wcs->cdelt[0]*p->wcs->pc[1];
      w[2]=p->wcs->cdelt[1]*p->wcs->pc[2];
      w[3]=p->wcs->cdelt[1]*p->wcs->pc[3];
    }
  if(p->wcs->altlin |= 2)
    {
      w[0]=p->wcs->cd[0];
      w[1]=p->wcs->cd[1];
      w[2]=p->wcs->cd[2];
      w[3]=p->wcs->cd[3];
    }
  else
    error(EXIT_FAILURE, 0, "currently the `--align' option only recognizes "
          "PCi_ja and CDi_ja keywords, not any others");

  /* Find the pixel scale along the two dimensions. Note that we will be
     using the scale along the image X axis for both values. */
  gal_wcs_pixel_scale_deg(p->wcs, &dx, &dy);

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

    Note that when the image is already aligned, a unity matrix should be
    output.
   */
  if( w[1]==0.0f && w[2]==0.0f )
    {
      amatrix[0]=1.0f;   amatrix[1]=0.0f;
      amatrix[2]=0.0f;   amatrix[3]=1.0f;
    }
  else
    {
      A = (w[3]/w[1]) - (w[2]/w[0]);
      amatrix[1] = dx / w[0] / A;
      amatrix[3] = dx / w[1] / A;
      amatrix[0] = -1 * amatrix[1] * w[3] / w[1];
      amatrix[2] = -1 * amatrix[3] * w[2] / w[0];
    }


  /* For a check:
  printf("dx: %e\n", dx);
  printf("w:\n");
  printf("  %.8e    %.8e\n", w[0], w[1]);
  printf("  %.8e    %.8e\n", w[2], w[3]);
  printf("x:\n");
  printf("  %.8e    %.8e\n", amatrix[0], amatrix[1]);
  printf("  %.8e    %.8e\n", amatrix[2], amatrix[3]);
  */


  /* Put the matrix elements into the output array: */
  tmatrix[0]=amatrix[0];  tmatrix[1]=amatrix[1]; tmatrix[2]=0.0f;
  tmatrix[3]=amatrix[2];  tmatrix[4]=amatrix[3]; tmatrix[5]=0.0f;
  tmatrix[6]=0.0f;        tmatrix[7]=0.0f;       tmatrix[8]=1.0f;
}



















/**************************************************************/
/***************       Prepare Matrix       *******************/
/**************************************************************/
/* This function is mainly for easy checking/debugging. */
void
printmatrix(double *matrix)
{
  printf("%-10.3f%-10.3f%-10.3f\n", matrix[0], matrix[1], matrix[2]);
  printf("%-10.3f%-10.3f%-10.3f\n", matrix[3], matrix[4], matrix[5]);
  printf("%-10.3f%-10.3f%-10.3f\n", matrix[6], matrix[7], matrix[8]);
}





void
inplace_matrix_multiply(double *in, double *with)
{
  /* `tin' will keep the values of the input array because we want to
     write the multiplication result in the input array. */
  double tin[9]={in[0],in[1],in[2],in[3],in[4],in[5],in[6],in[7],in[8]};

  /* For easy checking, here are the matrix/memory layouts:
          tin[0] tin[1] tin[2]     with[0] with[1] with[2]
          tin[3] tin[4] tin[5]  *  with[3] with[4] with[5]
          tin[6] tin[7] tin[8]     with[6] with[7] with[8]   */
  in[0] = tin[0]*with[0] + tin[1]*with[3] + tin[2]*with[6];
  in[1] = tin[0]*with[1] + tin[1]*with[4] + tin[2]*with[7];
  in[2] = tin[0]*with[2] + tin[1]*with[5] + tin[2]*with[8];

  in[3] = tin[3]*with[0] + tin[4]*with[3] + tin[5]*with[6];
  in[4] = tin[3]*with[1] + tin[4]*with[4] + tin[5]*with[7];
  in[5] = tin[3]*with[2] + tin[4]*with[5] + tin[5]*with[8];

  in[6] = tin[6]*with[0] + tin[7]*with[3] + tin[8]*with[6];
  in[7] = tin[6]*with[1] + tin[7]*with[4] + tin[8]*with[7];
  in[8] = tin[6]*with[2] + tin[7]*with[5] + tin[8]*with[8];
}





/* Fill in the warping matrix elements based on the options/arguments */
void
prepare_modular_matrix(struct imgwarpparams *p)
{
  int f1, f2;                   /* For flipping. */
  double s, c, tmatrix[9];
  struct uiparams *up=&p->up;
  struct optionwarpsll *tmp, *next;


  /* Allocate space for the matrix, then initialize it. */
  p->matrix[0]=1.0f;     p->matrix[1]=0.0f;     p->matrix[2]=0.0f;
  p->matrix[3]=0.0f;     p->matrix[4]=1.0f;     p->matrix[5]=0.0f;
  p->matrix[6]=0.0f;     p->matrix[7]=0.0f;     p->matrix[8]=1.0f;


  /* The linked list is last-in-first-out, so we need to reverse it to
     easily apply the changes in the same order that was read in. */
  prepare_optionwapsll(p);


  /* Do all the operations */
  tmp=up->owll;
  while(tmp!=NULL)
    {
      /* Fill `tmatrix' depending on the type of the warp. */
      switch(tmp->type)
        {
        case ALIGN_WARP:
          makealignmatrix(p, tmatrix);
          break;

        case ROTATE_WARP:
          s = sin( tmp->v1*M_PI/180 );
          c = cos( tmp->v1*M_PI/180 );
          tmatrix[0]=c;        tmatrix[1]=s;     tmatrix[2]=0.0f;
          tmatrix[3]=-1.0f*s;  tmatrix[4]=c;     tmatrix[5]=0.0f;
          tmatrix[6]=0.0f;     tmatrix[7]=0.0f;  tmatrix[8]=1.0f;
          break;

        case SCALE_WARP:
          if( isnan(tmp->v2) ) tmp->v2=tmp->v1;
          tmatrix[0]=tmp->v1;  tmatrix[1]=0.0f;     tmatrix[2]=0.0f;
          tmatrix[3]=0.0f;     tmatrix[4]=tmp->v2;  tmatrix[5]=0.0f;
          tmatrix[6]=0.0f;     tmatrix[7]=0.0f;     tmatrix[8]=1.0f;
          break;

        case FLIP_WARP:
          /* For the flip, the values dont really matter! As long as the
             value is non-zero, the flip in the respective axis will be
             made. Note that the second axis is optional (can be NaN), but
             the first axis is required.*/
          f1 = tmp->v1==0.0f ? 0 : 1;
          f2 = isnan(tmp->v2) ? 0 : ( tmp->v2==0.0f ? 0 : 1);
          if( f1 && !f2  )
            {
              tmatrix[0]=1.0f;   tmatrix[1]=0.0f;
              tmatrix[3]=0.0f;   tmatrix[4]=-1.0f;
            }
          else if ( !f1 && f2 )
            {
              tmatrix[0]=-1.0f;  tmatrix[1]=0.0f;
              tmatrix[3]=0.0f;   tmatrix[4]=1.0f;
            }
          else
            {
              tmatrix[0]=-1.0f;  tmatrix[1]=0.0f;
              tmatrix[3]=0.0f;   tmatrix[4]=-1.0f;
            }
                                                      tmatrix[2]=0.0f;
                                                      tmatrix[5]=0.0f;
          tmatrix[6]=0.0f;       tmatrix[7]=0.0f;     tmatrix[8]=1.0f;
          break;

        case SHEAR_WARP:
          if( isnan(tmp->v2) ) tmp->v2=tmp->v1;
          tmatrix[0]=1.0f;     tmatrix[1]=tmp->v1;    tmatrix[2]=0.0f;
          tmatrix[3]=tmp->v2;  tmatrix[4]=1.0f;       tmatrix[5]=0.0f;
          tmatrix[6]=0.0f;     tmatrix[7]=0.0f;       tmatrix[8]=1.0f;
          break;

        case TRANSLATE_WARP:
          if( isnan(tmp->v2) ) tmp->v2=tmp->v1;
          tmatrix[0]=1.0f;     tmatrix[1]=0.0f;       tmatrix[2]=tmp->v1;
          tmatrix[3]=0.0f;     tmatrix[4]=1.0f;       tmatrix[5]=tmp->v2;
          tmatrix[6]=0.0f;     tmatrix[7]=0.0f;       tmatrix[8]=1.0f;
          break;

        case PROJECT_WARP:
          if( isnan(tmp->v2) ) tmp->v2=tmp->v1;
          tmatrix[0]=1.0f;     tmatrix[1]=0.0f;       tmatrix[2]=0.0f;
          tmatrix[3]=0.0f;     tmatrix[4]=1.0f;       tmatrix[5]=0.0f;
          tmatrix[6]=tmp->v1;  tmatrix[7]=tmp->v2;    tmatrix[8]=1.0f;
          break;

        default:
          error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can "
                "address the problem. For some reason the value of tmp->type "
                "in `prepare_modular_matrix' of ui.c is not recognized. "
                "This is an internal, not a user issue. So please let us "
                "know.", PACKAGE_BUGREPORT);
        }

      /* Multiply this matrix with the main matrix in-place. */
      inplace_matrix_multiply(p->matrix, tmatrix);

      /* Keep the next element and free the node's allocated space. */
      next = tmp->next;
      free(tmp);
      tmp = next;

      /* For a check:
      printf("tmatrix:\n");
      printmatrix(tmatrix);
      printf("out:\n");
      printmatrix(p->matrix);
      */
    }
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

  /* If an actual matrix is given, then it will be used and all modular
     warpings will be ignored. */
  if(p->up.matrixstring || p->up.matrixname)
    read_matrix(p);
  else if (p->up.owll)
    prepare_modular_matrix(p);
  else
    error(EXIT_FAILURE, 0, "No input matrix specified.\n\nPlease either "
          "use the modular warp options like `--rotate' or `--scale', "
          "or directly specify the matrix on the command-line, or in the "
          "configuration files.\n\nRun with `--help' for the full list of "
          "modular warpings (among other options), or see the manual's "
          "`Warping basics' section for more on the matrix.");


  /* Check if there are any non-normal numbers in the matrix: */
  df=(d=p->matrix)+9;
  do
    if(!isfinite(*d++))
      {
        printmatrix(p->matrix);
        error(EXIT_FAILURE, 0, "%f is not a `normal' number in the "
              "input matrix shown above", *(d-1));
      }
  while(d<df);

  /* Check if the determinant is not zero: */
  if( m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7]
      - m[2]*m[4]*m[6] - m[1]*m[3]*m[8] - m[0]*m[5]*m[7] == 0 )
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
    error(EXIT_FAILURE, errno, "%zu bytes for the inverse array",
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
  p->up.owll        = NULL;

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

  /* Read the input image WCS structure. We are doing this here because
     some of the matrix operations might need it. */
  gal_fits_read_wcs(p->up.inputname, p->cp.hdu, p->hstartwcs,
                    p->hendwcs, &p->nwcs, &p->wcs);

  /* Do a sanity check. */
  sanitycheck(p);
  gal_checkset_check_remove_file(GAL_TXTARRAY_LOG, 0);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf(" Using %zu CPU thread%s\n", p->cp.numthreads,
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
  free(p->cp.hdu);
  free(p->inverse);
  free(p->cp.output);

  if(p->wcs)
    wcsvfree(&p->nwcs, &p->wcs);

  /* Print the final message. */
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in: ", 0);
}
