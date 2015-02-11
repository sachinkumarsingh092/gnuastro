/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <setjmp.h>
#include <string.h>
#include <jpeglib.h>

#include "main.h"
#include "jpeg.h"







/*************************************************************
 **************      Acceptable JPEG names      **************
 *************************************************************/
int
nameisjpeg(char *name)
{
  size_t len;
  len=strlen(name);
  if (strcmp(&name[len-4], ".jpg") == 0
      || strcmp(&name[len-4], ".JPG") == 0
      || strcmp(&name[len-4], ".jpeg") == 0
      || strcmp(&name[len-5], ".JPEG") == 0
      || strcmp(&name[len-5], ".jpe") == 0
      || strcmp(&name[len-5], ".jif") == 0
      || strcmp(&name[len-5], ".jfif") == 0
      || strcmp(&name[len-5], ".jfi") == 0)
    return 1;
  else
    return 0;
}




















/*************************************************************
 **************        Read a JPEG image        **************
 *************************************************************/
/* Read the example.c in libjpeg's source code to understand the
   details of what is going on here.  */
struct my_error_mgr
{
  struct jpeg_error_mgr pub;	/* "public" fields */
  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr *my_error_ptr;


METHODDEF(void)
jpeg_error_exit(j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}





void
makejsample(JSAMPLE **a, size_t size)
{
  JSAMPLE *jsarr;

  if(sizeof *jsarr!=1)
    {
      printf("\n\nJSAMPLE has to be unsigned char!\n\n");
      exit(EXIT_FAILURE);
    }

  errno=0;
  jsarr=malloc(size*sizeof *jsarr);
  if(jsarr==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for jsarr", size*sizeof *jsarr);

  *a=jsarr;
}





double **
readjpg(char *inname, size_t *outs0, size_t *outs1, size_t *numcolors)
{
  double **all;
  FILE * infile;
  JSAMPROW jrow;
  JSAMPLE *jsamp;
  int rowstride, c;
  JSAMPARRAY jsarr;
  struct my_error_mgr jerr;
  size_t i, j, size, nc, s0, s1;
  struct jpeg_decompress_struct cinfo;

  /* Open the input file */
  if ((infile = fopen(inname, "rb")) == NULL)
    {
      fprintf(stderr, "\n\nError: Can't open %s\n", inname);
      exit(EXIT_FAILURE);
    }

  /* Set up the error and decompressing (reading) functions. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = jpeg_error_exit;
  if (setjmp(jerr.setjmp_buffer))
    {
      printf("\n\nError in JPEG code.\n\n");
      jpeg_destroy_decompress(&cinfo);
      fclose(infile);
      exit(EXIT_FAILURE);
    }
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, infile);

  /* Read the JPEG header information and start de-compressing: */
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  /* Get the array width and height and number of color channels: */
  s0=*outs0=cinfo.output_height;
  s1=*outs1=cinfo.output_width;
  size=s0*s1;
  nc=*numcolors=cinfo.output_components;
  rowstride=s1*nc;
  makejsample(&jsamp, size*nc);

  /* Allocate all the arrays for each color: */
  errno=0;
  all=malloc(nc*sizeof *all);
  if(all==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for all", nc*sizeof *all);
  for(i=0;i<nc;++i)
    {
      errno=0;
      all[i]=malloc(s0*s1*sizeof *all[i]);
      if(all[i]==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for all[%lu]",
              s0*s1*sizeof *all[i], i);
    }

  /* Read the image line by line: */
  c=s0-1;
  while (cinfo.output_scanline < cinfo.output_height)
    {
      jrow=&jsamp[c-- * rowstride];
      jsarr=&jrow;
      jpeg_read_scanlines(&cinfo, jsarr, 1);
    }

  /* Put the different colors into the different arrays */
  for(i=0;i<size;++i)
    for(j=0;j<nc;++j)
      all[j][i]=jsamp[i*nc+j];

  /* Finish decompression, destroy it and close file: */
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);
  free(jsamp);

  return all;
}





void
preparejpeg(struct converttparams *p, char *filename)
{
  double **allcolors;
  size_t i, s0, s1, numcolors;

  /* Read the JPEG image into the all array. */
  allcolors=readjpg(filename, &s0, &s1, &numcolors);

  /* Check if the number of colors that will be added with this JPEG
     image does not exceed 4. */
  if(p->numch+numcolors>4)
    error(EXIT_FAILURE, 0, "The number of channels in %s added with the "
          "previous inputs will exceed 4 (the maximum number of color "
          "channels). Can't continue.", filename);

  /* Add each channel to the final set of channels. */
  for(i=0;i<numcolors;++i)
    {
      p->s0[p->numch]=s0;
      p->s1[p->numch]=s1;
      p->ch[p->numch]=allcolors[i];
      ++p->numch;
    }

  /* Free the array keeping the pointers to each channel. Note that
     each channel was allocated separately so we can safely remove it.*/
  free(allcolors);
}
