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

#include "fitsarrayvv.h"

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
  if ( ( len>=3 && strcmp(&name[len-3], "jpg") == 0 )
       || ( len>=3 && strcmp(&name[len-3], "JPG") == 0 )
       || ( len>=4 && strcmp(&name[len-4], "jpeg") == 0 )
       || ( len>=4 && strcmp(&name[len-4], "JPEG") == 0 )
       || ( len>=3 && strcmp(&name[len-3], "jpe") == 0 )
       || ( len>=3 && strcmp(&name[len-3], "jif") == 0 )
       || ( len>=4 && strcmp(&name[len-4], "jfif") == 0 )
       || ( len>=3 && strcmp(&name[len-3], "jfi") == 0 ) )
    return 1;
  else
    return 0;
}





int
nameisjpegsuffix(char *name)
{
  if (strcmp(name, "jpg") == 0   || strcmp(name, ".jpg") == 0
      || strcmp(name, "JPG") == 0 || strcmp(name, ".JPG") == 0
      || strcmp(name, "jpeg") == 0 || strcmp(name, ".jpeg") == 0
      || strcmp(name, "JPEG") == 0 || strcmp(name, ".JPEG") == 0
      || strcmp(name, "jpe") == 0 || strcmp(name, ".jpe") == 0
      || strcmp(name, "jif") == 0 || strcmp(name, ".jif") == 0
      || strcmp(name, "jfif") == 0 || strcmp(name, ".jfif") == 0
      || strcmp(name, "jfi") == 0 || strcmp(name, ".jfi") == 0)
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
  errno=0;
  if ((infile = fopen(inname, "rb")) == NULL)
    error(EXIT_FAILURE, errno, inname);

  /* Set up the error and decompressing (reading) functions. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = jpeg_error_exit;
  if (setjmp(jerr.setjmp_buffer))
    {
      jpeg_destroy_decompress(&cinfo);
      fclose(infile);
      error(EXIT_FAILURE, 0, "Problem in reading %s", inname);
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
      p->bitpixs[p->numch]=BYTE_IMG;
      ++p->numch;
    }

  /* Free the array keeping the pointers to each channel. Note that
     each channel was allocated separately so we can safely remove it.*/
  free(allcolors);
}




















/*************************************************************
 **************       Write a JPEG image        **************
 *************************************************************/
void
writejpeg(JSAMPLE *jsr, struct converttparams *p)
{
  JSAMPROW r[1];
  FILE * outfile;
  int row_stride, c;
  struct jpeg_error_mgr jerr;
  size_t s0=p->s0[0], s1=p->s1[0];
  struct jpeg_compress_struct cinfo;

  /* Begin the JPEG writing, following libjpeg's example.c  */
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  errno=0;
  if ((outfile = fopen(p->cp.output, "wb")) == NULL)
    error(EXIT_FAILURE, errno, p->cp.output);
  jpeg_stdio_dest(&cinfo, outfile);

  cinfo.image_width  = s1;
  cinfo.image_height = s0;
  switch(p->numch)
    {
    case 1:
      row_stride=s1;
      cinfo.input_components = 1;
      cinfo.in_color_space = JCS_GRAYSCALE;
      break;
    case 3:
      row_stride=3*s1;
      cinfo.input_components = 3;
      cinfo.in_color_space = JCS_RGB;
      break;
    case 4:
      row_stride=4*s1;
      cinfo.input_components = 4;
      cinfo.in_color_space = JCS_CMYK;
      break;
    default:
      error(EXIT_FAILURE, 0, "A bug! The number of channels in writejpeg "
            "is not 1, 3 or 4, but %lu. This should not happen. Please "
            "contact us so we can fix the problem.", p->numch);
    }

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, p->quality, TRUE);
  cinfo.density_unit=1;
  cinfo.Y_density=cinfo.X_density=s1/(p->widthincm/2.54);
  jpeg_start_compress(&cinfo, TRUE);

  /* cinfo.next_scanline is 'unsigned int' */
  c=s0-1; /* In FITS the first row is on the bottom!  */
  while (cinfo.next_scanline < cinfo.image_height)
    {
      r[0] = & jsr[c-- * row_stride];
      (void) jpeg_write_scanlines(&cinfo, r, 1);
    }

  jpeg_finish_compress(&cinfo);
  fclose(outfile);
  jpeg_destroy_compress(&cinfo);
}





void
savejpeg(struct converttparams *p)
{
  JSAMPLE *jsr;
  uint8_t **ech=p->ech;
  size_t pixel, color, size, numch=p->numch;

  /* Make sure the JSAMPLE is 8bits, then allocate the necessary space
     based on the number of channels. */
  if(sizeof *jsr!=1)
    error(EXIT_FAILURE, 0, "JSAMPLE has to be 8bit!");
  size=p->numch*p->s0[0]*p->s1[0];
  errno=0;
  jsr=malloc(size*sizeof *jsr);
  if(jsr==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for jsr", size*sizeof *jsr);

  /* Write the different colors into jsr. */
  size=p->s0[0]*p->s1[0];
  for(pixel=0;pixel<size;++pixel)
    for(color=0;color<numch;++color)
      {
        jsr[pixel*numch+color]=ech[color][pixel];
        /*
        printf("color: %lu, pixel: %lu, jsr: %d\n", color, pixel,
               (int)jsr[pixel*numch+color]);
        */
      }


  /* Write jsr to a JPEG image. */
  writejpeg(jsr, p);

  free(jsr);
}
