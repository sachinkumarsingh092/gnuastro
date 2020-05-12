/*********************************************************************
jpeg -- functions to read and write JPEG files.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <setjmp.h>
#include <string.h>
#ifdef HAVE_LIBJPEG
  #include <jpeglib.h>
#endif

#include <gnuastro/list.h>
#include <gnuastro/jpeg.h>

#include <gnuastro-internal/checkset.h>





/*************************************************************
 **************      Acceptable JPEG names      **************
 *************************************************************/
int
gal_jpeg_name_is_jpeg(char *name)
{
  size_t len;

  if(name)
    {
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
  else return 0;
}





int
gal_jpeg_suffix_is_jpeg(char *name)
{
  if(name)
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
  else return 0;
}




















/*************************************************************
 **************        Read a JPEG image        **************
 *************************************************************/
#ifdef HAVE_LIBJPEG
/* Read the example.c in libjpeg's source code to understand the
   details of what is going on here.  */
struct my_error_mgr
{
  struct jpeg_error_mgr pub;        /* "public" fields */
  jmp_buf setjmp_buffer;            /* for return to caller */
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





static void
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
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for jsarr",
          __func__, size*sizeof *jsarr);

  *a=jsarr;
}





static unsigned char **
readjpg(char *inname, size_t *outs0, size_t *outs1, size_t *numcolors)
{
  FILE * infile;
  JSAMPROW jrow;
  JSAMPLE *jsamp;
  int rowstride, c;
  JSAMPARRAY jsarr;
  unsigned char **all;
  struct my_error_mgr jerr;
  size_t i, j, size, nc, s0, s1;
  struct jpeg_decompress_struct cinfo;

  /* Open the input file */
  errno=0;
  if ((infile = fopen(inname, "rb")) == NULL)
    error(EXIT_FAILURE, errno, "%s", inname);

  /* Set up the error and decompressing (reading) functions. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = jpeg_error_exit;
  if (setjmp(jerr.setjmp_buffer))
    {
      jpeg_destroy_decompress(&cinfo);
      fclose(infile);
      error(EXIT_FAILURE, 0, "%s: problem in reading %s", __func__, inname);
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
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'all'",
          __func__, nc*sizeof *all);
  for(i=0;i<nc;++i)
    {
      errno=0;
      all[i]=malloc(s0*s1*sizeof *all[i]);
      if(all[i]==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'all[%zu]'",
              __func__, s0*s1*sizeof *all[i], i);
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
#endif  /* HAVE_LIBJPEG */





/* Read each color channel of a JPEG image as a separate array and put them
   in a linked list of data-structures. */
gal_data_t *
gal_jpeg_read(char *filename, size_t minmapsize, int quietmmap)
{
#ifdef HAVE_LIBJPEG
  char *name;
  gal_data_t *out=NULL;
  size_t ndim=2, dsize[2];
  unsigned char **allcolors;
  size_t i, s0, s1, numcolors;

  /* Read the JPEG image into the all array. */
  allcolors=readjpg(filename, &s0, &s1, &numcolors);

  /* Add the arrays to the linked list. */
  for(i=0;i<numcolors;++i)
    {
      dsize[0]=s0;
      dsize[1]=s1;
      if( asprintf(&name, "JPEG_CH_%zu", i+1)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_list_data_add_alloc(&out, allcolors[i], GAL_TYPE_UINT8, ndim,
                              dsize, NULL, 0, minmapsize, quietmmap,
                              name, NULL, NULL);
      free(name);
    }

  /* Free the array keeping the pointers to each channel. Note that each
     channel was allocated separately and goes out of this function with
     the data structure, so we just have to free the outer array that kept
     all the channels. */
  free(allcolors);

  /* Return the number of color channels. */
  return out;
#else
  error(EXIT_FAILURE, 0, "%s: libjpeg was not found during the "
        "configuration of %s on this system. To read from JPEG files, "
        "libjpeg is required. Please install libjpeg and configure, make "
        "and install %s again", __func__, PACKAGE_STRING, PACKAGE_STRING);
  return NULL;
#endif
}




















/*************************************************************
 **************       Write a JPEG image        **************
 *************************************************************/
#ifdef HAVE_LIBJPEG
static void
jpeg_write_array(JSAMPLE *jsr, gal_data_t *in, char *filename,
                 uint8_t quality, float widthincm)
{
  JSAMPROW r[1];
  FILE * outfile;
  int row_stride=0, c;
  size_t *dsize=in->dsize;
  struct jpeg_error_mgr jerr;
  struct jpeg_compress_struct cinfo;
  size_t numch=gal_list_data_number(in);

  /* A small sanity check. */
  if(quality > 100)
    error(EXIT_FAILURE, 0, "%s: quality value %u not acceptable. It must be "
          "a value between zero and 100 (inclusive)", __func__, quality);

  /* Begin the JPEG writing, following libjpeg's example.c  */
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  errno=0;
  if ((outfile = fopen(filename, "wb")) == NULL)
    error(EXIT_FAILURE, errno, "%s", filename);
  jpeg_stdio_dest(&cinfo, outfile);

  cinfo.image_width  = dsize[1];
  cinfo.image_height = dsize[0];
  switch(numch)
    {
    case 1:
      row_stride=dsize[1];
      cinfo.input_components = 1;
      cinfo.in_color_space = JCS_GRAYSCALE;
      break;
    case 3:
      row_stride=3*dsize[1];
      cinfo.input_components = 3;
      cinfo.in_color_space = JCS_RGB;
      break;
    case 4:
      row_stride=4*dsize[1];
      cinfo.input_components = 4;
      cinfo.in_color_space = JCS_CMYK;
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! The number of channels is not 1, 3 "
            "or 4, but %zu. This should not happen. Please contact us so we "
            "can fix the problem", __func__, numch);
    }

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, quality, TRUE);
  cinfo.density_unit=1;
  cinfo.Y_density=cinfo.X_density=dsize[1]/(widthincm/2.54);
  jpeg_start_compress(&cinfo, TRUE);

  /* cinfo.next_scanline is 'unsigned int' */
  c=dsize[0]-1; /* In JPEG the first row is on the bottom!  */
  while (cinfo.next_scanline < cinfo.image_height)
    {
      r[0] = & jsr[c-- * row_stride];
      (void) jpeg_write_scanlines(&cinfo, r, 1);
    }

  jpeg_finish_compress(&cinfo);
  fclose(outfile);
  jpeg_destroy_compress(&cinfo);
}
#endif  /* HAVE_LIBJPEG */





void
gal_jpeg_write(gal_data_t *in, char *filename, uint8_t quality,
               float widthincm)
{
#ifdef HAVE_LIBJPEG
  JSAMPLE *jsr;
  gal_data_t *channel;
  unsigned char *colors[4];
  size_t i, pixel, color;
  size_t numch=gal_list_data_number(in);

  /* Small sanity checks. */
  if(numch==2 || numch>4)
    error(EXIT_FAILURE, 0, "%s: only 1, 3, and 4 color channels are "
          "acceptable, input is a list of %zu data sets", __func__, numch);
  if(in->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "%s: input has a '%s' type, but JPEG images can "
          "only have a 'uint8' type", __func__, gal_type_name(in->type, 1));

  /* Make sure the file doesn't exist and that we have write
     permission. Note that the JPEG standard doesn't have multple
     extensions.*/
  if( gal_checkset_writable_notexist(filename)==0 )
    error(EXIT_FAILURE, 0, "%s: already exists or its directory doesn't "
          "write permssion. Note that the JPEG standard only allows one "
          "image per file", filename);

  /* Make sure the JSAMPLE is 8bits, then allocate the necessary space
     based on the number of channels. */
  if(sizeof *jsr!=1)
    error(EXIT_FAILURE, 0, "%s: JSAMPLE has to be 8bit", __func__);
  errno=0;
  jsr=malloc(numch * in->size * sizeof *jsr);
  if(jsr==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for jsr",
          __func__, numch * in->size * sizeof *jsr );

  /* Set the pointers to each color. */
  i=0;
  for(channel=in; channel!=NULL; channel=channel->next)
    colors[i++]=channel->array;

  /* Write the different colors into jsr. */
  for(pixel=0; pixel<in->size; ++pixel)
    for(color=0;color<numch;++color)
      {
        jsr[pixel*numch+color] = colors[color][pixel];
        /*
        printf("color: %zu, pixel: %zu, jsr: %d\n", color, pixel,
               (int)jsr[pixel*numch+color]);
        */
      }

  /* Write jsr to a JPEG image and clean up. */
  jpeg_write_array(jsr, in, filename, quality, widthincm);
  free(jsr);
#else
  error(EXIT_FAILURE, 0, "%s: libjpeg was not found during the "
        "configuration of %s on this system. To write JPEG files, libjpeg "
        "is required. Please install libjpeg, then configure, make and "
        "install %s again", __func__, PACKAGE_STRING, PACKAGE_STRING);
#endif  /* HAVE_LIBJPEG */
}
