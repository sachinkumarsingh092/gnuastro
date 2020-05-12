/*********************************************************************
tiff -- functions to read and write TIFF files.
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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_LIBTIFF
  #include <tiffio.h>
#endif

#include <gnuastro/fits.h>
#include <gnuastro/data.h>
#include <gnuastro/list.h>
#include <gnuastro/tiff.h>
#include <gnuastro/pointer.h>









/*************************************************************
 **************      Acceptable TIFF names      **************
 *************************************************************/
int
gal_tiff_name_is_tiff(char *name)
{
  size_t len;

  if(name)
    {
      len=strlen(name);
      if ( ( len>=3 && strcmp(&name[len-3], "tif") == 0 )
           || ( len>=3 && strcmp(&name[len-3], "TIF") == 0 )
           || ( len>=4 && strcmp(&name[len-4], "tiff") == 0 )
           || ( len>=4 && strcmp(&name[len-4], "TIFF") == 0 ) )
        return 1;
      else
        return 0;
    }
  else return 0;
}





int
gal_tiff_suffix_is_tiff(char *name)
{
  if(name)
    {
    if (strcmp(name, "tif") == 0   || strcmp(name, ".tif") == 0
        || strcmp(name, "TIF") == 0 || strcmp(name, ".TIF") == 0
        || strcmp(name, "tiff") == 0 || strcmp(name, ".tiff") == 0
        || strcmp(name, "TIFF") == 0 || strcmp(name, ".TIFF") == 0 )
      return 1;
    else
      return 0;
    }
  else return 0;
}





/* Users may define the TIFF directory to read as a string, in that case,
   this function can be used to convert it to a 'size_t' for use in
   'gal_tiff_read'.  */
size_t
gal_tiff_dir_string_read(char *string)
{
  long dir;
  char *tailptr=NULL;

  /* Read the given directory string. */
  errno=0;
  dir=strtol(string, &tailptr, 0);
  if(tailptr==string)
    error(EXIT_FAILURE, 0, "%s: '%s' couldn't be read as an integer",
          __func__, string);
  if(errno)
    error(EXIT_FAILURE, errno, "%s: reading %s", __func__, string);
  if(dir<0)
    error(EXIT_FAILURE, 0, "%s: %ld is a negative integer, it must be "
          "positive", __func__, dir);

  /* Return the result. */
  return dir;
}




















/*************************************************************
 **************        Read a JPEG image        **************
 *************************************************************/
#ifdef HAVE_LIBTIFF
static void
tiff_read_tag(TIFF *tif, ttag_t tag, void *out, char *filename, size_t dir)
{
  /* Read the tag */
  if( !TIFFGetField(tif, tag, out) )
    error(EXIT_FAILURE, 0, "%s: %s (dir %zu): tag %d couldn't be fetched",
          __func__, filename, dir, tag);
}





/* Convert the TIFF type code into Gnuastro's type code.*/
static uint8_t
tiff_type_read(TIFF *tif, char *filename, size_t dir)
{
  uint16_t bitspersample, sampleformat;

  /* Read the number of bits and the corresponding type. */
  tiff_read_tag(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample, filename, dir);

  /* Read the formatting of each pixel. If no such keyword exists, use the
     value of 'SAMPLEFORMAT_UINT'. */
  if( TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sampleformat) != 1 )
    sampleformat=SAMPLEFORMAT_UINT;

  /* Read the datatype. */
  switch(sampleformat)
    {
    /* Unsigned integer types. */
    case SAMPLEFORMAT_UINT:
      switch(bitspersample)
        {
        case 8:  return GAL_TYPE_UINT8;
        case 16: return GAL_TYPE_UINT16;
        case 32: return GAL_TYPE_UINT32;
        case 64: return GAL_TYPE_UINT64;
        default:
          error(EXIT_FAILURE, 0, "%s: %s (dir %zu): %u-bit samples not "
                "recognized for UNSIGNED-int format", __func__, filename,
                dir, bitspersample);
        }
      break;

    /* Signed integer types. */
    case SAMPLEFORMAT_INT:
      switch(bitspersample)
        {
        case 8:  return GAL_TYPE_INT8;
        case 16: return GAL_TYPE_INT16;
        case 32: return GAL_TYPE_INT32;
        case 64: return GAL_TYPE_INT64;
        default:
          error(EXIT_FAILURE, 0, "%s: %s (dir %zu): %u-bit samples not "
                "recognized for SIGNED-int format", __func__, filename,
                dir, bitspersample);
        }
      break;


    /* Floating point types. */
    case SAMPLEFORMAT_IEEEFP:
      switch(bitspersample)
        {
        case 32: return GAL_TYPE_FLOAT32;
        case 64: return GAL_TYPE_FLOAT64;
        default:
          error(EXIT_FAILURE, 0, "%s: %s (dir %zu): %u-bit samples not "
                "recognized for floating point format", __func__, filename,
                dir, bitspersample);
        }
      break;


    /* The reported value is not recognized. */
    default:
      error(EXIT_FAILURE, 0, "%s: %s (dir %zu): value %u not recognized "
            "for SAMPLEFORMAT tag", __func__, filename, dir, sampleformat);
    }

  /* Control should never reach here, but to avoid compiler warnings (and
     hard to find bugs when it does reach here), we'll just return an
     invalid type. */
  return GAL_TYPE_INVALID;
}





/* Get the basic TIFF image information. */
static void
tiff_img_info(TIFF *tif, uint8_t *type, size_t *ndim, size_t *dsize,
              size_t *numch, char *filename, size_t dir)
{
  size_t d=0;
  uint16_t u16;
  uint32_t u32;

  /* Based on if 'IMAGEDEPTH' is defined in the TIFF header, set the
     dimensions. */
  if( TIFFGetField(tif, TIFFTAG_IMAGEDEPTH, &u32) )
    dsize[ d++ ]=u32;

  /* Read the other sizes. Note that in the TIFF standard, IMAGELENGTH is
     the vertical length of the image and IMAGEWIDTH is the horizontal
     length. */
  tiff_read_tag(tif, TIFFTAG_IMAGELENGTH, &u32, filename, dir);
  dsize[ d++ ]=u32;
  tiff_read_tag(tif, TIFFTAG_IMAGEWIDTH, &u32, filename, dir);
  dsize[ d++ ]=u32;

  /* Write the dimensions. */
  *ndim=d;

  /* Read the type of the image. */
  *type=tiff_type_read(tif, filename, dir);

  /* Read the number of channels in the image. */
  *numch = TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &u16) ? u16 : 1;
}





/* Based on the 'TIFFReadContigStripData' function of 'tools/tiffinfo.c' of
   Libtiff's source. */
void
tiff_read_contig_strip_data(TIFF *tif, char *filename, size_t dir,
                            gal_data_t *out, size_t numch)
{
  tstrip_t strip;
  size_t ostart=0;
  unsigned char *buf;
  uint32_t row, rowsperstrip = (uint32_t)-1;
  size_t nrow=0, scanline=TIFFScanlineSize(tif);
  uint32 h=out->ndim==2?out->dsize[0]:out->dsize[1];

  /* Allocate the buffer. */
  errno=0;
  buf = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
  if(buf==NULL)
    error(EXIT_FAILURE, errno, "%s: %s (dir %zu): couldn't allocate "
          "necessary space to load image (%zu bytes)", __func__, filename,
          dir, scanline);

  /* Parse over the rows and read everything. */
  TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
  for(row=0; row<h; row+=rowsperstrip)
    {
      /* Read the part of the image into the buffer. */
      strip = TIFFComputeStrip(tif, row, 0);
      nrow  = (row+rowsperstrip > h ? h-row : rowsperstrip);
      if( TIFFReadEncodedStrip(tif, strip, buf, nrow*scanline) < 0 )
        error(EXIT_FAILURE, 0, "%s: %s (dir %zu): couldn't read data",
              __func__, filename, dir);

      /* Copy the contents of the buffer to the output array. Note that
         'ostart' is the byte count already, so the type is
         irrelevant. Thus, we can read 'out->array' as a 'char *'
         pointer.*/
      memcpy( (char *)(out->array)+ostart, buf, nrow*scanline);
      ostart+=nrow*scanline;
    }

  /* Clean up. */
  _TIFFfree(buf);
}





/* Based on the 'TIFFReadSeparateStripData' function of 'tools/tiffinfo.c'
   of Libtiff's source. */
static void
tiff_read_separate_strip_data(TIFF* tif, char *filename, size_t dir,
                              gal_data_t *out)
{
  tsample_t s;
  gal_data_t *ch;
  tstrip_t strip;
  unsigned char *buf;
  uint32_t rowsperstrip=(uint32_t)-1;
  size_t nrow=0, scanline=TIFFScanlineSize(tif);
  size_t ostart=0, numch=gal_list_data_number(out);
  uint32 row, h=out->ndim==2?out->dsize[0]:out->dsize[1];

  /* Allocate the buffer. */
  errno=0;
  buf = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
  if(buf==NULL)
    error(EXIT_FAILURE, errno, "%s: %s (dir %zu): couldn't allocate "
          "necessary space to load image (%zu bytes)", __func__, filename,
          dir, scanline);

  /* Parse over the dataset and read them into the output. */
  TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
  for(row=0; row<h; row+=rowsperstrip)
    {
      /* Restart the channel list. */
      ch=out;

      /* Go over each channel. */
      for(s=0; s<numch; s++)
        {
          /* Read this part of the data into the buffer. */
          strip = TIFFComputeStrip(tif, row, s);
          nrow = (row+rowsperstrip > h ? h-row : rowsperstrip);
          if( TIFFReadEncodedStrip(tif, strip, buf, nrow*scanline) < 0 )
            error(EXIT_FAILURE, 0, "%s: %s (dir %zu): couldn't read data",
                  __func__, filename, dir);

          /* Write the buffer into the output array. */
          memcpy( (char *)(ch->array)+ostart, buf, nrow*scanline);

          /* Go onto the next channel. */
          ch=ch->next;
        }

      /* Increment the starting pointer. */
      ostart+=nrow*scanline;
    }

  /* Clean up. */
  _TIFFfree(buf);
}





/* The data have been read contiguously (the pixels for each color are
   beside each other). We need to separate the color channels into
   different datasets. We will also use this chance to reverse the order of
   the rows */
static gal_data_t *
tiff_separate_channels_reverse(gal_data_t *out, size_t numch,
                               size_t minmapsize, int quietmmap)
{
  size_t k, l, i, j;
  gal_data_t *tch, *ch=NULL;
  size_t so=gal_type_sizeof(out->type);
  size_t width=out->dsize[1]*so/numch, lwidth=out->dsize[1]*so;

  /* A small sanity check. */
  if(out->ndim==3)
    error(EXIT_FAILURE, 0, "%s: currently only 2D datasets are supported, "
          "please get in touch with us at %s to add 3D support", __func__,
          PACKAGE_BUGREPORT);


  /* Make the separated datasets (temporarily fix the extra width). */
  out->dsize[1] /= numch;
  for(k=0; k<numch; ++k)
    gal_list_data_add_alloc(&ch, NULL, out->type, out->ndim, out->dsize,
                            NULL, 0, minmapsize, quietmmap, NULL, NULL, NULL);
  out->dsize[1] *= numch;


  /* Parse over the rows and write them in the output. */
  for(i=0;i<out->dsize[0];++i)
    {
      /* 'j' is the output row. */
      j=out->dsize[0]-1-i;

      /* 'k' is the element in each row and 'l' is the color/channel. */
      for(k=0;k<ch->dsize[1];++k)
        {
          /* Initialize the color/channel counter ocpy the elements. */
          l=0;
          for(tch=ch; tch!=NULL; tch=tch->next)
            {
              /* Elements of the 'j'th row into the 'i'th. */
              memcpy( (char *)(tch->array) + i*width  + k*so,
                      (char *)(out->array) + j*lwidth + k*numch*so + so*l,
                      so);

              /* Increment the color. */
              ++l;
            }
        }
    }

  /* Clean up and return. */
  return ch;
}





/* The standard TIFF format is up-side-down when viewed in FITS (which is
   the base of Gnuastro also). So we need to reverse the array for an
   identical orientation. */
static void
tiff_reverse_rows(gal_data_t *out)
{
  gal_data_t *ch=out;
  size_t c, i, j, numch=gal_list_data_number(out);
  size_t width=out->dsize[1]*gal_type_sizeof(out->type);
  void *tmp=gal_pointer_allocate(out->type, out->dsize[1], 0, __func__,
                                 "tmp");

  /* A small sanity check. */
  if(out->ndim==3)
    error(EXIT_FAILURE, 0, "%s: currently only 2D datasets are supported, "
          "please get in touch with us at %s to add 3D support", __func__,
          PACKAGE_BUGREPORT);

  /* Parse over the rows and reverse them. */
  for(c=0;c<numch;++c)
    {
      /* Initialize the parsing counters. */
      i=0;
      j=out->dsize[0]-1;

      /* Go over the channel. */
      while(j>i)
        {
          /* Copy the 'i'th row into a temporary array. */
          memcpy(tmp, (char *)(ch->array)+i*width, width);

          /* Put the 'j'th row into the 'i'th row. */
          memcpy( (char *)(ch->array)+i*width, (char *)(ch->array)+j*width,
                  width );

          /* Put the 'tmp' row into 'j'. */
          memcpy( (char *)(ch->array)+j*width, tmp, width);

          /* Increment the points. */
          ++i;
          --j;
        }

      /* Go to the next channel. */
      ch=ch->next;
    }

  /* Clean up. */
  free(tmp);
}





/* Read the data following the 'TIFFReadData' of Libtiff's
   'tools/tiffinfo.c' in the libtiff source code. */
static gal_data_t *
tiff_img_read(TIFF *tif, char *filename, size_t dir, size_t minmapsize,
              int quietmmap)
{
  uint8_t type;
  uint16_t config;
  gal_data_t *sep, *out=NULL;
  size_t i, ndim, numch, dsize[3];


  /* Get the basic image information. */
  tiff_img_info(tif, &type, &ndim, dsize, &numch, filename, dir);


  /* Find the planar state of the input: are the channels separate or
     contiguous? Based on that, allocate the output data structure. */
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
  if( config==PLANARCONFIG_CONTIG )
    {
      /* We'll multiply the last dimension's length by the number of
         channels to keep the contiguous color information in one array. */
      dsize[ndim-1] *= numch;
      out=gal_data_alloc(NULL, type, ndim, dsize, NULL, 0, minmapsize,
                         quietmmap, NULL, NULL, NULL);
    }
  else
    for(i=0; i<numch; ++i)
      gal_list_data_add_alloc(&out, NULL, type, ndim, dsize, NULL, 0,
                              minmapsize, quietmmap, NULL, NULL, NULL);


  /* The reading of the dataset depends on how it is organized, so first
     we'll look into the 'planarconfig' field. */
  if( TIFFIsTiled(tif) )
    {
      if(config==PLANARCONFIG_CONTIG)
        error(EXIT_FAILURE, 0, "%s: %s (dir %zu) is a contiguous tiled TIFF "
              "file which is not yet supported in Gnuastro, please get in "
              "touch with us at %s to add this feature", __func__, filename,
              dir, PACKAGE_BUGREPORT);
      else
        error(EXIT_FAILURE, 0, "%s: %s (dir %zu) is a non-contiguous tiled "
              "TIFF file which is not yet supported in Gnuastro, please "
              "get in touch with us at %s to add this feature", __func__,
              filename, dir, PACKAGE_BUGREPORT);
    }
  else
    {
      if(config==PLANARCONFIG_CONTIG)
        tiff_read_contig_strip_data(tif, filename, dir, out, numch);
      else
        tiff_read_separate_strip_data(tif, filename, dir, out);
    }


  /* When there are more than one channels and the colors are stored
     contiguously, we need to break up the array into multiple arrays. When
     any of these conditions don't hold, the channels are already
     separated, we just need to reverse them.*/
  if( numch>1 && config==PLANARCONFIG_CONTIG )
    {
      sep=tiff_separate_channels_reverse(out, numch, minmapsize, quietmmap);
      gal_data_free(out);
      out=sep;
    }
  else
    tiff_reverse_rows(out);


  /* Return the output. */
  return out;
}
#endif





gal_data_t *
gal_tiff_read(char *filename, size_t dir, size_t minmapsize, int quietmmap)
{
#ifdef HAVE_LIBTIFF
  TIFF *tif;
  gal_data_t *out;
  size_t dircount=0;

  /* Open the TIFF file. */
  tif=TIFFOpen(filename, "r");
  if(tif==NULL)
    error(EXIT_FAILURE, 0, "%s: '%s' couldn't be opened for reading",
          __func__, filename);

  /* If anything other than the first directory (value of zero) is
     requested, then change the directories. */
  if(dir)
    {
      if( TIFFSetDirectory(tif, dir)==0 )
        {
          /* For an informative error message, count how many directories
             are in the TIFF file. */
          do ++dircount; while( TIFFReadDirectory(tif) );

          /* Close the TIFF file and return. */
          TIFFClose(tif);
          error(EXIT_FAILURE, 0, "%s: '%s' has %zu director%s/extension%s, "
                "and directories are counted from 0. You have asked for "
                "directory %zu", __func__, filename, dircount,
                dircount==1?"y":"ies", dircount==1?"":"s", dir);
        }
    }

  /* Read the image. */
  out=tiff_img_read(tif, filename, dir, minmapsize, quietmmap);

  /* Close file, clean up and return. */
  TIFFClose(tif);
  return out;
#else
  error(EXIT_FAILURE, 0, "%s: libtiff was not found during the "
        "configuration of %s on this system. To read from TIFF files, "
        "libtiff is required. Please install libtiff, then configure, make "
        "and install %s again", __func__, PACKAGE_STRING, PACKAGE_STRING);
  return NULL;
#endif  /* HAVE_LIBTIFF */
}
