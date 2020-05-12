/*********************************************************************
eps -- functions to write EPS files.
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
#include <string.h>

#include <gnuastro/eps.h>
#include <gnuastro/list.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>






/*************************************************************
 **************      Acceptable EPS names      ***************
 *************************************************************/
int
gal_eps_name_is_eps(char *name)
{
  size_t len;
  if(name)
    {
      len=strlen(name);
      if ( ( len>=3 && strcmp(&name[len-3], "eps") == 0 )
           || ( len>=3 && strcmp(&name[len-3], "EPS") == 0 )
           || ( len>=4 && strcmp(&name[len-4], "epsf") == 0 )
           || ( len>=4 && strcmp(&name[len-4], "epsi") == 0 ) )
        return 1;
      else
        return 0;
    }
  else return 0;
}





int
gal_eps_suffix_is_eps(char *name)
{
  if(name)
    {
      if (strcmp(name, "eps") == 0 || strcmp(name, ".eps") == 0
          || strcmp(name, "EPS") == 0 || strcmp(name, ".EPS") == 0
          || strcmp(name, "epsf") == 0 || strcmp(name, ".epsf") == 0
          || strcmp(name, "epsi") == 0 || strcmp(name, ".epsi") == 0)
        return 1;
      else
        return 0;
    }
  else return 0;
}




















/*************************************************************
 **************       Write an EPS image        **************
 *************************************************************/
static int
eps_is_binary(gal_data_t *in, uint8_t *bitone)
{
  gal_data_t *channel;
  uint8_t second_set=0;
  unsigned char *i, *fi, first=0, second=0;

  /* Go through all the channels. */
  for(channel=in; channel!=NULL; channel=channel->next)
    {
      /* Go through all the values and see if there is more than two values
         in the array. */
      fi = (i=channel->array) + channel->size;
      first=*i;
      do
        if(*i!=first)
          {
            if(second_set)
              { if(*i!=second) break; }
            else
              { second=*i; second_set=1; }
          }
      while(++i<fi);

      /* If we didn't get to the end of the array, then we have a
         multi-valued (not binary) image. */
      if(i!=fi)
        return 0;
    }

  /* If we get to this point, then all the channels were binary, so return
     success. */
  *bitone = first>second ? first : second;
  return 1;
}





/* Convert the channels into into a 0 and 1 bit stream. This function is
   only called when the image is binary (has only two values). NOTE: each
   row has to have an integer number of bytes, so when the number of pixels
   in a row is not a multiple of 8, we'll add one. */
static gal_data_t *
eps_convert_to_bitstream(gal_data_t *in, size_t *numbytes, uint8_t bitone)
{
  size_t i, j, k, bytesinrow;
  gal_data_t *channel, *out=NULL;
  unsigned char *bits, byte, curbit, *arr;
  size_t s0=in->dsize[0], s1=in->dsize[1];

  /* Find the size values and allocate the array. */
  if( s1 % 8 ) bytesinrow = s1/8 + 1;
  else         bytesinrow = s1/8;
  *numbytes = bytesinrow*s0;

  /* Go over all the channels. */
  for(channel=in; channel!=NULL; channel=channel->next)
    {
      /* Allocate the array. Note that we currently don't have an
         allocation system for bits, so we'll allocate space in bytes, then
         convert  */
      gal_list_data_add_alloc(&out, NULL, GAL_TYPE_UINT8, 1, numbytes,
                              NULL, 0, -1, 1, NULL, NULL, NULL);
      out->type=GAL_TYPE_BIT;
      bits=out->array;

      /* Put the values in. */
      arr=channel->array;
      for(i=0;i<s0;++i)           /* i*s0+j is the byte, not bit position. */
        for(j=0;j<bytesinrow;++j)
          {
            /* Set the 8 bits to zero. */
            byte=0;

            /* Current bit position. */
            curbit=0x80;

            /* Write the next 8 values as bits. */
            for(k=0;k<8;++k)
              {
                if( j*8+k < s1 )
                  {
                    if(arr[i*s1+j*8+k]==bitone)
                      byte |= curbit;
                    curbit >>= 1;
                  }
                else break;
              }

            /* Write the byte into the array. */
            bits[i*bytesinrow+j]=byte;
          }
    }

  /* Reverse the list and return it. */
  gal_list_data_reverse(&out);
  return out;
}





static void
eps_write_hex(gal_data_t *write, FILE *fp, size_t numbytes)
{
  unsigned char *arr;
  gal_data_t *channel;
  size_t i=0, j, elem_for_newline=35;

  for(channel=write; channel!=NULL; channel=channel->next)
    {
      if(channel->status)       /* A blank channel has status==1. */
        fprintf(fp, "{<00>} %% Channel %zu is blank\n", i);
      else
        {
          arr=channel->array;
          fprintf(fp, "{<");
          for(j=0;j<numbytes;++j)
            {
              fprintf(fp, "%02X", arr[j]);
              if(j%elem_for_newline==0) fprintf(fp, "\n");
            }
          fprintf(fp, ">}\n");
        }
      ++i;
    }
}





static void
eps_write_ascii85(gal_data_t *write, FILE *fp, size_t numbytes)
{
  unsigned char *arr;
  gal_data_t *channel;
  uint32_t anint, base;
  size_t i=0, j, k, elem_for_newline=15;   /* 15*5=75 */

  for(channel=write; channel!=NULL; channel=channel->next)
    {
      if(channel->status)
        fprintf(fp, "{<00>} %% Channel %zu is blank\n", i);
      else
        {
          arr=channel->array;
          fprintf(fp, "{<~");
          for(j=0;j<numbytes;j+=4)
            {
              /* This is the last four bytes */
              if(numbytes-j<4)
                {
                  anint=arr[j]*256*256*256;
                  if(numbytes-j>1)  anint+=arr[j+1]*256*256;
                  if(numbytes-j==3) anint+=arr[j+2]*256;
                }
              else
                anint=( arr[j]*256*256*256 + arr[j+1]*256*256
                        + arr[j+2]*256     + arr[j+3]         );

              /* If all four bytes are zero, then just print 'z'. */
              if(anint==0) fprintf(fp, "z");
              else
                {
                  /* To check, just change the fprintf below to printf:
                     printf("\n\n");
                     printf("%u %u %u %u\n", in[i], in[i+1],
                            in[i+2], in[i+3]);
                  */
                  base=85*85*85*85;
                  /* Do the ASCII85 encoding: */
                  for(k=0;k<5;++k)
                    {
                      fprintf(fp, "%c", anint/base+33);
                      anint%=base;
                      base/=85;
                    }
                }
              /* Go to the next line if on the right place: */
              if(j%elem_for_newline==0) fprintf(fp, "\n");
            }
          fprintf(fp, "~>}\n");
        }
      ++i;
    }
}





static void
eps_write_image(gal_data_t *in, FILE *fp, int hex, int dontoptimize)
{
  int bpc=8;
  uint8_t bitone;
  gal_data_t *write;
  size_t i, numbytes, *dsize=in->dsize;
  size_t numch=gal_list_data_number(in);

  /* Set the number of bits per component. */
  if( numch==1 && dontoptimize==0 && eps_is_binary(in, &bitone) )
    {
      bpc=1;
      write=eps_convert_to_bitstream(in, &numbytes, bitone);
    }
  else
    {
      write=in;
      numbytes=in->size;
    }

  /* Write the basic meta data. */
  switch(numch)
    {
    case 1: fprintf(fp, "/DeviceGray setcolorspace\n"); break;
    case 3: fprintf(fp, "/DeviceRGB setcolorspace\n");  break;
    case 4: fprintf(fp, "/DeviceCMYK setcolorspace\n"); break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! The number of channels (%zu) is "
            "not 1, 3 or 4. Please contact us so we can find the issue and "
            "fix it", __func__, numch);
    }
  fprintf(fp, "<<\n");
  fprintf(fp, "  /ImageType 1\n");
  fprintf(fp, "  /Width %zu\n", dsize[1]);
  fprintf(fp, "  /Height %zu\n", dsize[0]);
  fprintf(fp, "  /ImageMatrix [ %zu 0 0 %zu 0 0 ]\n", dsize[1], dsize[0]);
  fprintf(fp, "  /MultipleDataSources true\n");
  fprintf(fp, "  /BitsPerComponent %d\n", bpc);
  fprintf(fp, "  /Decode[");
  for(i=0;i<numch;++i) {fprintf(fp, " 0 1");} fprintf(fp, " ]\n");
  fprintf(fp, "  /Interpolate false\n");
  fprintf(fp, "  /DataSource [\n");

  /* Based on the encoding, write the contents of the image. */
  if(hex)
    eps_write_hex(write, fp, numbytes);
  else
    eps_write_ascii85(write, fp, numbytes);

  /* Finish the file. */
  fprintf(fp, "  ]\n");
  fprintf(fp, ">>\n");
  fprintf(fp, "image\n\n");

  /* Clean up. */
  if(write!=in)
    gal_list_data_free(write);
}





void
gal_eps_to_pt(float widthincm, size_t *dsize, size_t *w_h_in_pt)
{
  w_h_in_pt[0] = widthincm*72.0f/2.54f;
  w_h_in_pt[1] = (float)( dsize[0] * w_h_in_pt[0] )/(float)(dsize[1]);
}





void
gal_eps_write(gal_data_t *in, char *filename, float widthincm,
              uint32_t borderwidth, int hex, int dontoptimize, int forpdf)
{
  FILE *fp;
  float hbw;
  time_t rawtime;
  size_t numch=gal_list_data_number(in);
  size_t w_h_in_pt[2], *dsize=in->dsize;

  /* Sanity checks. */
  if(numch==2 || numch>4)
    error(EXIT_FAILURE, 0, "%s: only 1, 3, and 4 color channels are "
          "acceptable, input is a list of %zu data sets", __func__, numch);
  if(in->type!=GAL_TYPE_UINT8)
    error(EXIT_FAILURE, 0, "%s: input has a '%s' type, but JPEG images can "
          "only have a 'uint8' type", __func__, gal_type_name(in->type, 1));


  /* Read the time to write in the output. */
  time(&rawtime);


  /* Find the bounding box  */
  hbw=(float)borderwidth/2.0;
  gal_eps_to_pt(widthincm, dsize, w_h_in_pt);


  /* Open the output file and write the top comments. */
  errno=0;
  fp=fopen(filename, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s", filename);
  fprintf(fp, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fp, "%%%%BoundingBox: 0 0 %zu %zu\n", w_h_in_pt[0]+2*borderwidth,
          w_h_in_pt[1]+2*borderwidth);
  fprintf(fp, "%%%%Creator: %s\n", PACKAGE_STRING);
  fprintf(fp, "%%%%CreationDate: %s", ctime(&rawtime));
  fprintf(fp, "%%%%LanuageLevel: 3\n");
  fprintf(fp, "%%%%EndComments\n\n");
  if(forpdf==0)
    fprintf(fp, "gsave\n\n");


  /* Commands to draw the border: */
  if(borderwidth)
    {
      fprintf(fp, "%% Draw the border:\n");
      fprintf(fp, "0 setgray\n");
      fprintf(fp, "%d setlinewidth\n", borderwidth);
      fprintf(fp, "%.1f %.1f moveto\n", hbw, hbw);
      fprintf(fp, "0 %zu rlineto\n", w_h_in_pt[1]+borderwidth);
      fprintf(fp, "%zu 0 rlineto\n", w_h_in_pt[0]+borderwidth);
      fprintf(fp, "0 -%zu rlineto\n", w_h_in_pt[1]+borderwidth);
      fprintf(fp, "closepath\n");
      fprintf(fp, "stroke\n\n");
    }


  /* Write the image: */
  fprintf(fp, "%% Draw the image:\n");
  fprintf(fp, "%d %d translate\n", borderwidth, borderwidth);
  fprintf(fp, "%zu %zu scale\n", w_h_in_pt[0], w_h_in_pt[1]);
  eps_write_image(in, fp, hex, dontoptimize);


  /* Ending of the EPS file: */
  if(forpdf) fprintf(fp, "showpage\n");
  else       fprintf(fp, "grestore\n");
  fprintf(fp, "%%%%EOF");
  fclose(fp);
}
