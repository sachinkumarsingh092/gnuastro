/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2018, Free Software Foundation, Inc.

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

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"






/*************************************************************
 **************      Acceptable EPS names      ***************
 *************************************************************/
int
nameiseps(char *name)
{
  size_t len;
  len=strlen(name);
  if ( ( len>=3 && strcmp(&name[len-3], "eps") == 0 )
       || ( len>=3 && strcmp(&name[len-3], "EPS") == 0 )
       || ( len>=4 && strcmp(&name[len-4], "epsf") == 0 )
       || ( len>=4 && strcmp(&name[len-4], "epsi") == 0 ) )
    return 1;
  else
    return 0;
}





int
nameisepssuffix(char *name)
{
  if (strcmp(name, "eps") == 0 || strcmp(name, ".eps") == 0
      || strcmp(name, "EPS") == 0 || strcmp(name, ".EPS") == 0
      || strcmp(name, "epsf") == 0 || strcmp(name, ".epsf") == 0
      || strcmp(name, "epsi") == 0 || strcmp(name, ".epsi") == 0)
    return 1;
  else
    return 0;
}





int
nameispdf(char *name)
{
  size_t len;
  len=strlen(name);
  if (strcmp(&name[len-3], "pdf") == 0
      || strcmp(&name[len-3], "PDF") == 0)
    return 1;
  else
    return 0;
}





int
nameispdfsuffix(char *name)
{
  if (strcmp(name, "pdf") == 0 || strcmp(name, ".pdf") == 0
      || strcmp(name, "PDF") == 0 || strcmp(name, ".PDF") == 0)
    return 1;
  else
    return 0;
}



















/*************************************************************
 **************       Write an EPS image        **************
 *************************************************************/
static int
eps_is_binary(struct converttparams *p)
{
  gal_data_t *channel;
  unsigned char *i, *fi;

  /* Go through all the channels. */
  for(channel=p->chll; channel!=NULL; channel=channel->next)
    {
      /* Go through all the values and see they are 0 and break out of the
         loop as soon as you get to a pixel that is not 0 or `maxbyte'. */
      fi = (i=p->chll->array) + p->chll->size;
      do
        if(*i!=p->maxbyte && *i!=0) break;
      while(++i<fi);

      /* If we didn't get to the end of the channel, then we have a
         non-binary image. */
      if(i!=fi)
        return 0;
    }

  /* If we get to this point, then all the channels were binary, so return
     success. */
  return 1;
}





/* Show the bit values in a uint8_t variable. It is included here as a
   test in debugging problems with blackandwhite. To test it use a
   very small valued input to make the outputs reasonable. For example
   I am now testing it with an input text array of 12 elements (while
   calling the --invert option):

   1 0 1 0 0 0
   0 0 0 1 0 1
*/
void
eps_show_bits(uint8_t x)
{
  int i;

  for(i=7;i>=0;--i)
    (x&(1<<i)) ? putchar('1') : putchar('0');
  putchar('\n');
}






/* Convert the channels into into a 0 and 1 bit stream. This function is
   only called when the image is binary (has only two values). NOTE: each
   row has to have an integer number of bytes, so when the number of pixels
   in a row is not a multiple of 8, we'll add one. */
size_t
eps_convert_to_bitstream(struct converttparams *p)
{
  gal_data_t *channel;
  size_t i, j, k, bytesinrow, bytesinimg;
  unsigned char *bits, byte, curbit, *in;
  size_t s0=p->chll->dsize[0], s1=p->chll->dsize[1];

  /* Find the size values and allocate the array. */
  if( s1 % 8 ) bytesinrow = s1/8 + 1;
  else         bytesinrow = s1/8;
  bytesinimg = bytesinrow*s0;

  /* Go over all the channels. */
  for(channel=p->chll; channel!=NULL; channel=channel->next)
    {
      /* Allocate the array. */
      bits=gal_data_malloc_array(GAL_TYPE_UINT8, bytesinimg, __func__,
                                 "bits");

      /* Put the values in. */
      in=channel->array;
      for(i=0;i<s0;++i)
        {
          for(j=0;j<bytesinrow;++j)
            {                  /* i*s0+j is the byte, not bit position. */
              byte=0;          /* Set the 8 bits to zero.               */
              curbit=0x80;     /* Current bit position, starting at:    */
              for(k=0;k<8;++k)
                {
                  if( j*8+k < s1 )
                    {
                      if(in[i*s1+j*8+k])
                        byte |= curbit;
                      curbit >>= 1;
                    }
                  else break;
                }
              /*eps_show_bits(byte);*/
              bits[i*bytesinrow+j]=byte;
            }
        }
      free(channel->array);
      channel->array=bits;
      channel->type=GAL_TYPE_BIT;
    }

  /* Return the total number of bytes in the image. */
  return bytesinimg;
}





void
eps_write_hex(struct converttparams *p, FILE *fp, size_t size)
{
  unsigned char *in;
  gal_data_t *channel;
  size_t i=0, j, elem_for_newline=35;

  for(channel=p->chll; channel!=NULL; channel=channel->next)
    {
      if(channel->status)       /* A blank channel has status==1. */
        fprintf(fp, "{<00>} %% Channel %zu is blank\n", i);
      else
        {
          in=channel->array;
          fprintf(fp, "{<");
          for(j=0;j<size;++j)
            {
              fprintf(fp, "%02X", in[j]);
              if(j%elem_for_newline==0) fprintf(fp, "\n");
            }
          fprintf(fp, ">}\n");
        }
      ++i;
    }
}





void
eps_write_ascii85(struct converttparams *p, FILE *fp, size_t size)
{
  unsigned char *in;
  gal_data_t *channel;
  uint32_t anint, base;
  size_t i=0, j, k, elem_for_newline=15;   /* 15*5=75 */

  for(channel=p->chll; channel!=NULL; channel=channel->next)
    {
      if(channel->status)
        fprintf(fp, "{<00>} %% Channel %zu is blank\n", i);
      else
        {
          in=channel->array;
          fprintf(fp, "{<~");
          for(j=0;j<size;j+=4)
            {
              /* This is the last four bytes */
              if(size-j<4)
                {
                  anint=in[j]*256*256*256;
                  if(size-j>1)  anint+=in[j+1]*256*256;
                  if(size-j==3) anint+=in[j+2]*256;
                }
              else
                anint=( in[j]*256*256*256 + in[j+1]*256*256
                        + in[j+2]*256     + in[j+3]         );

              /* If all four bytes are zero, then just print `z'. */
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
eps_write_image(struct converttparams *p, FILE *fp)
{
  int bpc=8;
  size_t i, size, *dsize=p->chll->dsize;

  /* Set the number of bits per component. */
  if( p->numch==1 && eps_is_binary(p) )
    {
      bpc=1;
      size=eps_convert_to_bitstream(p);
    }
  else size=p->chll->size;

  switch(p->numch)
    {
    case 1: fprintf(fp, "/DeviceGray setcolorspace\n"); break;
    case 3: fprintf(fp, "/DeviceRGB setcolorspace\n");  break;
    case 4: fprintf(fp, "/DeviceCMYK setcolorspace\n"); break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! The number of channels (%zu) is not "
            "1, 3 or 4. Please contact us so we can find the issue and fix it",
            __func__, p->numch);
    }
  fprintf(fp, "<<\n");
  fprintf(fp, "  /ImageType 1\n");
  fprintf(fp, "  /Width %zu\n", dsize[1]);
  fprintf(fp, "  /Height %zu\n", dsize[0]);
  fprintf(fp, "  /ImageMatrix [ %zu 0 0 %zu 0 0 ]\n", dsize[1], dsize[0]);
  fprintf(fp, "  /MultipleDataSources true\n");
  fprintf(fp, "  /BitsPerComponent %d\n", bpc);
  fprintf(fp, "  /Decode[");
  for(i=0;i<p->numch;++i) {fprintf(fp, " 0 1");} fprintf(fp, " ]\n");
  fprintf(fp, "  /Interpolate false\n");
  fprintf(fp, "  /DataSource [\n");
  if(p->hex) eps_write_hex(p, fp, size);
  else eps_write_ascii85(p, fp, size);
  fprintf(fp, "  ]\n");
  fprintf(fp, ">>\n");
  fprintf(fp, "image\n\n");
}




void
eps_write_eps_or_pdf(struct converttparams *p)
{
  FILE *fp;
  float hbw;
  char command[20000], *epsfilename=NULL;
  size_t winpt, hinpt, *dsize=p->chll->dsize;


  /* Find the bounding box  */
  winpt=p->widthincm*72.0f/2.54f;
  hinpt=(float)( dsize[0] * winpt )/(float)(dsize[1]);
  hbw=(float)p->borderwidth/2.0f;


  /* EPS filename */
  if(p->outformat==OUT_FORMAT_EPS)
    {
      epsfilename=p->cp.output;
      gal_checkset_writable_remove(epsfilename, 0, p->cp.dontdelete);
    }
  else if (p->outformat==OUT_FORMAT_PDF)
    {
      gal_checkset_writable_remove(p->cp.output, 0, p->cp.dontdelete);
      epsfilename=gal_checkset_automatic_output(&p->cp, p->cp.output, ".ps");
    }
  else
    error(EXIT_FAILURE, 0, "%s: a bug! code %d not recognized for "
          "`p->outformat'", __func__, p->outformat);


  /* Open the output file and write the top comments. */
  errno=0;
  fp=fopen(epsfilename, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s", p->cp.output);
  fprintf(fp, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fp, "%%%%BoundingBox: 0 0 %zu %zu\n", winpt+2*p->borderwidth,
          hinpt+2*p->borderwidth);
  fprintf(fp, "%%%%Creator: %s\n", PROGRAM_STRING);
  fprintf(fp, "%%%%CreationDate: %s", ctime(&p->rawtime));
  fprintf(fp, "%%%%LanuageLevel: 3\n");
  fprintf(fp, "%%%%EndComments\n\n");
  if(p->outformat==OUT_FORMAT_EPS)
    fprintf(fp, "gsave\n\n");


  /* Commands to draw the border: */
  if(p->borderwidth)
    {
      fprintf(fp, "%% Draw the border:\n");
      fprintf(fp, "0 setgray\n");
      fprintf(fp, "%d setlinewidth\n", p->borderwidth);
      fprintf(fp, "%.1f %.1f moveto\n", hbw, hbw);
      fprintf(fp, "0 %zu rlineto\n", hinpt+p->borderwidth);
      fprintf(fp, "%zu 0 rlineto\n", winpt+p->borderwidth);
      fprintf(fp, "0 -%zu rlineto\n", hinpt+p->borderwidth);
      fprintf(fp, "closepath\n");
      fprintf(fp, "stroke\n\n");
    }



  /* Write the image: */
  fprintf(fp, "%% Draw the image:\n");
  fprintf(fp, "%d %d translate\n", p->borderwidth, p->borderwidth);
  fprintf(fp, "%zu %zu scale\n", winpt, hinpt);
  eps_write_image(p, fp);



  /* Ending of the EPS file: */
  if(p->outformat==OUT_FORMAT_EPS)
    fprintf(fp, "grestore\n");
  else
    fprintf(fp, "showpage\n");
  fprintf(fp, "%%%%EOF");
  fclose(fp);



  if(p->outformat==OUT_FORMAT_PDF)
    {
      sprintf(command, "gs -q -o %s -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS=%zu"
              " -dDEVICEHEIGHTPOINTS=%zu -dPDFFitPage %s", p->cp.output,
              winpt+2*p->borderwidth, hinpt+2*p->borderwidth, epsfilename);
      if(system(command))
        error(EXIT_FAILURE, 0, "the command to convert a PostScript file to "
              "PDF (`%s') was not successful! The PostScript file (%s) is "
              "left if you want to convert or use it through any other "
              "means", command, epsfilename);
      sprintf(command, "rm %s", epsfilename);
      if(system(command))
        error(EXIT_FAILURE, 0, "The PDF output (%s) was created, but the "
              "PostScript file which was used to make it (%s) could not be"
              "removed", p->cp.output, epsfilename);
      free(epsfilename);
    }
}
