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
#include <string.h>

#include <timing.h>
#include <checkset.h>

#include "main.h"






/*************************************************************
 **************      Acceptable EPS names      ***************
 *************************************************************/
int
nameiseps(char *name)
{
  size_t len;
  len=strlen(name);
  if (strcmp(&name[len-3], "eps") == 0
      || strcmp(&name[len-3], "EPS") == 0
      || strcmp(&name[len-4], "epsf") == 0
      || strcmp(&name[len-4], "epsi") == 0)
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
void
channelsinhex(struct converttparams *p, FILE *fp)
{
  uint8_t *ech;
  size_t i, j, size, numelem=35;
  /*
  p->ech[0][0]=255;
  printf("\n\n%d: %02X\n\n", (int)p->ech[0][0], p->ech[0][0]);
  */
  size=p->s0[0]*p->s1[0];
  for(i=0;i<p->numch;++i)
    {
      if(p->isblank[i])
        fprintf(fp, "{<00>} %% Channel %lu is blank\n", i);
      else
        {
          ech=p->ech[i];
          fprintf(fp, "{<");
          for(j=0;j<size;++j)
            {
              fprintf(fp, "%02X", ech[j]);
              if(j%numelem==0) fprintf(fp, "\n");
            }
          fprintf(fp, ">}\n");
        }
    }
}





void
saveepsorpdf(struct converttparams *p)
{
  FILE *fp;
  float hbw;
  size_t i, winpt, hinpt;
  char command[20000], *epsfilename;



  /* EPS filename */
  if(p->outputtype==EPSFORMAT)
    epsfilename=p->cp.output;
  else if (p->outputtype==PDFFORMAT)
    {
      /* In ui.c we removed the output if it already existed, so it
         doesn't exist now. But automaticoutput is based on an input
         (which must exist), so temporarily make the file. */
      sprintf(command, "touch %s", p->cp.output);
      system(command);
      automaticoutput(p->cp.output, ".eps", 0, p->cp.dontdelete,
                      &epsfilename);
      sprintf(command, "rm %s", p->cp.output);
      system(command);
    }
  else
    error(EXIT_FAILURE, 0, "A bug! In `saveeps`, for outputtype is "
          "neither eps or pdf! Please contact us so we fix it.");



  /* Find the bounding box  */
  winpt=p->widthincm*72.0f/2.54f;
  hinpt=(float)(p->s0[0]*winpt)/(float)(p->s1[0]);
  hbw=(float)p->borderwidth/2.0f;



  /* Open the output file and write the top comments. */
  errno=0;
  fp=fopen(epsfilename, "w");
  if(fp==NULL)
    error(EXIT_FAILURE, errno, p->cp.output);
  fprintf(fp, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fp, "%%%%BoundingBox: 0 0 %lu %lu\n", winpt+2*p->borderwidth,
          hinpt+2*p->borderwidth);
  fprintf(fp, "%%%%Creator: %s\n", SPACK_STRING);
  fprintf(fp, "%%%%CreationDate: %s", ctime(&p->rawtime));
  fprintf(fp, "%%%%LanuageLevel: 3\n");
  fprintf(fp, "%%%%EndComments\n\n");
  if(p->outputtype==EPSFORMAT)
    fprintf(fp, "gsave\n\n");



  /* Commands to draw the border: */
  if(p->borderwidth)
    {
      fprintf(fp, "%% Draw the border:\n");
      fprintf(fp, "0 setgray\n");
      fprintf(fp, "%d setlinewidth\n", p->borderwidth);
      fprintf(fp, "%.1f %.1f moveto\n", hbw, hbw);
      fprintf(fp, "0 %lu rlineto\n", hinpt+p->borderwidth);
      fprintf(fp, "%lu 0 rlineto\n", winpt+p->borderwidth);
      fprintf(fp, "0 -%lu rlineto\n", hinpt+p->borderwidth);
      fprintf(fp, "closepath\n");
      fprintf(fp, "stroke\n\n");
    }



  /* Write the image: */
  fprintf(fp, "%% Draw the image:\n");
  if(p->numch==1)
    fprintf(fp, "/DeviceGray setcolorspace\n");
  else if(p->numch==3)
    fprintf(fp, "/DeviceRGB setcolorspace\n");
  else if(p->numch==4)
    fprintf(fp, "/DeviceCMYK setcolorspace\n");
  else
    error(EXIT_FAILURE, 0, "A bug! In saveepsorpdf the number of channels "
          "is not 1, 3 or 4. Please contact us so we can find the issue "
          "and fix it.");
  fprintf(fp, "%d %d translate\n", p->borderwidth, p->borderwidth);
  fprintf(fp, "%lu %lu scale\n", winpt, hinpt);
  fprintf(fp, "<<\n");
  fprintf(fp, "  /ImageType 1\n");
  fprintf(fp, "  /Width %lu\n", p->s1[0]);
  fprintf(fp, "  /Height %lu\n", p->s0[0]);
  fprintf(fp, "  /ImageMatrix [ %lu 0 0 %lu 0 0 ]\n", p->s1[0], p->s0[0]);
  fprintf(fp, "  /MultipleDataSources true\n");
  fprintf(fp, "  /BitsPerComponent 8\n");
  fprintf(fp, "  /Decode[");
  for(i=0;i<p->numch;++i) fprintf(fp, " 0 1"); fprintf(fp, " ]\n");
  fprintf(fp, "  /Interpolate false\n");
  fprintf(fp, "  /DataSource [\n");
  channelsinhex(p, fp);
  fprintf(fp, "  ]\n");
  fprintf(fp, ">>\n");
  fprintf(fp, "image\n\n");



  /* Ending of the EPS file: */
  if(p->outputtype==EPSFORMAT)
    fprintf(fp, "grestore\n");
  else
    fprintf(fp, "showpage\n");
  fprintf(fp, "%%%%EOF");
  fclose(fp);



  if(p->outputtype==PDFFORMAT)
    {
      sprintf(command, "gs -o %s -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS=%lu "
              "-dDEVICEHEIGHTPOINTS=%lu -dPDFFitPage %s", p->cp.output,
              winpt+2*p->borderwidth, hinpt+2*p->borderwidth, epsfilename);
      system(command);
      sprintf(command, "rm %s", epsfilename);
      system(command);
      free(epsfilename);
    }
}
