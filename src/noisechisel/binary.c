/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>










/****************************************************************
 *****************      Dilate and Erode     ********************
 ****************************************************************/
/* 4 connected dilation and erosion: b0_f1==0: Dilate the
   foreground. b0_f1==1: Erode the foreground.  */
void
dilate0_erode1_4con(unsigned char *byt, size_t nr, size_t nc,
                    unsigned char b0_f1)
{
  size_t i, j, ind;
  unsigned char f, b, *pt, *fpt;

  /* Do a sanity check: */
  if(b0_f1!=1 && b0_f1!=0)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix "
          "this problem. In dilate0_erode1_4con (binary.c), the value to "
          "b0_f1 is %u while it should be 0 or 1.", PACKAGE_BUGREPORT, b0_f1);

  /* Set the foreground and background values. */
  if(b0_f1==0) {f=1; b=0;}
  else         {f=0; b=1;}

  /* Check the 4 corners: */
  if(byt[0]==b && (byt[1]==f || byt[nc]==f) )
    byt[0]=2;

  if(byt[nc-1]==b && (byt[nc-2]==f || byt[2*nc-1]==f))
    byt[nc-1]=2;

  if(byt[(nr-1)*nc]==b
     && (byt[(nr-2)*nc]==f || byt[(nr-1)*nc+1]==f) )
    byt[(nr-1)*nc]=2;

  if(byt[nr*nc-1]==b
     && (byt[nr*nc-2]==f || byt[nr*nc-1-nc]==f) )
    byt[nr*nc-1]=2;

  /* Check the 4 sides: */
  for(j=1;j<nc-1;++j)
    if(byt[j]==b
       && (byt[j+1]==f || byt[j-1]==f || byt[j+nc]==f) )
      byt[j]=2;

  for(j=1;j<nc-1;++j)
    {
      ind=(nr-1)*nc+j;
      if(byt[ind]==b
	 && (byt[ind+1]==f || byt[ind-1]==f || byt[ind-nc]==f) )
	byt[ind]=2;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=i*nc;
      if(byt[ind]==b
	 && (byt[ind+1]==f || byt[ind+nc]==f || byt[ind-nc]==f) )
	byt[ind]=2;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=(i+1)*nc-1;
      if(byt[ind]==b
	 && (byt[ind-1]==f || byt[ind+nc]==f || byt[ind-nc]==f) )
	byt[ind]=2;
    }

  /* Check the body: */
  for(i=1;i<nr-1;++i)
    for(j=1;j<nc-1;++j)
      {
	ind=i*nc+j;
	if(byt[ind]==b
	   && (byt[ind-1]==f     || byt[ind+1]==f
	       || byt[ind+nc]==f || byt[ind-nc]==f) )
	  byt[ind]=2;
      }

  /* Set all the changed pixels to the proper values: */
  fpt=(pt=byt)+nr*nc;
  do *pt = *pt==2 ? f : *pt; while(++pt<fpt);
}





/* 8 connected dilation and erosion. b0_f1==0: Dilate the
   foreground. b0_f1==1: Erode the foreground. */
void
dilate0_erode1_8con(unsigned char *byt, size_t nr, size_t nc,
                    unsigned char b0_f1)
{
  size_t i, j, ind;
  unsigned char f, b, *pt, *fpt;

  /* Do a sanity check: */
  if(b0_f1!=1 && b0_f1!=0)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix "
          "this problem. In dilate0_erode1_4con (binary.c), the value to "
          "b0_f1 is %u while it should be 0 or 1.", PACKAGE_BUGREPORT, b0_f1);

  /* Set the foreground and background values: */
  if(b0_f1==0) {f=1; b=0;}
  else         {f=0; b=1;}

  /* Check the 4 corners: */
  if(byt[0]==b && (byt[1]==f
		   || byt[nc]==f || byt[nc+1]==f) )
    byt[0]=2;

  if(byt[nc-1]==b && (byt[nc-2]==f
		      || byt[2*nc-1]==f
		      || byt[2*nc-2]==f) )
    byt[nc-1]=2;

  if(byt[(nr-1)*nc]==b
     && ( byt[(nr-2)*nc]==f || byt[(nr-1)*nc+1]==f
	  || byt[(nr-2)*nc+1]==f) )
    byt[(nr-1)*nc]=2;

  if(byt[nr*nc-1]==b
     && ( byt[nr*nc-2]==f || byt[nr*nc-1-nc]==f
	  || byt[nr*nc-2-nc]==f) )
    byt[nr*nc-1]=2;

  /* Check the 4 sides: */
  for(j=1;j<nc-1;++j)
    if(byt[j]==b
       && ( byt[j+1]==f || byt[j-1]==f || byt[j+nc]==f
	    || byt[j-1+nc]==f || byt[j+1+nc]==f) )
      byt[j]=2;

  for(j=1;j<nc-1;++j)
    {
      ind=(nr-1)*nc+j;
      if(byt[ind]==b
	 && ( byt[ind+1]==f || byt[ind-1]==f || byt[ind-nc]==f
	      || byt[ind-1-nc]==f || byt[ind+1-nc]==f) )
	byt[ind]=2;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=i*nc;
      if(byt[ind]==b
	 && ( byt[ind+1]==f || byt[ind+nc]==f || byt[ind-nc]==f
	      || byt[ind+1-nc]==f || byt[ind+1+nc]==f) )
	byt[ind]=2;
    }

  for(i=1;i<nr-1;++i)
    {
      ind=(i+1)*nc-1;
      if(byt[ind]==b
	 && (byt[ind-1]==f || byt[ind+nc]==f || byt[ind-nc]==f
	     || byt[ind-1-nc]==f || byt[ind-1+nc]==f) )
	byt[ind]=2;
    }

  /* Check the body: */
  for(i=1;i<nr-1;++i)
    for(j=1;j<nc-1;++j)
      {
	ind=i*nc+j;
	if(byt[ind]==b
	   && (byt[ind-1]==f        || byt[ind+1]==f
	       || byt[ind+nc]==f    || byt[ind-nc]==f
	       || byt[ind-1-nc]==f  || byt[ind+1+nc]==f
	       || byt[ind-1+nc]==f  || byt[ind+1+nc]==f) )
	  byt[ind]=2;
      }

  /* Set all the changed pixels to the proper values: */
  fpt=(pt=byt)+nr*nc;
  do *pt = *pt==2 ? f : *pt; while(++pt<fpt);
}






/* Opening: erode n times then dilate n times. */
void
opening(unsigned char *byt, size_t s0, size_t s1,
        size_t depth, int con_type)
{
  size_t i;
  void (*func)(unsigned char *, size_t, size_t, unsigned char);

  /* Sanity check: */
  if(con_type!=4 && con_type!=8)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix "
          "this problem. For some reason, the value to con_type in opening "
          "(binary.c) is %d while it should be 4 or 8.", PACKAGE_BUGREPORT,
          con_type);

  /* Do the opening: */
  if(con_type==4)
    func=dilate0_erode1_4con;
  else func=dilate0_erode1_8con;

  for(i=0;i<depth;++i)
    func(byt, s0, s1, 1);

  for(i=0;i<depth;++i)
    func(byt, s0, s1, 0);
}
