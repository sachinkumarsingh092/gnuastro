/*********************************************************************
neighbors.h -- Find the neighbours around a pixel.
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
#ifndef NEIGHBORS_H
#define NEIGHBORS_H


/* These macros can be put any where you want to know the neighbors of
   an index, the important values that have to be declared and defined
   before this macro are:

   Inputs:
   ind        : Pointer to the index being considered.
   (is0, is1) : Input image number of columns.
   (x0,y0)    : The bottom left corner of image or mesh box.
   (x1,y1)    : The Top right corner of image or mesh box.

   Outputs:
   numngb     : Number of filled neighbors.
   *ngb       : Array keeping the neighbours indexs (allocated outside).
*/

#define FILL_NGB_4_REGION {			                   \
    numngb=0;							   \
    if (*ind/is1!=x0  ) ngb[numngb++]=*ind-is1;			   \
    if (*ind/is1!=x1-1) ngb[numngb++]=*ind+is1;			   \
    if (*ind%is1!=y0  ) ngb[numngb++]=*ind-1;			   \
    if (*ind%is1!=y1-1) ngb[numngb++]=*ind+1;			   \
  }


#define FILL_NGB_4_ALLIMG {					   \
    numngb=0;							   \
    if (*ind/is1!=0    ) ngb[numngb++]=*ind-is1;		   \
    if (*ind/is1!=is0-1) ngb[numngb++]=*ind+is1;		   \
    if (*ind%is1!=0    ) ngb[numngb++]=*ind-1;			   \
    if (*ind%is1!=is1-1) ngb[numngb++]=*ind+1;			   \
  }


#define FILL_NGB_8_REGION {					   \
    unsigned char bl=0, br=0, tl=0, tr=0;			   \
    numngb=0;						           \
    if (*ind/is1!=x0  ) {ngb[numngb++]=*ind-is1; ++bl; ++br;}	   \
    if (*ind/is1!=x1-1) {ngb[numngb++]=*ind+is1; ++tl; ++tr;}	   \
    if (*ind%is1!=y0  ) {ngb[numngb++]=*ind-1;   ++bl; ++tl;}	   \
    if (*ind%is1!=y1-1) {ngb[numngb++]=*ind+1;   ++tr; ++br;}	   \
    if (numngb==4)						   \
      {								   \
	numngb=8;						   \
	ngb[4]=*ind-is1-1; ngb[5]=*ind-is1+1;			   \
	ngb[6]=*ind+is1-1; ngb[7]=*ind+is1+1;			   \
      }								   \
    else							   \
      {								   \
	if(bl==2) ngb[numngb++]=*ind-is1-1;			   \
	if(br==2) ngb[numngb++]=*ind-is1+1;			   \
	if(tl==2) ngb[numngb++]=*ind+is1-1;			   \
	if(tr==2) ngb[numngb++]=*ind+is1+1;			   \
      }								   \
  }

/* Needs:
      size_t ngb[8], is1, numngb*/
#define FILL_NGB_8_ALLIMG {					   \
    unsigned char bl=0, br=0, tl=0, tr=0;			   \
    numngb=0;						           \
    if (*ind/is1!=0     ) {ngb[numngb++]=*ind-is1; ++bl; ++br;}	   \
    if (*ind/is1!=is0-1 ) {ngb[numngb++]=*ind+is1; ++tl; ++tr;}	   \
    if (*ind%is1!=0     ) {ngb[numngb++]=*ind-1;   ++bl; ++tl;}	   \
    if (*ind%is1!=is1-1 ) {ngb[numngb++]=*ind+1;   ++tr; ++br;}	   \
    if (numngb==4)						   \
      {								   \
	numngb=8;						   \
	ngb[4]=*ind-is1-1; ngb[5]=*ind-is1+1;			   \
	ngb[6]=*ind+is1-1; ngb[7]=*ind+is1+1;			   \
      }								   \
    else							   \
      {								   \
	if(bl==2) ngb[numngb++]=*ind-is1-1;			   \
	if(br==2) ngb[numngb++]=*ind-is1+1;			   \
	if(tl==2) ngb[numngb++]=*ind+is1-1;			   \
	if(tr==2) ngb[numngb++]=*ind+is1+1;			   \
      }								   \
  }


#define FILL_NGB_8_ALLIMG_IJ {					   \
    unsigned char bl=0, br=0, tl=0, tr=0;			   \
    numngb=0;						           \
    if (i!=0     ) {ngb[numngb++]=(i-1)*is1+j; ++bl; ++br;}	   \
    if (i!=is0-1 ) {ngb[numngb++]=(i+1)*is1+j; ++tl; ++tr;}	   \
    if (j!=0     ) {ngb[numngb++]=i*is1+j-1;   ++bl; ++tl;}	   \
    if (j!=is1-1 ) {ngb[numngb++]=i*is1+j+1;   ++tr; ++br;}	   \
    if (numngb==4)						   \
      {								   \
	numngb=8;						   \
	ngb[4]=(i-1)*is1+j-1; ngb[5]=(i-1)*is1+j+1;		   \
	ngb[6]=(i+1)*is1+j-1; ngb[7]=(i+1)*is1+j+1;		   \
      }								   \
    else							   \
      {								   \
	if(bl==2) ngb[numngb++]=(i-1)*is1+j-1;			   \
	if(br==2) ngb[numngb++]=(i-1)*is1+j+1;			   \
	if(tl==2) ngb[numngb++]=(i+1)*is1+j-1;			   \
	if(tr==2) ngb[numngb++]=(i+1)*is1+j+1;			   \
      }								   \
  }

#endif
