/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2019, Free Software Foundation, Inc.

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

#include <gnuastro/data.h>

#include "main.h"







/***********************************************************************/
/**************            From mono-channel           *****************/
/***********************************************************************/
/* This algorithm is a translation of the primary algorithm in this page:
https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both */
void
color_from_mono_hsv(struct converttparams *p)
{
  int i;
  gal_data_t *R, *G, *B, *channel;
  float *r, *g, *b, *f, *fp, min, max;
  float *params=p->colormap->next->array;
  float h, s=1, v, hh, ff, P, q, t, h_min, h_max;

  /* Set the input values. */
  h_min = params[0];
  h_max = params[1];

  /* Sanity checks. */
  if(h_min>h_max)
    error(EXIT_FAILURE, 0, "the minimum angle value (%g) is not smaller "
          "than the maximum (%f)", h_min, h_max);
  if(h_min<0)
    error(EXIT_FAILURE, 0, "the minimum angle (%g) must be larger than 0",
          h_min);
  if(h_max>360)
    error(EXIT_FAILURE, 0, "the maximum angle (%g) must be smaller than "
          "360", h_max);

  /* Convert the dataset to floating point, then change its range to the
     given angle values. */
  gal_type_min(GAL_TYPE_FLOAT32, &max);
  gal_type_max(GAL_TYPE_FLOAT32, &min);
  channel=gal_data_copy_to_new_type_free(p->chll, GAL_TYPE_FLOAT32);
  fp=(f=channel->array)+channel->size;
  do {if(*f<min) min=*f; if(*f>max) max=*f;} while(++f<fp);

  /* Allocate the three datasets to keep the RGB colors. */
  R=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   "RED", NULL, "Red color channel.");
  G=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   "GREEN", NULL, "Green color channel.");
  B=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   "BLUE", NULL, "Blue color channel.");

  /* Start the conversion. Note that the "Choroma" (`C') is fixed by our
     definition. */
  r=R->array;
  g=G->array;
  b=B->array;
  fp=(f=channel->array)+channel->size;
  do
    {
      if(isnan(*f))
        *r=*g=*b=0.0;
      else
        {
          /* Shift all the values to start from the desired angle. We'll
             set the "value" (brightness: 0: dark, 1: bright) using the
             pixel value, then scale it to fix the color. */
          v=(*f-min) / (max-min);
          h = v * (h_max-h_min) + h_min;
          if(h==360) h=0;

          /* Prepare the intermediate values. */
          hh=h/60;
          i  = (int)hh;
          ff = hh - i;
          P  = v * ( 1.0 -  s );
          q  = v * ( 1.0 - (s * ff) );
          t  = v * ( 1.0 - (s * (1.0 - ff) ));

          /* Based on the integer phase, set the r,g,b values. */
          switch(i)
            {
            case 0:          *r=v; *g=t; *b=P; break;
            case 1:          *r=q; *g=v; *b=P; break;
            case 2:          *r=P; *g=v; *b=t; break;
            case 3:          *r=P; *g=q; *b=v; break;
            case 4:          *r=t; *g=P; *b=v; break;
            case 5: default: *r=v; *g=P; *b=q; break;
            }
        }

      /* Convert the RGB values to 0 and 255.  With the steps above, they
         are between 0 and 1. */
      *r++ *= 255;
      *g++ *= 255;
      *b++ *= 255;
    }
  while(++f<fp);

  /* Convert the type to unsigned char. */
  R=gal_data_copy_to_new_type_free(R, GAL_TYPE_UINT8);
  G=gal_data_copy_to_new_type_free(G, GAL_TYPE_UINT8);
  B=gal_data_copy_to_new_type_free(B, GAL_TYPE_UINT8);
  p->chll=R;
  p->chll->next=G;
  p->chll->next->next=B;

  /* Clean up. */
  gal_data_free(channel);
}





/* From SAO DS9: */
void
color_from_mono_sls(struct converttparams *p)
{
  gal_data_t *R, *G, *B, *channel;
  float *r, *g, *b, *f, *fp, min, max;

  /* Convert the dataset to floating point, then find its minimum and
     maximum values. */
  gal_type_min(GAL_TYPE_FLOAT32, &max);
  gal_type_max(GAL_TYPE_FLOAT32, &min);
  channel=gal_data_copy_to_new_type_free(p->chll, GAL_TYPE_FLOAT32);
  fp=(f=channel->array)+channel->size;
  do {if(*f<min) min=*f; if(*f>max) max=*f;} while(++f<fp);

  /* Allocate the three datasets to keep the RGB colors. */
  R=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   "RED", NULL, "Red color channel.");
  G=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   "GREEN", NULL, "Green color channel.");
  B=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   "BLUE", NULL, "Blue color channel.");

  /* Start the conversion. Note that the "Choroma" (`C') is fixed by our
     definition. */
  r=R->array;
  g=G->array;
  b=B->array;
  fp=(f=channel->array)+channel->size;
  do
    {
      if(isnan(*f))
        *r=*g=*b=0.0;
      else
        switch( (int)((*f-min)/(max-min)*200) )
          {
          case 0:   *r=0.000000; *g=0.000000; *b=0.000000; break;
          case 1:   *r=0.043442; *g=0.000000; *b=0.052883; break;
          case 2:   *r=0.086883; *g=0.000000; *b=0.105767; break;
          case 3:   *r=0.130325; *g=0.000000; *b=0.158650; break;
          case 4:   *r=0.173767; *g=0.000000; *b=0.211533; break;
          case 5:   *r=0.217208; *g=0.000000; *b=0.264417; break;
          case 6:   *r=0.260650; *g=0.000000; *b=0.317300; break;
          case 7:   *r=0.304092; *g=0.000000; *b=0.370183; break;
          case 8:   *r=0.347533; *g=0.000000; *b=0.423067; break;
          case 9:   *r=0.390975; *g=0.000000; *b=0.475950; break;
          case 10:  *r=0.434417; *g=0.000000; *b=0.528833; break;
          case 11:  *r=0.477858; *g=0.000000; *b=0.581717; break;
          case 12:  *r=0.521300; *g=0.000000; *b=0.634600; break;
          case 13:  *r=0.506742; *g=0.000000; *b=0.640217; break;
          case 14:  *r=0.492183; *g=0.000000; *b=0.645833; break;
          case 15:  *r=0.477625; *g=0.000000; *b=0.651450; break;
          case 16:  *r=0.463067; *g=0.000000; *b=0.657067; break;
          case 17:  *r=0.448508; *g=0.000000; *b=0.662683; break;
          case 18:  *r=0.433950; *g=0.000000; *b=0.668300; break;
          case 19:  *r=0.419392; *g=0.000000; *b=0.673917; break;
          case 20:  *r=0.404833; *g=0.000000; *b=0.679533; break;
          case 21:  *r=0.390275; *g=0.000000; *b=0.685150; break;
          case 22:  *r=0.375717; *g=0.000000; *b=0.690767; break;
          case 23:  *r=0.361158; *g=0.000000; *b=0.696383; break;
          case 24:  *r=0.346600; *g=0.000000; *b=0.702000; break;
          case 25:  *r=0.317717; *g=0.000000; *b=0.712192; break;
          case 26:  *r=0.288833; *g=0.000000; *b=0.722383; break;
          case 27:  *r=0.259950; *g=0.000000; *b=0.732575; break;
          case 28:  *r=0.231067; *g=0.000000; *b=0.742767; break;
          case 29:  *r=0.202183; *g=0.000000; *b=0.752958; break;
          case 30:  *r=0.173300; *g=0.000000; *b=0.763150; break;
          case 31:  *r=0.144417; *g=0.000000; *b=0.773342; break;
          case 32:  *r=0.115533; *g=0.000000; *b=0.783533; break;
          case 33:  *r=0.086650; *g=0.000000; *b=0.793725; break;
          case 34:  *r=0.057767; *g=0.000000; *b=0.803917; break;
          case 35:  *r=0.028883; *g=0.000000; *b=0.814108; break;
          case 36:  *r=0.000000; *g=0.000000; *b=0.824300; break;
          case 37:  *r=0.000000; *g=0.019817; *b=0.838942; break;
          case 38:  *r=0.000000; *g=0.039633; *b=0.853583; break;
          case 39:  *r=0.000000; *g=0.059450; *b=0.868225; break;
          case 40:  *r=0.000000; *g=0.079267; *b=0.882867; break;
          case 41:  *r=0.000000; *g=0.099083; *b=0.897508; break;
          case 42:  *r=0.000000; *g=0.118900; *b=0.912150; break;
          case 43:  *r=0.000000; *g=0.138717; *b=0.926792; break;
          case 44:  *r=0.000000; *g=0.158533; *b=0.941433; break;
          case 45:  *r=0.000000; *g=0.178350; *b=0.956075; break;
          case 46:  *r=0.000000; *g=0.198167; *b=0.970717; break;
          case 47:  *r=0.000000; *g=0.217983; *b=0.985358; break;
          case 48:  *r=0.000000; *g=0.237800; *b=1.000000; break;
          case 49:  *r=0.000000; *g=0.268533; *b=1.000000; break;
          case 50:  *r=0.000000; *g=0.299267; *b=1.000000; break;
          case 51:  *r=0.000000; *g=0.330000; *b=1.000000; break;
          case 52:  *r=0.000000; *g=0.360733; *b=1.000000; break;
          case 53:  *r=0.000000; *g=0.391467; *b=1.000000; break;
          case 54:  *r=0.000000; *g=0.422200; *b=1.000000; break;
          case 55:  *r=0.000000; *g=0.452933; *b=1.000000; break;
          case 56:  *r=0.000000; *g=0.483667; *b=1.000000; break;
          case 57:  *r=0.000000; *g=0.514400; *b=1.000000; break;
          case 58:  *r=0.000000; *g=0.545133; *b=1.000000; break;
          case 59:  *r=0.000000; *g=0.575867; *b=1.000000; break;
          case 60:  *r=0.000000; *g=0.606600; *b=1.000000; break;
          case 61:  *r=0.000000; *g=0.631733; *b=0.975300; break;
          case 62:  *r=0.000000; *g=0.656867; *b=0.950600; break;
          case 63:  *r=0.000000; *g=0.682000; *b=0.925900; break;
          case 64:  *r=0.000000; *g=0.707133; *b=0.901200; break;
          case 65:  *r=0.000000; *g=0.732267; *b=0.876500; break;
          case 66:  *r=0.000000; *g=0.757400; *b=0.851800; break;
          case 67:  *r=0.000000; *g=0.782533; *b=0.827100; break;
          case 68:  *r=0.000000; *g=0.807667; *b=0.802400; break;
          case 69:  *r=0.000000; *g=0.832800; *b=0.777700; break;
          case 70:  *r=0.000000; *g=0.857933; *b=0.753000; break;
          case 71:  *r=0.000000; *g=0.883067; *b=0.728300; break;
          case 72:  *r=0.000000; *g=0.908200; *b=0.703600; break;
          case 73:  *r=0.000000; *g=0.901908; *b=0.676675; break;
          case 74:  *r=0.000000; *g=0.895617; *b=0.649750; break;
          case 75:  *r=0.000000; *g=0.889325; *b=0.622825; break;
          case 76:  *r=0.000000; *g=0.883033; *b=0.595900; break;
          case 77:  *r=0.000000; *g=0.876742; *b=0.568975; break;
          case 78:  *r=0.000000; *g=0.870450; *b=0.542050; break;
          case 79:  *r=0.000000; *g=0.864158; *b=0.515125; break;
          case 80:  *r=0.000000; *g=0.857867; *b=0.488200; break;
          case 81:  *r=0.000000; *g=0.851575; *b=0.461275; break;
          case 82:  *r=0.000000; *g=0.845283; *b=0.434350; break;
          case 83:  *r=0.000000; *g=0.838992; *b=0.407425; break;
          case 84:  *r=0.000000; *g=0.832700; *b=0.380500; break;
          case 85:  *r=0.000000; *g=0.832308; *b=0.354858; break;
          case 86:  *r=0.000000; *g=0.831917; *b=0.329217; break;
          case 87:  *r=0.000000; *g=0.831525; *b=0.303575; break;
          case 88:  *r=0.000000; *g=0.831133; *b=0.277933; break;
          case 89:  *r=0.000000; *g=0.830742; *b=0.252292; break;
          case 90:  *r=0.000000; *g=0.830350; *b=0.226650; break;
          case 91:  *r=0.000000; *g=0.829958; *b=0.201008; break;
          case 92:  *r=0.000000; *g=0.829567; *b=0.175367; break;
          case 93:  *r=0.000000; *g=0.829175; *b=0.149725; break;
          case 94:  *r=0.000000; *g=0.828783; *b=0.124083; break;
          case 95:  *r=0.000000; *g=0.828392; *b=0.098442; break;
          case 96:  *r=0.000000; *g=0.828000; *b=0.072800; break;
          case 97:  *r=0.033167; *g=0.834167; *b=0.066733; break;
          case 98:  *r=0.066333; *g=0.840333; *b=0.060667; break;
          case 99:  *r=0.099500; *g=0.846500; *b=0.054600; break;
          case 100: *r=0.132667; *g=0.852667; *b=0.048533; break;
          case 101: *r=0.165833; *g=0.858833; *b=0.042467; break;
          case 102: *r=0.199000; *g=0.865000; *b=0.036400; break;
          case 103: *r=0.232167; *g=0.871167; *b=0.030333; break;
          case 104: *r=0.265333; *g=0.877333; *b=0.024267; break;
          case 105: *r=0.298500; *g=0.883500; *b=0.018200; break;
          case 106: *r=0.331667; *g=0.889667; *b=0.012133; break;
          case 107: *r=0.364833; *g=0.895833; *b=0.006067; break;
          case 108: *r=0.398000; *g=0.902000; *b=0.000000; break;
          case 109: *r=0.430950; *g=0.902000; *b=0.000000; break;
          case 110: *r=0.463900; *g=0.902000; *b=0.000000; break;
          case 111: *r=0.496850; *g=0.902000; *b=0.000000; break;
          case 112: *r=0.529800; *g=0.902000; *b=0.000000; break;
          case 113: *r=0.562750; *g=0.902000; *b=0.000000; break;
          case 114: *r=0.595700; *g=0.902000; *b=0.000000; break;
          case 115: *r=0.628650; *g=0.902000; *b=0.000000; break;
          case 116: *r=0.661600; *g=0.902000; *b=0.000000; break;
          case 117: *r=0.694550; *g=0.902000; *b=0.000000; break;
          case 118: *r=0.727500; *g=0.902000; *b=0.000000; break;
          case 119: *r=0.760450; *g=0.902000; *b=0.000000; break;
          case 120: *r=0.793400; *g=0.902000; *b=0.000000; break;
          case 121: *r=0.810617; *g=0.897133; *b=0.003983; break;
          case 122: *r=0.827833; *g=0.892267; *b=0.007967; break;
          case 123: *r=0.845050; *g=0.887400; *b=0.011950; break;
          case 124: *r=0.862267; *g=0.882533; *b=0.015933; break;
          case 125: *r=0.879483; *g=0.877667; *b=0.019917; break;
          case 126: *r=0.896700; *g=0.872800; *b=0.023900; break;
          case 127: *r=0.913917; *g=0.867933; *b=0.027883; break;
          case 128: *r=0.931133; *g=0.863067; *b=0.031867; break;
          case 129: *r=0.948350; *g=0.858200; *b=0.035850; break;
          case 130: *r=0.965567; *g=0.853333; *b=0.039833; break;
          case 131: *r=0.982783; *g=0.848467; *b=0.043817; break;
          case 132: *r=1.000000; *g=0.843600; *b=0.047800; break;
          case 133: *r=0.995725; *g=0.824892; *b=0.051600; break;
          case 134: *r=0.991450; *g=0.806183; *b=0.055400; break;
          case 135: *r=0.987175; *g=0.787475; *b=0.059200; break;
          case 136: *r=0.982900; *g=0.768767; *b=0.063000; break;
          case 137: *r=0.978625; *g=0.750058; *b=0.066800; break;
          case 138: *r=0.974350; *g=0.731350; *b=0.070600; break;
          case 139: *r=0.970075; *g=0.712642; *b=0.074400; break;
          case 140: *r=0.965800; *g=0.693933; *b=0.078200; break;
          case 141: *r=0.961525; *g=0.675225; *b=0.082000; break;
          case 142: *r=0.957250; *g=0.656517; *b=0.085800; break;
          case 143: *r=0.952975; *g=0.637808; *b=0.089600; break;
          case 144: *r=0.948700; *g=0.619100; *b=0.093400; break;
          case 145: *r=0.952975; *g=0.600408; *b=0.085617; break;
          case 146: *r=0.957250; *g=0.581717; *b=0.077833; break;
          case 147: *r=0.961525; *g=0.563025; *b=0.070050; break;
          case 148: *r=0.965800; *g=0.544333; *b=0.062267; break;
          case 149: *r=0.970075; *g=0.525642; *b=0.054483; break;
          case 150: *r=0.974350; *g=0.506950; *b=0.046700; break;
          case 151: *r=0.978625; *g=0.488258; *b=0.038917; break;
          case 152: *r=0.982900; *g=0.469567; *b=0.031133; break;
          case 153: *r=0.987175; *g=0.450875; *b=0.023350; break;
          case 154: *r=0.991450; *g=0.432183; *b=0.015567; break;
          case 155: *r=0.995725; *g=0.413492; *b=0.007783; break;
          case 156: *r=1.000000; *g=0.394800; *b=0.000000; break;
          case 157: *r=0.998342; *g=0.361900; *b=0.000000; break;
          case 158: *r=0.996683; *g=0.329000; *b=0.000000; break;
          case 160: *r=0.995025; *g=0.296100; *b=0.000000; break;
          case 161: *r=0.993367; *g=0.263200; *b=0.000000; break;
          case 162: *r=0.991708; *g=0.230300; *b=0.000000; break;
          case 163: *r=0.990050; *g=0.197400; *b=0.000000; break;
          case 164: *r=0.988392; *g=0.164500; *b=0.000000; break;
          case 165: *r=0.986733; *g=0.131600; *b=0.000000; break;
          case 166: *r=0.985075; *g=0.098700; *b=0.000000; break;
          case 167: *r=0.983417; *g=0.065800; *b=0.000000; break;
          case 168: *r=0.981758; *g=0.032900; *b=0.000000; break;
          case 169: *r=0.980100; *g=0.000000; *b=0.000000; break;
          case 170: *r=0.955925; *g=0.000000; *b=0.000000; break;
          case 171: *r=0.931750; *g=0.000000; *b=0.000000; break;
          case 172: *r=0.907575; *g=0.000000; *b=0.000000; break;
          case 173: *r=0.883400; *g=0.000000; *b=0.000000; break;
          case 174: *r=0.859225; *g=0.000000; *b=0.000000; break;
          case 175: *r=0.835050; *g=0.000000; *b=0.000000; break;
          case 176: *r=0.810875; *g=0.000000; *b=0.000000; break;
          case 177: *r=0.786700; *g=0.000000; *b=0.000000; break;
          case 178: *r=0.762525; *g=0.000000; *b=0.000000; break;
          case 179: *r=0.738350; *g=0.000000; *b=0.000000; break;
          case 180: *r=0.714175; *g=0.000000; *b=0.000000; break;
          case 181: *r=0.690000; *g=0.000000; *b=0.000000; break;
          case 182: *r=0.715833; *g=0.083333; *b=0.083333; break;
          case 183: *r=0.741667; *g=0.166667; *b=0.166667; break;
          case 184: *r=0.767500; *g=0.250000; *b=0.250000; break;
          case 185: *r=0.793333; *g=0.333333; *b=0.333333; break;
          case 186: *r=0.819167; *g=0.416667; *b=0.416667; break;
          case 187: *r=0.845000; *g=0.500000; *b=0.500000; break;
          case 188: *r=0.870833; *g=0.583333; *b=0.583333; break;
          case 189: *r=0.896667; *g=0.666667; *b=0.666667; break;
          case 190: *r=0.922500; *g=0.750000; *b=0.750000; break;
          case 191: *r=0.948333; *g=0.833333; *b=0.833333; break;
          case 192: *r=0.974167; *g=0.916667; *b=0.916667; break;
          case 193: *r=1.000000; *g=1.000000; *b=1.000000; break;
          case 194: *r=1.000000; *g=1.000000; *b=1.000000; break;
          case 195: *r=1.000000; *g=1.000000; *b=1.000000; break;
          case 196: *r=1.000000; *g=1.000000; *b=1.000000; break;
          case 197: *r=1.000000; *g=1.000000; *b=1.000000; break;
          case 198: *r=1.000000; *g=1.000000; *b=1.000000; break;
          case 199: *r=1.000000; *g=1.000000; *b=1.000000; break;
          case 200: *r=1.000000; *g=1.000000; *b=1.000000; break;
          }

      /* Convert the RGB values to 0 and 255.  With the steps above, they
         are between 0 and 1. */
      *r++ *= 255;
      *g++ *= 255;
      *b++ *= 255;
    }
  while(++f<fp);

  /* Convert the type to unsigned char. */
  R=gal_data_copy_to_new_type_free(R, GAL_TYPE_UINT8);
  G=gal_data_copy_to_new_type_free(G, GAL_TYPE_UINT8);
  B=gal_data_copy_to_new_type_free(B, GAL_TYPE_UINT8);
  p->chll=R;
  p->chll->next=G;
  p->chll->next->next=B;

  /* Clean up. */
  gal_data_free(channel);
}




















/***********************************************************************/
/**************            From mono-channel           *****************/
/***********************************************************************/
/* This algorithm is a translation of the primary algorithm in this page:
https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both */
void
color_rgb_to_hsv(struct converttparams *p)
{
  float *h, *s, *v;
  gal_data_t *H, *S, *V;
  uint8_t *r, *g, *b, *rr, min, max, delta;

  /* Basic sanity checks. */
  if(gal_list_data_number(p->chll)!=3)
    error(EXIT_FAILURE, 0, "%s: three color channels must be input",
          __func__);
  if( p->chll->type!=GAL_TYPE_UINT8
      || p->chll->next->type!=GAL_TYPE_UINT8
      || p->chll->next->next->type!=GAL_TYPE_UINT8 )
    error(EXIT_FAILURE, 0, "when converting RGB to HSV, all three input "
          "color channels must have an 8-bit unsigned integer type");

  /* Allocate the three datasets to keep the RGB colors. */
  H=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, p->chll->ndim,
                   p->chll->dsize, p->chll->wcs, 0, p->cp.minmapsize,
                   "HUE", NULL, NULL);
  S=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, p->chll->ndim,
                   p->chll->dsize, p->chll->wcs, 0, p->cp.minmapsize,
                   "SATURATION", NULL, NULL);
  V=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, p->chll->ndim,
                   p->chll->dsize, p->chll->wcs, 0, p->cp.minmapsize,
                   "VALUE", NULL, NULL);

  /* Initiate the pointer arrays. */
  h=H->array;
  s=S->array;
  v=V->array;
  r=p->chll->array;
  g=p->chll->next->array;
  b=p->chll->next->next->array;

  /* Parse the dataset and do the conversion. */
  rr=r+p->chll->size;
  do
    {
      /* Get the minimum and maximum RGB values. */
      min = *r  < *g ? *r  : *g;
      min = min < *b ? min : *b;
      max = *r  > *g ? *r  : *g;
      max = max > *b ? max : *b;

      /* The "value" is the maximum. */
      *v = (float)(max)/255.0;

      /* See what the difference between the minimum and maximum are. */
      delta=max-min;
      if(delta)
        {
          if(max)
            {
              /* Set the Saturation and hue. */
              *s = (float)(delta)/(float)(max);
              *h = ( *r==max
                     /* Between yellow and magenta. */
                     ? ((float)(*g-*b)/(float)(delta))
                     : ( *g==max
                         /* Between cyan & yellow. */
                         ? (2.0+(float)(*b-*r)/(float)(delta))
                         /* Between magenta & cyan. */
                         : (4.0+(float)(*r-*g)/(float)(delta)) ) );

              /* Correct the hue: */
              *h *= 60.0;
              if( *h<0.0 ) *h += 360.0;
            }
          else
            /* When `max==0', then *r=*g=*b=0, so s=h=0. */
            *s=*h=0.0;
        }
      else
        /* When there is no difference, then its actually a grayscale
           dataset, so `*v' is the only parameter that matters. */
        *s=*h=0.0;


      /* Increment all the pointers. */
      ++g; ++b; ++h; ++s; ++v;
    }
  while(++r<rr);

  /* Free the old channels linked list and replace it with the new ones. */
  gal_list_data_free(p->chll);
  p->chll=H;
  p->chll->next=S;
  p->chll->next->next=V;
}
