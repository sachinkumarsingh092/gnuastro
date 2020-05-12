/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <gnuastro/data.h>
#include <gnuastro/statistics.h>

#include "main.h"
#include "convertt.h"





/***********************************************************************/
/**************            From mono-channel           *****************/
/***********************************************************************/
static void
color_min_max(struct converttparams *p, float *value, int min0max1)
{
  gal_data_t *tmp;
  gal_data_t *given  = min0max1 ? p->fluxhigh : p->fluxlow;
  uint8_t fixedlimit = min0max1 ? p->forcemax : p->forcemin;

  /* Find the value to write. */
  if(fixedlimit && given)
    tmp=gal_data_copy_to_new_type(given, GAL_TYPE_FLOAT32);
  else
    {
      tmp = ( min0max1
              ? gal_statistics_maximum(p->chll)
              : gal_statistics_minimum(p->chll) );
      tmp=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
    }
  *value=((float *)(tmp->array))[0];

  /* Clean up. */
  gal_data_free(tmp);
}





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

  /* Set the range of values and convert the dataset to float. */
  color_min_max(p, &min, 0);
  color_min_max(p, &max, 1);
  channel=gal_data_copy_to_new_type_free(p->chll, GAL_TYPE_FLOAT32);

  /* Allocate the three datasets to keep the RGB colors. */
  R=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "RED", NULL, "Red color channel.");
  G=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "GREEN", NULL, "Green color channel.");
  B=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "BLUE", NULL, "Blue color channel.");

  /* Start the conversion. Note that the "Choroma" ('C') is fixed by our
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

  /* Set the range, then convert the dataset to floating point. */
  color_min_max(p, &min, 0);
  color_min_max(p, &max, 1);
  channel=gal_data_copy_to_new_type_free(p->chll, GAL_TYPE_FLOAT32);

  /* Allocate the three datasets to keep the RGB colors. */
  R=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "RED", NULL, "Red color channel.");
  G=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "GREEN", NULL, "Green color channel.");
  B=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "BLUE", NULL, "Blue color channel.");

  /* Start the conversion. Note that the "Choroma" ('C') is fixed by our
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
          case 159: *r=0.995025; *g=0.296100; *b=0.000000; break;
          case 160: *r=0.993367; *g=0.263200; *b=0.000000; break;
          case 161: *r=0.991708; *g=0.230300; *b=0.000000; break;
          case 162: *r=0.990050; *g=0.197400; *b=0.000000; break;
          case 163: *r=0.988392; *g=0.164500; *b=0.000000; break;
          case 164: *r=0.986733; *g=0.131600; *b=0.000000; break;
          case 165: *r=0.985075; *g=0.098700; *b=0.000000; break;
          case 166: *r=0.983417; *g=0.065800; *b=0.000000; break;
          case 167: *r=0.981758; *g=0.032900; *b=0.000000; break;
          case 168: *r=0.980100; *g=0.000000; *b=0.000000; break;
          case 169: *r=0.955925; *g=0.000000; *b=0.000000; break;
          case 170: *r=0.931750; *g=0.000000; *b=0.000000; break;
          case 171: *r=0.907575; *g=0.000000; *b=0.000000; break;
          case 172: *r=0.883400; *g=0.000000; *b=0.000000; break;
          case 173: *r=0.859225; *g=0.000000; *b=0.000000; break;
          case 174: *r=0.835050; *g=0.000000; *b=0.000000; break;
          case 175: *r=0.810875; *g=0.000000; *b=0.000000; break;
          case 176: *r=0.786700; *g=0.000000; *b=0.000000; break;
          case 177: *r=0.762525; *g=0.000000; *b=0.000000; break;
          case 178: *r=0.738350; *g=0.000000; *b=0.000000; break;
          case 179: *r=0.714175; *g=0.000000; *b=0.000000; break;
          case 180: *r=0.690000; *g=0.000000; *b=0.000000; break;
          case 181: *r=0.715833; *g=0.083333; *b=0.083333; break;
          case 182: *r=0.741667; *g=0.166667; *b=0.166667; break;
          case 183: *r=0.767500; *g=0.250000; *b=0.250000; break;
          case 184: *r=0.793333; *g=0.333333; *b=0.333333; break;
          case 185: *r=0.819167; *g=0.416667; *b=0.416667; break;
          case 186: *r=0.845000; *g=0.500000; *b=0.500000; break;
          case 187: *r=0.870833; *g=0.583333; *b=0.583333; break;
          case 188: *r=0.896667; *g=0.666667; *b=0.666667; break;
          case 189: *r=0.922500; *g=0.750000; *b=0.750000; break;
          case 190: *r=0.948333; *g=0.833333; *b=0.833333; break;
          case 191: *r=0.974167; *g=0.916667; *b=0.916667; break;
          case 192: *r=1.000000; *g=1.000000; *b=1.000000; break;
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





/* Values here copied from PGFPlots 1.16 source. They are a little more
   precise than the matplotlib values:
     https://github.com/BIDS/colormap/blob/master/colormaps.py
   It was created by Stéfan van der Walt and Nathaniel Smith for
   Matplotlib, see here for a more complete explanation:
     https://bids.github.io/colormap */
void
color_from_mono_viridis(struct converttparams *p)
{
  gal_data_t *R, *G, *B, *channel;
  float *r, *g, *b, *f, *fp, min, max;

  /* Set the range, then convert the dataset to floating point. */
  color_min_max(p, &min, 0);
  color_min_max(p, &max, 1);
  channel=gal_data_copy_to_new_type_free(p->chll, GAL_TYPE_FLOAT32);

  /* Allocate the three datasets to keep the RGB colors. */
  R=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "RED", NULL, "Red color channel.");
  G=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "GREEN", NULL, "Green color channel.");
  B=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, channel->ndim,
                   channel->dsize, channel->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "BLUE", NULL, "Blue color channel.");

  /* Start the conversion. Note that the "Choroma" ('C') is fixed by our
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
        switch( (int)((*f-min)/(max-min)*255) )
          {
          case 0:   *r=0.26700401; *g=0.00487433; *b=0.32941519; break;
          case 1:   *r=0.26851048; *g=0.00960483; *b=0.33542652; break;
          case 2:   *r=0.26994384; *g=0.01462494; *b=0.34137895; break;
          case 3:   *r=0.27130489; *g=0.01994186; *b=0.34726862; break;
          case 4:   *r=0.27259384; *g=0.02556309; *b=0.35309303; break;
          case 5:   *r=0.27380934; *g=0.03149748; *b=0.35885256; break;
          case 6:   *r=0.27495242; *g=0.03775181; *b=0.36454323; break;
          case 7:   *r=0.27602238; *g=0.04416723; *b=0.37016418; break;
          case 8:   *r=0.2770184 ; *g=0.05034437; *b=0.37571452; break;
          case 9:   *r=0.27794143; *g=0.05632444; *b=0.38119074; break;
          case 10:  *r=0.27879067; *g=0.06214536; *b=0.38659204; break;
          case 11:  *r=0.2795655 ; *g=0.06783587; *b=0.39191723; break;
          case 12:  *r=0.28026658; *g=0.07341724; *b=0.39716349; break;
          case 13:  *r=0.28089358; *g=0.07890703; *b=0.40232944; break;
          case 14:  *r=0.28144581; *g=0.0843197 ; *b=0.40741404; break;
          case 15:  *r=0.28192358; *g=0.08966622; *b=0.41241521; break;
          case 16:  *r=0.28232739; *g=0.09495545; *b=0.41733086; break;
          case 17:  *r=0.28265633; *g=0.10019576; *b=0.42216032; break;
          case 18:  *r=0.28291049; *g=0.10539345; *b=0.42690202; break;
          case 19:  *r=0.28309095; *g=0.11055307; *b=0.43155375; break;
          case 20:  *r=0.28319704; *g=0.11567966; *b=0.43611482; break;
          case 21:  *r=0.28322882; *g=0.12077701; *b=0.44058404; break;
          case 22:  *r=0.28318684; *g=0.12584799; *b=0.44496   ; break;
          case 23:  *r=0.283072  ; *g=0.13089477; *b=0.44924127; break;
          case 24:  *r=0.28288389; *g=0.13592005; *b=0.45342734; break;
          case 25:  *r=0.28262297; *g=0.14092556; *b=0.45751726; break;
          case 26:  *r=0.28229037; *g=0.14591233; *b=0.46150995; break;
          case 27:  *r=0.28188676; *g=0.15088147; *b=0.46540474; break;
          case 28:  *r=0.28141228; *g=0.15583425; *b=0.46920128; break;
          case 29:  *r=0.28086773; *g=0.16077132; *b=0.47289909; break;
          case 30:  *r=0.28025468; *g=0.16569272; *b=0.47649762; break;
          case 31:  *r=0.27957399; *g=0.17059884; *b=0.47999675; break;
          case 32:  *r=0.27882618; *g=0.1754902 ; *b=0.48339654; break;
          case 33:  *r=0.27801236; *g=0.18036684; *b=0.48669702; break;
          case 34:  *r=0.27713437; *g=0.18522836; *b=0.48989831; break;
          case 35:  *r=0.27619376; *g=0.19007447; *b=0.49300074; break;
          case 36:  *r=0.27519116; *g=0.1949054 ; *b=0.49600488; break;
          case 37:  *r=0.27412802; *g=0.19972086; *b=0.49891131; break;
          case 38:  *r=0.27300596; *g=0.20452049; *b=0.50172076; break;
          case 39:  *r=0.27182812; *g=0.20930306; *b=0.50443413; break;
          case 40:  *r=0.27059473; *g=0.21406899; *b=0.50705243; break;
          case 41:  *r=0.26930756; *g=0.21881782; *b=0.50957678; break;
          case 42:  *r=0.26796846; *g=0.22354911; *b=0.5120084 ; break;
          case 43:  *r=0.26657984; *g=0.2282621 ; *b=0.5143487 ; break;
          case 44:  *r=0.2651445 ; *g=0.23295593; *b=0.5165993 ; break;
          case 45:  *r=0.2636632 ; *g=0.23763078; *b=0.51876163; break;
          case 46:  *r=0.26213801; *g=0.24228619; *b=0.52083736; break;
          case 47:  *r=0.26057103; *g=0.2469217 ; *b=0.52282822; break;
          case 48:  *r=0.25896451; *g=0.25153685; *b=0.52473609; break;
          case 49:  *r=0.25732244; *g=0.2561304 ; *b=0.52656332; break;
          case 50:  *r=0.25564519; *g=0.26070284; *b=0.52831152; break;
          case 51:  *r=0.25393498; *g=0.26525384; *b=0.52998273; break;
          case 52:  *r=0.25219404; *g=0.26978306; *b=0.53157905; break;
          case 53:  *r=0.25042462; *g=0.27429024; *b=0.53310261; break;
          case 54:  *r=0.24862899; *g=0.27877509; *b=0.53455561; break;
          case 55:  *r=0.2468114 ; *g=0.28323662; *b=0.53594093; break;
          case 56:  *r=0.24497208; *g=0.28767547; *b=0.53726018; break;
          case 57:  *r=0.24311324; *g=0.29209154; *b=0.53851561; break;
          case 58:  *r=0.24123708; *g=0.29648471; *b=0.53970946; break;
          case 59:  *r=0.23934575; *g=0.30085494; *b=0.54084398; break;
          case 60:  *r=0.23744138; *g=0.30520222; *b=0.5419214 ; break;
          case 61:  *r=0.23552606; *g=0.30952657; *b=0.54294396; break;
          case 62:  *r=0.23360277; *g=0.31382773; *b=0.54391424; break;
          case 63:  *r=0.2316735 ; *g=0.3181058 ; *b=0.54483444; break;
          case 64:  *r=0.22973926; *g=0.32236127; *b=0.54570633; break;
          case 65:  *r=0.22780192; *g=0.32659432; *b=0.546532  ; break;
          case 66:  *r=0.2258633 ; *g=0.33080515; *b=0.54731353; break;
          case 67:  *r=0.22392515; *g=0.334994  ; *b=0.54805291; break;
          case 68:  *r=0.22198915; *g=0.33916114; *b=0.54875211; break;
          case 69:  *r=0.22005691; *g=0.34330688; *b=0.54941304; break;
          case 70:  *r=0.21812995; *g=0.34743154; *b=0.55003755; break;
          case 71:  *r=0.21620971; *g=0.35153548; *b=0.55062743; break;
          case 72:  *r=0.21429757; *g=0.35561907; *b=0.5511844 ; break;
          case 73:  *r=0.21239477; *g=0.35968273; *b=0.55171011; break;
          case 74:  *r=0.2105031 ; *g=0.36372671; *b=0.55220646; break;
          case 75:  *r=0.20862342; *g=0.36775151; *b=0.55267486; break;
          case 76:  *r=0.20675628; *g=0.37175775; *b=0.55311653; break;
          case 77:  *r=0.20490257; *g=0.37574589; *b=0.55353282; break;
          case 78:  *r=0.20306309; *g=0.37971644; *b=0.55392505; break;
          case 79:  *r=0.20123854; *g=0.38366989; *b=0.55429441; break;
          case 80:  *r=0.1994295 ; *g=0.38760678; *b=0.55464205; break;
          case 81:  *r=0.1976365 ; *g=0.39152762; *b=0.55496905; break;
          case 82:  *r=0.19585993; *g=0.39543297; *b=0.55527637; break;
          case 83:  *r=0.19410009; *g=0.39932336; *b=0.55556494; break;
          case 84:  *r=0.19235719; *g=0.40319934; *b=0.55583559; break;
          case 85:  *r=0.19063135; *g=0.40706148; *b=0.55608907; break;
          case 86:  *r=0.18892259; *g=0.41091033; *b=0.55632606; break;
          case 87:  *r=0.18723083; *g=0.41474645; *b=0.55654717; break;
          case 88:  *r=0.18555593; *g=0.4185704 ; *b=0.55675292; break;
          case 89:  *r=0.18389763; *g=0.42238275; *b=0.55694377; break;
          case 90:  *r=0.18225561; *g=0.42618405; *b=0.5571201 ; break;
          case 91:  *r=0.18062949; *g=0.42997486; *b=0.55728221; break;
          case 92:  *r=0.17901879; *g=0.43375572; *b=0.55743035; break;
          case 93:  *r=0.17742298; *g=0.4375272 ; *b=0.55756466; break;
          case 94:  *r=0.17584148; *g=0.44128981; *b=0.55768526; break;
          case 95:  *r=0.17427363; *g=0.4450441 ; *b=0.55779216; break;
          case 96:  *r=0.17271876; *g=0.4487906 ; *b=0.55788532; break;
          case 97:  *r=0.17117615; *g=0.4525298 ; *b=0.55796464; break;
          case 98:  *r=0.16964573; *g=0.45626209; *b=0.55803034; break;
          case 99:  *r=0.16812641; *g=0.45998802; *b=0.55808199; break;
          case 100: *r=0.1666171 ; *g=0.46370813; *b=0.55811913; break;
          case 101: *r=0.16511703; *g=0.4674229 ; *b=0.55814141; break;
          case 102: *r=0.16362543; *g=0.47113278; *b=0.55814842; break;
          case 103: *r=0.16214155; *g=0.47483821; *b=0.55813967; break;
          case 104: *r=0.16066467; *g=0.47853961; *b=0.55811466; break;
          case 105: *r=0.15919413; *g=0.4822374 ; *b=0.5580728 ; break;
          case 106: *r=0.15772933; *g=0.48593197; *b=0.55801347; break;
          case 107: *r=0.15626973; *g=0.4896237 ; *b=0.557936  ; break;
          case 108: *r=0.15481488; *g=0.49331293; *b=0.55783967; break;
          case 109: *r=0.15336445; *g=0.49700003; *b=0.55772371; break;
          case 110: *r=0.1519182 ; *g=0.50068529; *b=0.55758733; break;
          case 111: *r=0.15047605; *g=0.50436904; *b=0.55742968; break;
          case 112: *r=0.14903918; *g=0.50805136; *b=0.5572505 ; break;
          case 113: *r=0.14760731; *g=0.51173263; *b=0.55704861; break;
          case 114: *r=0.14618026; *g=0.51541316; *b=0.55682271; break;
          case 115: *r=0.14475863; *g=0.51909319; *b=0.55657181; break;
          case 116: *r=0.14334327; *g=0.52277292; *b=0.55629491; break;
          case 117: *r=0.14193527; *g=0.52645254; *b=0.55599097; break;
          case 118: *r=0.14053599; *g=0.53013219; *b=0.55565893; break;
          case 119: *r=0.13914708; *g=0.53381201; *b=0.55529773; break;
          case 120: *r=0.13777048; *g=0.53749213; *b=0.55490625; break;
          case 121: *r=0.1364085 ; *g=0.54117264; *b=0.55448339; break;
          case 122: *r=0.13506561; *g=0.54485335; *b=0.55402906; break;
          case 123: *r=0.13374299; *g=0.54853458; *b=0.55354108; break;
          case 124: *r=0.13244401; *g=0.55221637; *b=0.55301828; break;
          case 125: *r=0.13117249; *g=0.55589872; *b=0.55245948; break;
          case 126: *r=0.1299327 ; *g=0.55958162; *b=0.55186354; break;
          case 127: *r=0.12872938; *g=0.56326503; *b=0.55122927; break;
          case 128: *r=0.12756771; *g=0.56694891; *b=0.55055551; break;
          case 129: *r=0.12645338; *g=0.57063316; *b=0.5498411 ; break;
          case 130: *r=0.12539383; *g=0.57431754; *b=0.54908564; break;
          case 131: *r=0.12439474; *g=0.57800205; *b=0.5482874 ; break;
          case 132: *r=0.12346281; *g=0.58168661; *b=0.54744498; break;
          case 133: *r=0.12260562; *g=0.58537105; *b=0.54655722; break;
          case 134: *r=0.12183122; *g=0.58905521; *b=0.54562298; break;
          case 135: *r=0.12114807; *g=0.59273889; *b=0.54464114; break;
          case 136: *r=0.12056501; *g=0.59642187; *b=0.54361058; break;
          case 137: *r=0.12009154; *g=0.60010387; *b=0.54253043; break;
          case 138: *r=0.11973756; *g=0.60378459; *b=0.54139999; break;
          case 139: *r=0.11951163; *g=0.60746388; *b=0.54021751; break;
          case 140: *r=0.11942341; *g=0.61114146; *b=0.53898192; break;
          case 141: *r=0.11948255; *g=0.61481702; *b=0.53769219; break;
          case 142: *r=0.11969858; *g=0.61849025; *b=0.53634733; break;
          case 143: *r=0.12008079; *g=0.62216081; *b=0.53494633; break;
          case 144: *r=0.12063824; *g=0.62582833; *b=0.53348834; break;
          case 145: *r=0.12137972; *g=0.62949242; *b=0.53197275; break;
          case 146: *r=0.12231244; *g=0.63315277; *b=0.53039808; break;
          case 147: *r=0.12344358; *g=0.63680899; *b=0.52876343; break;
          case 148: *r=0.12477953; *g=0.64046069; *b=0.52706792; break;
          case 149: *r=0.12632581; *g=0.64410744; *b=0.52531069; break;
          case 150: *r=0.12808703; *g=0.64774881; *b=0.52349092; break;
          case 151: *r=0.13006688; *g=0.65138436; *b=0.52160791; break;
          case 152: *r=0.13226797; *g=0.65501363; *b=0.51966086; break;
          case 153: *r=0.13469183; *g=0.65863619; *b=0.5176488 ; break;
          case 154: *r=0.13733921; *g=0.66225157; *b=0.51557101; break;
          case 155: *r=0.14020991; *g=0.66585927; *b=0.5134268 ; break;
          case 156: *r=0.14330291; *g=0.66945881; *b=0.51121549; break;
          case 157: *r=0.1466164 ; *g=0.67304968; *b=0.50893644; break;
          case 158: *r=0.15014782; *g=0.67663139; *b=0.5065889 ; break;
          case 159: *r=0.15389405; *g=0.68020343; *b=0.50417217; break;
          case 160: *r=0.15785146; *g=0.68376525; *b=0.50168574; break;
          case 161: *r=0.16201598; *g=0.68731632; *b=0.49912906; break;
          case 162: *r=0.1663832 ; *g=0.69085611; *b=0.49650163; break;
          case 163: *r=0.1709484 ; *g=0.69438405; *b=0.49380294; break;
          case 164: *r=0.17570671; *g=0.6978996 ; *b=0.49103252; break;
          case 165: *r=0.18065314; *g=0.70140222; *b=0.48818938; break;
          case 166: *r=0.18578266; *g=0.70489133; *b=0.48527326; break;
          case 167: *r=0.19109018; *g=0.70836635; *b=0.48228395; break;
          case 168: *r=0.19657063; *g=0.71182668; *b=0.47922108; break;
          case 169: *r=0.20221902; *g=0.71527175; *b=0.47608431; break;
          case 170: *r=0.20803045; *g=0.71870095; *b=0.4728733 ; break;
          case 171: *r=0.21400015; *g=0.72211371; *b=0.46958774; break;
          case 172: *r=0.22012381; *g=0.72550945; *b=0.46622638; break;
          case 173: *r=0.2263969 ; *g=0.72888753; *b=0.46278934; break;
          case 174: *r=0.23281498; *g=0.73224735; *b=0.45927675; break;
          case 175: *r=0.2393739 ; *g=0.73558828; *b=0.45568838; break;
          case 176: *r=0.24606968; *g=0.73890972; *b=0.45202405; break;
          case 177: *r=0.25289851; *g=0.74221104; *b=0.44828355; break;
          case 178: *r=0.25985676; *g=0.74549162; *b=0.44446673; break;
          case 179: *r=0.26694127; *g=0.74875084; *b=0.44057284; break;
          case 180: *r=0.27414922; *g=0.75198807; *b=0.4366009 ; break;
          case 181: *r=0.28147681; *g=0.75520266; *b=0.43255207; break;
          case 182: *r=0.28892102; *g=0.75839399; *b=0.42842626; break;
          case 183: *r=0.29647899; *g=0.76156142; *b=0.42422341; break;
          case 184: *r=0.30414796; *g=0.76470433; *b=0.41994346; break;
          case 185: *r=0.31192534; *g=0.76782207; *b=0.41558638; break;
          case 186: *r=0.3198086 ; *g=0.77091403; *b=0.41115215; break;
          case 187: *r=0.3277958 ; *g=0.77397953; *b=0.40664011; break;
          case 188: *r=0.33588539; *g=0.7770179 ; *b=0.40204917; break;
          case 189: *r=0.34407411; *g=0.78002855; *b=0.39738103; break;
          case 190: *r=0.35235985; *g=0.78301086; *b=0.39263579; break;
          case 191: *r=0.36074053; *g=0.78596419; *b=0.38781353; break;
          case 192: *r=0.3692142 ; *g=0.78888793; *b=0.38291438; break;
          case 193: *r=0.37777892; *g=0.79178146; *b=0.3779385 ; break;
          case 194: *r=0.38643282; *g=0.79464415; *b=0.37288606; break;
          case 195: *r=0.39517408; *g=0.79747541; *b=0.36775726; break;
          case 196: *r=0.40400101; *g=0.80027461; *b=0.36255223; break;
          case 197: *r=0.4129135 ; *g=0.80304099; *b=0.35726893; break;
          case 198: *r=0.42190813; *g=0.80577412; *b=0.35191009; break;
          case 199: *r=0.43098317; *g=0.80847343; *b=0.34647607; break;
          case 200: *r=0.44013691; *g=0.81113836; *b=0.3409673 ; break;
          case 201: *r=0.44936763; *g=0.81376835; *b=0.33538426; break;
          case 202: *r=0.45867362; *g=0.81636288; *b=0.32972749; break;
          case 203: *r=0.46805314; *g=0.81892143; *b=0.32399761; break;
          case 204: *r=0.47750446; *g=0.82144351; *b=0.31819529; break;
          case 205: *r=0.4870258 ; *g=0.82392862; *b=0.31232133; break;
          case 206: *r=0.49661536; *g=0.82637633; *b=0.30637661; break;
          case 207: *r=0.5062713 ; *g=0.82878621; *b=0.30036211; break;
          case 208: *r=0.51599182; *g=0.83115784; *b=0.29427888; break;
          case 209: *r=0.52577622; *g=0.83349064; *b=0.2881265 ; break;
          case 210: *r=0.5356211 ; *g=0.83578452; *b=0.28190832; break;
          case 211: *r=0.5455244 ; *g=0.83803918; *b=0.27562602; break;
          case 212: *r=0.55548397; *g=0.84025437; *b=0.26928147; break;
          case 213: *r=0.5654976 ; *g=0.8424299 ; *b=0.26287683; break;
          case 214: *r=0.57556297; *g=0.84456561; *b=0.25641457; break;
          case 215: *r=0.58567772; *g=0.84666139; *b=0.24989748; break;
          case 216: *r=0.59583934; *g=0.84871722; *b=0.24332878; break;
          case 217: *r=0.60604528; *g=0.8507331 ; *b=0.23671214; break;
          case 218: *r=0.61629283; *g=0.85270912; *b=0.23005179; break;
          case 219: *r=0.62657923; *g=0.85464543; *b=0.22335258; break;
          case 220: *r=0.63690157; *g=0.85654226; *b=0.21662012; break;
          case 221: *r=0.64725685; *g=0.85839991; *b=0.20986086; break;
          case 222: *r=0.65764197; *g=0.86021878; *b=0.20308229; break;
          case 223: *r=0.66805369; *g=0.86199932; *b=0.19629307; break;
          case 224: *r=0.67848868; *g=0.86374211; *b=0.18950326; break;
          case 225: *r=0.68894351; *g=0.86544779; *b=0.18272455; break;
          case 226: *r=0.69941463; *g=0.86711711; *b=0.17597055; break;
          case 227: *r=0.70989842; *g=0.86875092; *b=0.16925712; break;
          case 228: *r=0.72039115; *g=0.87035015; *b=0.16260273; break;
          case 229: *r=0.73088902; *g=0.87191584; *b=0.15602894; break;
          case 230: *r=0.74138803; *g=0.87344918; *b=0.14956101; break;
          case 231: *r=0.75188414; *g=0.87495143; *b=0.14322828; break;
          case 232: *r=0.76237342; *g=0.87642392; *b=0.13706449; break;
          case 233: *r=0.77285183; *g=0.87786808; *b=0.13110864; break;
          case 234: *r=0.78331535; *g=0.87928545; *b=0.12540538; break;
          case 235: *r=0.79375994; *g=0.88067763; *b=0.12000532; break;
          case 236: *r=0.80418159; *g=0.88204632; *b=0.11496505; break;
          case 237: *r=0.81457634; *g=0.88339329; *b=0.11034678; break;
          case 238: *r=0.82494028; *g=0.88472036; *b=0.10621724; break;
          case 239: *r=0.83526959; *g=0.88602943; *b=0.1026459 ; break;
          case 240: *r=0.84556056; *g=0.88732243; *b=0.09970219; break;
          case 241: *r=0.8558096 ; *g=0.88860134; *b=0.09745186; break;
          case 242: *r=0.86601325; *g=0.88986815; *b=0.09595277; break;
          case 243: *r=0.87616824; *g=0.89112487; *b=0.09525046; break;
          case 244: *r=0.88627146; *g=0.89237353; *b=0.09537439; break;
          case 245: *r=0.89632002; *g=0.89361614; *b=0.09633538; break;
          case 246: *r=0.90631121; *g=0.89485467; *b=0.09812496; break;
          case 247: *r=0.91624212; *g=0.89609127; *b=0.1007168 ; break;
          case 248: *r=0.92610579; *g=0.89732977; *b=0.10407067; break;
          case 249: *r=0.93590444; *g=0.8985704 ; *b=0.10813094; break;
          case 250: *r=0.94563626; *g=0.899815  ; *b=0.11283773; break;
          case 251: *r=0.95529972; *g=0.90106534; *b=0.11812832; break;
          case 252: *r=0.96489353; *g=0.90232311; *b=0.12394051; break;
          case 253: *r=0.97441665; *g=0.90358991; *b=0.13021494; break;
          case 254: *r=0.98386829; *g=0.90486726; *b=0.13689671; break;
          case 255: *r=0.99324789; *g=0.90615657; *b=0.1439362 ; break;
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





void
color_map_prepare(struct converttparams *p)
{
  switch(p->colormap->status)
    {
    case COLOR_HSV:          color_from_mono_hsv(p); break;
    case COLOR_SLS:          color_from_mono_sls(p); break;
    case COLOR_VIRIDIS:      color_from_mono_viridis(p); break;
    case COLOR_GRAY:         convertt_scale_to_uchar(p); break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
            "the problem. The value %d is not a recognized color-space "
            "code", __func__, PACKAGE_BUGREPORT, p->colormap->status);
    }
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
                   p->cp.quietmmap, "HUE", NULL, NULL);
  S=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, p->chll->ndim,
                   p->chll->dsize, p->chll->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "SATURATION", NULL, NULL);
  V=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, p->chll->ndim,
                   p->chll->dsize, p->chll->wcs, 0, p->cp.minmapsize,
                   p->cp.quietmmap, "VALUE", NULL, NULL);

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
            /* When 'max==0', then *r=*g=*b=0, so s=h=0. */
            *s=*h=0.0;
        }
      else
        /* When there is no difference, then its actually a grayscale
           dataset, so '*v' is the only parameter that matters. */
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
