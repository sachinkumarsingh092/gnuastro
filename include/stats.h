/*********************************************************************
Functions for simple statistical analysis.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef STATS_H
#define STATS_H

#define MINFD -1e30
#define MAXFD 1e30


void
floatmin(float *in, size_t size, float *min);

void
floatmax(float *in, size_t size, float *max);

void
floatmaxmasked(float *in, unsigned char *mask, size_t size, float *max);

void
floatsecondmax(float *in, size_t size, float *secondmax);

void
floatsecondmin(float *in, size_t size, float *secondmin);

void
fminmax(float *in, size_t size, float *min, float *max);

void
dmax_withindex(double *in, size_t size, double *max, size_t *index);

void
fmax_withindex(float *in, size_t size, float *max, size_t *index);

void
dmin_withindex(double *in, size_t size, double *min, size_t *index);

void
fmin_withindex(float *in, size_t size, float *min, size_t *index);






float
floatsum(float *in, size_t size);

float
floatsumsquared(float *in, size_t size);

float
floatsummask(float *in, unsigned char *mask, size_t size, size_t *nsize);

float
floatsummaskl(float *in, long *mask,
        size_t size, size_t *nsize);

float
floatsumsquaredmask(float *in, unsigned char *mask, size_t size,
		    size_t *nsize);

float
floatsumsquaredmaskl(float *in, long *mask, size_t size, size_t *nsize);






void
fave(float *in, size_t size, float *ave, unsigned char *mask);

void
favel(float *in, size_t size, float *ave, long *mask);

void
favestd(float *in, size_t size, float *ave, float *std,
    unsigned char *mask);

void
favestdl(float *in, size_t size, float *ave, float *std,
    long *mask);
void
floatavestdmaskbyt0inregion(float *in, unsigned char *byt,
			    unsigned char *mask, size_t startind,
			    size_t s0, size_t s1, size_t is1,
			    float *ave, float *std);
void
floatavestdmaskbyt0inregionsclip(float *in, unsigned char *byt,
				 unsigned char *mask, size_t startind,
				 size_t s0, size_t s1, size_t is1,
				 size_t numback, float *ave, float *std);



void
setbins(float *sorted, size_t size, size_t numbins, float min,
	float max, int binonzero, float quant, float **obins);

void
histogram(float *sorted, size_t size, float *bins, size_t numbins,
	  int normhist, int maxhistone);

void
cumulativefp(float *sorted, size_t size, float *bins, size_t numbins,
	     int normcfp);

void
savehist(float *sorted, size_t size, size_t numbins,
	 char *filename, char *histname, size_t id);






void
removeoutliers_flatcdf(float *arr, size_t *outsize);






int
sigmaclip_converge(float *array, int o1_n0, size_t num_elem,
		   float sigma_multiple, float accuracy,
		   float *outave, float *outmed, float *outstd);

int
sigmaclip_certainnum(float *array, int o1_n0, size_t num_elem,
		     float sigma_multiple, size_t numtimes,
		     float *outave, float *outmed, float *outstd);












size_t
indexfromquantile(size_t size, float quant);

void
valuefromquantile(float *data, size_t size, float quant, float *quantflux,
		  unsigned char *mask);

void
multivaluefromquantile(float *data, size_t size, float *quants,
		       float *quantfluxs, size_t numquants,
		       unsigned char *mask);

void
valuefromquantile_nocopy(float *data, size_t size, float quant,
			 float *quantflux, unsigned char *mask);

void
valuefromquantile_sortednocopy(float *data, size_t size, float quant,
			       float *quantflux, unsigned char *mask);

void
quantilefromvalue(float *data, size_t size, float *quant, float quantflux,
		  unsigned char *mask);


void
quantilefromvalue_nocopy(float *data, size_t size, float *quant,
			 float quantflux, unsigned char *mask);

void
quantilefromvalue_sorted(float *sorted, size_t size,
			 float *quant, float quantflux);

#endif
