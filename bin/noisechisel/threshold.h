/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef THRESHOLD_H
#define THRESHOLD_H

#define THRESHOLD_NO_ERODE_VALUE 2


enum threshold_type
  {
    THRESHOLD_QUANTILES,
    THRESHOLD_SKY_STD,
  };


void
threshold_apply(struct noisechiselparams *p, float *value1, float *value2,
                int type);

void
threshold_write_sn_table(struct noisechiselparams *p, gal_data_t *sntable,
                         gal_data_t *snind, char *filename,
                         gal_list_str_t *comments, char *extname);

void
threshold_interp_smooth(struct noisechiselparams *p, gal_data_t **first,
                        gal_data_t **second, gal_data_t **third,
                        char *filename);

void
threshold_no_outlier(struct noisechiselparams *p, gal_data_t *first,
                     gal_data_t *second, gal_data_t *third,
                     char *filename);

void
threshold_quantile_find_apply(struct noisechiselparams *p);


#endif
