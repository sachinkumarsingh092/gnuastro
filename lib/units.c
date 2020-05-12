/*********************************************************************
Units -- Convert data from one unit to other.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Kartik Ohri <kartikohri13@gmail.com>
Contributing author(s):
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Copyright (C) 2020, Free Software Foundation, Inc.

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

#include <math.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/type.h>
#include <gnuastro/pointer.h>



/**********************************************************************/
/****************      Functions to parse strings     *****************/
/**********************************************************************/
/* Parse the input string consisting of numbers separated by given
   delimiter into an array. */
int
gal_units_extract_decimal(char *convert, const char *delimiter,
                          double *args, size_t n)
{
  size_t i = 0;
  char *copy, *token, *end;

  /* Create a copy of the string to be parsed and parse it. This is because
     it will be modified during the parsing. */
  copy=strdup(convert);
  do
    {
      /* Check if the required number of arguments are passed */
      if(i==n+1)
        {
          free(copy);
          error(0, 0, "%s: input '%s' exceeds maximum number of arguments "
                "(%zu)", __func__, convert, n);
          return 0;
        }

      /* Extract the substring till the next delimiter */
      token=strtok(i==0?copy:NULL, delimiter);
      if(token)
        {
          /* Parse extracted string as a number, and check if it worked. */
          args[i++] = strtod (token, &end);
          if (*end && *end != *delimiter)
            {
              free(copy);
              error(0, 0, "%s: unable to parse element %zu in '%s'\n",
                    __func__, i, convert);
              return 0;
            }
        }
    }
  while(token && *token);
  free (copy);

  /* Check if the number of elements parsed. */
  if (i != n)
    {
      error (0, 0, "%s: input '%s' must contain %lu numbers, but has "
             "%lu numbers\n", __func__, convert, n, i);
      return 0;
    }

  /* Numbers are written, return successfully. */
  return 1;
}


















/**********************************************************************/
/****************      Convert string to decimal      *****************/
/**********************************************************************/

/* Parse the right ascension input as a string in form of hh:mm:ss to a
 * single decimal value calculated by (hh + mm / 60 + ss / 3600 ) * 15. */
double
gal_units_ra_to_degree(char *convert)
{
  double val[3];
  double decimal=0.0;

  /* Check whether the string is successfully parsed */
  if(gal_units_extract_decimal (convert, ":", val, 3))
    {
      /* Check whether the first value is in within limits, and add it. */
      if(val[0]<0.0 || val[0]>24.0)
        {
          error(0, 0, "%s: value of first decimal (%g) in '%s' should be "
                "between 0 and 24", __func__, val[0], convert);
          return NAN;
        }
      decimal += val[0];

      /* Check whether value of minutes is in within limits, and add it. */
      if(val[1]<0.0 || val[1]>60.0)
        {
          error(0, 0, "%s: value of second decimal (%g) in '%s' should be "
                "between 0 and 60", __func__, val[0], convert);
          return NAN;
        }
      decimal += val[1] / 60;

      /* Check whether value of seconds is in within limits, and add it. */
      if(val[2]<0.0 || val[2]>60.0)
        {
          error(0, 0, "%s: value of third decimal (%g) in '%s' should be "
                "between 0 and 60", __func__, val[0], convert);
          return NAN;
        }
      decimal += val[2] / 3600;

      /* Convert value to degrees and return. */
      decimal *= 15.0;
      return decimal;
    }
  else
    {
      error(0, 0, "%s: input '%s' couldn't be parsed", __func__, convert);
      return NAN;
    }

  /* Control shouldn't reach this point. If it does, its a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NAN;
}





/* Parse the declination input as a string in form of dd:mm:ss to a decimal
 * calculated by (dd + mm / 60 + ss / 3600 ). */
double
gal_units_dec_to_degree (char *convert)
{
  int sign;
  double val[3], decimal=0.0;

  /* Parse the values in the input string. */
  if(gal_units_extract_decimal (convert, ":", val, 3))
    {
      /* Check whether the first value is in within limits. */
      if(val[0]<-90.0 || val[0]>90.0)
        {
          error(0, 0, "%s: value of first decimal (%g) in '%s' should be "
                "between -90 and 90", __func__, val[0], convert);
          return NAN;
        }

      /* If declination is negative, the first value in the array will be
         negative and all other values will be positive. In that case, we
         set sign equal to -1. Therefore, we multiply the first value by
         sign to make it positive. The final answer is again multiplied by
         sign to make its sign same as original. */
      sign = val[0]<0.0 ? -1 : 1;
      decimal += val[0] * sign;

      /* Check whether value of arc-minutes is in within limits. */
      if(val[1]<0.0 || val[1]>60.0)
        {
          error(0, 0, "%s: value of second decimal (%g) in '%s' should be "
                "between 0 and 60", __func__, val[1], convert);
          return NAN;
        }
      /* Convert arc-minutes to decimal and add to the decimal value */
      decimal += val[1] / 60;

      /* Check whether value of arc-seconds is in within limits */
      if (val[2] < 0.0 || val[2] > 60.0)
        {
          error(0, 0, "%s: value of third decimal (%g) in '%s' should be "
                "between 0 and 60", __func__, val[2], convert);
          return NAN;
        }

      /* Convert arc-seconds to decimal and add to the decimal value */
      decimal += val[2] / 3600;

      /* Make the sign of the decimal value same as input and return. */
      decimal *= sign;
      return decimal;
    }
  else
    {
      error(0, 0, "%s: input '%s' couldn't be parsed", __func__, convert);
      return NAN;
    }

  /* Control shouldn't reach this point. If it does, its a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NAN;
}




















/**********************************************************************/
/****************      Convert decimal to string      *****************/
/**********************************************************************/

/* Max-length of output string. */
#define UNITS_RADECSTR_MAXLENGTH 50

/* Parse the right ascension input as a decimal to a string in form of
   hh:mm:ss.ss . */
char *
gal_units_degree_to_ra(double decimal)
{
  size_t nchars;
  int hours=0, minutes=0;
  float seconds=0.0; /* For sub-second accuracy */

  /* Allocate string of length 15 which is large enough for string of
     format hh:mm:ss.ss and sign */
  char *ra=gal_pointer_allocate(GAL_TYPE_UINT8, UNITS_RADECSTR_MAXLENGTH,
                                0, __func__, "ra");

  /* Check if decimal value is within bounds otherwise return error */
  if (decimal<0 || decimal>360)
    {
      error (0, 0, "%s: value of decimal should be between be 0 and 360, "
             "but is %g\n", __func__, decimal);
      return NULL;
    }

  /* Divide decimal value by 15 and extract integer part of decimal value
     to obtain hours */
  decimal /= 15.0;
  hours = (int)decimal;

  /* Subtract hours from decimal and multiply remaining value by 60 to
     obtain minutes. */
  minutes = (int)((decimal - hours) * 60);

  /* Subtract hours and minutes from decimal and multiply remaining value
     by 3600 to obtain seconds. */
  seconds = (decimal - hours - minutes / 60.0) * 3600;

  /* Format the extracted hours, minutes and seconds as a string with
     leading zeros if required, in hh:mm:ss format */
  nchars = snprintf(ra, UNITS_RADECSTR_MAXLENGTH-1, "%02d:%02d:%g",
                    hours, minutes, seconds);
  if(nchars>UNITS_RADECSTR_MAXLENGTH)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to address "
          "the problem. The output string has an unreasonable length of "
          "%zu characters", __func__, PACKAGE_BUGREPORT, nchars);

  /* Return the final string. */
  return ra;
}





/* Parse the declination input as a decimal to a string in form of dd:mm:ss*/
char *
gal_units_degree_to_dec(double decimal)
{
  size_t nchars;
  float arc_seconds=0.0;
  int sign, degrees=0, arc_minutes=0;

  /* Allocate string of fixed length which is large enough for string of
   * format hh:mm:ss.ss and sign */
  char *dec=gal_pointer_allocate(GAL_TYPE_UINT8, UNITS_RADECSTR_MAXLENGTH,
                                 0, __func__, "ra");

  /* Check if decimal value is within bounds otherwise return error */
  if(decimal<-90 || decimal>90)
    {
      error (0, 0, "%s: value of decimal should be between be -90 and 90, "
             "but is %g\n", __func__, decimal);
      return NULL;
    }

  /* If declination is negative, we set 'sign' equal to -1. We multiply the
     decimal by to make sure it is positive. We then extract degrees,
     arc-minutes and arc-seconds from the decimal. Finally, we add a minus
     sign in beginning of string if input was negative. */
  sign = decimal<0.0 ? -1 : 1;
  decimal *= sign;

  /* Extract integer part of decimal value to obtain degrees. */
  degrees=(int)decimal;

  /* Subtract degrees from decimal and multiply remaining value by 60 to
     obtain arc-minutes. */
  arc_minutes=(int)( (decimal - degrees) * 60 );

  /* Subtract degrees and arc-minutes from decimal and multiply remaining
     value by 3600 to obtain arc-seconds. */
  arc_seconds = (decimal - degrees - arc_minutes / 60.0) * 3600;

  /* Format the extracted degrees, arc-minutes and arc-seconds as a string
     with leading zeros if required, in hh:mm:ss format with correct
     sign. */
  nchars = snprintf(dec, UNITS_RADECSTR_MAXLENGTH-1, "%s%02d:%02d:%g",
                    sign<0?"-":"+", degrees, arc_minutes, arc_seconds);
  if(nchars>UNITS_RADECSTR_MAXLENGTH)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to address "
          "the problem. The output string has an unreasonable length of "
          "%zu characters", __func__, PACKAGE_BUGREPORT, nchars);

  /* Return the final string. */
  return dec;
}
