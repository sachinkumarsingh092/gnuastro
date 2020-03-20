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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <error.h>





/**********************************************************************/
/****************      Functions to parse strings     *****************/
/**********************************************************************/


/* Parse the input string consisting of numbers separated by given
 * delimiter into an array. */
int
gal_units_extract_decimal (char *convert, const char *delimiter,
                           double *args, size_t n)
{
  char *copy, *token, *end;
  size_t i = 0;
  /* Create a copy of the string to be parsed */
  copy = strdup (convert);

  do
    {
      /* Check if the required number of arguments are passed */
      if (i == n + 1)
        {
          free (copy);
          error (0, 0, "'%s' exceeds maximum number of "
                       "arguments\n", convert);
          return 0;
        }


      /* Extract the substring till the next delimiter */
      token = strtok (i == 0 ? copy : NULL, delimiter);

      if (token)
        {
          /* Parse extracted string as a number */
          args[i++] = strtod (token, &end);
          /* Check if strtod fails to parse any character in substring */
          if (*end && *end != delimiter)
            {
              free (copy);
              error (0, 0, "Unable to parse element %lu in "
                           "'%s'\n", i, convert);
              return 0;
            }
        }

    }
  while (token && *token);

  free (copy);
  /* Check if the number of elements parsed is unequal to numbers of
   * elements requested */
  if (i != n)
    {
      error (0, 0, "`%s' must contain %lu numbers, but has "
                   "%lu numbers\n", convert, n, i);
      return 0;
    }

  return 1;
}


















/**********************************************************************/
/****************      Convert string to decimal      *****************/
/**********************************************************************/


/* Parse the right ascension input as a string in form of hh:mm:ss to a
 * single decimal value calculated by (hh + mm / 60 + ss / 3600 ) * 15. */
double
gal_units_ra_to_decimal (char *convert)
{
  double val[3];
  double decimal = 0.0;

  /* Check whether the string is successfully parsed */
  if (gal_units_extract_decimal (convert, ":", val, 3))
    {
      /* Check whether value of hours is in within limits */
      if (val[0] < 0.0 || val[0] > 24.0)
        {
          error (0, 0, "Error: Value of hours should be "
                       "between be 0 and 24. %s\n", convert);
          return NAN;
        }
      /* Add hours to the decimal value */
      decimal += val[0];

      /* Check whether value of minutes is in within limits */
      if (val[1] < 0.0 || val[1] > 60.0)
        {
          error (0, 0, "Error: Value of minutes should be "
                       "between be 0 and 60. %s\n", convert);
          return NAN;
        }
      /* Convert minutes to hours and add to the decimal value */
      decimal += val[1] / 60;

      /* Check whether value of seconds is in within limits */
      if (val[2] < 0.0 || val[2] > 60.0)
        {
          error (0, 0, "Error: Value of seconds should be "
                       "between be 0 and 60. %s\n", convert);
          return NAN;
        }
      /* Convert seconds to hours and add to the decimal value */
      decimal += val[2] / 3600;

      /* Convert value to decimal */
      decimal *= 15.0;

      return decimal;
    }
  else
    return NAN;
}

/* Parse the declination input as a string in form of dd:mm:ss to a decimal
 * calculated by (dd + mm / 60 + ss / 3600 ). */
double
gal_units_dec_to_decimal (char *convert)
{
  double val[3];
  double decimal = 0.0;
  int sign = 1;

  if (gal_units_extract_decimal (convert, ":", val, 3))
    {
      /* Check whether value of decimal is in within limits */
      if (val[0] < -90.0 || val[0] > 90.0)
        {
          error (0, 0, "Error: Value of decimal should be "
                       "between be -90 and 90. %s\n", convert);
          return NAN;
        }

      /* If declination is negative, the first value in the array will be
       * negative and all other values will be positive. In that case, we
       * set sign equal to -1. Therefore, we multiply the first value by
       * sign to make it positive. The final answer is again multiplied by
       * sign to make its sign same as original. */
      if (val[0] < 0.0)
        sign = -1;

      /* Add decimal to the decimal value after making it positive */
      decimal += val[0] * sign;

      /* Check whether value of arc-minutes is in within limits */
      if (val[1] < 0.0 || val[1] > 60.0)
        {
          error (0, 0, "Error: Value of arc-minutes should be "
                       "between be 0 and 60. %s\n", convert);
          return NAN;
        }
      /* Convert arc-minutes to decimal and add to the decimal value */
      decimal += val[1] / 60;

      /* Check whether value of arc-seconds is in within limits */
      if (val[2] < 0.0 || val[2] > 60.0)
        {
          error (0, 0, "Error: Value of arc-seconds should be "
                       "between be 0 and 60. %s\n", convert);
          return NAN;
        }
      /* Convert arc-seconds to decimal and add to the decimal value */
      decimal += val[2] / 3600;

      /* Make the sign of the decimal value same as input */
      decimal *= sign;

      return decimal;
    }
  else
    return NAN;
}

















/**********************************************************************/
/****************      Convert decimal to string      *****************/
/**********************************************************************/


/* Parse the right ascension input as a decimal to a string in form of
 * hh:mm:ss.ss . */
char *
gal_units_decimal_to_ra (double decimal)
{
  int hours = 0, minutes = 0;
  float seconds = 0.0; /* For sub-second accuracy */
  /* Allocate string of length 15 which is large enough for string of
   * format hh:mm:ss.ss and sign */
  char *ra = ( char * ) malloc (sizeof (char) * 15);

  /* Check if decimal value is within bounds otherwise return error */
  if (decimal < 0 || decimal > 360)
    {
      error (0, 0, "Error: Value of decimal should be "
                   "between be 0 and 360. %.10f\n", decimal);
      return 0;
    }

  /* Divide decimal value by 15 */
  decimal /= 15.0;

  /* Extract integer part of decimal value to obtain hours */
  hours = ( int ) (decimal);

  /* Subtract hours from decimal and multiply remaining value by 60 to
   * obtain minutes. */
  minutes = ( int ) ((decimal - hours) * 60);

  /* Subtract hours and minutes from decimal and multiply remaining value
   * by 3600 to obtain seconds. */
  seconds = (decimal - hours - minutes / 60.0) * 3600;

  /* Format the extracted hours, minutes and seconds as a string with
   * leading zeros if required, in hh:mm:ss format */
  snprintf (ra,
            sizeof (char) * 15, "%02d:%02d:%05.2f", hours, minutes, seconds);

  return ra;
}

/* Parse the declination input as a decimal to a string in form of dd:mm:ss*/
char *
gal_units_decimal_to_dec (double decimal)
{
  int degrees = 0, arc_minutes = 0;
  float arc_seconds = 0.0;
  /* Allocate string of length 15 which is large enough for string of
   * format hh:mm:ss.ss and sign */
  char *dec = ( char * ) malloc (sizeof (char) * 15);
  int sign = 1;

  /* Check if decimal value is within bounds otherwise return error */
  if (decimal < -90 || decimal > 90)
    {
      error (0, 0, "Error: Value of decimal should be "
                   "between be -90 and 90. %.10f\n", decimal);
      return 0;
    }

  /* If declination is negative, we set sign equal to -1. We multiply the
   * decimal by to make sure it is positive. We then extract degrees,
   * arc-minutes and arc-seconds from the decimal. Finally, we add a minus
   * sign in beginning of string if input was negative. */
  if (decimal < 0)
    sign = -1;

  /* Multiply decimal by sign to make its positive. */
  decimal *= sign;

  /* Extract integer part of decimal value to obtain degrees. */
  degrees = ( int ) (decimal);

  /* Subtract degrees from decimal and multiply remaining value by 60 to
   * obtain arc-minutes. */
  arc_minutes = ( int ) ((decimal - degrees) * 60);

  /* Subtract degrees and arc-minutes from decimal and multiply remaining
   * value by 3600 to obtain arc-seconds. */
  arc_seconds = (decimal - degrees - arc_minutes / 60.0) * 3600;

  /* Format the extracted degrees, arc-minutes and arc-seconds as a string
   * with leading zeros if required, in hh:mm:ss format with correct sign */
  if (sign < 0)
    snprintf (dec, sizeof (char)
                   * 15, "-%02d:%02d:%05.2f", degrees, arc_minutes, arc_seconds);
  else
    snprintf (dec, sizeof (char)
                   * 15, "%02d:%02d:%05.2f", degrees, arc_minutes, arc_seconds);

  return dec;
}