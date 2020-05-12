/*********************************************************************
Spectral lines.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2019-2020, Free Software Foundation, Inc.

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
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/speclines.h>


/*********************************************************************/
/*************        Internal names and codes         ***************/
/*********************************************************************/
/* Return line's name as literal string. */
char *
gal_speclines_line_name(int linecode)
{
  switch(linecode)
    {
    case GAL_SPECLINES_SIIRED:           return GAL_SPECLINES_NAME_SIIRED;
    case GAL_SPECLINES_SII:              return GAL_SPECLINES_NAME_SII;
    case GAL_SPECLINES_SIIBLUE:          return GAL_SPECLINES_NAME_SIIBLUE;
    case GAL_SPECLINES_NIIRED:           return GAL_SPECLINES_NAME_NIIRED;
    case GAL_SPECLINES_NII:              return GAL_SPECLINES_NAME_NII;
    case GAL_SPECLINES_HALPHA:           return GAL_SPECLINES_NAME_HALPHA;
    case GAL_SPECLINES_NIIBLUE:          return GAL_SPECLINES_NAME_NIIBLUE;
    case GAL_SPECLINES_OIIIRED:          return GAL_SPECLINES_NAME_OIIIRED;
    case GAL_SPECLINES_OIII:             return GAL_SPECLINES_NAME_OIII;
    case GAL_SPECLINES_OIIIBLUE:         return GAL_SPECLINES_NAME_OIIIBLUE;
    case GAL_SPECLINES_HBETA:            return GAL_SPECLINES_NAME_HBETA;
    case GAL_SPECLINES_HEIIRED:          return GAL_SPECLINES_NAME_HEIIRED;
    case GAL_SPECLINES_HGAMMA:           return GAL_SPECLINES_NAME_HGAMMA;
    case GAL_SPECLINES_HDELTA:           return GAL_SPECLINES_NAME_HDELTA;
    case GAL_SPECLINES_HEPSILON:         return GAL_SPECLINES_NAME_HEPSILON;
    case GAL_SPECLINES_NEIII:            return GAL_SPECLINES_NAME_NEIII;
    case GAL_SPECLINES_OIIRED:           return GAL_SPECLINES_NAME_OIIRED;
    case GAL_SPECLINES_OII:              return GAL_SPECLINES_NAME_OII;
    case GAL_SPECLINES_OIIBLUE:          return GAL_SPECLINES_NAME_OIIBLUE;
    case GAL_SPECLINES_BLIMIT:           return GAL_SPECLINES_NAME_BLIMIT;
    case GAL_SPECLINES_MGIIRED:          return GAL_SPECLINES_NAME_MGIIRED;
    case GAL_SPECLINES_MGII:             return GAL_SPECLINES_NAME_MGII;
    case GAL_SPECLINES_MGIIBLUE:         return GAL_SPECLINES_NAME_MGIIBLUE;
    case GAL_SPECLINES_CIIIRED:          return GAL_SPECLINES_NAME_CIIIRED;
    case GAL_SPECLINES_CIII:             return GAL_SPECLINES_NAME_CIII;
    case GAL_SPECLINES_CIIIBLUE:         return GAL_SPECLINES_NAME_CIIIBLUE;
    case GAL_SPECLINES_HEIIBLUE:         return GAL_SPECLINES_NAME_HEIIBLUE;
    case GAL_SPECLINES_LYALPHA:          return GAL_SPECLINES_NAME_LYALPHA;
    case GAL_SPECLINES_LYLIMIT:          return GAL_SPECLINES_NAME_LYLIMIT;
    default:                           return NULL;
    }
  return NULL;
}





/* Return the code of the given line name. */
int
gal_speclines_line_code(char *name)
{
  if( !strcmp(name, GAL_SPECLINES_NAME_SIIRED) )
    return GAL_SPECLINES_SIIRED;
  else if( !strcmp(name, GAL_SPECLINES_NAME_SII) )
    return GAL_SPECLINES_SII;
  else if( !strcmp(name, GAL_SPECLINES_NAME_SIIBLUE) )
    return GAL_SPECLINES_SIIBLUE;
  if( !strcmp(name, GAL_SPECLINES_NAME_NIIRED) )
    return GAL_SPECLINES_NIIRED;
  else if( !strcmp(name, GAL_SPECLINES_NAME_NII) )
    return GAL_SPECLINES_NII;
  else if( !strcmp(name, GAL_SPECLINES_NAME_HALPHA) )
    return GAL_SPECLINES_HALPHA;
  else if( !strcmp(name, GAL_SPECLINES_NAME_NIIBLUE) )
    return GAL_SPECLINES_NIIBLUE;
  else if( !strcmp(name, GAL_SPECLINES_NAME_OIIIRED) )
    return GAL_SPECLINES_OIIIRED;
  else if( !strcmp(name, GAL_SPECLINES_NAME_OIII) )
    return GAL_SPECLINES_OIII;
  else if( !strcmp(name, GAL_SPECLINES_NAME_OIIIBLUE) )
    return GAL_SPECLINES_OIIIBLUE;
  else if( !strcmp(name, GAL_SPECLINES_NAME_HBETA) )
    return GAL_SPECLINES_HBETA;
  else if( !strcmp(name, GAL_SPECLINES_NAME_HEIIRED) )
    return GAL_SPECLINES_HEIIRED;
  else if( !strcmp(name, GAL_SPECLINES_NAME_HGAMMA) )
    return GAL_SPECLINES_HGAMMA;
  else if( !strcmp(name, GAL_SPECLINES_NAME_HDELTA) )
    return GAL_SPECLINES_HDELTA;
  else if( !strcmp(name, GAL_SPECLINES_NAME_HEPSILON) )
    return GAL_SPECLINES_HEPSILON;
  else if( !strcmp(name, GAL_SPECLINES_NAME_NEIII) )
    return GAL_SPECLINES_NEIII;
  else if( !strcmp(name, GAL_SPECLINES_NAME_OIIRED) )
    return GAL_SPECLINES_OIIRED;
  else if( !strcmp(name, GAL_SPECLINES_NAME_OII) )
    return GAL_SPECLINES_OII;
  else if( !strcmp(name, GAL_SPECLINES_NAME_OIIBLUE) )
    return GAL_SPECLINES_OIIBLUE;
  else if( !strcmp(name, GAL_SPECLINES_NAME_BLIMIT) )
    return GAL_SPECLINES_BLIMIT;
  else if( !strcmp(name, GAL_SPECLINES_NAME_MGIIRED) )
    return GAL_SPECLINES_MGIIRED;
  else if( !strcmp(name, GAL_SPECLINES_NAME_MGII) )
    return GAL_SPECLINES_MGII;
  else if( !strcmp(name, GAL_SPECLINES_NAME_MGIIBLUE) )
    return GAL_SPECLINES_MGIIBLUE;
  else if( !strcmp(name, GAL_SPECLINES_NAME_CIIIRED) )
    return GAL_SPECLINES_CIIIRED;
  else if( !strcmp(name, GAL_SPECLINES_NAME_CIII) )
    return GAL_SPECLINES_CIII;
  else if( !strcmp(name, GAL_SPECLINES_NAME_CIIIBLUE) )
    return GAL_SPECLINES_CIIIBLUE;
  else if( !strcmp(name, GAL_SPECLINES_NAME_HEIIBLUE) )
    return GAL_SPECLINES_HEIIBLUE;
  else if( !strcmp(name, GAL_SPECLINES_NAME_LYALPHA) )
    return GAL_SPECLINES_LYALPHA;
  else if( !strcmp(name, GAL_SPECLINES_NAME_LYLIMIT) )
    return GAL_SPECLINES_LYLIMIT;
  else return GAL_SPECLINES_INVALID;
  return GAL_SPECLINES_INVALID;
}





/* Return the wavelength (in Angstroms) of given line. */
double
gal_speclines_line_angstrom(int linecode)
{
  switch(linecode)
    {
    case GAL_SPECLINES_SIIRED:        return GAL_SPECLINES_ANGSTROM_SIIRED;
    case GAL_SPECLINES_SII:           return GAL_SPECLINES_ANGSTROM_SII;
    case GAL_SPECLINES_SIIBLUE:       return GAL_SPECLINES_ANGSTROM_SIIBLUE;
    case GAL_SPECLINES_NIIRED:        return GAL_SPECLINES_ANGSTROM_NIIRED;
    case GAL_SPECLINES_NII:           return GAL_SPECLINES_ANGSTROM_NII;
    case GAL_SPECLINES_HALPHA:        return GAL_SPECLINES_ANGSTROM_HALPHA;
    case GAL_SPECLINES_NIIBLUE:       return GAL_SPECLINES_ANGSTROM_NIIBLUE;
    case GAL_SPECLINES_OIIIRED:       return GAL_SPECLINES_ANGSTROM_OIIIRED;
    case GAL_SPECLINES_OIII:          return GAL_SPECLINES_ANGSTROM_OIII;
    case GAL_SPECLINES_OIIIBLUE:      return GAL_SPECLINES_ANGSTROM_OIIIBLUE;
    case GAL_SPECLINES_HBETA:         return GAL_SPECLINES_ANGSTROM_HBETA;
    case GAL_SPECLINES_HEIIRED:       return GAL_SPECLINES_ANGSTROM_HEIIRED;
    case GAL_SPECLINES_HGAMMA:        return GAL_SPECLINES_ANGSTROM_HGAMMA;
    case GAL_SPECLINES_HDELTA:        return GAL_SPECLINES_ANGSTROM_HDELTA;
    case GAL_SPECLINES_HEPSILON:      return GAL_SPECLINES_ANGSTROM_HEPSILON;
    case GAL_SPECLINES_NEIII:         return GAL_SPECLINES_ANGSTROM_NEIII;
    case GAL_SPECLINES_OIIRED:        return GAL_SPECLINES_ANGSTROM_OIIRED;
    case GAL_SPECLINES_OII:           return GAL_SPECLINES_ANGSTROM_OII;
    case GAL_SPECLINES_OIIBLUE:       return GAL_SPECLINES_ANGSTROM_OIIBLUE;
    case GAL_SPECLINES_BLIMIT:        return GAL_SPECLINES_ANGSTROM_BLIMIT;
    case GAL_SPECLINES_MGIIRED:       return GAL_SPECLINES_ANGSTROM_MGIIRED;
    case GAL_SPECLINES_MGII:          return GAL_SPECLINES_ANGSTROM_MGII;
    case GAL_SPECLINES_MGIIBLUE:      return GAL_SPECLINES_ANGSTROM_MGIIBLUE;
    case GAL_SPECLINES_CIIIRED:       return GAL_SPECLINES_ANGSTROM_CIIIRED;
    case GAL_SPECLINES_CIII:          return GAL_SPECLINES_ANGSTROM_CIII;
    case GAL_SPECLINES_CIIIBLUE:      return GAL_SPECLINES_ANGSTROM_CIIIBLUE;
    case GAL_SPECLINES_HEIIBLUE:      return GAL_SPECLINES_ANGSTROM_HEIIBLUE;
    case GAL_SPECLINES_LYALPHA:       return GAL_SPECLINES_ANGSTROM_LYALPHA;
    case GAL_SPECLINES_LYLIMIT:       return GAL_SPECLINES_ANGSTROM_LYLIMIT;
    default:
      error(EXIT_FAILURE, 0, "%s: '%d' not recognized line identifier",
            __func__, linecode);
    }
  return NAN;
}




















/*********************************************************************/
/*************             Redshifted lines            ***************/
/*********************************************************************/
double
gal_speclines_line_redshift(double obsline, double restline)
{
  return (obsline/restline)-1;
}





double
gal_speclines_line_redshift_code(double obsline, int linecode)
{
  double restline=gal_speclines_line_angstrom(linecode);
  return (obsline/restline)-1;
}
