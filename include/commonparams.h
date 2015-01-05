/*********************************************************************
Common parameters for all the utilities.
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
#ifndef COMMONPARAMS_H
#define COMMONPARAMS_H


/* The structure keeping all the common parameters in AstrUtils. This
   could be part of commonargs.h, but since main.h in all those
   programs will need this structure, there will be problems about the
   unused `commonargp` structure.*/
struct commonparams
{
  char        *spack;  /* Subpackage name.                              */

  /* Input/Output: */
  char       *output;  /* Directory containg output.                    */
  char          *hdu;  /* Image extension.                              */
  int     dontdelete;  /* ==1: Don't delete existing.                   */
  int  removedirinfo;  /* ==1: Remove directory information.            */

  /* Operating modes: */
  int           verb;  /* ==1: report steps. ==0 don't.                 */
  int    printparams;  /* Only print the used values.                   */
  int     setdirconf;  /* ==1: Set the current directory default values.*/
  int     setusrconf;  /* ==1: Set the user default values.             */
  size_t  numthreads;  /* Number of threads to use.                     */

  /* Check: */
  int  numthreadsset;  /* If the number of threads are set.             */
  int         hduset;  /* If the input image extension is set.          */
  int      outputset;  /* If the output is set.                         */
};

#endif
