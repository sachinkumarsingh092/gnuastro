/*********************************************************************
wcsdistortion -- Manipulation of distortions in WCS.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Sachin Kumar Singh <sachinkumarsingh092@gmail.com>
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

#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <wcslib/wcslib.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>

#include <gnuastro-internal/wcsdistortion.h>





/* Internally used macro(s) to help in the processing */
#define wcsdistortion_max(a,b) \
   ({ __typeof__ (a) _a = (a);  \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/* Declarations to avoid ordering problems. */
static void
wcsdistortion_calc_tpveq(struct wcsprm *wcs, double cd[2][2],
                         double tpvu[8][8], double tpvv[8][8]);
static double
wcsdistortion_calcsip(size_t axis, size_t m, size_t n, double tpvu[8][8],
                      double tpvv[8][8]);







/**************************************************************/
/**********          Reading utilities             ************/
/**************************************************************/
/* Extract the required pv parameters from wcs. */
static void
wcsdistortion_get_tpvparams(struct wcsprm *wcs, double cd[2][2],
                            double *pv1, double *pv2)
{
  size_t pv_m=0;
  size_t i, j, k,index=0;

  /* Make sure a WCS structure is actually given. */
  if(wcs==NULL)
    error(EXIT_FAILURE, 0, "%s: input WCS structure is NULL", __func__);

  for(i=0,k=0; i<2; ++i)
    for(j=0; j<2; ++j)
      {
        /* If a value is present store it, else store 0.*/
        if((wcs->cd)[i] != 0)   cd[i][j]=(wcs->cd)[k++];
        else                    cd[i][j]=0;

        /* For a check:
        printf("CD%ld_%ld\t%.10lf\n", i+1, j+1, cd[i][j]);
        */
      }

  for(j=0; j < wcs->npvmax; ++j)
    {
      if(wcs->pv[j].i == 1)
        {
          /*pv_m is used to check the index in the header.*/
          pv_m=wcs->pv[j].m;

          /* `index` is the index of the pv* array.*/
          index = pv_m;

          if( wcs->pv[pv_m].value != 0 && index == pv_m )
            pv1[index]=wcs->pv[j].value;
          else
            pv1[index]=0;

          /* For a check
          printf("PV1_%d\t%.10f\n", pv_m, pv1[k]);
          */
        }
      else if(wcs->pv[j].i == 2)
        {
          /*pv_m is used to check the index in the header.*/
          pv_m=wcs->pv[j].m;

          /* `index` is the index of the pv* array.*/
          index = pv_m;

          if( wcs->pv[pv_m].value != 0 && index == pv_m )
            pv2[index]=wcs->pv[j].value;
          else
            pv2[index]=0;

          /* For a check
          printf("PV2_%d\t%.10f\n", pv_m, pv2[k]);
          */
        }
      else
        {
          error(EXIT_FAILURE, 0, "%s: No such axis present!", __func__);
          break;
        }
    }

  /* To check a full array:
     for(i = 0; i < 40; ++i)
     printf("PV2_%ld\t%.10f\n", i, pv2[i]);
  */

}






/* Extract the required sip parameters from wcs. */
static void
wcsdistortion_get_sipparams(struct wcsprm *wcs, double cd[2][2],
                            double a_coeff[5][5], double b_coeff[5][5])
{
  const char *cp;
  double *temp_cd=NULL;
  struct dpkey  *keyp=NULL;
  struct disprm *dispre=NULL;
  size_t i=0, j=0, m=0, n=0, k=0;

  /* Make sure a WCS structure is actually given. Note that 'wcs->lin' is
     not a pointer, but the full structure. So it is always present. */
  if(wcs==NULL)
    error(EXIT_FAILURE, 0, "%s: input WCS structure is NULL", __func__);
  if(wcs->lin.dispre == NULL)
    error(EXIT_FAILURE, 0, "%s: input WCS structure's 'lin.dispre' is NULL",
          __func__);

  /* For easy reading. */
  dispre=wcs->lin.dispre;
  keyp=dispre->dp;

  /* Fill the 2-times allocated CD array (cd[][]). Note that the required
     CD matix is extracted using the `gal_wcs_wrap_matrix` as a single
     allocated array (temp_cd[]), that is then used to fill cd[][]. */
  temp_cd=gal_wcs_warp_matrix(wcs);
  for(i=0,k=0; i<2; ++i)
    for(j=0; j<2; ++j)
      {
        /* If a value is present store it, else store 0.*/
        if(temp_cd[k] != 0)   cd[i][j]=temp_cd[k++];
        else                  cd[i][j]=0;

        /* For a check:
        printf("CD%ld_%ld\t%.10E\n", i+1, j+1, cd[i][j]);
        */
      }

  /* Extract the SIP values from the strings and numbers inside the
     wcsprm. */
  for(i=0; i<dispre->ndp; ++i, ++keyp)
    {
      /* For axis1. */
      if (keyp->j == 1)
        {
          /* Ignore any missing keyvalues. */
          if ( keyp->field == NULL ) continue;

          cp = strchr(keyp->field, '.') + 1;
          if (strncmp(cp, "SIP.FWD.", 8) != 0) continue;
          cp += 8;

          sscanf(cp, "%ld_%ld", &m, &n);
          a_coeff[m][n]=dispre->dp[i].value.f;

          /*For a check.
          printf("A_%ld_%ld =\t %.10E\n", m, n, dispre->dp[i].value.f);
          */
        }

      /* For axis2. */
      else if (keyp->j == 2)
        {
          /* Ignore any missing keyvalues. */
          if ( keyp->field == NULL ) continue;

          cp = strchr(keyp->field, '.') + 1;
          if (strncmp(cp, "SIP.FWD.", 8) != 0) continue;
          cp += 8;

          sscanf(cp, "%ld_%ld", &m, &n);

          b_coeff[m][n]=dispre->dp[i].value.f;

          /*For a check.
          printf("B_%ld_%ld =\t %.10E\n", m, n, dispre->dp[i].value.f);
          */
        }
      else
        error(EXIT_FAILURE, 0, "%s: No such axis present!", __func__);
    }
}






/* Extract sip coefficients from the polynomial equation. */
static void
wcsdistortion_get_sipcoeff(struct wcsprm *wcs, size_t *a_order,
                           size_t *b_order, double a_coeff[5][5],
                           double b_coeff[5][5])
{
  size_t m, n;
  double val=0;
  size_t i, j, k;
  double cd[2][2];
  double tpvu[8][8], tpvv[8][8];

  /* Initialise the 2d matrices. */
  for(i=0; i<2; ++i) for(j=0; j<2; ++j) cd[i][j]=0;
  for(i=0; i<8; ++i) for(j=0; j<8; ++j) { tpvu[i][j]=0; tpvv[i][j]=0; }

  /* Calculate the TPV elements. */
  wcsdistortion_calc_tpveq(wcs, cd, tpvu, tpvv);

  for(i=0,k=0; i<2; ++i)
    for(j=0; j<2; ++j)
      {
        /* If a value is present store it, else store 0.*/
        if((wcs->cd)[i] != 0)   cd[i][j]=(wcs->cd)[k++];
        else                    cd[i][j]=0;

        /* For a check:
        printf("CD%ld_%ld\t%.10lf\n", i+1, j+1, cd[i][j]);
        */
      }

  for(m=0; m<=4; ++m)
    for(n=0; n<=4; ++n)
      {
        /*For axis = 1*/
        val=wcsdistortion_calcsip(1, m, n, tpvu, tpvv);
        if(val != 0)
          {
            a_coeff[m][n]=val;
            *a_order=wcsdistortion_max(*a_order, wcsdistortion_max(m,n));
          }
        else
          a_coeff[m][n]=0;

        /*For axis = 2*/
        val=wcsdistortion_calcsip(2, m, n, tpvu, tpvv);
        if(val != 0)
          {
            b_coeff[m][n]=val;
            *b_order=wcsdistortion_max(*b_order, wcsdistortion_max(m,n));
          }
        else
          b_coeff[m][n]=0;
      }

  /*For a check.
  for(m=0;m<=4;++m)
  for(n=0;n<=4;++n)
  printf("A_%ld_%ld = %.10E \t B_%ld_%ld = %.10E\n",
         m, n, a_coeff[m][n], m, n, b_coeff[m][n]);
  */
}




















/**************************************************************/
/**********       Intermidiate Equations           ************/
/**************************************************************/
/* Intermidiate equations formed during pv->sip conversions. */
static void
wcsdistortion_intermidate_tpveq(double cd[2][2], double *pv1, double *pv2,
                                double k[5][5], double l[5][5])
{
  /*Intermidate polynomials, k[i][j] and l[i][j].
    See Appendix A of the SPIE proceedings paper at
    http://web.ipac.caltech.edu/staff/shupe/reprints/SIP_to_PV_SPIE2012.pdf

    The work described in that paper is extended to 7th order.

     References and citations:
     "More flexibility in representing geometric distortion in
     astronomical images,"
     Shupe, David L.; Laher, Russ R.; Storrie-Lombardi, Lisa; Surace, Jason;
     Grillmair, Carl; Levitan, David; Sesar, Branimir, 2012, in Software and
     Cyberinfrastructure for Astronomy II.
     Proceedings of the SPIE, Volume 8451, article id. 84511M.
     */

  k[0][0] = pv1[0];
  l[0][0] = pv2[0];

  k[0][1] = ( cd[0][1]   * pv1[1]
             + cd[1][1] * pv1[2] );
  l[0][1] = ( cd[0][1]   * pv2[2]
             + cd[1][1] * pv2[1] );


  k[1][0]= ( cd[0][0]   * pv1[1]
             + cd[1][0] * pv1[2] );
  l[1][0]= ( cd[0][0]   * pv2[2]
             + cd[1][0] * pv2[1] );


  k[0][2] = ( cd[0][1]   * cd[0][1] * pv1[4]
              + cd[0][1] * cd[1][1] * pv1[5]
              + cd[1][1] * cd[1][1] * pv1[6] );
  l[0][2] = ( cd[0][1]   * cd[0][1] * pv2[6]
              + cd[0][1] * cd[1][1] * pv2[5]
              + cd[1][1] * cd[1][1] * pv2[4] );


  k[1][1] = ( 2*cd[0][0]   * cd[0][1] * pv1[4]
              + cd[0][0]   * cd[1][1] * pv1[5]
              + cd[0][1]   * cd[1][0] * pv1[5]
              + 2*cd[1][0] * cd[1][1] * pv1[6] );
  l[1][1] = ( 2*cd[0][0]   * cd[0][1] * pv2[6]
              + cd[0][0]   * cd[1][1] * pv2[5]
              + cd[0][1]   * cd[1][0] * pv2[5]
              + 2*cd[1][0] * cd[1][1] * pv2[4] );


  k[2][0] = ( cd[0][0]   * cd[0][0] * pv1[4]
              + cd[0][0] * cd[1][0] * pv1[5]
              + cd[1][0] * cd[1][0] * pv1[6] );
  l[2][0] = ( cd[0][0]   * cd[0][0] * pv2[6]
              + cd[0][0] * cd[1][0] * pv2[5]
              + cd[1][0] * cd[1][0] * pv2[4] );


  k[0][3] = ( cd[0][1]   * cd[0][1] * cd[0][1] * pv1[7]
              + cd[0][1] * cd[0][1] * cd[1][1] * pv1[8]
              + cd[0][1] * cd[1][1] * cd[1][1] * pv1[9]
              + cd[1][1] * cd[1][1] * cd[1][1] * pv1[10] );
  l[0][3] = ( cd[0][1]   * cd[0][1] * cd[0][1] * pv2[10]
              + cd[0][1] * cd[0][1] * cd[1][1] * pv2[9]
              + cd[0][1] * cd[1][1] * cd[1][1] * pv2[8]
              + cd[1][1] * cd[1][1] * cd[1][1] * pv2[7] );


  k[1][2] = ( 3*cd[0][0]   * cd[0][1] * cd[0][1] * pv1[7]
              + 2*cd[0][0] * cd[0][1] * cd[1][1] * pv1[8]
              +   cd[0][0] * cd[1][1] * cd[1][1] * pv1[9]
              +   cd[0][1] * cd[0][1] * cd[1][0] * pv1[8]
              + 2*cd[0][1] * cd[1][0] * cd[1][1] * pv1[9]
              + 3*cd[1][0] * cd[1][1] * cd[1][1] * pv1[10] );
  l[1][2] = ( 3*cd[0][0]   * cd[0][1] * cd[0][1] * pv2[10]
              + 2*cd[0][0] * cd[0][1] * cd[1][1] * pv2[9]
              +   cd[0][0] * cd[1][1] * cd[1][1] * pv2[8]
              +   cd[0][1] * cd[0][1] * cd[1][0] * pv2[9]
              + 2*cd[0][1] * cd[1][0] * cd[1][1] * pv2[8]
              + 3*cd[1][0] * cd[1][1] * cd[1][1] * pv2[7] );


  k[2][1] = ( 3*cd[0][0]   * cd[0][0] * cd[0][1] * pv1[7]
              + 2*cd[0][0] * cd[0][1] * cd[1][0] * pv1[8]
              +   cd[0][0] * cd[0][0] * cd[1][1] * pv1[8]
              +   cd[0][1] * cd[1][0] * cd[1][0] * pv1[9]
              + 2*cd[0][0] * cd[1][0] * cd[1][1] * pv1[9]
              + 3*cd[1][0] * cd[1][0] * cd[1][1] * pv1[10] );
  l[2][1] = ( 3*cd[0][0]   * cd[0][0] * cd[0][1] * pv2[10]
              + 2*cd[0][0] * cd[0][1] * cd[1][0] * pv2[9]
              +   cd[0][0] * cd[0][0] * cd[1][1] * pv2[9]
              +   cd[0][1] * cd[1][0] * cd[1][0] * pv2[8]
              + 2*cd[0][0] * cd[1][0] * cd[1][1] * pv2[8]
              + 3*cd[1][0] * cd[1][0] * cd[1][1] * pv2[7] );


  k[3][0] = ( cd[0][0]   * cd[0][0] * cd[0][0] * pv1[7]
              + cd[0][0] * cd[0][0] * cd[1][0] * pv1[8]
              + cd[0][0] * cd[1][0] * cd[1][0] * pv1[9]
              + cd[1][0] * cd[1][0] * cd[1][0] * pv1[10] );
  l[3][0] = ( cd[0][0]   * cd[0][0] * cd[0][0] * pv2[10]
              + cd[0][0] * cd[0][0] * cd[1][0] * pv2[9]
              + cd[0][0] * cd[1][0] * cd[1][0] * pv2[8]
              + cd[1][0] * cd[1][0] * cd[1][0] * pv2[7] );


  k[0][4] = ( cd[0][1]   * cd[0][1] * cd[0][1] * cd[0][1] * pv1[12]
              + cd[0][1] * cd[0][1] * cd[1][1] * cd[1][1] * pv1[13]
              + cd[0][1] * cd[0][1] * cd[1][1] * cd[1][1] * pv1[14]
              + cd[0][1] * cd[1][1] * cd[1][1] * cd[1][1] * pv1[15]
              + cd[1][1] * cd[1][1] * cd[1][1] * cd[1][1] * pv1[16] );
  l[0][4] = ( cd[0][1]   * cd[0][1] * cd[0][1] * cd[0][1] * pv2[16]
              + cd[0][1] * cd[0][1] * cd[1][1] * cd[1][1] * pv2[15]
              + cd[0][1] * cd[0][1] * cd[1][1] * cd[1][1] * pv2[14]
              + cd[0][1] * cd[1][1] * cd[1][1] * cd[1][1] * pv2[13]
              + cd[1][1] * cd[1][1] * cd[1][1] * cd[1][1] * pv2[12] );


  k[1][3] = ( 4*cd[0][0]   * cd[0][1] * cd[0][1] * cd[0][1] * pv1[12]
              + 3*cd[0][0] * cd[0][1] * cd[0][1] * cd[1][1] * pv1[13]
              + 2*cd[0][0] * cd[0][1] * cd[1][1] * cd[1][1] * pv1[14]
              +   cd[0][0] * cd[1][1] * cd[1][1] * cd[1][1] * pv1[15]
              +   cd[0][1] * cd[0][1] * cd[0][1] * cd[1][0] * pv1[13]
              + 2*cd[0][1] * cd[0][1] * cd[1][0] * cd[1][1] * pv1[14]
              + 3*cd[0][1] * cd[1][0] * cd[1][1] * cd[1][1] * pv1[15]
              + 4*cd[1][0] * cd[1][1] * cd[1][1] * cd[1][1] * pv1[16] );
  l[1][3] = ( 4*cd[0][0]   * cd[0][1] * cd[0][1] * cd[0][1] * pv2[16]
              + 3*cd[0][0] * cd[0][1] * cd[0][1] * cd[1][1] * pv2[15]
              + 2*cd[0][0] * cd[0][1] * cd[1][1] * cd[1][1] * pv2[14]
              +   cd[0][0] * cd[1][1] * cd[1][1] * cd[1][1] * pv2[13]
              +   cd[0][1] * cd[0][1] * cd[0][1] * cd[1][0] * pv2[15]
              + 2*cd[0][1] * cd[0][1] * cd[1][0] * cd[1][1] * pv2[14]
              + 3*cd[0][1] * cd[1][0] * cd[1][1] * cd[1][1] * pv2[13]
              + 4*cd[1][0] * cd[1][1] * cd[1][1] * cd[1][1] * pv2[12] );


  k[2][2] = ( 6*cd[0][0]   * cd[0][0] * cd[0][1] * cd[0][1] * pv1[12]
              + 3*cd[0][0] * cd[0][0] * cd[0][1] * cd[1][1] * pv1[13]
              +   cd[0][0] * cd[0][0] * cd[1][1] * cd[1][1] * pv1[14]
              + 3*cd[0][0] * cd[0][1] * cd[0][1] * cd[1][0] * pv1[13]
              + 4*cd[0][0] * cd[0][1] * cd[1][0] * cd[1][1] * pv1[14]
              + 3*cd[0][0] * cd[1][0] * cd[1][1] * cd[1][1] * pv1[15]
              +   cd[0][1] * cd[0][1] * cd[1][0] * cd[1][0] * pv1[14]
              + 3*cd[0][1] * cd[1][0] * cd[1][0] * cd[1][1] * pv1[15]
              + 6*cd[1][0] * cd[1][0] * cd[1][1] * cd[1][1] * pv1[16] );
  l[2][2] = ( 6*cd[0][0]   * cd[0][0] * cd[0][1] * cd[0][1] * pv2[16]
              + 3*cd[0][0] * cd[0][0] * cd[0][1] * cd[1][1] * pv2[15]
              +   cd[0][0] * cd[0][0] * cd[1][1] * cd[1][1] * pv2[14]
              + 3*cd[0][0] * cd[0][1] * cd[0][1] * cd[1][0] * pv2[15]
              + 4*cd[0][0] * cd[0][1] * cd[1][0] * cd[1][1] * pv2[14]
              + 3*cd[0][0] * cd[1][0] * cd[1][1] * cd[1][1] * pv2[15]
              +   cd[0][1] * cd[0][1] * cd[1][0] * cd[1][0] * pv2[14]
              + 3*cd[0][1] * cd[1][0] * cd[1][0] * cd[1][1] * pv2[13]
              + 6*cd[1][0] * cd[1][0] * cd[1][1] * cd[1][1] * pv2[12] );


  k[3][1] = ( 4*cd[0][0]   * cd[0][0] * cd[0][0] * cd[0][1] * pv1[12]
              +   cd[0][0] * cd[0][0] * cd[0][0] * cd[1][1] * pv1[13]
              + 3*cd[0][0] * cd[0][0] * cd[0][1] * cd[1][0] * pv1[13]
              + 2*cd[0][0] * cd[0][0] * cd[1][0] * cd[1][1] * pv1[14]
              + 2*cd[0][0] * cd[0][1] * cd[1][0] * cd[1][0] * pv1[14]
              + 3*cd[0][0] * cd[1][0] * cd[1][0] * cd[1][1] * pv1[15]
              +   cd[0][1] * cd[1][0] * cd[1][0] * cd[1][0] * pv1[15]
              + 4*cd[1][0] * cd[1][0] * cd[1][0] * cd[1][1] * pv1[16] );
  l[3][1] = ( 4*cd[0][0]   * cd[0][0] * cd[0][0] * cd[0][1] * pv2[16]
              +   cd[0][0] * cd[0][0] * cd[0][0] * cd[1][1] * pv2[15]
              + 3*cd[0][0] * cd[0][0] * cd[0][1] * cd[1][0] * pv2[15]
              + 2*cd[0][0] * cd[0][0] * cd[1][0] * cd[1][1] * pv2[14]
              + 2*cd[0][0] * cd[0][1] * cd[1][0] * cd[1][0] * pv2[14]
              + 3*cd[0][0] * cd[1][0] * cd[1][0] * cd[1][1] * pv2[13]
              +   cd[0][1] * cd[1][0] * cd[1][0] * cd[1][0] * pv2[13]
              + 4*cd[1][0] * cd[1][0] * cd[1][0] * cd[1][1] * pv2[12] );


  k[4][0] = ( cd[0][0]   * cd[0][0] * cd[0][0] * cd[0][0] * pv1[12]
              + cd[0][0] * cd[0][0] * cd[0][0] * cd[1][0] * pv1[13]
              + cd[0][0] * cd[0][0] * cd[1][0] * cd[1][0] * pv1[14]
              + cd[0][0] * cd[1][0] * cd[1][0] * cd[1][0] * pv1[15]
              + cd[1][0] * cd[1][0] * cd[1][0] * cd[1][0] * pv1[16] );
  l[4][0] = ( cd[0][0]   * cd[0][0] * cd[0][0] * cd[0][0] * pv2[16]
              + cd[0][0] * cd[0][0] * cd[0][0] * cd[1][0] * pv2[15]
              + cd[0][0] * cd[0][0] * cd[1][0] * cd[1][0] * pv2[14]
              + cd[0][0] * cd[1][0] * cd[1][0] * cd[1][0] * pv2[13]
              + cd[1][0] * cd[1][0] * cd[1][0] * cd[1][0] * pv2[12] );


  /* For a check:
  {
    size i, j;
    for(i=0; i<=4; ++i)
      for(j=0;j<=4;++j)
	{
	  printf("k%ld_%ld \t %.10E\n",   i, j, k[i][j]);
	  printf("l%ld_%ld \t %.10E\n\n", i, j, l[i][j]);
	}
  }
  */
}






/* Intermidiate equations formed during sip->pv conversions. */
static void
wcsdistortion_intermidate_sipeq(double cd[2][2], double cd_inv[2][2],
                                double a_coeff[5][5], double b_coeff[5][5],
                                double *pv1, double *pv2)
{
  /*Intemidiate equations which give the value of PVi_j
    excluding radial terms PV[i,3], PV[i,11]

    For the intermidiate PVi_j calculations, the Eqs. 1 and 2 from
    section 2 of the SPIE proceedings paper at
    http://web.ipac.caltech.edu/staff/shupe/reprints/SIP_to_PV_SPIE2012.pdf
    are used which were further modified as shown in
    https://github.com/stargaser/sip_tpv/blob/master/sip_tpv/pvsiputils.py
    to generate these equtions.

    The exact script used for generation of these equations are given here:
    https://gitlab.com/sachinkumarsingh092/gnuastro-test-files/-/blob/
    master/scripts/equations.py

    References and citations:
    "More flexibility in representing geometric distortion in
    astronomical images,"
    Shupe, David L.; Laher, Russ R.; Storrie-Lombardi, Lisa; Surace, Jason;
    Grillmair, Carl; Levitan, David; Sesar, Branimir, 2012, in Software and
    Cyberinfrastructure for Astronomy II.
    Proceedings of the SPIE, Volume 8451, article id. 84511M.
  */

  /* pvi_0 */
  pv1[0]=a_coeff[0][0]*cd[0][0] + b_coeff[0][0]*cd[0][1];
  pv2[0]=a_coeff[0][0]*cd[1][0] + b_coeff[0][0]*cd[1][1];


  /* pvi_1 */
  pv1[1] = ( a_coeff[0][1]   * cd[0][0] * cd_inv[1][0]
             + a_coeff[1][0] * cd[0][0] * cd_inv[0][0]
             + b_coeff[0][1] * cd[0][1] * cd_inv[1][0]
             + b_coeff[1][0] * cd[0][1] * cd_inv[0][0]
             +      cd[0][0] * cd_inv[0][0]
             +      cd[0][1] * cd_inv[1][0]  );

  pv2[1] = ( a_coeff[0][1]   * cd[1][0] * cd_inv[1][1]
             + a_coeff[1][0] * cd[1][0] * cd_inv[0][1]
             + b_coeff[0][1] * cd[1][1] * cd_inv[1][1]
             + b_coeff[1][0] * cd[1][1] * cd_inv[0][1]
             +      cd[1][0] * cd_inv[0][1]
             +      cd[1][1] * cd_inv[1][1]  );


  /* pvi_2 */
  pv1[2] = ( a_coeff[0][1]   * cd[0][0] * cd_inv[1][1]
             + a_coeff[1][0] * cd[0][0] * cd_inv[0][1]
             + b_coeff[0][1] * cd[0][1] * cd_inv[1][1]
             + b_coeff[1][0] * cd[0][1] * cd_inv[0][1]
             +      cd[0][0] * cd_inv[0][1]
             +      cd[0][1] * cd_inv[1][1]  );

  pv2[2] = ( a_coeff[0][1]   * cd[1][0] * cd_inv[1][0]
             + a_coeff[1][0] * cd[1][0] * cd_inv[0][0]
             + b_coeff[0][1] * cd[1][1] * cd_inv[1][0]
             + b_coeff[1][0] * cd[1][1] * cd_inv[0][0]
             +      cd[1][0] * cd_inv[0][0]
             +      cd[1][1] * cd_inv[1][0]  );


  /* pvi_4 */
  pv1[4] = ( a_coeff[0][2]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][0]
             + a_coeff[1][1] * cd[0][0] * cd_inv[0][0] * cd_inv[1][0]
             + a_coeff[2][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0]
             + b_coeff[0][2] * cd[0][1] * cd_inv[1][0] * cd_inv[1][0]
             + b_coeff[1][1] * cd[0][1] * cd_inv[0][0] * cd_inv[1][0]
             + b_coeff[2][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] );

  pv2[4] = ( a_coeff[0][2]   * cd[1][0] * cd_inv[1][1] * cd_inv[1][1]
             + a_coeff[1][1] * cd[1][0] * cd_inv[0][1] * cd_inv[1][1]
             + a_coeff[2][0] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1]
             + b_coeff[0][2] * cd[1][1] * cd_inv[1][1] * cd_inv[1][1]
             + b_coeff[1][1] * cd[1][1] * cd_inv[0][1] * cd_inv[1][1]
             + b_coeff[2][0] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] );


  /* pvi_5 */
  pv1[5] = ( 2*a_coeff[0][2]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][1]
             +   a_coeff[1][1] * cd[0][0] * cd_inv[0][0] * cd_inv[1][1]
             +   a_coeff[1][1] * cd[0][0] * cd_inv[0][1] * cd_inv[1][0]
             + 2*a_coeff[2][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1]
             + 2*b_coeff[0][2] * cd[0][1] * cd_inv[1][0] * cd_inv[1][1]
             +   b_coeff[1][1] * cd[0][1] * cd_inv[0][0] * cd_inv[1][1]
             +   b_coeff[1][1] * cd[0][1] * cd_inv[0][1] * cd_inv[1][0]
             + 2*b_coeff[2][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] );

  pv2[5] = ( 2*a_coeff[0][2]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][1]
             +   a_coeff[1][1] * cd[1][0] * cd_inv[0][0] * cd_inv[1][1]
             +   a_coeff[1][1] * cd[1][0] * cd_inv[0][1] * cd_inv[1][0]
             + 2*a_coeff[2][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1]
             + 2*b_coeff[0][2] * cd[1][1] * cd_inv[1][0] * cd_inv[1][1]
             +   b_coeff[1][1] * cd[1][1] * cd_inv[0][0] * cd_inv[1][1]
             +   b_coeff[1][1] * cd[1][1] * cd_inv[0][1] * cd_inv[1][0]
             + 2*b_coeff[2][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] );


  /* pvi_6 */
  pv1[6]= ( a_coeff[0][2]   * cd[0][0] * cd_inv[1][1] * cd_inv[1][1]
            + a_coeff[1][1] * cd[0][0] * cd_inv[0][1] * cd_inv[1][1]
            + a_coeff[2][0] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1]
            + b_coeff[0][2] * cd[0][1] * cd_inv[1][1] * cd_inv[1][1]
            + b_coeff[1][1] * cd[0][1] * cd_inv[0][1] * cd_inv[1][1]
            + b_coeff[2][0] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] );

  pv2[6]= ( a_coeff[0][2]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][0]
            + a_coeff[1][1] * cd[1][0] * cd_inv[0][0] * cd_inv[1][0]
            + a_coeff[2][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0]
            + b_coeff[0][2] * cd[1][1] * cd_inv[1][0] * cd_inv[1][0]
            + b_coeff[1][1] * cd[1][1] * cd_inv[0][0] * cd_inv[1][0]
            + b_coeff[2][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] );


  /* pvi_7 */
  pv1[7] = ( a_coeff[0][3]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
             + a_coeff[1][2] * cd[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
             + a_coeff[2][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
             + a_coeff[3][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0]
             + b_coeff[0][3] * cd[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
             + b_coeff[1][2] * cd[0][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
             + b_coeff[2][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
             + b_coeff[3][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] );

  pv2[7] = ( a_coeff[0][3]   * cd[1][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
             + a_coeff[1][2] * cd[1][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
             + a_coeff[2][1] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
             + a_coeff[3][0] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1]
             + b_coeff[0][3] * cd[1][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
             + b_coeff[1][2] * cd[1][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
             + b_coeff[2][1] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
             + b_coeff[3][0] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] );


  /* pvi_8 */
  pv1[8] = ( 3*a_coeff[0][3]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
             + 2*a_coeff[1][2] * cd[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
             +   a_coeff[1][2] * cd[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
             +   a_coeff[2][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
             + 2*a_coeff[2][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
             + 3*a_coeff[3][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1]
             + 3*b_coeff[0][3] * cd[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
             + 2*b_coeff[1][2] * cd[0][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
             +   b_coeff[1][2] * cd[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
             +   b_coeff[2][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
             + 2*b_coeff[2][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
             + 3*b_coeff[3][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] );

  pv2[8] = ( 3*a_coeff[0][3]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
             +   a_coeff[1][2] * cd[1][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
             + 2*a_coeff[1][2] * cd[1][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
             + 2*a_coeff[2][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
             +   a_coeff[2][1] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
             + 3*a_coeff[3][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1]
             + 3*b_coeff[0][3] * cd[1][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
             +   b_coeff[1][2] * cd[1][1] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
             + 2*b_coeff[1][2] * cd[1][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
             + 2*b_coeff[2][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
             +   b_coeff[2][1] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
             + 3*b_coeff[3][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] );


  /* pvi_9 */
  pv1[9] = ( 3*a_coeff[0][3]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
             +   a_coeff[1][2] * cd[0][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
             + 2*a_coeff[1][2] * cd[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
             + 2*a_coeff[2][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
             +   a_coeff[2][1] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
             + 3*a_coeff[3][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1]
             + 3*b_coeff[0][3] * cd[0][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
             +   b_coeff[1][2] * cd[0][1] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
             + 2*b_coeff[1][2] * cd[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
             + 2*b_coeff[2][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
             +   b_coeff[2][1] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
             + 3*b_coeff[3][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] );

  pv2[9] = ( 3*a_coeff[0][3]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
             + 2*a_coeff[1][2] * cd[1][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
             +   a_coeff[1][2] * cd[1][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
             +   a_coeff[2][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
             + 2*a_coeff[2][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
             + 3*a_coeff[3][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1]
             + 3*b_coeff[0][3] * cd[1][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
             + 2*b_coeff[1][2] * cd[1][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
             +   b_coeff[1][2] * cd[1][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
             +   b_coeff[2][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
             + 2*b_coeff[2][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
             + 3*b_coeff[3][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] );


  /* pvi_10 */
  pv1[10] = ( a_coeff[0][3]   * cd[0][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[1][2] * cd[0][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[2][1] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              + a_coeff[3][0] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1]
              + b_coeff[0][3] * cd[0][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[1][2] * cd[0][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[2][1] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              + b_coeff[3][0] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] );


  pv2[10] = ( a_coeff[0][3]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[1][2] * cd[1][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[2][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
              + a_coeff[3][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0]
              + b_coeff[0][3] * cd[1][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[1][2] * cd[1][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[2][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
              + b_coeff[3][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] );


  /* pvi_12 */
  pv1[12] = ( a_coeff[0][4]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[1][3] * cd[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[2][2] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[3][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
              + a_coeff[4][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0]
              + b_coeff[0][4] * cd[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[1][3] * cd[0][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[2][2] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[3][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
              + b_coeff[4][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] );

  pv2[12] = ( a_coeff[0][4]   * cd[1][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[1][3] * cd[1][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[2][2] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[3][1] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              + a_coeff[4][0] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1]
              + b_coeff[0][4] * cd[1][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[1][3] * cd[1][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[2][2] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[3][1] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              + b_coeff[4][0] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] );


  /* pvi_13 */
  pv1[13] = ( 4*a_coeff[0][4]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   a_coeff[1][3] * cd[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + 2*a_coeff[2][2] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
              + 2*a_coeff[2][2] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              +   a_coeff[3][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
              + 3*a_coeff[3][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
              + 4*a_coeff[4][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1]
              + 4*b_coeff[0][4] * cd[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[0][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   b_coeff[1][3] * cd[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + 2*b_coeff[2][2] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
              + 2*b_coeff[2][2] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              +   b_coeff[3][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
              + 3*b_coeff[3][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
              + 4*b_coeff[4][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] );


  pv2[13] = ( 4*a_coeff[0][4]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              +   a_coeff[1][3] * cd[1][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[1][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 2*a_coeff[2][2] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + 2*a_coeff[2][2] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              + 3*a_coeff[3][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              +   a_coeff[3][1] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 4*a_coeff[4][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1]
              + 4*b_coeff[0][4] * cd[1][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              +   b_coeff[1][3] * cd[1][1] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[1][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 2*b_coeff[2][2] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + 2*b_coeff[2][2] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              + 3*b_coeff[3][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              +   b_coeff[3][1] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 4*b_coeff[4][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] );


  /* pvi_14 */
  pv1[14] = ( 6*a_coeff[0][4]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   a_coeff[2][2] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
              + 4*a_coeff[2][2] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              +   a_coeff[2][2] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              + 3*a_coeff[3][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
              + 3*a_coeff[3][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 6*a_coeff[4][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1]
              + 6*b_coeff[0][4] * cd[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[0][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   b_coeff[2][2] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
              + 4*b_coeff[2][2] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              +   b_coeff[2][2] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              + 3*b_coeff[3][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
              + 3*b_coeff[3][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 6*b_coeff[4][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] );

  pv2[14] = ( 6*a_coeff[0][4]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[1][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[1][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   a_coeff[2][2] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
              + 4*a_coeff[2][2] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              +   a_coeff[2][2] * cd[1][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              + 3*a_coeff[3][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
              + 3*a_coeff[3][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 6*a_coeff[4][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1]
              + 6*b_coeff[0][4] * cd[1][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[1][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[1][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   b_coeff[2][2] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1]
              + 4*b_coeff[2][2] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              +   b_coeff[2][2] * cd[1][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              + 3*b_coeff[3][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1]
              + 3*b_coeff[3][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 6*b_coeff[4][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] );


  /* pvi_15 */
  pv1[15] = ( 4*a_coeff[0][4]   * cd[0][0] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              +  a_coeff[1][3]  * cd[0][0] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 2*a_coeff[2][2] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + 2*a_coeff[2][2] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              + 3*a_coeff[3][1] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              +   a_coeff[3][1] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 4*a_coeff[4][0] * cd[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1]
              + 4*b_coeff[0][4] * cd[0][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              +   b_coeff[1][3] * cd[0][1] * cd_inv[0][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1] * cd_inv[1][1]
              + 2*b_coeff[2][2] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + 2*b_coeff[2][2] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][1]
              + 3*b_coeff[3][1] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              +   b_coeff[3][1] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][0]
              + 4*b_coeff[4][0] * cd[0][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] );

  pv2[15] = ( 4*a_coeff[0][4]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              + 3*a_coeff[1][3] * cd[1][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   a_coeff[1][3] * cd[1][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + 2*a_coeff[2][2] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
              + 2*a_coeff[2][2] * cd[1][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              +   a_coeff[3][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
              + 3*a_coeff[3][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
              + 4*a_coeff[4][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1]
              + 4*b_coeff[0][4] * cd[1][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              + 3*b_coeff[1][3] * cd[1][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][1]
              +   b_coeff[1][3] * cd[1][1] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + 2*b_coeff[2][2] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][1]
              + 2*b_coeff[2][2] * cd[1][1] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0] * cd_inv[1][0]
              +   b_coeff[3][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][1]
              + 3*b_coeff[3][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] * cd_inv[1][0]
              + 4*b_coeff[4][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][1] );


  /* pvi_16 */
  pv1[16] = ( a_coeff[0][4]   * cd[0][0] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[1][3] * cd[0][0] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[2][2] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + a_coeff[3][1] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              + a_coeff[4][0] * cd[0][0] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1]
              + b_coeff[0][4] * cd[0][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[1][3] * cd[0][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[2][2] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1] * cd_inv[1][1]
              + b_coeff[3][1] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[1][1]
              + b_coeff[4][0] * cd[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] * cd_inv[0][1] );

  pv2[16] = ( a_coeff[0][4]   * cd[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[1][3] * cd[1][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[2][2] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
              + a_coeff[3][1] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
              + a_coeff[4][0] * cd[1][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0]
              + b_coeff[0][4] * cd[1][1] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[1][3] * cd[1][1] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[2][2] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0] * cd_inv[1][0]
              + b_coeff[3][1] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[1][0]
              + b_coeff[4][0] * cd[1][1] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] * cd_inv[0][0] );


  /* For a check:
  {
    size_t j;
    for(j=0; j<=16 ;++j)
      {
        if (j == 3 || j == 11) continue;
        printf("pv1_%ld \t %.12E\n", j, pv1[j]);
        printf("pv2_%ld \t %.12E\n", j, pv2[j]);
      }
  }
  */
}




















/**************************************************************/
/**********             Calculations               ************/
/**************************************************************/
/* Calcualte the final PV equation from intermidiate equations. */
static void
wcsdistortion_calc_tpveq(struct wcsprm *wcs, double cd[2][2],
                         double tpvu[8][8], double tpvv[8][8])
{
  /* tpvu and tpvv are u-v distortion equations in TPV convention. */
  size_t i=0,j=0;
  double cd_inv[2][2];
  double determinant=0;
  double k[5][5], l[5][5];
  double pv1[17]={0}, pv2[17]={0};

  /* Initialise the 2d matrices. */
  for(i=0; i<2; ++i) for(j=0; j<2; ++j) cd_inv[i][j]=0;
  for(i=0; i<5; ++i) for(j=0; j<5; ++j) { k[i][j]=0; l[i][j]=0; }

  /* Estimate the parameters. */
  wcsdistortion_get_tpvparams(wcs, cd, pv1, pv2);
  wcsdistortion_intermidate_tpveq(cd, pv1, pv2, k, l);

  /* We will find matrix tpvu and tpvv by finding inverse of
     CD matrix and multiplying with tpv* matrix.
     For inverse of a 2x2 matrix we use the below trasformations:
              inverse(|a  b|) =  1 *|d  -b|
                      |c  d|    |A| |-c  a|
      where |A| is the determinant of the matrix which is calculated by:
                          |A| = a*d-b*c.
    */
  determinant = cd[0][0]*cd[1][1] - cd[0][1]*cd[1][0];

  /* Inverse matrix */
  cd_inv[0][0]=cd[1][1]/determinant;      /* a */
  cd_inv[0][1]=(-1*cd[0][1])/determinant; /* b */
  cd_inv[1][0]=(-1*cd[1][0])/determinant; /* c */
  cd_inv[1][1]=cd[0][0]/determinant;      /* d */

  /*For a check.
  printf("%.10lf\t%.10lf\t%.10lf\t%.10lf\n", cd_inv[0][0], cd_inv[0][1], \
                                             cd_inv[1][0], cd_inv[1][1]);
  printf("%.10lf\n", determinant);
  */

  /* For matrix tpvv and tpvu, we have to use the following
     matrix equation:

                  |tpvu| = cd_inv*|k[i][j]|
                  |tpvv|          |l[i][j]|
    though intermidate polynomial equations have to be calculated prior
    to multiplycation with cd_inv.
  */

  for(i=0; i<=4; ++i)
    for(j=0; j<=4; ++j)
      {
        tpvu[i][j]=cd_inv[0][0]*k[i][j]+cd_inv[0][1]*l[i][j];
        tpvv[i][j]=cd_inv[1][0]*k[i][j]+cd_inv[1][1]*l[i][j];

        /*For a check:
        printf("%.10E, %.10E\n", tpvu[i][j], tpvv[i][j]);
        */
      }
}





/* Calcualte the final SIP equation from intermidiate equations
   and extract pv coefficients from it. */
static void
wcsdistortion_calc_sipeq(struct wcsprm *wcs, double cd[2][2],
                         double *pv1, double *pv2)
{
  /* tpvu and tpvv are u-v distortion equations in TPV convention. */
  size_t i=0, j=0;
  double cd_inv[2][2];
  double determinant=0;
  double a_coeff[5][5], b_coeff[5][5];

  /* Initialise the 2d matrices. */
  for(i=0;i<2;++i) for(j=0;j<2;++j) cd_inv[i][j]=0;
  for(i=0;i<5;++i) for(j=0;j<5;++j) {a_coeff[i][j]=0; b_coeff[i][j]=0;}

  /* We will find matrix tpvu and tpvv by finding inverse of
     CD matrix and multiplying with tpv* matrix.
     For inverse of a 2x2 matrix we use the below trasformations:
              inverse(|a  b|) =  1 *|d  -b|
                      |c  d|    |A| |-c  a|
      where |A| is the determinant of the matrix which is calculated by:
                          |A| = a*d-b*c.
    */
  wcsdistortion_get_sipparams(wcs, cd, a_coeff, b_coeff);
  determinant = cd[0][0]*cd[1][1] - cd[0][1]*cd[1][0];

  /* Inverse matrix */
  cd_inv[0][0]=cd[1][1]/determinant;      /* a */
  cd_inv[0][1]=(-1*cd[0][1])/determinant; /* b */
  cd_inv[1][0]=(-1*cd[1][0])/determinant; /* c */
  cd_inv[1][1]=cd[0][0]/determinant;      /* d */

  /* Find the parameters. */
  wcsdistortion_intermidate_sipeq( cd, cd_inv, a_coeff, b_coeff, pv1, pv2);

  /*For a check.
  printf("%.10lf\t%.10lf\t%.10lf\t%.10lf\n", cd_inv[0][0], cd_inv[0][1], \
                                             cd_inv[1][0], cd_inv[1][1]);
  printf("%.10lf\n", determinant);
  */
}





/* Calculate the SIP coefficients from CD matrix parameters and
   PV coefficients. */
static double
wcsdistortion_calcsip(size_t axis, size_t m, size_t n, double tpvu[8][8],
                      double tpvv[8][8])
{
  double sip_coeff;
  if(     axis == 1) sip_coeff=tpvu[m][n];
  else if(axis == 2) sip_coeff=tpvv[m][n];
  else error(EXIT_FAILURE, 0, "%s: axis does not exists! ", __func__);

  if(      (axis == 1) && (m == 1) && (n == 0) ) sip_coeff = sip_coeff - 1.0;
  else if( (axis == 2) && (m == 0) && (n == 1) ) sip_coeff = sip_coeff - 1.0;
  else error(EXIT_FAILURE, 0, "%s: axis does not exists! ", __func__);

  return sip_coeff;

}





/* To calculate reverse sip coefficients, we make a grid and polpulate
   it with forward coefficients. Then we fit a polynomial through these
   points and finally extract coefficeint from these polynomial which are
   the reverse sip cefficents. */
static void
wcsdistortion_fitreverse(double *u, double *v, size_t a_order, size_t b_order,
                         size_t naxis1, size_t naxis2, double a_coeff[5][5],
                         double b_coeff[5][5], double ap_coeff[5][5],
                         double bp_coeff[5][5])
{
  double chisq_ap, chisq_bp;
  double *udiff=NULL, *vdiff=NULL;
  double *uprime=NULL, *vprime=NULL;
  double **udict=NULL, **vdict=NULL;
  size_t tsize=(naxis1/4)*(naxis2/4);
  double **updict=NULL, **vpdict=NULL;
  gsl_vector *y_ap, *y_bp, *c_ap, *c_bp;
  size_t ap_order=a_order, bp_order=b_order;
  gsl_matrix *X_ap, *X_bp, *cov_ap, *cov_bp;
  size_t i=0, j=0, k=0, p_ap=0, p_bp=0, ij=0;
  gsl_multifit_linear_workspace *work_ap, *work_bp;

  /* Storage structures used:

     updict, vpdict - dictionaries (key-value pairs) with key being an
                      integer starting from 0 to max(aporder,
                      bporder)(inclusive) and its values are the uprime,
                      vprime raised to the powers of corresponding keys for
                      each key.

     udict, vdict   - dictionaries (key-value pairs) with key being an
                      integer starting from 0 to max(aorder,
                      border)(inclusive) and its values are the u, v raised
                      to the powers of corresponding keys for each key.

     u, v           - The 1d representation of 2d grid of all points
                      strating from -CRPIXi to NAXISi - CRPIXi (with a
                      stride of 4). CRPIXi is subtracted to bring pixels
                      in world coordinate system (wcs).

     uprime, vprime - The grid (represented internally as a 1d array)
                      with forward coefficients fitted on them
                      using the relevant udict, vdict values.

     udiff, vdiff   - 1d array with the values of uprime, vprim subtracted
                      from u, v arrays.


     In matrix equation AX=B,
     For axis 1 - A = transpose of a matrix of udict*vdict
                  B = udiff
     For axis 2 - A = transpose of a matrix of updict*vpdict
                  B = vdiff
    */

  /* Allocate updict and vpdict and initialize them. */
  updict=malloc((ap_order+1)*sizeof(*updict));
  vpdict=malloc((bp_order+1)*sizeof(*vpdict));
  for(i=0; i<=wcsdistortion_max(ap_order, bp_order); ++i)
    {
      updict[i]=malloc(tsize*sizeof(updict));
      vpdict[i]=malloc(tsize*sizeof(vpdict));
    }

  /* Allocate and initialize uprime and vprime. */
  uprime=malloc(tsize*sizeof(*uprime));
  vprime=malloc(tsize*sizeof(*vprime));
  for(i=0; i<tsize; ++i)
    {
      uprime[i]=u[i];
      vprime[i]=v[i];
    }

  /* Allocate and initialize udict and vdict.*/
  udict=malloc((a_order+1)*sizeof(*udict));
  vdict=malloc((b_order+1)*sizeof(*vdict));
  for(i=0; i<=wcsdistortion_max(a_order, b_order); ++i)
    {
      udict[i]=malloc(tsize*sizeof(udict));
      vdict[i]=malloc(tsize*sizeof(vdict));
    }
  for(i=0; i<tsize; ++i)
    {
      udict[0][i]=1;
      vdict[0][i]=1;
    }


  /* Fill the values from the in the dicts. The rows of the
      dicts act as a key to achieve a key-value functionality. */
  for(i=1; i<=a_order; ++i)
    for(j=0; j<tsize; ++j)
      udict[i][j]=udict[i-1][j]*u[j];

  for(i=1; i<=b_order; ++i)
    for(j=0; j<tsize; ++j)
      vdict[i][j]=vdict[i-1][j]*v[j];


  /* Populating forward coefficients on a grid. */
  for(i=0; i<=a_order; ++i)
    for(j=0; j<=a_order-i; ++j)
      {
        for(k=0; k<tsize; ++k)
          uprime[k]+=a_coeff[i][j]*udict[i][k]*vdict[j][k];

        /* The number of parameters for AP_* coefficients. */
        ++p_ap;
      }
  for(i=0; i<=b_order; ++i)
    for(j=0; j<=b_order-i; ++j)
      {
        for(k=0; k<tsize; ++k)
          vprime[k]+=b_coeff[i][j]*udict[i][k]*vdict[j][k];

        /* The number of parameters for BP_* coefficients. */
        ++p_bp;
      }

  /*For a check.
  for(i=0; i<=a_order; ++i)
    for(j=0; j<tsize; ++j)
    {
      printf("udict[%ld][%ld] = %.8lf\n", i, j, uprime[j]);
      printf("u%ld = %.10E\n", i, u[i]);
    }
  */

  /* Now we have a grid populated with forward coeffiecients.  Now we fit a
     reverse polynomial through points using multiparameter linear least
     square fittings.*/

  /* Initialize dicts. */
  for(i=0; i<tsize; ++i)
    {
      updict[0][i]=1;
      vpdict[0][i]=1;
    }

  /* Fill the values from the in the dicts. The rows of the
      dicts act as a key to achieve a key-value functionality. */
  for(i=1; i<=ap_order; ++i)
    for(j=0; j<tsize; ++j)
      updict[i][j]=updict[i-1][j]*uprime[j];

  for(i=1; i<=bp_order; ++i)
    for(j=0; j<tsize; ++j)
      vpdict[i][j]=vpdict[i-1][j]*vprime[j];

  /* Allocate memory for Multi-parameter Linear Regressions. */
  X_ap = gsl_matrix_alloc (tsize, p_ap);
  X_bp = gsl_matrix_alloc (tsize, p_bp);

  y_ap = gsl_vector_alloc (tsize);
  y_bp = gsl_vector_alloc (tsize);

  c_ap = gsl_vector_alloc (p_ap);
  c_bp = gsl_vector_alloc (p_bp);

  cov_ap = gsl_matrix_alloc (p_ap, p_ap);
  cov_bp = gsl_matrix_alloc (p_bp, p_bp);

  /* Allocate and initialize udiff and vdiff. */
  udiff=malloc(tsize*sizeof(*udiff));
  vdiff=malloc(tsize*sizeof(*vdiff));
  for(i=0; i<tsize; ++i)
    {
      udiff[i]=u[i]-uprime[i];
      vdiff[i]=v[i]-vprime[i];

      gsl_vector_set (y_ap, i, udiff[i]);
      gsl_vector_set (y_bp, i, vdiff[i]);
    }

  /* Filling up he X_ap matrix in column-wise order. */
  for(i=0; i<=ap_order; ++i)
    for(j=0; j<=ap_order-i; ++j)
      {
        for(k=0; k<tsize; ++k)
          {
            gsl_matrix_set(X_ap, k, ij, updict[i][k]*vpdict[j][k]);

            /*For a check.
            printf("x_ap[%ld] = %.8lf\n", ij, updict[i][k]*vpdict[j][k]);
            */
          }
        ++ij;
      }

  ij=0;
  /* Filling up he X_bp matrix in column-wise order. */
  for(i=0; i<=bp_order; ++i)
    for(j=0; j<=bp_order-i; ++j)
      {
        for(k=0; k<tsize; ++k)
          {
            gsl_matrix_set(X_bp, k, ij, updict[i][k]*vpdict[j][k]);

            /*For a check.
            printf("x_bp[%ld] = %.8lf\n", ij, updict[i][k]*vpdict[j][k]);
            */
          }
        ++ij;
      }

  /* Initialize and do the linear least square fitting. */
  work_ap = gsl_multifit_linear_alloc (tsize, p_ap);
  work_bp = gsl_multifit_linear_alloc (tsize, p_bp);
  gsl_multifit_linear(X_ap, y_ap, c_ap, cov_ap, &chisq_ap, work_ap);
  gsl_multifit_linear(X_bp, y_bp, c_bp, cov_bp, &chisq_bp, work_bp);

  /* Extract the reverse coefficients. */
  p_ap=0;
  for(i=0; i<=ap_order; ++i)
    for(j=0; j<=ap_order-i; ++j)
    {
        ap_coeff[i][j]=gsl_vector_get(c_ap, p_ap);

        /*For a check.
        printf("AP_%ld_%ld = %.8E\n", i, j, ap_coeff[i][j]);
        */

        ++p_ap;
    }
  p_bp=0;
  for(i=0; i<=bp_order; ++i)
    for(j=0; j<=bp_order-i; ++j)
    {
        bp_coeff[i][j]=gsl_vector_get(c_bp, p_bp);

        /*For a check.
        printf("BP_%ld_%ld = %.8E\n", i, j, bp_coeff[i][j]);
        */

        ++p_bp;
    }

  /*For a check.
  for(i=0; i<p_ap; ++i)
    for(j=0; j<tsize; ++j)
      printf("X[%ld][%ld] = %.8lf\n", i, j, gsl_matrix_get(X_ap, j, i));
  */

  /* Free the memory allocations.*/
  gsl_multifit_linear_free(work_bp);
  gsl_multifit_linear_free(work_ap);
  gsl_matrix_free(cov_bp);
  gsl_matrix_free(cov_ap);
  gsl_vector_free(c_bp);
  gsl_vector_free(c_ap);
  gsl_vector_free(y_bp);
  gsl_vector_free(y_ap);
  gsl_matrix_free(X_bp);
  gsl_matrix_free(X_ap);
  free(vdiff);
  free(udiff);
  free(vprime);
  free(uprime);

  for(i=0; i<ap_order; ++i) { free(vpdict[i]); free(updict[i]); }
  for(i=0; i<a_order; ++i) { free(vdict[i]); free(udict[i]); }
  free(vpdict);
  free(updict);
  free(vdict);
  free(udict);
}






/* Calculate the reverse sip coefficients. */
static void
wcsdistortion_get_revkeyvalues(struct wcsprm *wcs, size_t *fitsize,
                               double ap_coeff[5][5], double bp_coeff[5][5])
{
  size_t tsize;
  size_t i, j, k;
  double *u=NULL, *v=NULL;
  size_t a_order=0, b_order=0;
  double a_coeff[5][5], b_coeff[5][5];
  size_t naxis1=fitsize[1], naxis2=fitsize[0];
  double crpix1=wcs->crpix[0], crpix2=wcs->crpix[1];

  /* Initialise the 2d matrices. */
  tsize=(naxis1/4)*(naxis2/4);
  for(i=0;i<5;++i) for(j=0;j<5;++j) {a_coeff[i][j]=0; b_coeff[i][j]=0;}

  /* Allocate the size of u,v arrays. */
  u=malloc(tsize*sizeof(*u));
  v=malloc(tsize*sizeof(*v));

  if(u==NULL)
    error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for `u'",
          __func__, tsize*sizeof(*u));
  if(v==NULL)
    error(EXIT_FAILURE, 0, "%s: allocating %zu bytes for `v'",
          __func__, tsize*sizeof(*v));


  /* Make the grid and bring it's origin to the world's origin*/
  k=0; for(i=0; i<naxis2; i+=4) for(j=0; j<naxis1; j+=4) u[k++]=j-crpix1;
  k=0; for(i=0; i<naxis2; i+=4) for(j=0; j<naxis1; j+=4) v[k++]=i-crpix2;
  /*For a check.
    for(i=0; i<tsize; ++i)
    printf("u%ld = %.10E\n", i, u[i]);
  */

  wcsdistortion_get_sipcoeff(wcs, &a_order, &b_order, a_coeff, b_coeff);

  wcsdistortion_fitreverse(u, v, a_order, b_order, naxis1, naxis2,
                           a_coeff, b_coeff, ap_coeff, bp_coeff);

  /* Free the memory allocations. */
  free(v);
  free(u);
}




















/**************************************************************/
/**********          Writing utilities             ************/
/**************************************************************/
/* Make the sip key-cards and add them to fullheader.

   Return:
   char *fullheader - string of keycards.*/
static char *
wcsdistortion_add_sipkeywords(struct wcsprm *wcs, size_t *fitsize,
                              double tpvu[8][8], double tpvv[8][8],
                              int add_reverse, int *nkeys)
{
  double val=0;
  uint8_t i, j, k=0;
  int size = wcs->naxis;
  size_t a_order=0, b_order=0;
  size_t ap_order=0, bp_order=0;
  size_t m, n, num=0, numkey=100;
  double ap_coeff[5][5], bp_coeff[5][5];
  char *fullheader, fmt[50], sipkey[8], keyaxis[9], pcaxis[10];

  /* Initialise the 2d matrices. */
  *nkeys = 0;
  for(i=0;i<5;++i) for(j=0;j<5;++j) {ap_coeff[i][j]=0; bp_coeff[i][j]=0;}

  /* The format for each card. */
  sprintf(fmt, "%%-8s= %%20.12E%%50s");

  /* Allcate memory for cards. */
  fullheader = malloc(numkey*80);
  if(fullheader==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `fullheader'",
          __func__, sizeof *fullheader);

  /* Add other necessary cards. */
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20d%50s", "WCSAXES",
          wcs->naxis, "");

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CRPIX%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.8lf%50s", keyaxis,
              wcs->crpix[i-1], "");
    }

  for(i=1; i<=size; ++i)
    for(j=1; j<=size; ++j)
      {
        sprintf(pcaxis, "PC%d_%d", i, j);
        sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", pcaxis,
                wcs->pc[k++], "");
      }

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CDELT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", keyaxis,
              wcs->cdelt[i-1], "");
    }

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CUNIT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", keyaxis,
              wcs->cunit[i-1]);
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE1",
          "'RA---TAN-SIP'");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE2",
          "'DEC--TAN-SIP'");

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CRVAL%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.10lf%50s", keyaxis,
              wcs->crval[i-1], "");
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LONPOLE",
          wcs->lonpole, "");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LATPOLE",
          wcs->latpole, "");

#if GAL_CONFIG_HAVE_WCSLIB_MJDREF == 1
  for(i=1; i<=size; ++i)
    sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "MJDREFI",
            wcs->mjdref[i-1], "");
#endif

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "RADESYS",
          wcs->radesys);

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "EQUINOX",
          wcs->equinox, "");


  for(m=0; m<=4; ++m)
    for(n=0; n<=4; ++n)
      {
        /*For axis = 1*/
        val=wcsdistortion_calcsip(1, m, n, tpvu, tpvv);
        if(val != 0)
          {
            /* Make keywords */
            sprintf(sipkey, "A_%ld_%ld", m, n);
            sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey, val, "");
            a_order=wcsdistortion_max(a_order, wcsdistortion_max(m,n));
          }

        /*For axis = 2*/
        val=wcsdistortion_calcsip(2, m, n, tpvu, tpvv);
        if(val != 0)
          {
            /* Make keywords */
            sprintf(sipkey, "B_%ld_%ld", m, n);
            sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey, val, "");
            b_order=wcsdistortion_max(b_order, wcsdistortion_max(m,n));
          }

      }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "A_ORDER",
          a_order, "");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "B_ORDER",
          b_order, "");

  /* If reverse coefficients are required. */
  if( add_reverse )
    {
      ap_order=a_order;
      bp_order=b_order;

      wcsdistortion_get_revkeyvalues(wcs, fitsize, ap_coeff, bp_coeff);

      for(m=0; m<=ap_order; ++m)
        for(n=0; n<=ap_order-m; ++n)
          {
            /*For axis = 1*/
            val=ap_coeff[m][n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(sipkey, "AP_%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey,
                        val, "");
              }

            /*For axis = 2*/
            val=bp_coeff[m][n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(sipkey, "BP_%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey,
                        val, "");
              }
          }

      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "AP_ORDER",
              ap_order, "");
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "BP_ORDER",
              bp_order, "");

    }

  /*For a check.
    printf("%s\n", fullheader);
  */
  *nkeys = num;
  return fullheader;
}






/* Make the pv key-cards and add them to fullheader.

   Return:
   char *fullheader - string of keycards. */
static char *
wcsdistortion_add_pvkeywords(struct wcsprm *wcs, double *pv1,
                             double *pv2, int *nkeys)
{

  double val=0;
  uint8_t i, j, k=0;
  int size = wcs->naxis;
  size_t m, n, num=0, numkey=100;
  char *fullheader, fmt[50], pvkey[8], keyaxis[9], pcaxis[10];

  /* Initialize values. */
  *nkeys = 0;

  /* The format for each card. */
  sprintf(fmt, "%%-8s= %%20.12E%%50s");

  /* Allcate memory for cards. */
  fullheader = malloc(numkey*80);
  if(fullheader==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for `fullheader'",
              __func__, sizeof *fullheader);

  /* Add other necessary cards. */
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20d%50s", "WCSAXES",
          wcs->naxis, "");

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CRPIX%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.8lf%50s", keyaxis,
              wcs->crpix[i-1], "");
    }

  for(i=1; i<=size; ++i)
    for(j=1; j<=size; ++j)
      {
        sprintf(pcaxis, "PC%d_%d", i, j);
        sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", pcaxis,
                wcs->pc[k++], "");
      }

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CDELT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", keyaxis,
              wcs->cdelt[i-1], "");
    }

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CUNIT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", keyaxis,
              wcs->cunit[i-1]);
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE1",
          "'RA---TPV'");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE2",
          "'DEC--TPV'");

  for(i=1; i<=size; ++i)
    {
      sprintf(keyaxis, "CRVAL%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.10lf%50s", keyaxis,
              wcs->crval[i-1], "");
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LONPOLE",
          wcs->lonpole, "");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LATPOLE",
          wcs->latpole, "");

#if GAL_CONFIG_HAVE_WCSLIB_MJDREF == 1
  for(i=1; i<=size; ++i)
    sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "MJDREFI",
            wcs->mjdref[i-1], "");
#endif

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "RADESYS",
          wcs->radesys);

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "EQUINOX",
          wcs->equinox, "");

  for(m=1; m<=2; ++m)
    for(n=0; n<=16; ++n)
      {
        /*For axis = 1*/
        if(m == 1)
          {
            val=pv1[n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(pvkey, "PV%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, pvkey,
                        val, "");
              }
          }

        /*For axis = 2*/
        if(m == 2)
          {
            val=pv2[n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(pvkey, "PV%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, pvkey,
                        val, "");
              }
          }

      }

  /*For a check.
  printf("%s\n", fullheader);
  */
  *nkeys = num;
  return fullheader;

}





 /* Set the internal structure and do sanity checks. */
static void
wcsdistortion_set_internalstruct(struct wcsprm *wcs, char *fullheader,
                                 size_t fulllen, int status, int nwcs,
                                 int sumcheck)
{
  size_t i;

  if( wcs == NULL )
  {

    fprintf(stderr, "\n##################\n"
            "WCSLIB Warning: wcspih ERROR %d: %s.\n"
            "##################\n",
            status, wcs_errmsg[status]);
    wcs=NULL; nwcs=0;
  }

  if(wcs)
    {
      /* It may happen that the WCS-related keyword values are stored as
         strings (they have single-quotes around them). In this case,
         WCSLIB will read the CRPIX and CRVAL values as zero. When this
         happens do a small check and abort, while informing the user about
         the problem. */
      sumcheck=0;
      for(i=0;i<wcs->naxis;++i)
        {sumcheck += (wcs->crval[i]==0.0f) + (wcs->crpix[i]==0.0f);}
      if(sumcheck==wcs->naxis*2)
        {
          /* We only care about the first set of characters in each
             80-character row, so we don't need to parse the last few
             characters anyway. */
          fulllen=strlen(fullheader)-12;
          for(i=0;i<fulllen;++i)
            if( strncmp(fullheader+i, "CRVAL1  = '", 11) == 0 )
              fprintf(stderr, "WARNING: WCS Keyword values are not "
                      "numbers.\n\n"
                      "WARNING: The values to the WCS-related keywords are "
                      "enclosed in single-quotes. In the FITS standard "
                      "this is how string values are stored, therefore "
                      "WCSLIB is unable to read them AND WILL PUT ZERO IN "
                      "THEIR PLACE (creating a wrong WCS in the output). "
                      "Please update the respective keywords of the input "
                      "to be numbers (see next line).\n\n"
                      "WARNING: You can do this with Gnuastro's `astfits' "
                      "program and the `--update' option. The minimal WCS "
                      "keywords that need a numerical value are: `CRVAL1', "
                      "`CRVAL2', `CRPIX1', `CRPIX2', `EQUINOX' and "
                      "`CD%%_%%' (or `PC%%_%%', where the %% are integers), "
                      "please see the FITS standard, and inspect your FITS "
                      "file to identify the full set of keywords that you "
                      "need correct (for example PV%%_%% keywords).\n\n");
        }

      /* CTYPE is a mandatory WCS keyword, so if it hasn't been given (its
         '\0'), then the headers didn't have a WCS structure. However,
         WCSLIB still fills in the basic information (for example the
         dimensionality of the dataset). */
      if(wcs->ctype[0][0]=='\0')
        {
          wcsfree(wcs);
          wcs=NULL;
          nwcs=0;
        }
      else
        {
          /* For a check.
          printf("flag: %d\n", wcs->flag);
          printf("naxis: %d\n", wcs->naxis);
          printf("crpix: %f, %f\n", wcs->crpix[0], wcs->crpix[1]);
          printf("pc: %f, %f, %f, %f\n", wcs->pc[0], wcs->pc[1],
                 wcs->pc[2], wcs->pc[3]);
          printf("cdelt: %f, %f\n", wcs->cdelt[0], wcs->cdelt[1]);
          printf("crval: %f, %f\n", wcs->crval[0], wcs->crval[1]);
          printf("cunit: %s, %s\n", wcs->cunit[0], wcs->cunit[1]);
          printf("ctype: %s, %s\n", wcs->ctype[0], wcs->ctype[1]);
          printf("lonpole: %f\n", wcs->lonpole);
          printf("latpole: %f\n", wcs->latpole);
          */

          /* Set the WCS structure. */
          status=wcsset(wcs);
          if(status)
            {
              fprintf(stderr, "\n##################\n"
                      "WCSLIB Warning: wcsset ERROR %d: %s.\n"
                      "##################\n",
                      status, wcs_errmsg[status]);
              wcsfree(wcs);
              wcs=NULL;
              nwcs=0;
            }
          else
            /* A correctly useful WCS is present. When no PC matrix
               elements were present in the header, the default PC matrix
               (a unity matrix) is used. In this case WCSLIB doesn't set
               `altlin' (and gives it a value of 0). In Gnuastro, later on,
               we might need to know the type of the matrix used, so in
               such a case, we will set `altlin' to 1. */
            if(wcs->altlin==0) wcs->altlin=1;
        }
    }

  /* For a check.
    wcsprt(wcs);
  */

  /* Clean up. */
  status=0;
  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in freeing the memory used to "
                      "keep all the headers");
}




















/************************************************************************/
/***************         High-level functions       *********************/
/************************************************************************/
/* After the actual TPV->SIP conversions are done and the keycards are
   created using `wcsdistortion_add_sipkeywords`, we then use `wcspih`
   to generate and return the new headers having SIP coefficients.

   Parameters:
     struct wcsprm *inwcs - The wcs parameters of the input fits file.
     size_t *fitsize      - The size of the array along each dimension.

   Return:
     struct wcsprm *outwcs - The transformed wcs parameters in the
                             sip distortion type.*/
struct wcsprm *
gal_wcsdistortion_tpv_to_sip(struct wcsprm *inwcs,
                             size_t *fitsize)
{
  int ctrl=0;                /* Don't report why a keyword wasn't used. */
  int nreject=0;             /* Number of keywords rejected for syntax. */
  double cd[2][2];
  size_t i=0, j=0;
  size_t fulllen=0;
  char *fullheader;
  int relax=WCSHDR_all;      /* Macro: use all informal WCS extensions. */
  int nwcs, sumcheck=0;
  int nkeys=0, status=0;
  struct wcsprm *outwcs=NULL;
  double tpvu[8][8], tpvv[8][8];

  /* Initialise the 2d matrices. */
  for(i=0; i<2; ++i) for(j=0; j<2; ++j) cd[i][j]=0;
  for(i=0; i<8; ++i) for(j=0; j<8; ++j) { tpvu[i][j]=0; tpvv[i][j]=0; }

  /* Calculate the pv equations and extract sip coefficients from it. */
  wcsdistortion_calc_tpveq(inwcs, cd, tpvu, tpvv);

  /* Add the sip keywords. */
  fullheader=wcsdistortion_add_sipkeywords(inwcs, fitsize, tpvu, tpvv,
                                           1, &nkeys);

  /* WCSlib function to parse the FITS headers. */
  status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, &nwcs, &outwcs);

  /* Set the internal structure. */
  wcsdistortion_set_internalstruct(outwcs, fullheader, fulllen, status,
                                   nwcs, sumcheck);

  /* Return the output WCS. */
  return outwcs;

}





/* After the actual SIP->TPV conversions are done and the keycards are
   created using `wcsdistortion_add_pvkeywords`, we then use `wcspih`
   to generate and return the new headers having PV coefficients.

   Parameters:
     struct wcsprm *inwcs - The wcs parameters of the input fits file.

   Return:
     struct wcsprm *outwcs - The transformed wcs parameters in the
                             pv distortion type.*/
struct wcsprm *
gal_wcsdistortion_sip_to_tpv(struct wcsprm *inwcs)
{
  int ctrl=0;                /* Don't report why a keyword wasn't used. */
  int nreject=0;             /* Number of keywords rejected for syntax. */
  double cd[2][2];
  size_t i=0, j=0;
  size_t fulllen=0;
  char *fullheader;
  int nwcs, sumcheck=0;
  int nkeys=0, status=0;
  int relax=WCSHDR_all;      /* Macro: use all informal WCS extensions. */
  struct wcsprm *outwcs=NULL;
  double pv1[17]={0}, pv2[17]={0};

  /* Initialise the 2d matrices. */
  for(i=0; i<2; ++i) for(j=0; j<2; ++j) cd[i][j]=0;

  /* Calculate the sip equations and extract pv coefficients from it. */
  wcsdistortion_calc_sipeq(inwcs, cd, pv1, pv2);

  /* Add the pv keywords. */
  fullheader=wcsdistortion_add_pvkeywords(inwcs, pv1, pv2, &nkeys);


  /* WCSlib function to parse the FITS headers. */
  status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, &nwcs, &outwcs);

  /* Set the internal structure. */
  wcsdistortion_set_internalstruct(outwcs, fullheader, fulllen, status,
                                   nwcs, sumcheck);

  return outwcs;
}
