/*=========================================================================
 * * Program :   $Insight Segmentation & Registration Toolkit $ * Module  :   $DMDBiomarkerTool: DMDBase.h $ * Purpose :   $The base class of Duchenne Muscle Dystrophy biomarker tools $
 * Language:   $C++ $
 * Date    :   $Date: 2010-06-04 12:36:34 $
 * Version :   $Revision: 0.10 $
 * Authors :   $Jiahui Wang, Martin Styner $
 * Update  :   1. Created the data class (06-04-10)
 * Copyright (c) Nero Image Research and Analysis Lab.@UNC All rights reserved.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even 
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.
 *
=========================================================================*/
#define false 0
#define true !false
#define progtitle    "POLYFIT"
#define progversion  "1.0.0"
#define progdate     "08.03.2003"
#define progauther   "K. Hough"
//#define MAXPAIRS 50
#define MAXPAIRS 4       // if 10 echos, then use 9 / if 5 echos then use 4
#define MAXORDER 8      /* Can do higher, but who cares?   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>      /* need functions 'toupper' and tolower'   */
#include <float.h>      /* defines LDBL_MAX*/
#include <math.h>
#include <time.h>       /* for date marking of output data files   */

//using namespace std;

#ifndef DMDCURVEFIT    // prevent for redefine
#define DMDCURVEFIT    
extern "C" { 
#include "f2c.h" 
void sgels_(char *trans, integer *m, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, real *work, integer *lwork, integer *info);
}

class DMDCurveFit
{
    private:
        int data_pairs, order;
        float correl_coef;
        long double xmin, xmax, ymin, ymax;
        long double polycoefs[MAXPAIRS];
        int debug;   /* enable display of results on screen -- reset via command line to true or false */

        int polynomfit( void );
        void datasort ( void );
        float poly_correlcoef( void );
        void get_xy_max_min ( void );
        void reporterror( char erstr[] );
        long double ldpow( long double n, unsigned p );
 
        double data_array[2][MAXPAIRS]; 

    public:
        DMDCurveFit( ) {
            
        }

        double polycurvefit( double data[2][MAXPAIRS] );
        double lslcurvefit( double data[2][MAXPAIRS] );

};
/*-------------------------------------------------------------------------*/
int DMDCurveFit::polynomfit( void )
{
     /* Polynomial Interpolation
    ~~~~~~~~~~~~~~~~~~~~~~~~*/
     /* NB Data array uses elements [x][1] to [x][data_pairs]  */
     /*    All data values MUST be positive                    */

     /*--------------------------------------------------------*/
     /*    This algorithm has been translated and developed    */
     /*    from an original Acorn BBC/BBC Basic programme      */
     /*    (ref. D.G.K.Guy, Practical Computing, May 1985,     */
     /*    page 135)                                           */
     /*    I still don't understand how it works!              */
     /*--------------------------------------------------------*/

     /* dimension arrays/variables
    ~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    long double a[MAXPAIRS+1][MAXPAIRS+1], B[MAXPAIRS+1], C[MAXPAIRS+1];
    long double S[MAXPAIRS+1];
    long double A1, A2, Y1, m, S1, x1;
    long double xscale, yscale;
    long double xlow = xmax;
    long double ylow = ymax;

    int i, j, k, L, R;

    /* set scaling factors
    ~~~~~~~~~~~~~~~~~~~~~~*/
    /* -- Avoids having to manage very large/very small numbers in polynomial
          calculation
              See rescaling ruitine at end of this function
    */

    for( i = 0; i < data_pairs; i++ )
        {
        if( data_array[0][i] < xlow && data_array[0][i]    != 0 ) xlow = data_array[0][i];
        if( data_array[1][i] < ylow && data_array[1][i]    != 0 ) ylow = data_array[1][i];
        }

    if( xlow < .001 && xmax < 1000 ) xscale = 1 / xlow;
        else if( xmax > 1000 && xlow > .001 ) xscale = 1 / xmax;
        else xscale = 1;

    if( ylow < .001 && ymax < 1000 ) yscale = 1 / ylow;
        else if( ymax > 1000 && ylow > .001 ) yscale = 1 / ymax;
        else yscale = 1;

        /* initialise array variables
    ~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    for(i = 0; i <= MAXPAIRS; i++ )
        {
        B[i] = 0; C[i] = 0; S[i] = 0;
        for( j = 0; j < MAXPAIRS; j++ )
            a[i][j] = 0;
        }

    for(i = 0; i <= MAXORDER; i++ )
        polycoefs[i] = 0;

     /* ensure all data is in ascending order wrt x
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    datasort();

     /* start of polynomial fit
    ~~~~~~~~~~~~~~~~~~~~~~~*/
    Y1 = 0;
    for( j = 1; j <= data_pairs; j++ ) {
        for( i = 1; i <= order; i++ ) {
            B[i] = B[i] + data_array[1][j] * yscale * ldpow( data_array[0][j] * xscale, i );
            if( B[i] == LDBL_MAX ) 
                return false;
            for( k = 1; k <= order; k++ ) {
                a[i][k] = a[i][k] + ldpow( data_array[0][j] * xscale, (i + k) );
                if( a[i][k] == LDBL_MAX ) 
                    return false;
            }
            S[i] = S[i] + ldpow( data_array[0][j] * xscale, i );
            if( S[i] == LDBL_MAX ) 
                return false;
        }
        Y1 = Y1 + data_array[1][j] * yscale;
        if( Y1 == LDBL_MAX ) 
            return false;
    }
    /*-----------------------------------------------------------------*/
    for( i = 1; i<= order; i++ ) {
        for( j = 1; j <= order; j++ ) {
            a[i][j] = a[i][j] - S[i] * S[j] / (long double) data_pairs;
            if( a[i][j] == LDBL_MAX ) return false;
        }
        B[i] = B[i] - Y1 * S[i] / (long double) data_pairs;
        if( B[i] == LDBL_MAX ) return false;
    }
    /*-----------------------------------------------------------------*/
    for( k = 1; k <= order; k++ )
        {
        R = k; A1 = 0;
        for( L = k; L <= order; L++ )
            {
            A2 = fabsl( a[L][k] );
            if( A2 > A1 ){ A1 = A2; R = L;}
            }
        if( A1 == 0 ) return false;
        if( R == k ) goto polfit1;
        for( j = k; j <= order; j++ )
            {
            x1 = a[R][j]; a[R][j] = a[k][j]; a[k][j] = x1;
            }
        x1 = B[R]; B[R] = B[k]; B[k] = x1;
polfit1:        for( i = k; i <= order; i++ )
            {
            m = a[i][k];
            for( j = k; j <= order; j++ )
                {
                if( i == k )
                    a[i][j] = a[i][j] / m;
                     else
                    a[i][j] = a[i][j] - m * a[k][j];
                }
            if( i == k )
                B[i] = B[i] / m;
                 else
                B[i] = B[i] - m * B[k];
            }
        }
    /*-----------------------------------------------------------------*/
    polycoefs[order] = B[order];
    for( k = 1; k <= order - 1; k++ ) {
        i = order - k; 
        S1 = 0;
        for( j = 1; j <= order; j++ ) {
            S1 = S1 + a[i][j] * polycoefs[j];
            if( S1 == LDBL_MAX ) 
                return false;
        }
        polycoefs[i] = B[i] - S1;
    }
    /*-----------------------------------------------------------------*/
    S1 = 0;
    for( i = 1; i <= order; i++ ) {
        S1 = S1 + polycoefs[i] * S[i] / (long double)data_pairs;
        if( S1 == LDBL_MAX ) 
            return false;
    }
    polycoefs[0] = (Y1 / (long double)data_pairs - S1);
    /*-----------------------------------------------------------------*/
    /* zero all coeficient values smaller than +/- .000001 (avoids -0)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    for( i = 0; i <= order; i++ )
        if ( fabsl( polycoefs[i] * 1000000 ) < 1 ) 
            polycoefs[i] = 0;

    /* rescale parameters
    ~~~~~~~~~~~~~~~~~~~~~*/
    for( i = 0; i <= order; i++ )
        polycoefs[i] = ( 1 / yscale ) * polycoefs[i] * ldpow( xscale, i );

    /* calc correlation coeficient
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    correl_coef = poly_correlcoef();

    /* No errors so return 'true'
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    return true;
}
/**************************************************************************/
void DMDCurveFit::datasort( void )
{
/* sort data in array data_array[i][j]
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	-- in ascending order for x
	~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	/*       i = 0 -- xdata
	 *       i = 1 -- ydata
	 *       j = number of data pair
*/

    int last, inner, outer, swapflag;
    long double k;

    last = data_pairs - 1;
    for( outer = 0; outer <= data_pairs; outer++ ) {
        swapflag = 0;
        for( inner = 0; inner <= last; inner++ ) {
            if( data_array[0][inner] > data_array[0][inner + 1] ) {
                k = data_array[0][inner];
                data_array[0][inner] = data_array[0][inner + 1];
                data_array[0][inner + 1] = k;

                k = data_array[1][inner];
                data_array[1][inner] = data_array[1][inner + 1];
                data_array[1][inner + 1] = k;
                swapflag = 1;
            }
        }
        last = last - 1;
        if( swapflag == 0 ) 
            outer = data_pairs;
    }
}
/**************************************************************************/
float DMDCurveFit::poly_correlcoef( void )
{
    int i, j;
    long double E[MAXORDER+1];
    long double ycalc[MAXPAIRS+1];
    long double my, mycalc, yc, SAA, SBB, SAB, A, B;
    long double SY, SYCALC, SYYCALC;

    for( i = 0; i <= order; i++ )
        E[i] = polycoefs[i];

    my = 0; mycalc = 0;
    for( i = 1; i <= data_pairs; i++ ){
        my = my + data_array[1][i];
        yc = 0;
        for( j = 0; j <= order; j++ )
             yc = yc + E[j] * ldpow( data_array[0][i], j );
        ycalc[i] = yc;
        mycalc = mycalc + ycalc[i];
    }

    my = my / data_pairs;
    mycalc = mycalc / data_pairs;
    SAA = 0; SBB = 0; SAB = 0;
    for( i = 1; i <= data_pairs; i++ ) {
        A = data_array[1][i] - my;
        B = ycalc[i] - mycalc;
        SAA = SAA + A * A;
        SBB = SBB + B * B;
        SAB = SAB + A * B;
    }

    SY = sqrtl( SAA / data_pairs );
    SYCALC = sqrtl( SBB / data_pairs );
    SYYCALC = SAB / data_pairs;
    
    return (float) ( SYYCALC / ( SY * SYCALC ) );

}
/**************************************************************************/
void DMDCurveFit::get_xy_max_min( void )
{
    int i;

    xmin = data_array[0][1]; xmax = data_array[0][1];
    ymin = data_array[1][1]; ymax = data_array[1][1];

    for( i = 1; i <= data_pairs; i++ ) {
        if( data_array[0][i] < xmin ) 
            xmin = data_array[0][i];
        if( data_array[0][i] > xmax ) 
            xmax = data_array[0][i];
        if( data_array[1][i] < ymin ) 
            ymin = data_array[1][i];
        if( data_array[1][i] > ymax ) 
            ymax = data_array[1][i];
    }
}
/**************************************************************************/
void DMDCurveFit::reporterror( char erstr[] )
{
    char key;

    if( debug == true ) {
        printf( "\n\n%s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
	printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "%s\n\n", erstr );
	printf( "Press <RETURN> to continue\n" );
	key = getc( stdin );
	return;
    }
    else {
	/*	if( write_results_file( false, erstr ) == false );
			{
			reporterror( "ERROR in writing results file" );
			exit(0);
			}*/
    }
}
/**************************************************************************/
long double DMDCurveFit::ldpow( long double n, unsigned p )
{
   /* Calculates n to the power (unsigned integer) p
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   /*--------------------------------------------------------*/
   /* This function returns a long double value whereas the  */
   /* standard function 'pow' returns a double value.        */
   /*                                                        */
   /* Also, this function gives smaller errors than 'pow'    */
   /* where only INTEGER powers are required.                */
   /*--------------------------------------------------------*/
     long double x = 1;
     unsigned i;

     for(i = 0; i < p; i++ ) x = x * n;

     return x;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double DMDCurveFit::polycurvefit( double data[2][MAXPAIRS] )
{
    //FILE *fp;
   // int i;
   // char s1[256];
   // char key;
    order = 1;//atoi( s1 );
    data_pairs = 10;
    xmin = xmax = ymin = ymax = 0;
    for (int i = 0; i < MAXPAIRS; i++) {
        data_array[0][i] = data[0][i];
        data_array[1][i] = data[1][i];
    }
     /* check min / max values for x and y
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    get_xy_max_min();

     /* check order v number of data pairs
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if( order > (data_pairs-1) ) {
        reporterror( "ERROR -- Specified order > (data_pairs - 1)." );
        exit(0);
    }

     /* do polynomial fit
    ~~~~~~~~~~~~~~~~~*/
    if( polynomfit() == false ) {
        reporterror( "ERROR in 'polynomfit' --  Data Error or Overflow ?" );
        exit(0);
    }
    else {
     /* write results file
        ~~~~~~~~~~~~~~~~~~*/
    /*    if( write_results_file( true, "" ) == false ) {
            reporterror( "ERROR in writing results file" );
            exit(0);
        }*/
    }

    return (double)polycoefs[1];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double DMDCurveFit::lslcurvefit( double data[2][MAXPAIRS] )
{
    real mA[MAXPAIRS*2];
    real v[MAXPAIRS];
/*
    mA[ 0]= 1;   mA[ 9]= 0.014;   v[0]= 6.1026;
    mA[ 1]= 1;   mA[10]= 0.021;   v[1]= 5.9216;
    mA[ 2]= 1;   mA[11]= 0.028;   v[2]= 5.8260;
    mA[ 3]= 1;   mA[12]= 0.035;   v[3]= 5.6348;
    mA[ 4]= 1;   mA[13]= 0.042;   v[4]= 5.5294;
    mA[ 5]= 1;   mA[14]= 0.049;   v[5]= 5.6348;
    mA[ 6]= 1;   mA[15]= 0.056;   v[6]= 5.5491;
    mA[ 7]= 1;   mA[16]= 0.063;   v[7]= 5.2149;
    mA[ 8]= 1;   mA[17]= 0.070;   v[8]= 5.3132;
*/
    for (int i = 0; i < MAXPAIRS; i++) {
        mA[i] = 1;
        mA[i + MAXPAIRS] = data[0][i];
        v[i] = data[1][i];
    }

    char trans = 'N';
    integer m=MAXPAIRS;
    integer n=2;
    integer nrhs = 1;
    integer lwork;
    real *mWork;
    integer info = 0;
  
    lwork=-1;
    mWork=(float *)malloc(1*sizeof(float));
    sgels_(&trans, &m, &n, &nrhs, mA, &m, v, &m, mWork, &lwork, &info);
    lwork = (int) mWork[0];
    free(mWork);
    mWork=(float *)malloc(lwork*sizeof(float));
    sgels_(&trans, &m, &n, &nrhs, mA, &m, v, &m, mWork, &lwork, &info);
    free(mWork);
    return (double)v[1];
}

#endif
