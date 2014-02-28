/*****************************************************
 *                 polyfit19.c                       *
 *                                                   *
 *                 (for Linux)                       *
 *                                                   *
 *              K.Hough 08.03.2003                   *
 *                                                   *
 *****************************************************/
/*
 Mods/debug:

 -- Added description of data file format in HELP
    (09.01.2003)


 Notes:

 1. For compilation using GNU gcc or g++

    gcc produces smaller (possibly faster?) compiled
    code, but MUST be used with the -lm flag to include
    maths libary.

 2. This version is clear of warnings with the -Wall
    option set

 3. Can also use the -O option to optimise code
    -O3 produces maximum optimisation for speed

    Possible optimal command lines for compilation are:

           gcc -Wall -O3 polyfitxx.c -o polyfit -lm

       or: g++ -Wall -O3 polyfitxx.c -o polyfit

    These command lines result in an executable file of about
    33kB. Adding the -static flag results in a standalone
    executable file which does not rely on external libraries,
    but will be approx 1.6MB.

 4. This version applies scaling to very large/very small data
    sets before polynomial fits are done. After fit is done,
    coefficients are rescaled to give correct results. This
    reduces the magnitudes of intermediate calculated values,
    and hence the chance for over/under flow.

 5. Function 'ldpow' is used instead of the standard 'powl'
    Refer function 'ldpow' in this programme for information

 6. For programme details, refer to function 'helpscreen'

*/
/*****************************************************/

        #include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <ctype.h>      /* need functions 'toupper' and tolower'   */
        #include <float.h>      /* defines LDBL_MAX                        */
	#include <math.h>
	#include <time.h>       /* for date marking of output data files   */

     /* set 'defines'
	~~~~~~~~~~~~*/
	#define false 0
	#define true !false
	#define progtitle    "POLYFIT"
	#define progversion  "1.0.0"
	#define progdate     "08.03.2003"
	#define progauther   "K. Hough"

	#define MAXPAIRS 50
	#define MAXORDER 8      /* Can do higher, but who cares?           */

     /* declare functions
	~~~~~~~~~~~~~~~~~*/
	int read_data_file( void );
	int write_results_file( int result, char error_report[] );
	int polynomfit( void );
	void datasort ( void );
	float poly_correlcoef( void );
	void get_xy_max_min ( void );
	void reporterror( char erstr[] );
	char *stripleadingspaces( char *s );
	char *striptrailingspaces( char *s );
	char *stripallspaces( char *s );
	char *strupr( char s[] );
	char *strlwr( char s[] );
	int instr( char *sin, char *stest );
	void clearscr( void );
	char *datestring( void );
	char *timestring( void );
	long double ldpow( long double n, unsigned p );
	void helpscreen( void );

     /* declare array for data
	~~~~~~~~~~~~~~~~~~~~~~*/
	double data_array[2][MAXPAIRS+1];

     /* other global variables
	~~~~~~~~~~~~~~~~~~~~~~*/
	char data_in_file[256];
	char results_out_file[256];

	int data_pairs, order;
	float correl_coef;
	long double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
	long double polycoefs[MAXORDER+1];

	int debug = true;   /* enable display of results on screen
			       -- reset via command line to true or false */
	char strbuf[256];

/*-------------------------------------------------------------------------*/
int main( int argc, char *argv[] )
{
	FILE *fp;
	int i;
	char s1[256];
	char key;

     /* check for correct number of command line arguments
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if( argc == 1 || strcmp( argv[1], "--help" ) == 0 )
		{
		helpscreen();
		clearscr();
		exit(0);
		}

	if( argc < 3 )
		{
		reporterror( "COMMAND LINE ERROR -- Incorrect command line." );
		printf( "\nCommand line arguments:\n" );
		for( i = 0; i < argc; i++ )
			printf( "\t\t\targv[%d] : %s\n" , i, argv[i] );
		printf( "\n\n" );
		exit(0);
		}

     /* place argv's into variables
	~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	strcpy( data_in_file, argv[1] );
	strcpy( s1, argv[2] ); order = atoi( s1 );

        strcpy( s1, argv[argc - 1] );
	if( strcmp( stripallspaces(strupr( s1 )), "D" ) == 0 || argc == 3 )
	        debug = true;
             else
	        debug = false;

        if( argc >= 4 )
		strcpy( results_out_file, argv[3] );
	     else
		strcpy( results_out_file, "" );

     /* check for existance of data file
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if( (fp = fopen( data_in_file, "r" )) == NULL )
		{
		sprintf( s1,"COMMAND LINE ERROR -- Data file '%s' does not exist\n", data_in_file );
		reporterror( s1 );
		helpscreen();
		exit(0);
		}
	     else
		fclose( fp );

     /* check that specified 'order' is within range
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if( order < 1 || order > MAXORDER )
		{
		sprintf( s1,"COMMAND LINE ERROR -- Specified order '%d' for polnomial fit is out of range.\n\n", order );
		reporterror( s1 );
		helpscreen();
		exit(0);
		}

     /* read datafile
	~~~~~~~~~~~~~*/
	if( read_data_file() == false )
		{
		reporterror( "ERROR in reading data file" );
		exit(0);
		}

     /* check min / max values for x and y
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	get_xy_max_min();

     /* check order v number of data pairs
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        if( order > (data_pairs-1) )
	        {
		reporterror( "ERROR -- Specified order > (data_pairs - 1)." );
		exit(0);
		}

     /* do polynomial fit
        ~~~~~~~~~~~~~~~~~*/
        if( polynomfit() == false )
	     {
	     reporterror( "ERROR in 'polynomfit' --  Data Error or Overflow ?" );
	     exit(0);
	     }
          else
	     {
	     /* write results file
	        ~~~~~~~~~~~~~~~~~~*/
	     if( write_results_file( true, "" ) == false )
		     {
		     reporterror( "ERROR in writing results file" );
		     exit(0);
		     }
	     }

     /* if set for debug
	~~~~~~~~~~~~~~~~*/
	if( debug == true )
		{
	/*	clearscr();
		printf( "%s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
		printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
		printf( "Information read from command line:\n" );
		printf( "     Name/path of Data file: %s\n", data_in_file );
		printf( "     Order of fit requested: %d\n", order );
		if( strcmp( results_out_file, "" ) != 0 )
			printf( " Name/path for results file: %s\n", results_out_file );
		printf( "\n Data pairs read from data file : %d", data_pairs );
		if( data_pairs == MAXPAIRS ) printf( "\tie max no. allowed\n" ); else printf( "\n" );
		printf( "\nResults:\n" );
		printf( "     Correlation coeficient = %+7.6f\n\n", correl_coef );
		printf( "     Polynomial coefficients:\n" );
		for(i = 0; i <= order; i++ )
			{
			sprintf( s1, "%e", (float)polycoefs[i] );
		        if( instr( s1, "-0.000000" ) == 0 ) s1[0] = '+';
			printf( "           coef[%d] = %s\n", i, s1 );
			}
		printf( "\nPress <RETURN> to continue\n" );
		key = getchar();*/
		}
	return 0;
}
/*-------------------------------------------------------------------------*/
int polynomfit( void )
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
		if( data_array[0][i] < xlow && data_array[0][i]	!= 0 ) xlow = data_array[0][i];
		if( data_array[1][i] < ylow && data_array[1][i]	!= 0 ) ylow = data_array[1][i];
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
	for( j = 1; j <= data_pairs; j++ )
		{
		for( i = 1; i <= order; i++ )
			{
			B[i] = B[i] + data_array[1][j] * yscale * ldpow( data_array[0][j] * xscale, i );
			if( B[i] == LDBL_MAX ) return false;
			for( k = 1; k <= order; k++ )
			        {
			        a[i][k] = a[i][k] + ldpow( data_array[0][j] * xscale, (i + k) );
			        if( a[i][k] == LDBL_MAX ) return false;
			        }
			S[i] = S[i] + ldpow( data_array[0][j] * xscale, i );
			if( S[i] == LDBL_MAX ) return false;
			}
			Y1 = Y1 + data_array[1][j] * yscale;
			if( Y1 == LDBL_MAX ) return false;
		}
	/*-----------------------------------------------------------------*/
	for( i = 1; i<= order; i++ )
		{
		for( j = 1; j <= order; j++ )
			{
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
	for( k = 1; k <= order - 1; k++ )
		{
		i = order - k; S1 = 0;
		for( j = 1; j <= order; j++ )
			{
			S1 = S1 + a[i][j] * polycoefs[j];
			if( S1 == LDBL_MAX ) return false;
			}
		polycoefs[i] = B[i] - S1;
		}
	/*-----------------------------------------------------------------*/
	S1 = 0;
	for( i = 1; i <= order; i++ )
                {
                S1 = S1 + polycoefs[i] * S[i] / (long double)data_pairs;
	        if( S1 == LDBL_MAX ) return false;
		}
	polycoefs[0] = (Y1 / (long double)data_pairs - S1);
	/*-----------------------------------------------------------------*/
	/* zero all coeficient values smaller than +/- .000001 (avoids -0)
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	for( i = 0; i <= order; i++ )
		if( fabsl( polycoefs[i] * 1000000 ) < 1 ) polycoefs[i] = 0;

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
void datasort( void )
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
	for( outer = 0; outer <= data_pairs; outer++ )
		{
		swapflag = 0;
		for( inner = 0; inner <= last; inner++ )
			{
			if( data_array[0][inner] > data_array[0][inner + 1] )
				{
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
			if( swapflag == 0 ) outer = data_pairs;
		}
}
/**************************************************************************/
float poly_correlcoef( void )
{
	int i, j;
	long double E[MAXORDER+1];
	long double ycalc[MAXPAIRS+1];
	long double my, mycalc, yc, SAA, SBB, SAB, A, B;
	long double SY, SYCALC, SYYCALC;

	for( i = 0; i <= order; i++ )
		E[i] = polycoefs[i];

	my = 0; mycalc = 0;
	for( i = 1; i <= data_pairs; i++ )
		{
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
	for( i = 1; i <= data_pairs; i++ )
		{
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
void get_xy_max_min( void )
{
	int i;

	xmin = data_array[0][1]; xmax = data_array[0][1];
	ymin = data_array[1][1]; ymax = data_array[1][1];

	for( i = 1; i <= data_pairs; i++ )
		{
		if( data_array[0][i] < xmin ) xmin = data_array[0][i];
		if( data_array[0][i] > xmax ) xmax = data_array[0][i];
		if( data_array[1][i] < ymin ) ymin = data_array[1][i];
		if( data_array[1][i] > ymax ) ymax = data_array[1][i];
		}
}
/**************************************************************************/
void reporterror( char erstr[] )
{
	char key;

	if( debug == true )
		{
		printf( "\n\n%s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
		printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
		printf( "%s\n\n", erstr );
		printf( "Press <RETURN> to continue\n" );
		key = getc( stdin );
		return;
		}
	     else
		{
		if( write_results_file( false, erstr ) == false );
			{
			reporterror( "ERROR in writing results file" );
			exit(0);
			}
		}
}
/**************************************************************************/
int read_data_file( void )
{
	FILE *fp;
	char s1[256], s2[256];
	int i, j;

	data_pairs = 0;

	if( (fp = fopen( data_in_file, "r" )) == NULL )
		return false;
	    else
		{
		/* read data file */
		do
			{
			/* read a line from file
			~~~~~~~~~~~~~~~~~~~~~~~~*/
rdloop1:                fgets( s1, 256, fp);

			/* catch end of file marker
			~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			if( feof( fp )) break;

			/* remove last character from line( chr 0A ?? )
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			s1[ strlen( s1 )-1] = '\0';

			/* jump over all blank lines and REMs
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			if( strlen( s1 ) == 0 || s1[0] == ';' ) goto rdloop1;

			/* assume a valid line
			~~~~~~~~~~~~~~~~~~~~~~*/
			data_pairs++;

			stripallspaces( s1 );

			/* get x value
			~~~~~~~~~~~~~~*/
			for( i = 0; i < instr( s1, "," ); i++ )
				s2[i] = s1[i];
			s2[i] = '\0';
                        data_array[0][data_pairs] = (long double)atof( s2 );

                        /*#############################################################*/
			/* NB function 'atol' doesn't seem to accept fractional values */
			/*    hence use of function 'atof' and cast of '(long double)' */
                        /*#############################################################*/

			/* get y value
			~~~~~~~~~~~~~~*/
			j = 0;
			for( i = instr( s1, "," ) + 1; i < (int)(strlen( s1 )); i++ )
				{ s2[j] = s1[i]; j++; }
			s2[j] = '\0';
			stripleadingspaces( s2 );
			data_array[1][data_pairs] = (long double)atof( s2 );
			}
		while ( data_pairs < MAXPAIRS );

		/* close file
		~~~~~~~~~~~~~*/
		fclose( fp );
		}
return true;
}
/**************************************************************************/
int write_results_file( int result, char error_report[] )
{
     FILE *fp;
     char s1[256], s2[256];
     int i;

     if( strcmp( results_out_file, "" ) == 0 ) return true;

     if( (fp = fopen( results_out_file,"w" )) != NULL )
	 {
	/* strcpy( s1, ";Programme : " );
	 strcat( s1, progtitle );
	 strcat( s1, "  Version " );
	 strcat( s1, progversion );
	 fputs( s1, fp ); fputc( '\n', fp );
	 sprintf( s1,";Run Date : %s (dd/mm/yy)", datestring() ); fputs( s1, fp ); fputc( '\n', fp );
         sprintf( s1,";Run Time : %s (hh:mm:ss)", timestring() ); fputs( s1, fp ); fputc( '\n', fp );
 	 if( result == true ) strcpy( s2, "SUCCESSFUL" ); else strcpy( s2, "FAILED" );
	 sprintf( s1,";RESULT : %s", s2 ); fputs( s1, fp ); fputc( '\n', fp );
	 if( result == false )
		{
		sprintf( s1,";%s", error_report ); fputs( s1, fp ); fputc( '\n', fp );
		fclose( fp );*/
	     /* terminate here
		~~~~~~~~~~~~~~*/
/*		exit(0);
		}
	 fputs( ";Data Pairs", fp ); fputc( '\n', fp );
	 sprintf( s1, "%d", data_pairs ); fputs( s1, fp ); fputc( '\n', fp );
	 fputs( ";Order", fp ); fputc( '\n', fp );
	 sprintf( s1, "%d", order ); fputs( s1, fp ); fputc( '\n', fp );
	 fputs( ";Correlation coefficient", fp ); fputc( '\n', fp );
	 sprintf( s1, "%+7.6f", correl_coef ); fputs( s1, fp ); fputc( '\n', fp );

	 fputs( ";Polynomial coefficients", fp ); fputc( '\n', fp );*/
	 for( i = 0; i <= order; i++ ) {
		sprintf( s1, "%f", (float)polycoefs[i] );
	//	if( instr( s1, "-0.000000" ) == 0 ) s1[0] = '+';
		fputs( s1, fp ); fputc( '\n', fp );
	 }
/*	 fputs( ";xmin, xmax, ymin, ymax", fp ); fputc( '\n', fp );
	 sprintf( s1, "%+LE\n", xmin ); fputs( s1, fp );
	 sprintf( s1, "%+LE\n", xmax ); fputs( s1, fp );
	 sprintf( s1, "%+LE\n", ymin ); fputs( s1, fp );
	 sprintf( s1, "%+LE\n", ymax ); fputs( s1, fp );
*/
	 fclose( fp );
	 return true;
	 }
      else
	 return false;
}
/**************************************************************************/
char *stripleadingspaces( char *s )
{
    /* strips ONLY leading spaces from character strings
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       char s2[256];
       int c = 0, x, l = strlen( s );

       while( s[c] == ' ' ) c++;
       for( x = c; x <= ( l + 1 ); x++ ) s2[x-c] = s[x];
       strcpy( s, s2 );

       return &s[0];
}
/**************************************************************************/
char *striptrailingspaces( char *s )
{
    /* strips ONLY trailing spaces from character strings
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       int c = strlen( s );

       while( c > 0 && s[c-1] == ' ' ) c--;
       s[c] = '\0';

       return &s[0];
}
/**************************************************************************/
char *stripallspaces( char s[] )
{
    /* strips ALL spaces from character strings
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       char s1[256];
       int c0 = 0, c1 = 0;

       while ( s[c0] != '\0' )
           {
           while( s[c0] == ' ' ) c0++;
           s1[c1] = s[c0];
           c0++;
           c1++;
           }
       s1[c1] = '\0';
       strcpy( s, s1 );

       return &s[0];
}
/**************************************************************************/
char *strupr( char s[] )
{
    /* convert string to upper case
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       int c = 0;

       while(  s[c] != '\0' ) { s[c] = toupper(s[c]); c++; }

       return &s[0];
}
/**************************************************************************/
char *strlwr( char s[] )
{
    /* convert string to lower case
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       int c = 0;

       while(  s[c] != '\0' ) { s[c] = tolower(s[c]); c++; }

       return &s[0];
}
/**************************************************************************/
int instr( char *sin, char *stest )
{
    /* find occurrence of a substring within a string
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       char buff[256];
       int lin = strlen(sin), ltest = strlen( stest );
       int c, n;

       for( c = 0; c <= ( lin - ltest ); c++ )
	  {
	  for(n = 0; n < ltest; n++ )
	      {
	      buff[n] = sin[c + n];
	      buff[n + 1] = '\0';
	      }
	  if( strcmp(stest, buff) == 0 ) return c;
	  }

       return -1;
}
/**************************************************************************/
void clearscr( void )
{
    /* clear screen using standard ANSI escape codes
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       printf("\033[2J");
       printf("\033[H");
}
/**************************************************************************/
char *datestring( void )
{
    /* find date and build formatted date string as dd/mm/yy
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
       int day, month, year;
       char daystr[3] = "", monthstr[3] = "", yearstr[5] = "";
       time_t timer;
       struct tm *my_time;
       time( &timer );
       my_time = localtime( &timer );
       day = my_time->tm_mday;
       if(day < 10) sprintf( daystr,"0%d", day); else sprintf( daystr, "%d", day );
       month = my_time->tm_mon + 1;
       if(month < 10) sprintf( monthstr,"0%d", month); else sprintf( monthstr, "%d", month );
       year = my_time->tm_year + 1900;
       sprintf( yearstr,"%d", year); yearstr[0] = yearstr[2]; yearstr[1] = yearstr[3]; yearstr[2]='\0';
       sprintf( strbuf,"%s/%s/%s", daystr, monthstr, yearstr );

       return &strbuf[0];
}
/**************************************************************************/
char *timestring( void )
{
   /* get time and build formatted time string as hh:mm:ss
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      char hours[3], minutes[3], seconds[3];
      time_t timer;
      struct tm *my_time;
      time( &timer );
      my_time = localtime( &timer );
      if( my_time->tm_hour < 10 ) sprintf( hours, "0%d", my_time->tm_hour ); else sprintf( hours, "%d", my_time->tm_hour );
      if( my_time->tm_min < 10 ) sprintf( minutes, "0%d", my_time->tm_min ); else sprintf( minutes, "%d", my_time->tm_min );
      if( my_time->tm_sec < 10 ) sprintf( seconds, "0%d", my_time->tm_sec ); else sprintf( seconds, "%d", my_time->tm_sec );
      sprintf( strbuf,"%s:%s:%s", hours, minutes, seconds );

      return &strbuf[0];
}
/**************************************************************************/
long double ldpow( long double n, unsigned p )
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
/**************************************************************************/
void helpscreen( void )
{
   /* help/description pages
      ~~~~~~~~~~~~~~~~~~~~~~*/
        char key, s1[256];
        int i;

        clearscr();
	printf( "                %s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
	printf( "                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "%s provides a means of deriving polynomial curve fits for data sets\n", progtitle );
	printf( "comprising up to %d data pairs. Polynomial fits of up to %dth order can be done.\n", MAXPAIRS, MAXORDER );
	printf( "THE PROGRAMME IS COMMAND LINE DRIVEN and there are TWO NECESSARY parts,\n" );
	printf( "separated by at least one space:\n" );
	printf( "\n  1. Full file path/name of data file    2. Order of polynomial fit required\n" );
	strcpy( s1, progtitle );
	printf( "\n                  eg.     %s  datafile.dat  3\n", strlwr( s1 ) );
        printf( "\nThe programme will then carry out the specified order of fit and display the\n" );
	printf( "results on screen.\n" );
	printf( "\nIn addition, the full file path/name for a results file can be included. By\n" );
	printf( "default, results will then be written to this file and will not be shown on\n" );
	printf( "screen. This mode of operation allows %s to be used from within other\n" , progtitle );
	printf( "programmes.\n" );
	printf( "            eg.     %s  datafile.dat  3  resultsfile.dat\n", strlwr( s1 ) );
	printf( "\n\n\n\n             Press <RETURN> to continue or <E><RETURN> to exit\n" );

	key = getc( stdin );
	if( key == 'e' || key == 'E' ) return;

        clearscr();
	printf( "                %s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
	printf( "                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "Adding a final 'd' to the command line, causes results again to be displayed\n" );
	printf( "on screen (and written to file).\n" );
	printf( "\n            eg.     %s  datafile.dat  3  resultsfile.dat  d\n", strlwr( s1 ) );
	printf( "\nThis mode of operation is intended for use while developing/debugging\n" );
	printf( "programmes which use %s.\n", progtitle );
 	printf( "\n  NB: %s does not know about 'ncurses' so screen output will not\n", progtitle );
	printf( "      be correctly formated if called from within a programme running\n" );
	printf( "      'ncurses'. Do not use the 'd' flag if running 'ncurses'.\n" );
        printf( "\nIf %s is called without a command line, or with '--help', these help\n", progtitle );
	printf( "pages will be shown.\n" );
	printf( "\nPolynomial curve fitting entails the use of a number of data arrays and\n" );
	printf( "possibly the calculation of very large numbers, especially for higher order\n" );
	printf( "fits. The maximum order of fit has therefore been limited to 8.\n" );
	printf( "If the data file contains more than %d data pairs, ONLY THE FIRST %d\n", MAXPAIRS, MAXPAIRS );
	printf( "PAIRS WILL BE READ.\n\n" );

	printf( "             Press <RETURN> to continue or <E><RETURN> to exit\n" );

	key = getc( stdin );
	if( key == 'e' || key == 'E' ) return;

	clearscr();
	printf( "                %s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
	printf( "                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "%s has been written using the 'C' language and makes use of 'long double'\n", progtitle );
	printf( "variables so as accommodate large numbers that can be generated during\n" );
	printf( "polynomial fitting. When appropriate, the magnitude of these numbers is\n" );
	printf( "managed by a scaling/rescaling process before/after fitting. Even so, errors\n" );
	printf( "may occur in the calculated coefficients for higher order polynomial fits\n" );
	printf( "and/or large data values. At this stage, some judgement and testing of\n" );
	printf( "results may be needed.\n" );
	printf( "\nFor guidance, the following limits were found to be OK prior to including\n" );
	printf( "the scaling process:\t (Note to myself: Further testing needed here!)\n" );
	printf( "\n                 1st order fit  --  data values +/- ???????\n" );
	printf( "     2nd, 3rd or 4th order fit  --  data values +/-100,000\n" );
	printf( "                 5th order fit  --  data values +/-10,000\n" );
	printf( "                 6th order fit  --  data values +/-1,000\n" );
	printf( "                 7th order fit  --  data values +/-100\n" );
	printf( "                 8th order fit  --  data values +/-100\n" );
        for( i = 0; i < 4; i++ ) printf( "\n" );
	printf( "             Press <RETURN> to continue or <E><RETURN> to exit\n" );

	key = getc( stdin );
	if( key == 'e' || key == 'E' ) return;

	clearscr();
	printf( "                %s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
	printf( "                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "Format for data file:\n" );
	printf( "\nData must be provided as comma separated data pairs as in the example below. No other\n" );
	printf( "information is needed and should not be included.\n ");
	printf( "\n" );
	printf( "\t\t0,\t2\n" );
	printf( "\t\t1,\t3\n" );
	printf( "\t\t2,\t5\n" );
	printf( "\t\t3,\t7\n" );
	printf( "\t\t4,\t9\n" );
	printf( "\t\t5,\t11\n" );
	printf( "\t\t6,\t13\n" );

	for( i = 0; i < 9; i++ ) printf( "\n" );
	printf( "             Press <RETURN> to continue or <E><RETURN> to exit\n" );

	key = getc( stdin );
	if( key == 'e' || key == 'E' ) return;



	clearscr();
	printf( "                %s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
	printf( "                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "Example of (SUCCESSFUL) results file:\n" );
	printf( "\n       ;Programme : %s  Version %s\n", progtitle, progversion );
	printf( "       ;Run Date : %s (dd/mm/yy)\n", datestring() );
	printf( "       ;Run Time : %s (hh:mm:ss)\n", timestring() );
	printf( "       ;RESULT : SUCCESSFUL           --------------------------------------\n" );
	printf( "       ;Data Pairs                   | Notes:                               |\n" );
	printf( "       10                            | 1. All non data lines start with a   |\n" );
	printf( "       ;Order                        |    semi-colon                        |\n" );
	printf( "       3                             | 2. These results show that a 3rd     |\n" );
	printf( "       ;Correlation coefficient      |    order fit has been done on data   |\n" );
	printf( "       +1.000000                     |    that is clearly 2nd order         |\n" );
	printf( "       ;Polynomial coefficients      |                                      |\n" );
	printf( "       +0.000000E+00                 | <-- coef(0)                          |\n" );
	printf( "       +0.000000E+00                 | <-- coef(1)                          |\n" );
	printf( "       +1.000000E+00                 | <-- coef(2)   2nd order coefficient  |\n" );
	printf( "       +0.000000E+00                 | <-- coef(3)                          |\n" );
	printf( "       ;xmin, xmax, ymin, ymax        --------------------------------------\n" );
	printf( "       +0.000000E+00\n" );
        printf( "       +1.000000E+01\n" );
	printf( "       +0.000000E+00\n" );
	printf( "       +1.000000E+02" );

	printf( "        Press <RETURN> to continue or <E><RETURN> to exit\n" );

	key = getc( stdin );
	if( key == 'e' || key == 'E' ) return;

	clearscr();
	printf( "                %s  version %s  %s  %s\n", progtitle, progversion, progdate, progauther );
	printf( "                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "Example of (FAILED) results file:\n" );
	printf( "\n       ;Programme : %s  Version %s\n", progtitle, progversion );
	printf( "       ;Run Date : %s (dd/mm/yy)\n", datestring() );
	printf( "       ;Run Time : %s (hh:mm:ss)\n", timestring() );
	printf( "       ;RESULT : FAILED\n" );
	printf( "       ;ERROR -- Specified order > (data_pairs - 1)" );
	for( i = 0; i < 14; i++ ) printf( "\n" );
	printf( "\n                           Press <RETURN> to end\n" );

	key = getc( stdin );

	return;
}
/**************************************************************************/
