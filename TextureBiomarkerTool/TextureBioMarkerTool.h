#define DARK_CALIB 100
#define BRIGHT_CALIB 1000
#define EPSILON 0.1
/*
//6 muscles
#define REC_FEM 1
#define SEMIT 4
#define CRAN_SART 5
#define BIC_FEM 2
#define ADDUCTOR 3
#define GRACILIS 6
#define SEMI_MEM 7
#define VAS_LAT 8
#define VAS_MED 10
#define VAS_INT 9
#define CAUD_SART 11
*/

//11 muscles
#define CRAN_SART 1
#define REC_FEM 2
#define BIC_FEM 3
#define GRACILIS 4
#define SEMI_MEM 5
#define ADDUCTOR 6
#define SEMIT 7
#define VAS_LAT 8
#define VAS_INT 9
#define VAS_MED 10
#define CAUD_SART 11

#define NUMBER_OF_MUSCLE 11
#define TRUE 1
#define FALSE 0
#define T2FIT
#define SEG_INTERPOLATE
#define false 0
#define true !false
#define progtitle    "POLYFIT"
#define progversion  "1.0.0"
#define progdate     "08.03.2003"
#define progauther   "K. Hough"
//#define MAXPAIRS 50
#define MAXPAIRS 9
#define MAXORDER 8      /* Can do higher, but who cares?   */
#define PI 3.1415926
#define ECHO_TIME 7    // unit ms 7 or 20
#define DATA_PAIRS 10    // unit ms 10 or 5
#define _POMPE
#define HISTOGRAM_TEXTURE
#define RUN_LENGTH_MATRIX
#define E_CONSTANT 2.71828
#define COOCURRENCE_MATRIX

//#include "DMDT2Fitting.h"
//#include "DMDCalib.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkTransformFileReader.h"
#include "itkAffineTransform.h"
#include "itkTransformFactory.h"
#include "itkResampleImageFilter.h"
#include "itkDivideImageFilter.h"
#include "TextureFeatureCal.h"
//#include "DMDCurveFit.h"
#include <itksys/Process.h>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>      /* need functions 'toupper' and tolower'   */
#include <float.h>      /* defines LDBL_MAX*/
#include <math.h>
#include <time.h>       /* for date marking of output data files   */


