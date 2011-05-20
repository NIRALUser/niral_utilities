#include "ImageRotation.h"
  
ImageRotation::ImageRotation()
{
	m_progressbar = 0;
}


ImageRotation::~ImageRotation()
{
}


void ImageRotation::SetInput(ImageType::Pointer inputimage)
{
	m_inputimage = inputimage;
}


void ImageRotation::SetDirection(int direction)
{
	m_direction = direction;
}


void ImageRotation::SetRotation(float angle)
{
	m_angle = -angle;
}

void ImageRotation::SetProgressBar(QProgressBar* progressbar)
{
	m_progressbar = progressbar;
}



void ImageRotation::Update()
{

	ImageType::SizeType m_size = m_inputimage->GetLargestPossibleRegion().GetSize();
	ImageType::SpacingType m_spacing = m_inputimage->GetSpacing();

	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(m_inputimage);
    duplicator->Update();
	m_outputimage = duplicator->GetOutput();

	int m_i,m_j,m_k;
	m_i = m_j = m_k = 0;

	if (m_direction == 0)
	{
		m_i = 1;
		m_j = 2;
		m_k = 0;
	}

	if (m_direction == 1)
	{
		m_i = 0;
		m_j = 2;
		m_k = 1;
	}

	if (m_direction == 2)
	{
		m_i = 0;
		m_j = 1;
		m_k = 2;
	}

	int xc = m_size[m_i]/2;
	int yc = m_size[m_j]/2;
	float xf,yf;
	int xi,yi;

	IteratorType m_itimage(m_outputimage,m_outputimage->GetLargestPossibleRegion());
	m_itimage.GoToBegin();
    unsigned short* cWinImData = new unsigned short[ m_size[m_i] * m_size[m_j] ];


	itk::Index<3> ind;

	if (m_progressbar)
	{
		m_progressbar->setProgress(0);
		m_progressbar->setTotalSteps(m_size[m_k]);
	}

	for (unsigned int m_slice=0;m_slice<m_size[m_k];m_slice++)
	{
		ind[m_k] = m_slice;
		for (unsigned int y=0;y<m_size[m_j];y++)
		{
			ind[m_j] = y;
			for (unsigned int x=0;x<m_size[m_i];x++)
			{	
				ind[m_i] = x;
				cWinImData[x+y*m_size[m_i]] = (unsigned short)m_inputimage->GetPixel(ind);
			}	
		}

		for (unsigned int y=0;y<m_size[m_j];y++)
		{
			ind[m_j] = y;
			for (unsigned int x=0;x<m_size[m_i];x++)
			{
				ind[m_i] = x;
				xf = cos(m_angle)*(x*m_spacing[m_i]-xc*m_spacing[m_i]) + sin(m_angle)*(y*m_spacing[m_j]-yc*m_spacing[m_j]) + xc*m_spacing[m_i];
				yf = -sin(m_angle)*(x*m_spacing[m_i]-xc*m_spacing[m_i]) + cos(m_angle)*(y*m_spacing[m_j]-yc*m_spacing[m_j]) + yc*m_spacing[m_j];	
				xf = xf/m_spacing[m_i];
				yf = yf/m_spacing[m_j];
				if ((xf<0) || (yf<0) || (xf>m_size[m_i]) || (yf>m_size[m_j]))
					m_outputimage->SetPixel(ind,0);
				else
				{
					xi = (int)xf;
					yi = (int)yf;
					m_outputimage->SetPixel(ind,(unsigned short)Bicubic_Interpol(cWinImData, xi, yi, m_size[m_i], m_size[m_j], xf-xi, yf-yi, 0));
				}
		}
		}
		if (m_progressbar)	m_progressbar->setProgress(m_slice);
	}
	delete [] cWinImData;
}

ImageRotation::ImageType::Pointer ImageRotation::GetOutput()
{
	return m_outputimage;
}


int ImageRotation::Bicubic_Interpol(unsigned short* m_image, int newx, int newy, int xsize, int ysize, float a, float b, int f)
{

  /* Perform the bicubic interpolation method */

   int m, n;
   float sum, cubsum;
   int index;
   int imageval;

   sum = 0;
 
   /* Find pixel value */ 
   index = ((newx + newy*xsize));   

   for (m = -1; m <=2; m++) {
       for (n = -1; n <= 2; n++) {
	 if ((newx+m >= 1) && (newx+m < xsize-1)) {
              if ((newy+n >=1) && (newy+n < ysize-1)) {  
		/* Find surrounding pixels values */ 
			imageval = m_image[index + (m + xsize*n)];
              } else {
		 imageval = 0;
              }
         } else {
              imageval=0;
         }
         /* Weight... */
         cubsum = CubicFunc(f, (float)(-a+m))*CubicFunc(f, (float)(-b+n));
         /* Add it to do a complete weighted average */
         sum = sum + imageval*cubsum;
      
       } 
   }

   return((int) sum);
}


/* Interpolation functions */

float ImageRotation::CubicBSpline(float x) {

  /* Return the value of the cubic B spline
     between -2 .. 2          */

   float y, y2, y3;
   float result;

   result=0;
 
   if (x < 0) {
      y = -x; 
   } else {
      y = x;
   }
   
   y2=y*y;
   y3=y2*y;

   if (y < 1.0) {
      result=(2.0/3.0 + 0.5 * y3 - y2);
   } else {
     if (y <= 2.0) {
         result=(1.0/6.0)*(2-y)*(2-y)*(2-y);
     }
   }
   return(result);
}

float ImageRotation::CubicInterpolation(float x) {

  /* Returns the value of an interpolated
     sine between -2 .. 2  */

   float a, y, y2, y3, result;
  
   a=0.1;  /* arbitrary value */
   result=0;
 
   if (x < 0) {
      y = -x;
   } else {
      y=x;
   }

   y2=y*y;
   y3=y2*y;

   if (y <= 1.0) {
      result=((a+2)*y3 - (a+3)*y2 + 1);
   } else {
     if (y <= 2.0) {
        result=(a*y3 - 5*a*y2 + 8*a*y - 4*a);
     }
   }
   
   return(result);
}

float ImageRotation::SquareFunction(float x) {

  /* Returns 1 between -0.5 .. 0.5 */

  float result;

  if ((x >= -0.5) && (x <= 0.5)) {
     result = 1.0;
  } else {
     result = 0;
  }
  return(result);
}

float ImageRotation::TriangleFunction(float x) {

  /* Returns the triangle function 
     between -1.0 .. 1.0            */ 

   float result;

   result = 0;

   if ((x >= -1) && (x <= 0)){ 
      result = x + 1;
   } else {
      if ((x > 0) && (x <= 1)) result = 1 - x;
   }

   return(result);
}

float ImageRotation::CubicFunc(int f, float x) {

  /* Returns the different cubic functions */
  /* 0 = B-Spline, 1 = Interpolated sine   */

   float result = 0;

   switch (f) {
        case 0: result = CubicBSpline(x);
                break;
        case 1: result = CubicInterpolation(x);
                break;
   }
   
   return(result);
}

float ImageRotation::BilinearFunc(int f, float x) {
    
  /* Returns the different linear functions */
  /* 0 = Square ; 1 = Triangle              */

   float result = 0;

   switch (f) {
        case 0: result = SquareFunction(x);
                break;
        case 1: result = TriangleFunction(x);
                break;
   }
   
   return(result);
}
