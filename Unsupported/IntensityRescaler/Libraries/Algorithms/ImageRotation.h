#ifndef __ImageRotation_h_
#define __ImageRotation_h_

#include "ItkFastIO.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"

#include <qprogressbar.h>

class ImageRotation
{
public:
	ImageRotation();
	~ImageRotation();
	
	typedef unsigned short						PixelType;
	typedef itk::Image< PixelType,  3 >			ImageType;
	typedef itk::ImageFileReader< ImageType  >  ReaderType;
	typedef itk::ImageFileWriter< ImageType  >  WriterType;
	typedef ImageType::Pointer					ImagePointer;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	typedef itk::ImageDuplicator<ImageType>		DuplicatorType;

	void SetProgressBar(QProgressBar* progressbar);
	void SetInput(ImageType::Pointer inputimage);
	void SetDirection(int direction);
	void SetRotation(float angle);
	void Update();
	ImageType::Pointer GetOutput();

protected:
	float BilinearFunc(int f, float x);
	float CubicFunc(int f, float x);
	float TriangleFunction(float x) ;
	float SquareFunction(float x);
	float CubicInterpolation(float x);
	float CubicBSpline(float x);
	int Bicubic_Interpol(unsigned short* m_image, int newx, int newy, int xsize, int ysize, float a, float b, int f);


private:
	ImageType::Pointer m_inputimage;
	ImageType::Pointer m_outputimage;
	float m_angle;
	int m_direction;
	QProgressBar* m_progressbar;

};


#endif // __ImageRotation_h_
