#ifndef __ImageCaster_h_
#define __ImageCaster_h_

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"

template< class TPixel >
class ImageCaster
{
public:
	ImageCaster();
	~ImageCaster();
	
	typedef TPixel PixelType;
	typedef itk::Image<PixelType,3> ImageType;

	typedef itk::Image<float,3> ImageFloatType;
	typedef itk::Image<unsigned short,3> ImageUShortType;

	typedef ImageType::Pointer  ImagePointer;
	typedef itk::ImageDuplicator<ImageType>		DuplicatorType;
    typedef itk::ImageRegionIterator<ImageType> IteratorType;

	ImageFloatType::Pointer ToFloat(ImagePointer image);
 	ImageUShortType::Pointer ToUShort(ImagePointer image);

};

#include "ImageCaster.txx"

#endif // __ImageCaster_h_
