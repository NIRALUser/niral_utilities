#ifndef ImageCaster_TXX
#define ImageCaster_TXX

#include "ImageCaster.h"


template< class TPixel >
ImageCaster<TPixel>
::ImageCaster()
{
}

template< class TPixel >
ImageCaster<TPixel>
::~ImageCaster()
{
}

template< class TPixel >
ImageCaster<TPixel>::ImageFloatType::Pointer
ImageCaster<TPixel>
::ToFloat(ImageType::Pointer image)
{
	ImageFloatType::Pointer m_newimage = ImageFloatType::New();

    m_newimage->SetRegions( image->GetLargestPossibleRegion() );
    m_newimage->Allocate();
    m_newimage->SetOrigin(image->GetOrigin());
    m_newimage->SetSpacing(image->GetSpacing());

	typedef itk::ImageRegionIterator<ImageFloatType> IteratorFloatType;
	IteratorFloatType itS(m_newimage,m_newimage->GetLargestPossibleRegion());
	IteratorType itSorig(image,image->GetLargestPossibleRegion());


	itS.GoToBegin();
	itSorig.GoToBegin();

	while (!itS.IsAtEnd())
	{
		itS.Set(itSorig.Get());
		++itS;
		++itSorig;
	}

	return m_newimage;
}

template< class TPixel >
ImageCaster<TPixel>::ImageUShortType::Pointer
ImageCaster<TPixel>
::ToUShort(ImageType::Pointer image)
{
	ImageUShortType::Pointer m_newimage = ImageUShortType::New();

    m_newimage->SetRegions( image->GetLargestPossibleRegion() );
    m_newimage->Allocate();
    m_newimage->SetOrigin(image->GetOrigin());
    m_newimage->SetSpacing(image->GetSpacing());

	typedef itk::ImageRegionIterator<ImageUShortType> IteratorFloatType;
	IteratorFloatType itS(m_newimage,m_newimage->GetLargestPossibleRegion());
	IteratorType itSorig(image,image->GetLargestPossibleRegion());


	itS.GoToBegin();
	itSorig.GoToBegin();

	while (!itS.IsAtEnd())
	{
		itS.Set((unsigned short)itSorig.Get());
		++itS;
		++itSorig;
	}

	return m_newimage;
}

#endif
