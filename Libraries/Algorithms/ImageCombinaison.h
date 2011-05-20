#ifndef __ImageCombinaison_h_
#define __ImageCombinaison_h_


#include <itkImage.h>
#include <itkImageDuplicator.h>
#include "itkImageRegionIterator.h"

class ImageCombinaison
{
public:
  ImageCombinaison();
  ~ImageCombinaison();
  
  typedef float								PixelType;
  typedef itk::Image< PixelType,  3 >			ImageType;
  typedef ImageType::Pointer					ImagePointer;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;
  typedef itk::ImageDuplicator<ImageType>		DuplicatorType;

  ImagePointer Multiply(ImagePointer image1, ImagePointer image2);
  ImagePointer Minus(float val,ImagePointer image1);
  ImagePointer Minus(ImagePointer image1,float val);
  ImagePointer Divide(ImagePointer image1,float val);
  ImagePointer DuplicateImage(ImagePointer image);


};

#endif
