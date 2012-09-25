#ifndef FTIO_H
#define FTIO_H

#include "itkImage.h"
#include "itkGiplImageIOFactory.h"
//#include "itkMetaImageIOFactory.h"
#include "itkAnalyzeImageIOFactory.h"
#include "itkGiplImageIO.h"
//#include "itkMetaImageIO.h"
#include "itkAnalyzeImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template< class TPixel >
class FTIO
{

  public:

  typedef TPixel PixelType;
  typedef itk::Image<PixelType,3> ImageType;

  typedef itk::GiplImageIO    GiplIOType;
 // typedef itk::MetaImageIO    MetaIOType;
  typedef itk::AnalyzeImageIO AnalyzeIOType;
  typedef typename ImageType::Pointer ImagePointer;

  
  void SetFilename(const char * name);
  void SetImage2Write(ImagePointer image);
  int Write( void );
  int Update( void );
  ImagePointer GetITKImage( void );

  FTIO();
  ~FTIO(); 

  protected:


  ImagePointer m_Image;
  const char * m_Filename;

};



#include "ItkFastIO.txx"


#endif
