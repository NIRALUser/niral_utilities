#ifndef _QtImagerLoader_H
#define _QtImagerLoader_H

#include "itkImage.h"
#include "ItkFastIO.h"
#include <qfiledialog.h>

template< class TPixel >
class QtImageLoader
{

  public:
	typedef TPixel PixelType;
	typedef itk::Image<PixelType,3> ImageType;
  	typedef itk::ImageFileReader<ImageType>   ReaderType;
  	typedef typename ImageType::Pointer ImagePointer;

	ImagePointer GetOutput( void );
	QString GetFileName();
	QString GetFilePath();

	QtImageLoader();
	~QtImageLoader(); 

  protected:
	ImagePointer m_Image;
	QString m_filename;

};


#include "QtImageLoader.txx"


#endif
