#ifndef _QtImageWriter_H
#define _QtImageWriter_H

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include "ItkFastIO.h"
#include <qfiledialog.h>

template< class TPixel >
class QtImageWriter
{
  public:
   typedef TPixel PixelType;
   typedef itk::Image<PixelType,3> ImageType;
   typedef itk::ImageFileWriter<ImageType> WriterType;
   
   void Write();
   void SetInput(typename ImageType::Pointer image);

   QtImageWriter();
   QtImageWriter(QString filename);
   ~QtImageWriter(); 

  protected:
   typename ImageType::Pointer m_image;
   QString m_filename;

};


#include "QtImageWriter.txx"


#endif
