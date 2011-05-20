#ifndef _QtImagerWriter_TXX
#define _QtImagerWriter_TXX

template< class TPixel >
QtImageWriter<TPixel>
::QtImageWriter()
{
}

template< class TPixel >
QtImageWriter<TPixel>
::QtImageWriter(QString filename)
{
   m_filename = filename;
}


template< class TPixel >
QtImageWriter<TPixel>
::~QtImageWriter()
{
}

template< class TPixel >
void
QtImageWriter<TPixel>
::Write()
{

  typename WriterType::Pointer writer = WriterType::New();
  QString s = QFileDialog::getSaveFileName(m_filename,"Images (*.*)", 0, "Open file dialog","Chose an image filename" );
  if (s.isEmpty())
   return;

  m_filename = s.latin1();
  writer->SetFileName( (const char *) m_filename);
  writer->SetInput(m_image);

  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
   return;
    }
}

template< class TPixel >
void
QtImageWriter<TPixel>
::SetInput(typename QtImageWriter<TPixel>::ImageType::Pointer image)
{
   m_image = image;
}

#endif
