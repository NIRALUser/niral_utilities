#ifndef _QtImagerLoader_TXX
#define _QtImagerLoader_TXX

template< class TPixel >
QtImageLoader<TPixel>
::QtImageLoader()
{
}

template< class TPixel >
QtImageLoader<TPixel>
::~QtImageLoader()
{
}

template< class TPixel >
typename QtImageLoader<TPixel>::ImagePointer
QtImageLoader<TPixel>
::GetOutput( void )
{

  typename ReaderType::Pointer reader = ReaderType::New();
  QString s = QFileDialog::getOpenFileName(QString::null,"Images (*.*)", 0, "Open file dialog","Chose an image filename" );
  if (s.isEmpty())
	return NULL;

  m_filename = s.latin1();
  reader->SetFileName((const char *) m_filename );

  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
	return 0;
    }
 
	QDir::setCurrent(s.mid(0,s.findRev("/")));

	return reader->GetOutput();
}

template< class TPixel >
QString
QtImageLoader<TPixel>
::GetFileName()
{
	int offset = m_filename.findRev("/");
	return m_filename.mid(offset+1,m_filename.findRev(".")-offset-1);
}

template< class TPixel >
QString
QtImageLoader<TPixel>
::GetFilePath()
{
	return m_filename;
}


#endif
