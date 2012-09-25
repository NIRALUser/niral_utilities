#ifndef ItkFastIO_TXX
#define ItkFastIO_TXX

#include "ItkFastIO.h"

template< class TPixel >
FTIO<TPixel>
::FTIO()
{
}

template< class TPixel >
FTIO<TPixel>
::~FTIO()
{
}

template< class TPixel >
void
FTIO<TPixel>
::SetFilename(const char * name)
{
	m_Filename=name;
}

template< class TPixel >
void
FTIO<TPixel>
::SetImage2Write(FTIO<TPixel>::ImagePointer image)
{
	m_Image=image;
}

template< class TPixel >
typename FTIO<TPixel>::ImagePointer
FTIO<TPixel>
::GetITKImage( void )
{
	return m_Image;
}

template< class TPixel >
int
FTIO<TPixel>
::Update( void )
{

	if(!m_Filename)
	{
		return -1;
	}

	const char * basename = strrchr(m_Filename,'/');
    if(!basename) basename=m_Filename;

    if(strrchr(basename,'.')==NULL)
    {
      return -1;
    }

	
	std::string m_string(basename);
	std::string extension = m_string.substr(m_string.rfind("."));
	if (extension == ".gz")
		extension = m_string.substr(m_string.rfind(".",m_string.length()-4),m_string.rfind(".")-m_string.rfind(".",m_string.length()-4));

	char *ext = (char*) extension.c_str();


    typename itk::ImageFileReader<ImageType>::Pointer reader 
                                  = itk::ImageFileReader<ImageType>::New();

    reader->SetFileName(m_Filename);
    std::cout << "Reading " << m_Filename << std::endl;
    
    if(strcmp(ext,".hdr")==0)
    {
      AnalyzeIOType::Pointer  readerType  = AnalyzeIOType::New();
      readerType->SetByteOrderToLittleEndian();
      reader->SetImageIO(readerType);
      reader->SetFileName(m_Filename);

      if(!readerType->CanReadFile(m_Filename))
      {
        std::cout << "Switchind to big endian..." << std::endl;
        readerType->SetByteOrderToBigEndian();
      }
      
    }
    if(strcmp(ext,".gipl")==0)
    {
      GiplIOType::Pointer  readerType  = GiplIOType::New();
      readerType->SetByteOrderToLittleEndian();
      reader->SetImageIO(readerType);
      reader->SetFileName(m_Filename);
      if(!readerType->CanReadFile(m_Filename))
      {
        std::cout << "Switching to big endian..." << std::endl;
        readerType->SetByteOrderToBigEndian();
      }
    
    }
   /* if(strcmp(ext,".mha")==0)
    {
    
      MetaIOType::Pointer  readerType  = MetaIOType::New();
      readerType->SetByteOrderToLittleEndian();
      reader->SetImageIO(readerType);
      reader->SetFileName(m_Filename);
      if(!readerType->CanReadFile(m_Filename))
      {
        std::cout << "Switching to big endian..." << std::endl;
        readerType->SetByteOrderToBigEndian();
      }

    }*/

    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      
      std::cerr << "exception in file reader " << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
       return -1;  //return EXIT_FAILURE
    }
  
    m_Image = reader->GetOutput();
    m_Image->SetSpacing(reader->GetOutput()->GetSpacing());

    return 0;

}


template< class TPixel >
int
FTIO<TPixel>
::Write( void )
{

if(!m_Image || !m_Filename)
{

	return 0;

}
else
{

	typename itk::ImageFileWriter<ImageType>::Pointer writer 
                                  = itk::ImageFileWriter<ImageType>::New();

  char * myFile = (char *)(m_Filename);
  char * basename = strrchr(myFile,'/');
  if (!basename) basename = myFile;//no directories, so start at beginning

  
  if(strrchr(basename, '.')==NULL)
  {
    std::cout << "No (valid) extension found, default is Analyze" << std::endl;
    myFile = new char[strlen(myFile)+4];
    sprintf(myFile, "%s%s", m_Filename, ".hdr");
    basename = strrchr(myFile, '/');
    if (!basename) basename = myFile;
          
  }

	std::string m_string(basename);
	std::string extension = m_string.substr(m_string.rfind("."));
	if (extension == ".gz")
		extension = m_string.substr(m_string.rfind(".",m_string.length()-4),m_string.rfind(".")-m_string.rfind(".",m_string.length()-4));

	char *ext = (char*) extension.c_str();

  
  if(0 == strncmp(ext,".gipl",5))
  {
    GiplIOType::Pointer  writerType  = GiplIOType::New();
    writerType->SetByteOrderToBigEndian();
    writer->SetImageIO(writerType);
    writer->SetFileName(myFile);
    writer->SetInput(m_Image);
  }
  /*else if(0 == strncmp(ext,".mha",4))
  {
    MetaIOType::Pointer writerType = MetaIOType::New();
    writerType->SetByteOrderToBigEndian();
    writer->SetImageIO(writerType);
    writer->SetFileName(myFile);
    writer->SetInput(m_Image);
  }*/
  else if(0 == strncmp(ext,".hdr",4))
  {
    AnalyzeIOType::Pointer writerType = AnalyzeIOType::New();
    writerType->SetByteOrderToBigEndian();
    writer->SetImageIO(writerType);
    writer->SetFileName(myFile);
    writer->SetInput(m_Image);
  }
  else
  {
    return -1;
  }

  
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "exception in file writer " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return -1;
  }

  return 0;
}
}


#endif
