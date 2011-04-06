/*=========================================================================

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>

#include <itkExtractImageFilter.h>

#include "argio.h"

using namespace itk;

int main(int argc, const char* argv[])
{


  if (argc <= 1 || ipExistsArgument(argv,"-h ") || ipExistsArgument(argv,"-help") || 
      ipExistsArgument(argv,"-usage"))
  {
    std::cout << "usage: " << argv[0] << std::endl
	      << "       infile [-o outfile] [-a axis] [-v] " << std::endl
	      << " -a axis : axis along which the file is split (default 3: time axis)"
	      << std::endl;
    exit(0);
  }
  // argument processing
  std::string fileName(argv[1]);
  
  bool debug =  ipExistsArgument(argv,"-v"); 
  bool axisid = ipExistsArgument(argv,"-a");
  std::string outfileName(ipGetStringArgument(argv, "-o", ""));

  //Set the axis along which the file is split, default is 3
  int splitaxis = 3;
  if(axisid)
    splitaxis = ipGetIntArgument(argv,"-a",0);

  if (outfileName.empty()) 
  {
    outfileName = "output.nrrd";
    std::cout << "no outputname specified using " << outfileName << std::endl;
  } 

  // Definitions
  const int inputDimension = 4;
  const int outputDimension = 3;

  typedef float ImagePixelType;

  typedef Image<ImagePixelType,inputDimension>  inputImageType;
  typedef Image<ImagePixelType,outputDimension>  outputImageType;

  typedef ImageFileReader< inputImageType > VolumeReaderType;
  typedef ImageFileWriter< outputImageType > VolumeWriterType;

  typedef ExtractImageFilter<inputImageType, outputImageType> ExtractFilterType;
  
  //Read the input file
  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  if (debug) std::cout << "Loading file " << std::endl;
  try
    {
      imageReader->SetFileName(fileName.c_str()) ;
      imageReader->Update() ;
    }
  catch (ExceptionObject e)
    {
      e.Print(std::cout) ;
      exit(0) ;
    }

  inputImageType::Pointer image4d = inputImageType::New();  
  image4d = imageReader->GetOutput();
  if(debug) std::cout << image4d << std::endl;

  //Get the image parameters
  inputImageType::SizeType size = image4d->GetLargestPossibleRegion().GetSize();
  inputImageType::SpacingType spacing = image4d->GetSpacing();

  unsigned int nb_volume = size[splitaxis];

  if(debug) std::cout << "Size: " << size << " | Spacing: " << spacing << std::endl;

  //Set the output file name
  //Get the first occurence of "/" to get the file name
  std::string::size_type slashlocation = outfileName.find_last_of("/");
  std::string filenameonly("");
  std::string path("./");
  if(slashlocation > outfileName.size())
    filenameonly = outfileName;
  //It means that only the file name was given, no path
  else
  {
    filenameonly = outfileName.substr(slashlocation+1,outfileName.size());
    path = outfileName.substr(0,slashlocation+1);
  }

  if(debug) std::cout << "Path: " << path << " File: " << filenameonly << std::endl;
  std::string::size_type dotlocation = filenameonly.find(".");
  std::string extension = filenameonly.substr(dotlocation+1,filenameonly.size());
  std::string filewoext = filenameonly.substr(0,dotlocation);
  if(debug) std::cout << "FilewoExt: " << filewoext << " ext: " << extension << std::endl;

  inputImageType::RegionType extractregion;
  //For the extraction filter to know which dimension is colapsed, we set the size to 0
  size[splitaxis] = 0;
  //Loop over each volume to extract
  for(unsigned int volid = 0 ; volid < nb_volume ; volid++)
  {

    ExtractFilterType::Pointer extract = ExtractFilterType::New();
    extract->SetInput(image4d);

    inputImageType::IndexType extractindex;
    extractindex.Fill(0);
    extractindex[splitaxis] = volid;

    extractregion.SetSize(size);
    extractregion.SetIndex(extractindex);

    extract->SetExtractionRegion(extractregion);

    extract->Update();

    //Set the outputfilename

    VolumeWriterType::Pointer writer = VolumeWriterType::New();

    std::ostringstream numext;
    numext << volid;
    std::string zeropadding("");
    if(volid < 10)
      zeropadding = "00";
    else if(volid < 100)
      zeropadding = "0";
    else
      zeropadding = "";

    std::string finaloutfile = path + filewoext + zeropadding + numext.str() + "." + extension;
    if(debug) std::cout << "Final out: " << finaloutfile << std::endl;
    writer->SetFileName(finaloutfile);
    writer->SetInput(extract->GetOutput());
    writer->Update();
  }
}
