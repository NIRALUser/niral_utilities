/**

 **/
#include <vector>
#include <string>

#include <itkImage.h>
#include <itkVector.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h> 
#include <itkImageSeriesReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>
#include <itkDICOMSeriesFileNames.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <DICOMAppHelper.h>
#include <DICOMParser.h>

#include "argio.h"

using namespace itk;
using namespace std;

static int debug = 1;

/*  main() that will instantiate the application  */
int main(int argc, const char* argv[])
{
  // argument processing
  if (argc <= 1 || ipExistsArgument(argv,"-help") || ipExistsArgument(argv,"-usage"))
    {
      cout << "usage: " << argv[0] << endl
      << "       infile1 [infile2 .....] outfile [-2Dfiles | -dicom | -dic2 | -dic3] [-float factor] [-int factor] [-v]" << endl
      << " converts a file in one ITK readable format into another, "
      << " the format is detected using the suffix of the files. "
      << " Data needs to be short. " << endl;
      cout << " if -2Dfiles is used then infile is the set of 2D files" << endl;
      cout << " if -dicom is used then infile the directory containing all the files [uses PatientPosition]" << endl;
      cout << " if -dic2 is used then infile the directory containing all the files [uses SlicePosition]" << endl;
      cout << " if -dic3 is used then infile is the set of 2D files, but DICOM is enforced" << endl;
      cout << endl;
      cout << " use -float if the input image is a double or float image (it will be converted to short)" << endl;
      cout << " -float has to supply a conversion multiplication factor" << endl;
      cout << " use -int if the input image is a int or long int image (it will be converted to short)" << endl;
      cout << " -int has to supply a conversion multiplication factor (can be a non-integer number)" << endl;
      exit(0) ;
    }
  string inputFileName(argv[1]); // input image

  int numFilename = 1;
  debug = ipExistsArgument(argv, "-v");

  while (argv[numFilename] && argv[numFilename][0] != '\0' && argv[numFilename][0] != '-')
    {
      numFilename++;
    }
  numFilename-- ;
  numFilename-- ;

  string outfileName;
  if (numFilename != 0)
    outfileName = argv[numFilename + 1];
//  if (numFilename == 0) outfileName.clear();

  int Files2DOn = 0, dicomOn = 0, dicom2On = 0, dicom3On = 0;
  Files2DOn = ipExistsArgument(argv, "-2Dfiles");
  dicomOn = ipExistsArgument(argv, "-dicom");
  dicom2On = ipExistsArgument(argv, "-dic2");
  dicom3On = ipExistsArgument(argv, "-dic3");
  bool doubleOn = ipExistsArgument(argv, "-float");
  double doubleConvFactor = ipGetDoubleArgument(argv, "-float",1.0);
  if (doubleConvFactor == 0.0 ) doubleConvFactor = 1.0;
  bool intOn = ipExistsArgument(argv, "-int");
  double intConvFactor = ipGetDoubleArgument(argv, "-int",1.0);
  if (intConvFactor == 0.0 ) intConvFactor = 1.0;

  vector<string> names;
  if (Files2DOn || dicom3On) {
    cout << numFilename << " inputfiles" << endl;
    for (int index = 1; index <= numFilename; index ++ ) {
      string inputName(argv[index]);
      names.push_back(inputName);
    }
  }
  bool writeFlag = true;
  if (outfileName.empty()) {
    writeFlag = false;
  }

  enum { ImageDimension = 4 };
  typedef   short PixelType;
  typedef   Image<PixelType,ImageDimension>  ImageType;

  typedef   int                                         IntPixelType;
  typedef   Image<IntPixelType,ImageDimension>          IntImageType;
  typedef   ImageFileReader< IntImageType >             VolumeIntReaderType;
  typedef   RescaleIntensityImageFilter< IntImageType,  IntImageType> RescaleIntFilterType;
  typedef   MinimumMaximumImageCalculator<IntImageType> MinMaxIntCalcType;
  typedef   CastImageFilter< IntImageType , ImageType>  CastIntFilterType;
  typedef   ImageRegionIterator< IntImageType >         IntIterator;

  typedef   double                                         DoublePixelType;
  typedef   Image<DoublePixelType,ImageDimension>          DoubleImageType;
  typedef   ImageFileReader< DoubleImageType >             VolumeDoubleReaderType;
  typedef   RescaleIntensityImageFilter< DoubleImageType,  DoubleImageType> RescaleDoubleFilterType;
  typedef   MinimumMaximumImageCalculator<DoubleImageType> MinMaxDoubleCalcType;
  typedef   CastImageFilter< DoubleImageType , ImageType>  CastDoubleFilterType;
  typedef   ImageRegionIterator< DoubleImageType >         DoubleIterator;

  typedef ImageFileReader< ImageType > VolumeReaderType;
  typedef ImageSeriesReader< ImageType > VolumeSeriesReaderType;
  typedef ImageFileWriter< ImageType > VolumeWriterType;
  try 
    {
      // load image
      if (debug) cout << "Loading file " << inputFileName << endl;
      if (Files2DOn || dicom3On) {
   VolumeSeriesReaderType::Pointer imageReader = VolumeSeriesReaderType::New();
   imageReader->SetFileNames(names) ;
   imageReader->Update() ;
   ImageType::Pointer output = imageReader->GetOutput();
   
   // Problems with size of slices -> take difference in position
   if (output->GetSpacing()[2] == 1.0) {
     DICOMPARSER_NAMESPACE::DICOMAppHelper m_AppHelper;
     DICOMPARSER_NAMESPACE::DICOMParser m_Parser;
     float pos1, pos2, *allpos;
     m_AppHelper.Clear();
     m_Parser.OpenFile((names[0]).c_str());
     m_Parser.ClearAllDICOMTagCallbacks();
     m_AppHelper.RegisterCallbacks(&m_Parser);
     m_Parser.ReadHeader();
     allpos = m_AppHelper.GetImagePositionPatient();
     pos1 = allpos[2];
     m_AppHelper.Clear();
     m_Parser.OpenFile((names[1]).c_str());
     m_Parser.ClearAllDICOMTagCallbacks();
     m_AppHelper.RegisterCallbacks(&m_Parser);
     m_Parser.ReadHeader();
     allpos = m_AppHelper.GetImagePositionPatient();
     pos2 = allpos[2];
     
     double slice_thickness = fabs(pos2 - pos1);
     if (slice_thickness == 0) {
       slice_thickness= 1.0;
     }
     if (debug) cout << "adjusting spacing to " << slice_thickness << "," << pos2 <<  endl;
     const itk::Vector<double, 4> origpixdims = output->GetSpacing();
     double newpixdims[4];
     newpixdims[0] = origpixdims[0];
     newpixdims[1] = origpixdims[1];
     newpixdims[2] = slice_thickness;
     newpixdims[3] = origpixdims[3];
     
     output->SetSpacing(newpixdims);
   }
     ImageType::Pointer image = imageReader->GetOutput();
     if (debug) {
		cout << "Input Image Info: ";
       ImageType::RegionType region = image->GetLargestPossibleRegion();
       std::cout << "dims: " << region.GetSize(0) << " , "  << region.GetSize(1) << " , "
            << region.GetSize(2) << " , " << region.GetSize(3) << "," << std::endl;
       const itk::Vector<double,4> inSpacing = image->GetSpacing();
       std::cout << "voxdims: " << inSpacing[0] << " , "  << inSpacing[1] << " , "
            << inSpacing[2] << " , " << inSpacing[3] << "," << std::endl;
     }

   
   // write output
   if (writeFlag) {
     if (debug) cout << "writing output data" << endl;
     VolumeWriterType::Pointer writer = VolumeWriterType::New();
	 writer->UseCompressionOn();
     writer->SetFileName(outfileName.c_str()); 
     
     writer->SetInput(imageReader->GetOutput());
     writer->Write();
   }
      } else if (dicomOn || dicom2On) {
     DICOMSeriesFileNames::Pointer sortFilename = DICOMSeriesFileNames::New();
     sortFilename->SetDirectory(inputFileName);
     if (dicomOn) {
       sortFilename->SetFileNameSortingOrder(DICOMSeriesFileNames::SortBySliceLocation);
     } else {
        sortFilename->SetFileNameSortingOrder(DICOMSeriesFileNames::SortByImagePositionPatient);
     }
     names = sortFilename->GetFileNames();
     
     VolumeSeriesReaderType::Pointer imageReader = VolumeSeriesReaderType::New();
     imageReader->SetFileNames(names) ;
     imageReader->Update() ;
     ImageType::Pointer output = imageReader->GetOutput();
     
     // Problems with size of slices -> take difference in position
     if (output->GetSpacing()[2] == 1.0) {
       DICOMPARSER_NAMESPACE::DICOMAppHelper m_AppHelper;
       DICOMPARSER_NAMESPACE::DICOMParser m_Parser;
       float pos1, pos2, *allpos;
       m_AppHelper.Clear();
       m_Parser.OpenFile((names[0]).c_str());
       m_Parser.ClearAllDICOMTagCallbacks();
       m_AppHelper.RegisterCallbacks(&m_Parser);
       m_Parser.ReadHeader();
       allpos = m_AppHelper.GetImagePositionPatient();
       pos1 = allpos[2];
       m_AppHelper.Clear();
       m_Parser.OpenFile((names[1]).c_str());
       m_Parser.ClearAllDICOMTagCallbacks();
       m_AppHelper.RegisterCallbacks(&m_Parser);
       m_Parser.ReadHeader();
       allpos = m_AppHelper.GetImagePositionPatient();
       pos2 = allpos[2];
       
       double slice_thickness = fabs(pos2 - pos1);
       if (debug) cout << "adjusting spacing to " << slice_thickness << endl;
       const itk::Vector<double,4> origpixdims = output->GetSpacing();
       double newpixdims[4];
       newpixdims[0] = origpixdims[0];
       newpixdims[1] = origpixdims[1];
       newpixdims[2] = slice_thickness;
       newpixdims[3] = origpixdims[3];
       
       output->SetSpacing(newpixdims);
     }
     ImageType::Pointer image = output;
     if (debug) {
	   std::cout <<"DICOM Output Image: ";
       ImageType::RegionType region = image->GetLargestPossibleRegion();
       std::cout << "dims: " << region.GetSize(0) << " , "  << region.GetSize(1) << " , "
            << region.GetSize(2) << " , " << region.GetSize(3) << std::endl;
       const itk::Vector<double,4> inSpacing = image->GetSpacing();
       std::cout << "voxdims: " << inSpacing[0] << " , "  << inSpacing[1] << " , "
            << inSpacing[2] << " , " << inSpacing[3] << std::endl;
     }

     // write output
     if (writeFlag) {
       if (debug) cout << "writing output data" << endl;
       VolumeWriterType::Pointer writer = VolumeWriterType::New();
       writer->SetFileName(outfileName.c_str()); 
	   writer->UseCompressionOn();
       
       writer->SetInput(output);
       writer->Write();
     }
      } else if (doubleOn){
   VolumeDoubleReaderType::Pointer imageReader = VolumeDoubleReaderType::New();
   imageReader->SetFileName(inputFileName.c_str()) ;
   imageReader->Update() ;
   DoubleImageType::Pointer inputImage = imageReader->GetOutput(); 
   MinMaxDoubleCalcType::Pointer minmaxCalc = MinMaxDoubleCalcType::New();
   minmaxCalc->SetImage(inputImage);
   minmaxCalc->Compute();
   cout << minmaxCalc->GetMinimum() << "," << minmaxCalc->GetMaximum() << "," << doubleConvFactor << endl;

   DoubleIterator iterImage (inputImage, inputImage->GetBufferedRegion());
   while ( !iterImage.IsAtEnd() )  {
     double value =  ( double) iterImage.Get() * doubleConvFactor;
     if (value < 0 ) cout <<"-";
     iterImage.Set(value);
     ++iterImage;
   }


   //RescaleDoubleFilterType::Pointer scaleFilter  = RescaleDoubleFilterType::New();
   //scaleFilter->SetInput(imageReader->GetOutput());
   
   //scaleFilter->SetOutputMinimum(minmaxCalc->GetMinimum() * doubleConvFactor);
   //scaleFilter->SetOutputMaximum(minmaxCalc->GetMaximum() * doubleConvFactor);
   
   CastDoubleFilterType::Pointer castFilter = CastDoubleFilterType::New();
   castFilter->SetInput(inputImage);
   castFilter->Update();
   ImageType::Pointer image = castFilter->GetOutput();
   if (debug) {
     ImageType::RegionType region = image->GetLargestPossibleRegion();
     std::cout << "dims: " << region.GetSize(0) << " , "  << region.GetSize(1) << " , "
          << region.GetSize(2) << " , " << region.GetSize(3) << std::endl;
          const itk::Vector<double,4> inSpacing = image->GetSpacing();
     std::cout << "voxdims: " << inSpacing[0] << " , "  << inSpacing[1] << " , "
          << inSpacing[2] << " , " << inSpacing[3] << std::endl;
   }
   
   // write output
   if (writeFlag) {
     if (debug) cout << "writing output data" << endl;
     VolumeWriterType::Pointer writer = VolumeWriterType::New();
     writer->SetFileName(outfileName.c_str()); 
     writer->UseCompressionOn();
     
     writer->SetInput(castFilter->GetOutput());
     writer->Write();
   }
      } else if (intOn){
			   VolumeIntReaderType::Pointer imageReader = VolumeIntReaderType::New();
			   imageReader->SetFileName(inputFileName.c_str()) ;
			   imageReader->Update() ;
			   MinMaxIntCalcType::Pointer minmaxCalc = MinMaxIntCalcType::New();
			   minmaxCalc->SetImage(imageReader->GetOutput());
			   minmaxCalc->Compute();
			   cout << minmaxCalc->GetMinimum() << "," << minmaxCalc->GetMaximum() << "," << intConvFactor << endl;
			   RescaleIntFilterType::Pointer scaleFilter  = RescaleIntFilterType::New();
			   scaleFilter->SetInput(imageReader->GetOutput());
			   
			   scaleFilter->SetOutputMinimum((int) ((double) minmaxCalc->GetMinimum() * intConvFactor));
			   scaleFilter->SetOutputMaximum((int) ((double) minmaxCalc->GetMaximum() * intConvFactor));
			   
			   CastIntFilterType::Pointer castFilter = CastIntFilterType::New();
			   castFilter->SetInput(scaleFilter->GetOutput());
			   castFilter->Update();
			   ImageType::Pointer image = castFilter->GetOutput();
			   if (debug) {
				image->Print(cout);
				 std::cout << "INT Type Output Image: ";
				 ImageType::RegionType region = image->GetLargestPossibleRegion();
				 std::cout << "dims: " << region.GetSize(0) << " , "  << region.GetSize(1) << " , "
				  << region.GetSize(2) << " , " << region.GetSize(3) << std::endl;
				  const itk::Vector<double,4> inSpacing = image->GetSpacing();
				 std::cout << "voxdims: " << inSpacing[0] << " , "  << inSpacing[1] << " , "
				  << inSpacing[2] << " , " << inSpacing[3] << std::endl;
			   }
			   
			   // write output
			   if (writeFlag) {
				 if (debug) cout << "writing output data" << endl;
				 VolumeWriterType::Pointer writer = VolumeWriterType::New();
				 writer->SetFileName(outfileName.c_str()); 
				 writer->UseCompressionOn();
				 
				 writer->SetInput(castFilter->GetOutput());
				 writer->Write();
			   }
      } else {
	   VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
	   imageReader->SetFileName(inputFileName.c_str()) ;
	   imageReader->Update() ;
	   ImageType::Pointer image = imageReader->GetOutput();
	   if (debug) {
		 image->Print(cout);
		 cout << image->GetImageDimension() << endl;
		 std::cout << "Input Image: ";
		 ImageType::RegionType region = image->GetLargestPossibleRegion();
		 std::cout << "dims: " << region.GetSize(0) << " , "  << region.GetSize(1) << " , "
			  << region.GetSize(2) << " , " << region.GetSize(3) << std::endl;
			  const itk::Vector<double,4> inSpacing = image->GetSpacing();
		 std::cout << "voxdims: " << inSpacing[0] << " , "  << inSpacing[1] << " , "
			  << inSpacing[2] << " , " << inSpacing[3] << std::endl;
	   }
	   
	   // write output
	   if (writeFlag) {
		 if (debug) cout << "writing output data " << outfileName << endl;
		 VolumeWriterType::Pointer writer = VolumeWriterType::New();
		 writer->SetFileName(outfileName.c_str()); 
		 writer->UseCompressionOn();
		 
		 writer->SetInput(imageReader->GetOutput());
		 writer->Write();
	   }
      }
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "ITK exception caught in main" << std::endl;
    std::cerr << e << std::endl;
    }
  catch( std::exception & e )
    {
    std::cerr << "STD exception caught in main" << std::endl;
    std::cerr << e.what() << std::endl;
    }
  catch( ... )
    {
    std::cerr << "unknown exception caught in main" << std::endl;
    }


  return 0;

}
