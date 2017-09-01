/**

 **/
#include <vector>
#include <string>

#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkVector.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include "argio.h"

using namespace itk;
using namespace std;

/*  main() that will instantiate the application  */
int main(int argc, const char* argv[])
{
  // argument processing
  if (argc <= 1 || ipExistsArgument(argv,"-help") || ipExistsArgument(argv,"-usage") || !ipExistsArgument(argv,"-dim"))
    {
      cout << "usage: " << argv[0] << endl
	   << "       infile outfile -dim <dimx,dimy> [options] " << endl
	   << " converts a text file in an ITK readable image " << endl;
      cout << " currently only work with 2D images" << endl;
      cout << " -dim <dimx>,<dimy>  dimension of input image" << endl;
      cout << " -pixsize <sizex>,<sizey>  dimension of input image" << endl;
      cout << " -v                  verbose" << endl;
      cout << " -help               create this message" << endl;
      exit(0) ;
    }
  string inputFileName(argv[1]); // input image
  string outfileName(argv[2]);

  const int numDimParam = 2;
  int dimParam[numDimParam];
  int dimOn = ipExistsArgument(argv,"-dim");
  if (dimOn)  {
    char *tmp_str    = ipGetStringArgument(argv, "-dim", NULL);
    int numDim       = ipExtractIntTokens(dimParam, tmp_str, numDimParam);
    if (numDim != numDimParam) {
      cerr << argv[0] << ": -dim needs "<< numDimParam << " parameters.\n";
      exit(-1);
    }
    free(tmp_str);
  }
  int pixdimParam[numDimParam];
  int pixdimOn = ipExistsArgument(argv,"-pixsize");
  if (pixdimOn)  {
    char *tmp_str    = ipGetStringArgument(argv, "-pixsize", NULL);
    int numDim       = ipExtractIntTokens(pixdimParam, tmp_str, numDimParam);
    if (numDim != numDimParam) {
      cerr << argv[0] << ": -pixSize needs "<< numDimParam << " parameters.\n";
      exit(-1);
    }
    free(tmp_str);
  } else {
    pixdimParam[0] = 1;
    pixdimParam[1] = 1;
  }

  int debug = ipExistsArgument(argv,"-v");

  bool writeFlag = true;
  if (outfileName.empty()) {
    writeFlag = false;
  }

  enum { ImageDimension = 2 };
  typedef   float PixelType;
  typedef   Image<PixelType,ImageDimension>  ImageType;

  typedef ImageFileWriter< ImageType > VolumeWriterType;
  typedef ImageRegionIterator< ImageType > IteratorType;

  try
    {

      ImageType::Pointer image = ImageType::New();

      ImageType::SizeType ImDim;
      ImDim[0] = static_cast<int>(dimParam[0]);
      ImDim[1] = static_cast<int>(dimParam[1]);

      ImageType::SpacingType spacing;
      spacing[0] = pixdimParam[0];
      spacing[1] = pixdimParam[1];

      image->SetSpacing(spacing);
      image->SetRegions(ImDim);
      image->Allocate();
      IteratorType iterImage (image, image->GetBufferedRegion());

      if (debug)
	cout << "reading " << inputFileName << endl;
      std::ifstream input(inputFileName.c_str(), std::ios::in);

      while (!input.eof())
	{
	  PixelType value;
	  input>>value;
	  iterImage.Set(value);
	  ++iterImage;
	}

		
      if (writeFlag) {
	VolumeWriterType::Pointer writer = VolumeWriterType::New();
	if (debug) cout << "writing output data " << outfileName << endl;
	writer->UseCompressionOn();
	writer->SetFileName(outfileName.c_str()); 	
	writer->SetInput(image);
	writer->Write();
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
