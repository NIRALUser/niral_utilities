/*
 * compute image math and combinations
 *
 * author:  Martin Styner
 *
 * changes: Yundi Shi, Ipek Oguz, and many many more
 *
 */


/****** Needed to build test ******/
#ifdef _WIN32
#define Module_EXPORT __declspec(dllexport)
#else
#define Module_EXPORT
#endif

#if defined(main) && !defined(REGISTER_TEST)
// If main defined as a preprocessor symbol, redefine it to the expected entry point.
#undef main
#define main ModuleEntryPoint

extern "C" {
  Module_EXPORT int ModuleEntryPoint(const int, const char**);
}
#endif
/****** Needed to build test ******/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4503 )
#endif

#ifndef isnan
  #define isnan _isnan // std::isnan -> isnan by Adrien Kaiser 01/22/2013 for windows compilation // windows knows _isnan and not isnan
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <sstream>

using namespace std;

#include <string.h>
#include <sys/types.h>
#include <stdlib.h>    // for exit, system
#include <math.h>
#include <errno.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageToImageFilter.h>

#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkCastImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkSquareImageFilter.h>
#include <itkSqrtImageFilter.h>

#include <itkCurvatureFlowImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkMeanImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkMinimumImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkRescaleIntensityImageFilter.h>

#include <itkDiffusionTensor3D.h>
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageRegionIterator.h>
#include <itkMetaDataObject.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkOtsuMultipleThresholdsImageFilter.h>
#include <itkMaskedImageToHistogramFilter.h>
#include <itkImageToHistogramFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkBinaryThinningImageFilter.h>

#include <itkConnectedComponentImageFilter.h>

#include "argio.h"
#include "ImageMath.h"

#define IMAGEMATH_VERSION "1.3.3"
#define DEFAULT_SAMP 2
// number of samples by default

using namespace std;
using namespace itk;

typedef float PixelType;
typedef short ShortPixelType;
typedef unsigned char BinaryPixelType;
enum { ImageDimension = 3 };
typedef Image<PixelType,ImageDimension>       ImageType;
typedef Image<ShortPixelType,ImageDimension>  ShortImageType;
typedef Image<BinaryPixelType,ImageDimension> BinaryImageType;
typedef ImageType::RegionType                 ImageRegionType;
typedef ImageRegionIterator< ImageType >      IteratorType;
typedef ImageRegionIteratorWithIndex< ImageType > IteratorWithIndexType;
typedef ImageRegionIterator< ShortImageType> ShortIteratorType;
typedef ImageRegionConstIterator<ImageType>   ConstIteratorType;
typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;
typedef ImageType::Pointer                    ImagePointer;
typedef ImageToImageFilter<ImageType,BinaryImageType> Conv_float2binType;
typedef ImageToImageFilter<BinaryImageType,ImageType> Conv_bin2floatType;

typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileWriter< ImageType >          VolumeWriterType;
typedef ImageFileWriter< BinaryImageType >    BinaryVolumeWriterType;
typedef ImageFileWriter< ShortImageType >     ShortVolumeWriterType;

typedef CastImageFilter< ImageType, BinaryImageType > castBinaryFilterType;
typedef CastImageFilter< ImageType,  ShortImageType > castShortFilterType;

typedef itk::OtsuThresholdImageFilter< ImageType, ImageType > OtsuFilterType;
typedef itk::OtsuMultipleThresholdsImageFilter< ImageType, ImageType > OtsuMultipleThresholdsFilterType;
typedef ThresholdImageFilter< ImageType > maskThreshFilterType;
typedef BinaryThresholdImageFilter< ImageType , ImageType > threshFilterType;
typedef BinaryThresholdImageFilter< ShortImageType , ImageType > ShortthreshFilterType;
typedef MaskImageFilter< ImageType, ImageType, ImageType >  maskFilterType;
typedef AddImageFilter< ImageType, ImageType,  ImageType > addFilterType;
typedef SquareImageFilter< ImageType, ImageType > squareFilterType;
typedef SubtractImageFilter< ImageType, ImageType,  ImageType > subFilterType;
typedef MultiplyImageFilter< ImageType, ImageType,  ImageType > mulFilterType;
typedef DivideImageFilter< ImageType, ImageType,  ImageType > divFilterType;
typedef SqrtImageFilter< ImageType, ImageType > sqrtFilterType;

typedef BinaryBallStructuringElement<PixelType,ImageDimension> StructuringElementType;

typedef BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> dilateFilterType;
typedef BinaryErodeImageFilter<ImageType, ImageType, StructuringElementType> erodeFilterType;
typedef HistogramMatchingImageFilter< ImageType , ImageType > matchHistogramFilterType;

typedef MeanImageFilter<ImageType, ImageType> meanFilterType;
typedef CurvatureFlowImageFilter<ImageType, ImageType> curvFilterType;
typedef DiscreteGaussianImageFilter<ImageType, ImageType> gaussFilterType;
typedef GrayscaleErodeImageFilter<ImageType, ImageType, StructuringElementType> erodegrayFilterType;
typedef GrayscaleDilateImageFilter<ImageType, ImageType, StructuringElementType> dilategrayFilterType;
typedef GrayscaleFillholeImageFilter<ImageType, ImageType> fillholegrayFilterType;
typedef GradientAnisotropicDiffusionImageFilter< ImageType, ImageType > anisoDiffFilterType;

typedef ConnectedComponentImageFilter<ImageType,ShortImageType> ConnectiveFilterType;
typedef RelabelComponentImageFilter<ShortImageType,ImageType> RelabelFilterType;

typedef ExtractImageFilter<ImageType, ImageType> CropFilterType;

typedef MaximumImageFilter<ImageType, ImageType> MaximumImageFilterType;
typedef MinimumImageFilter<ImageType, ImageType> MinimumImageFilterType;

typedef itk::MinimumMaximumImageCalculator<ImageType> MaxFilterType;
typedef itk::Statistics::MaskedImageToHistogramFilter< ImageType , ImageType > MaskedHistogramFilterType;
typedef itk::Statistics::ImageToHistogramFilter< ImageType > HistogramFilterType;
typedef itk::ShiftScaleImageFilter< ImageType, ImageType > ShiftScaleFilterType;

typedef itk::BinaryThinningImageFilter<ImageType,ImageType> SkeletonFilterType;

static int debug = 0;
static const int BGVAL = 0;
static const int FGVAL = 1;

double Average( vector< double > vec )
{
  double Av = 0 ;
  for(unsigned int i = 0 ; i < vec.size() ; i++ )
  {
    Av += vec[ i ] ;
  }
  Av /= (double)vec.size() ;
  return Av ;
}

double DistanceToMean( std::vector< double > vec1 , std::vector< double > vec2 )
{
    double Av1 = Average( vec1 ) ;
    double Av2 = Average( vec2 ) ;
    double var = 0 ;
    for(unsigned int i = 0 ; i < vec1.size() ; i++ )
    {
      var += ( vec1[ i ] - Av1 ) * ( vec2[ i ] - Av2 ) ;
    }
    return var ;
}

double ImagesCorrelation( ImageType::Pointer image1 , ImageType::Pointer image2 , ShortImageType::Pointer mask )
{
  typedef ImageType::IndexType IndexType ;
  IndexType index ;
  ShortIteratorType it( mask , mask->GetLargestPossibleRegion() ) ;
  std::vector< double > val1 ;
  std::vector< double > val2 ;
  for( it.GoToBegin() ; !it.IsAtEnd() ; ++it )
  {
    if( it.Get() )
    {
      index = it.GetIndex() ;
      val1.push_back( image1->GetPixel( index ) ) ;
      val2.push_back( image2->GetPixel( index ) ) ;
    }
  }
  double s1 = DistanceToMean( val1 , val1 ) ;
  double s2 = DistanceToMean( val2, val2 ) ;
  double s12 = DistanceToMean( val1 , val2 ) ;
  return s12/sqrt(s1*s2) ;
}

//Added by Yundi Shi
//Calculate the average intensity of an image with a mask file
double ImagesAverage( ImageType::Pointer image ,  ShortImageType::Pointer mask )
{
  typedef ImageType::IndexType IndexType ;
  IndexType index ;
  ShortIteratorType it( mask , mask->GetLargestPossibleRegion() ) ;
  std::vector< double > val ;
  for( it.GoToBegin() ; !it.IsAtEnd() ; ++it )
  {
    if( it.Get() )
    {
      index = it.GetIndex() ;
      val.push_back( image->GetPixel( index ) ) ;
    }
  }
  double s = Average( val ) ;
  return s ;
}

//Added by Yundi Shi
//Calculate the standard deviation of intensities of an image with a mask file
double ImagesStd( ImageType::Pointer image ,  ShortImageType::Pointer mask )
{
  typedef ImageType::IndexType IndexType ;
  IndexType index ;
  ShortIteratorType it( mask , mask->GetLargestPossibleRegion() ) ;
  std::vector< double > val ;
  for( it.GoToBegin() ; !it.IsAtEnd() ; ++it )
  {
    if( it.Get() )
    {
      index = it.GetIndex() ;
      val.push_back( image->GetPixel( index ) ) ;
    }
  }
  double image_sum = DistanceToMean( val , val ) ;
  double image_std = sqrt(image_sum/(val.size()-1.0));

  return image_std ;
}


//What pixeltype is the image
void GetImageType( char* fileName ,
                   itk::ImageIOBase::IOPixelType &pixelType ,
                   itk::ImageIOBase::IOComponentType &componentType )
{
  itk::ImageFileReader< BinaryImageType >::Pointer imageReader =
    itk::ImageFileReader< BinaryImageType >::New();
  imageReader->SetFileName( fileName ) ;
  imageReader->UpdateOutputInformation() ;
  pixelType = imageReader->GetImageIO()->GetPixelType() ;
  componentType = imageReader->GetImageIO()->GetComponentType() ;
}


ImagePointer RunCleanComponent(ImagePointer image, int cleanLabel, float cleanPercent)
{
  // 0. Get max on the label image
  MaxFilterType::Pointer maxFilter = MaxFilterType::New();
  maxFilter->SetImage(image);
  maxFilter->ComputeMaximum();
  int MaxLabel = (int) maxFilter->GetMaximum();

  // 1. first extract the label
  threshFilterType::Pointer threshFilter = threshFilterType::New();
  threshFilter->SetInput(image);
  threshFilter->SetLowerThreshold(cleanLabel);
  threshFilter->SetUpperThreshold(cleanLabel);
  threshFilter->SetOutsideValue (BGVAL);
  threshFilter->SetInsideValue (FGVAL);
  try {
    threshFilter->Update();
  }
  catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    exit(EXIT_FAILURE);
  }

  // 2. get the connected components of the label
  ConnectiveFilterType::Pointer Connective = ConnectiveFilterType::New();
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  Connective->SetInput(threshFilter->GetOutput());
  try {
    Connective->Update();
  }
  catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    exit(EXIT_FAILURE);
  }

  //3. Sort the components according to their size
  relabelFilter->SetInput(Connective->GetOutput());
  try {
    relabelFilter->Update();
  }
  catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    exit(EXIT_FAILURE);
  }
  // 4. get largest component that is smaller than the cutoff
  RelabelFilterType::ObjectSizeInPixelsContainerType compSizes = relabelFilter->GetSizeOfObjectsInPixels();
  // note, sizes are in pixels -> int
  int numComp = compSizes.size();
  int totalSize = 0;
  for (int comp = 0; comp < numComp; comp++) { totalSize += compSizes[comp];}
  int cutoff = numComp;
  float minSize = (totalSize / 100 * cleanPercent);
  for (int comp = 0; comp < cutoff; comp++) {
    if (compSizes[comp] < minSize) {
      cutoff = comp;
    }
  }
  if (debug) cout << "numComps: " << numComp << "; total Size: " << totalSize << "; cutOff size " << minSize << "; cutOff component " << cutoff + 1 << endl;
  if (cutoff == numComp) {
      // nothing to do, return
      if (debug) cout << "No components smaller than cutoff, doing nothing" << endl;
      return image;
    }

  // relabel all components that are too small
  for (int comp = cutoff; comp < numComp; comp++) {
    PixelType compLabel = comp + 1; // component indexing is 1 off from label number

    if (debug) cout << "relabeling sorted component " << comp + 1 << " with size " << compSizes[comp] << endl;

    // 5. relabel the selected component that is too small
    // 5a) find neighboring labels (with number of neighboring voxels)
    IteratorType iterLabel (image, image->GetBufferedRegion());
    IteratorType iterComp (relabelFilter->GetOutput(), image->GetBufferedRegion());

    //  Relabeling: considering neighborhood, 6 connectedness
    NeighborhoodIteratorType::RadiusType Radius;
    Radius.Fill(1);
    NeighborhoodIteratorType NeighborhoodIterator(Radius,image,image->GetBufferedRegion());
    NeighborhoodIteratorType::OffsetType offset1 = {{-1,0,0}};
    NeighborhoodIteratorType::OffsetType offset2 = {{1,0,0}};
    NeighborhoodIteratorType::OffsetType offset3 = {{0,-1,0}};
    NeighborhoodIteratorType::OffsetType offset4 = {{0,1,0 }};
    NeighborhoodIteratorType::OffsetType offset5 = {{0,0,-1}};
    NeighborhoodIteratorType::OffsetType offset6 = {{0,0,1}};

    PixelType *LabelVotes = new PixelType[MaxLabel + 1];
    for (int Label = 0; Label <= MaxLabel; Label++)
      LabelVotes[Label] = 0;

    while ( !iterComp.IsAtEnd() )  {
      PixelType compVal =  iterComp.Get();

      if (compLabel == compVal) {
    // get votes from neighbors
    LabelVotes[(ShortPixelType)NeighborhoodIterator.GetPixel(offset1)]++;
    LabelVotes[(ShortPixelType)NeighborhoodIterator.GetPixel(offset2)]++;
    LabelVotes[(ShortPixelType)NeighborhoodIterator.GetPixel(offset3)]++;
    LabelVotes[(ShortPixelType)NeighborhoodIterator.GetPixel(offset4)]++;
    LabelVotes[(ShortPixelType)NeighborhoodIterator.GetPixel(offset5)]++;
    LabelVotes[(ShortPixelType)NeighborhoodIterator.GetPixel(offset6)]++;
      }

      ++iterComp;
      ++NeighborhoodIterator;
    }
    // reset current label to 0, as we don't care about it
    LabelVotes[cleanLabel] = 0;
    LabelVotes[0] = 0;

    PixelType MaxVoteLabel = 0;
    PixelType MaxVoteValue = 0;
    // now get the label with the max vote
    for (int Label = 0; Label <= MaxLabel; Label++)
      {
    if (LabelVotes[Label] > MaxVoteValue)
      {
        MaxVoteValue = LabelVotes[Label];
        MaxVoteLabel = Label;
      }
      }

    delete[] LabelVotes;
    iterComp.GoToBegin();
    if (debug) cout << "relabeling to  " << MaxVoteLabel << endl;

    // 5b. now reassign that component to the majority vote
    while ( !iterComp.IsAtEnd() )  {
      PixelType compVal =  iterComp.Get();

      if (compLabel == compVal) {
    iterLabel.Set(MaxVoteLabel);
      }

      ++iterComp;
      ++iterLabel;
    }

  }

  return image;

}

int main(const int argc, const char **argv)
{
  if (argc <=1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help")) {
    cout << "ImageMath version:"<<IMAGEMATH_VERSION<<" ("<<__DATE__<<")" << endl;
    cout << " computes base operations on an image" << endl;
    cout << "usage: ImageMath -version    Prints current version" << endl ;
    cout << "       ImageMath infile <options>" << endl;
    cout << endl;
    cout << "infile                 input dataset" << endl;;
    cout << "-outbase outbase       base-outputfilename, if omitted then the same as base of input" << endl;
    cout << "-outfile outfile       outfilename, will only be applied to main output (if multiple output are given)" << endl;
    cout << "-verbose               verbose mode " << endl;
    cout << "-extension ext         Extension of the output (determines output FORMAT)" << endl;
    cout << "-nocomp                don't automatically compress the output image" << endl;
    cout << "-combine infile2       combine the inputfile by interpreting them as labelfiles. " << endl
       << "   -relabel              labels in infile2 will be relabeled to succeed labels in infile (no label overlap)" << endl;
    cout << "       labels in infile2 overwrite only the background label in infile1" << endl;
    cout << "-extractLabel label    extract the mentioned label from the file" << endl;
    cout << "-otsu                  apply a standard 2 class otsu threshold" << endl;
    cout << "-otsuMultipleThresholds        apply an n class otsu threshold. Parameters can be set with -otsuPara nbThresholds,labelOffset,nbHistogramBins" << endl;
    cout << "-otsuPara nbThresholds,labelOffset,nbHistogramBins    optional parameters for otsuMultipleThresholds" <<endl;
    cout << "-threshold min,max     threshold: everything I < min , I > max is set to 0, otherwise 1" << endl;
    cout << "-threshMask min,max    mask threshold: everything I < min , I > max is set to 0, otherwise is left as is" << endl;
    cout << "-mask infile2          use infile2 as a binary mask (intensities > 0 ) and combine with inputfile (DTI compatible)" << endl ;
    cout << "-constOper opID,val    apply the following operation to the image: I op val, op = +/0, -/1, */2, //3" << endl;
    cout << "-add infile2           apply the following operation to the image: I1 + I2" << endl;
    cout << "-sub infile2           apply the following operation to the image: I1 - I2" << endl;
    cout << "-mul infile2           apply the following operation to the image: I1 * I2" << endl;
    cout << "-div infile2           apply the following operation to the image: I1 / I2" << endl;
    cout << "-pwr arg               compute each voxels to power arg" << endl;
    cout << "-normalizeEMS count -EMSfile prob_1,prob_2,...,prob_count(infile should be grayscale template) normalizes the EMS prob maps" << endl;
    cout << "-Normalize prob_1 prob_2 ...   normalizes input probability maps using input image as a mask (sum of probability maps equals to NormValue)" << endl;
    cout << "  -NormValue Value         normalization value (default:255)" << endl;
    cout << "-editPixdims px,py,pz   simply change the pixdims (without reslicing) of the image (DTI compatible)" << endl;
    cout << "-dilate radius,val      apply isotropic dilation with ball element of given radius and value" << endl;
    cout << "-erode radius,val       apply isotropic erosion with ball element of given radius and value" << endl;
    cout << "-matchHistogram infile2  match the image histogram to the one in infile2" << endl;
    cout << "-matchHistoPara bins,points,thresh  optional parameters for matchHistogram (bins= number of histogram bins" << endl;
    cout << "                        points = number of control points, thresh = boolean for mean intensity threshol [0/1])" << endl;
    cout << "-smooth -gauss -curveEvol -grayOpen -grayClose -grayDilate -grayErode -grayFillHole -meanFilter -anisoDiff [-size val] [-iter num] [-siz3 val1,val2,val3]" << endl;
    cout << "                       smoothing of image using any of the mentioned smoothing filters" << endl;
    cout << "                       size is stuctural element size, variance, or timestep depending on the filter choice, alternatively can specify 3 elements via siz3" << endl;
    cout << "                       iter is number of iterations" << endl;
    cout << "-type byte|short|float Type of processing image (computations are always done with float images), default is short" << endl;
    cout << "-conComp n           For a binary image, rank all the connected components according to their sizes and create a binary image with the 'n' biggest ones" << endl;
    cout << "-cleanComp Label,perc  Extract Label 'label' from image, run connected components and clean up any components of size smaller tthn 'perc' percent size via majority vote over neighboring label" << endl;
    cout << "                       if Lbl=0 outputs the labeled image with all the components labeled according to their size" << endl;
    cout << "-changeOrig px,py,pz [or filename]   Change the origin of the image (DTI compatible)" << endl;
    cout << "-createIm X,Y,Z,px,py,pz Create an empty image with the specified parameters: image dimensions X,Y,Z, image resolution px,py,pz" << endl;
    cout << " -crop px,py,pz,w,h,d   cropimage: origin px,py,pz (startindex is 0) dimensions width(w),height(h), depth(d) (DTI compatible)" << endl;
    cout << "-changeSp spx,spy,spz [or filename] Change the spacing of the image (DTI compatible)" << endl;
    cout << "-max infile2           Compute the maximum between the corresponding pixels of the two input images" << endl;
    cout << "-min infile2           Compute the minimum between the corresponding pixels of the two input images" << endl;
    cout << "-avg infile2 infile3...       Compute the average image" << endl;
    cout << "-majorityVoting infile2 infile3...       Compute an accurate parcellation map considering a majority voting process" << endl;
    cout << "-weightedMajorityVoting infile2 infile3... -weights w1,w2,...     Compute an accurate parcellation map considering a weighted majority voting process" << endl;
    cout << "-center                Center image (DTI compatible)" << endl;
    cout << "-flip [x,][y,][z]      Flip image" << endl;
    cout << "-NaNCor                Removes NaN and set the value of those voxels the average of the (non-NaN) neighbor values" << endl;
    cout << "-std infile2 infile3...             Compute the standard deviation from a set of images (one image has to be specified for the input, but it's not included in the process)" << endl;
    cout << "-rescale min,max       Applies a linear transformation to the intensity levels of the input Image" << endl;
    cout << "-rescalePerc lowPerc,lowVal,upPerc,upVal       Applies a linear transformation to the intensity levels of the input Image, setting the provided lower and upper percentiles to the provided values. Does not use background in percentile computation." << endl;
    cout << "   -rescaleMask mask   Uses the provided mask for the computation of the percentiles" << endl;
    cout << "-correl images2,mask   Computes the correlation between 2 images voxelwize" << endl;
    cout << "-RegionAvg mask        Computes the average intensities over a masked region" << endl;
    cout << "-RegionStd mask        Computes the standard deviation of image intensities over a masked region" << endl;
    cout << "-pixelLookup x,y,z     Outputs the pixel intensity for a given coordinate" << endl ;
    cout << "-danDistanceMap        Outputs the Danielsson Distance Map" << endl ;
    cout << "-mosaic                Creates a mosaic image from 2 images" << endl ;
    cout << "  -mosaicStep size     Size of each window in the mosaic (default:10 voxels)" << endl;
    cout << "-setLocationTolerance  Sets the coordinate and direction tolerance allowed for the ITK filters" << endl ;
    cout << "-skeleton              Copmute Skeleton via 3D thinning, assumes image is binary" << endl ;
    cout << endl << endl;
    return EXIT_SUCCESS ;
  }
  if( ipExistsArgument(argv, "-version") )
  {
    cout << "ImageMath version:"<<IMAGEMATH_VERSION<<" ("<<__DATE__<<")" << endl;
    return EXIT_SUCCESS ;
  }

  const unsigned int MaxNumFiles = 1000;
  unsigned int NbFiles = 0;
  char *files[MaxNumFiles];
  float wparameters[MaxNumFiles];
  vector<string> InputFiles;

  char *inputFileName = strdup(argv[1]);
  char *outputFileName = ipGetStringArgument(argv, "-outfile", NULL);
  char *outbase    = ipGetStringArgument(argv, "-outbase", NULL);
  char * base_string;
  if (!outbase) {
    if (!outputFileName) {
      base_string = strdup(ipGetBaseName(inputFileName));
    } else {
      base_string = strdup(ipGetBaseName(outputFileName));
    }
  } else {
    base_string = outbase;
  }
  string outFileName ("dummy");
  char *typeChat       = ipGetStringArgument(argv, "-type", NULL);
  bool correlOn = ipExistsArgument(argv, "-correl");
  char *correl_str = ipGetStringArgument( argv, "-correl" , NULL ) ;
  vector< string > correlFiles ;
  if( correlOn )
  {
    char *correlFiles_tmp[ 2 ] ;
    int num = ipExtractStringTokens(correlFiles_tmp, correl_str, 2);
    if ( num != 2 )
    {
      cerr << "correl needs 2 comma separated entries" << endl;
      return EXIT_FAILURE ;
    }
    // read in the names of the ems files
    else{
      for(int i = 0 ; i < 2 ; i++){
    if(debug) cout << "reading file name " << correlFiles_tmp[ i ] << endl ;
    correlFiles.push_back(correlFiles_tmp[i]) ; }
    }
  }

  bool writeFloat= false;
  bool writeByte = false;
  if (typeChat && !strcmp(typeChat,"byte")) writeByte = true;
  if (typeChat && !strcmp(typeChat,"float"))writeFloat = true;

  char * formatChar = ipGetStringArgument(argv, "-extension", ".nrrd");
  string format;
  if (! strchr(formatChar, '.')) {
    format = string(".") + string(formatChar);
  } else {
    format = string(formatChar);
  }

  debug      = ipExistsArgument(argv, "-verbose");

  bool nocompOn = false;
  nocompOn = ipExistsArgument(argv, "-nocomp");
  bool center = false;
  center = ipExistsArgument(argv, "-center");

  char *combineFile    = ipGetStringArgument(argv, "-combine", NULL);
  bool relabelOn = ipExistsArgument(argv, "-relabel");

  bool extractLabelOn   = ipExistsArgument(argv, "-extractLabel");
  int extractLabel   = ipGetIntArgument(argv, "-extractLabel", 0);

  bool maskOn   = ipExistsArgument(argv, "-mask");
  char *maskFile    = ipGetStringArgument(argv, "-mask", NULL);

  bool addOn   = ipExistsArgument(argv, "-add");
  char *addFile    = ipGetStringArgument(argv, "-add", NULL);
  bool subOn   = ipExistsArgument(argv, "-sub");
  char *subFile    = ipGetStringArgument(argv, "-sub", NULL);
  bool mulOn   = ipExistsArgument(argv, "-mul");
  char *mulFile    = ipGetStringArgument(argv, "-mul", NULL);
  bool divOn   = ipExistsArgument(argv, "-div");
  char *divFile    = ipGetStringArgument(argv, "-div", NULL);

  bool skeletonOn = ipExistsArgument(argv, "-skeleton");

  bool OtsuOn   = ipExistsArgument(argv, "-otsu");

  bool OtsuMultipleThresholdsOn = ipExistsArgument(argv, "-otsuMultipleThresholds");
  bool thresholdOn    = ipExistsArgument(argv, "-threshold");
  char * tmp_str      = ipGetStringArgument(argv, "-threshold", NULL);
  PixelType tmin = 0;
  PixelType tmax = 0;
  float textend[3];
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "threshold needs 2 comma separated entries: min,max" << endl;
      return EXIT_FAILURE ;
    } else {
      tmin = (PixelType) textend[0];
      tmax = (PixelType) textend[1];
    }
  }

  tmp_str      = ipGetStringArgument(argv, "-otsuPara", NULL);
  int otsuPara[3];
  int nbThresholds, labelOffset, nbHistogramBins;
  if (tmp_str) {
    int num = ipExtractIntTokens(otsuPara, tmp_str, 3);
    if (3 != num) {
      cerr << "otsuPara needx 3 comma separated entries: nbThresholds,labelOffset,nbHistogramBins" << endl;
      return EXIT_FAILURE ;
    } else {
      nbThresholds = otsuPara[0];
      labelOffset = otsuPara[1];
      nbHistogramBins = otsuPara[2];
    }
  } else {
      nbThresholds = 2;
      labelOffset = 0;
      nbHistogramBins = 128;
  }

  bool threshMaskOn    = ipExistsArgument(argv, "-threshMask");
  tmp_str      = ipGetStringArgument(argv, "-threshMask", NULL);
  PixelType tmaskmin = 0;
  PixelType tmaskmax = 0;
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "mask threshold needs 2 comma separated entries: min,max" << endl;
      return EXIT_FAILURE ;
    } else {
      tmaskmin = (PixelType) textend[0];
      tmaskmax = (PixelType) textend[1];
    }
  }

  int erodeRadius, dilateRadius;
  erodeRadius = dilateRadius = 1;
  PixelType erodeVal=1, dilateVal=1;
  bool dilateOn   = ipExistsArgument(argv, "-dilate");
  tmp_str      = ipGetStringArgument(argv, "-dilate", NULL);
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "dilate needs 2 comma separated entries: radius,value" << endl;
      return EXIT_FAILURE ;
    } else {
      dilateRadius = (int) textend[0];
      dilateVal = (PixelType) textend[1];
    }
  }
  bool erodeOn   = ipExistsArgument(argv, "-erode");
  tmp_str      = ipGetStringArgument(argv, "-erode", NULL);
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "erode needs 2 comma separated entries: radius,value" << endl;
      return EXIT_FAILURE ;
    } else {
      erodeRadius = (int) textend[0];
      erodeVal = (PixelType) textend[1];
    }
  }

  bool editPixdimsOn    = ipExistsArgument(argv, "-editPixdims");
  tmp_str      = ipGetStringArgument(argv, "-editPixdims", NULL);
  float pixdims[3];
  if (tmp_str) {
    int num = ipExtractFloatTokens(pixdims, tmp_str, 3);
    if (3 != num) {
      cerr << "editPixdims needs 3 comma separated entries: px,py,pz " << endl;
      return EXIT_FAILURE ;
    }
  }

  bool constOperOn    = ipExistsArgument(argv, "-constOper");
  tmp_str      = ipGetStringArgument(argv, "-constOper", NULL);
  int operID = 0;
  PixelType operVal = 0;
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "oper needs 2 comma separated entries: opID,val" << endl;
      return EXIT_FAILURE ;
    } else {
      operID = (int) textend[0];
      operVal = (PixelType) textend[1];
    }
  }

  bool connectiveCompOn  =  ipExistsArgument(argv, "-conComp");
  tmp_str =     ipGetStringArgument(argv, "-conComp", NULL);
  int Lbl = 1;
  if(tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 1);
    if(1 != num){
      cerr << "conCompt only needs 1 value" << endl;
      return EXIT_FAILURE ;
    } else {
      Lbl = (int) textend[0];
    }
  }

  //-cleanComp Label,perc  Extract Label 'label' from image, run connected components and clean up any components of size smaller tthn 'perc' percent size via majority vote over neighboring label
  bool cleanCompOn  =  ipExistsArgument(argv, "-cleanComp");
  tmp_str =     ipGetStringArgument(argv, "-cleanComp", NULL);
  int cleanLabel = 1;
  float cleanPercent = 5;
  if(tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if(2 != num){
      cerr << "cleanComp needs exactlt 2 values" << endl;
      return EXIT_FAILURE ;
    } else {
      cleanLabel = (int) textend[0];
      cleanPercent = textend[1];
    }
  }


  bool normalizeEMSOn    = ipExistsArgument(argv, "-normalizeEMS");
  int EMScount = ipGetIntArgument(argv,"-normalizeEMS",1);
  tmp_str      = ipGetStringArgument(argv, "-EMSfile", NULL);
  char ** probFiles = new char*[EMScount]; // char * probFiles[EMScount]; replaced by 'new' by Adrien Kaiser 01/22/2013 for windows compilation
  vector<string>  EMSFiles;
  if (tmp_str) {
    if(debug) cout<<"there are "<<EMScount<<" files"<<endl;
    int num = ipExtractStringTokens(probFiles, tmp_str, EMScount);
    if (EMScount != num)
    {
      cerr << "normalizeEMS needs "<<EMScount<<" comma separated entries" << endl;
      return EXIT_FAILURE ;
    }
    // read in the names of the ems files
    else{
      for(int i = 0 ; i < EMScount ; i++){
    if(debug) cout<<"reading file name "<<probFiles[i]<<endl;
    EMSFiles.push_back(probFiles[i]);}
    }
  }
delete []probFiles ; // Added because 'new' by Adrien Kaiser 01/22/2013 for windows compilation
  bool NormalizeOn = ipExistsArgument(argv, "-Normalize");
  vector<string>  NormalizeFiles;
  if (NormalizeOn)
    {
      NbFiles = ipGetStringMultipArgument(argv, "-Normalize", files, MaxNumFiles);
      for(unsigned int i = 0 ; i < NbFiles ; i++)
    {
      if(debug) cout<<"reading file name "<<files[i]<<endl;
      NormalizeFiles.push_back(files[i]);
    }
    }
  double NormValue = ipGetDoubleArgument(argv,"-NormValue",255);

  bool changeOrigOn    = ipExistsArgument(argv, "-changeOrig");
  tmp_str      = ipGetStringArgument(argv, "-changeOrig", NULL);
  float origCoor[3];
  if (tmp_str) {
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(tmp_str) ;
    try
    {
      imageReader->UpdateOutputInformation();
      ImageType::PointType origin = imageReader->GetOutput()->GetOrigin() ;
      for( unsigned i = 0 ; i < 3 ; i++ )
      {
        origCoor[ i ] = origin[ i ] ;
      }
    }
    catch (ExceptionObject & err)
    {
      int num = ipExtractFloatTokens(textend, tmp_str, 3);
      if (3 == num)
      {
        origCoor[0] = textend[0];
        origCoor[1] = textend[1];
        origCoor[2] = textend[2];
      }
      else
      {
        cerr << "changeOrig needs 3 comma separated entries: px,py,pz " << endl;
        return EXIT_FAILURE ;
      }
    }
  }
  bool changeSpOn    = ipExistsArgument(argv, "-changeSp");
  tmp_str      = ipGetStringArgument(argv, "-changeSp", NULL);
  float spacingval[3];
  if (tmp_str) {
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(tmp_str) ;
    try
    {
      imageReader->UpdateOutputInformation();
      ImageType::SpacingType spacing = imageReader->GetOutput()->GetSpacing() ;
      for( unsigned i = 0 ; i < 3 ; i++ )
      {
        spacingval[ i ] = spacing[ i ] ;
      }
    }
    catch (ExceptionObject & err)
    {
      int num = ipExtractFloatTokens(textend, tmp_str, 3);
      if (3 != num)
      {
        cerr << "changeSp needs 3 comma separated entries: spx,spy,spz " << endl;
        return EXIT_FAILURE ;
      }
      else
      {
        spacingval[0] = static_cast<float>(textend[0]);
        spacingval[1] = static_cast<float>(textend[1]);
        spacingval[2] = static_cast<float>(textend[2]);
      }
    }
  }
  bool matchHistogramOn = ipExistsArgument(argv, "-matchHistogram");
  char *matchHistogramFile    = ipGetStringArgument(argv, "-matchHistogram", NULL);
  tmp_str      = ipGetStringArgument(argv, "-matchHistoPara", NULL);
  int matchHistoPara[3];
  int matchHistoNumBins, matchHistoNumPoints, matchHistoThresh;
  if (tmp_str) {
    int num = ipExtractIntTokens(matchHistoPara, tmp_str, 3);
    if (3 != num) {
      cerr << "matchHistoPara needx 3 comma separated entries: bins,points,threshBool" << endl;
      return EXIT_FAILURE ;
    } else {
      matchHistoNumBins = matchHistoPara[0];
      matchHistoNumPoints = matchHistoPara[1];
      matchHistoThresh = matchHistoPara[2];
    }
  } else {
      matchHistoNumBins = 1024;
      matchHistoNumPoints = 50;
      matchHistoThresh = 1;
  }

  bool imageCreationOn = ipExistsArgument(argv, "-createIm");
  tmp_str = ipGetStringArgument(argv, "-createIm", NULL);
  float Dims[6];
  if(tmp_str) {
    int num = ipExtractFloatTokens(Dims, tmp_str, 6);
    if(6 != num){
      cerr << "createIm needs 6 parameters separated by commas" << endl;
    } else {

      std::cout << "Val: " << Dims[0] << " | " << Dims[1] << " | " << Dims[2] << " | " << Dims[3] << " | " << Dims[4] << " | " << Dims[5] << std::endl;
    }
  }

  const int numCropParam = 6;
  int cropParam[numCropParam];
  bool cropOn = ipExistsArgument(argv,"-crop");
  if (cropOn) {
    tmp_str    = ipGetStringArgument(argv, "-crop", NULL);
    int numDim       = ipExtractIntTokens(cropParam, tmp_str, numCropParam);
    if (numDim != numCropParam) {
      cerr << argv[0] << ": crop needs "<< numCropParam << " parameters.\n";
      return EXIT_FAILURE ;
    }
    free(tmp_str);
  }

  bool smoothOn =  ipExistsArgument(argv,"-smooth");
  bool gaussianOn =  ipExistsArgument(argv,"-gauss");
  bool curveEvolOn = ipExistsArgument(argv,"-curveEvol");
  bool grayOpenOn = ipExistsArgument(argv,"-grayOpen");
  bool grayCloseOn = ipExistsArgument(argv,"-grayClose");
  bool grayDilateOn = ipExistsArgument(argv,"-grayDilate");
  bool grayErodeOn = ipExistsArgument(argv,"-grayErode");
  bool grayFillHoleOn = ipExistsArgument(argv, "-grayFillHole");
  bool meanFilterOn = ipExistsArgument(argv,"-meanFilter");
  bool anisoDiffOn = ipExistsArgument(argv,"-anisoDiff");
  bool pixelLookupOn = ipExistsArgument(argv,"-pixelLookup");

  tmp_str      = ipGetStringArgument(argv, "-pixelLookup", NULL);
  float pixelLookupIndex[3];
  for( unsigned int i = 0 ; i < 3 ; i++ )
  {
    pixelLookupIndex[ i ] = 0 ;
  }
  if (tmp_str) {
      int num = ipExtractFloatTokens(textend, tmp_str, 3);
      if (3 != num)
      {
        cerr << "pixelLookup needs 3 comma separated entries: x,y,z " << endl;
        return EXIT_FAILURE ;
      }
      else
      {
        pixelLookupIndex[0] = static_cast<float>(textend[0]);
        pixelLookupIndex[1] = static_cast<float>(textend[1]);
        pixelLookupIndex[2] = static_cast<float>(textend[2]);
      }
    }


  if (smoothOn && !gaussianOn && !curveEvolOn
        && !grayOpenOn && !grayCloseOn && !grayDilateOn && !grayErodeOn && !grayFillHoleOn
        && !meanFilterOn && !anisoDiffOn) {
    curveEvolOn = true;
  };
  double smoothSize = ipGetDoubleArgument(argv,"-size",-1);
  tmp_str      = ipGetStringArgument(argv, "-siz3", NULL);
  float smoothSizes[3];
  for( unsigned int i = 0 ; i < 3 ; i++ )
  {
    smoothSizes[ i ] = -1.0 ;
  }
  if (tmp_str) {
      int num = ipExtractFloatTokens(textend, tmp_str, 3);
      if (3 != num)
      {
        cerr << "siz3 needs 3 comma separated entries: Sx,Sy,Sz " << endl;
        return EXIT_FAILURE ;
      }
      else
      {
        smoothSizes[0] = static_cast<float>(textend[0]);
        smoothSizes[1] = static_cast<float>(textend[1]);
        smoothSizes[2] = static_cast<float>(textend[2]);
      }
    }


  unsigned int numIter = ipGetIntArgument(argv,"-iter",1);

  bool MaxOn   = ipExistsArgument(argv, "-max");
  char *MaxFile    = ipGetStringArgument(argv, "-max", NULL);

  bool MinOn   = ipExistsArgument(argv, "-min");
  char *MinFile    = ipGetStringArgument(argv, "-min", NULL);

  bool pwrOn   = ipExistsArgument(argv, "-pwr");
  tmp_str = ipGetStringArgument(argv, "-pwr", NULL);
  float pwrval = 1;
  if (tmp_str) {
    int num = ipExtractFloatTokens(&pwrval, tmp_str, 1);
    if (1 != num) {
      cerr << "pwr option requires one entry: the value of the power" << endl;
      return EXIT_FAILURE ;
    }
  }
  bool RegionAvgOn = ipExistsArgument(argv, "-RegionAvg");
  bool RegionStdOn = ipExistsArgument(argv, "-RegionStd");
  char *RegionMaskFile  = ipGetStringArgument(argv, "-RegionAvg", NULL);
  if(RegionStdOn) {
    RegionMaskFile  = ipGetStringArgument(argv, "-RegionStd", NULL);
  }

  bool AvgOn = ipExistsArgument(argv, "-avg");
  if (AvgOn)
    {
      NbFiles = ipGetStringMultipArgument(argv, "-avg", files, MaxNumFiles);
      for(unsigned int i = 0 ; i < NbFiles ; i++)
    InputFiles.push_back(files[i]);
    }

  bool MajorityVotingOn = ipExistsArgument(argv, "-majorityVoting");
  if (MajorityVotingOn)
    {
      NbFiles = ipGetStringMultipArgument(argv, "-majorityVoting", files, MaxNumFiles);
      for(unsigned int i = 0 ; i < NbFiles ; i++)
    InputFiles.push_back(files[i]);
    }

  bool WeightedMajorityVotingOn = ipExistsArgument(argv, "-weightedMajorityVoting");
  if (WeightedMajorityVotingOn)
    {
      NbFiles = ipGetStringMultipArgument(argv, "-weightedMajorityVoting", files, MaxNumFiles);
      char *mv_tmp_str   = ipGetStringArgument(argv, "-weights", NULL);
      unsigned int num = ipExtractFloatTokens(wparameters, mv_tmp_str, MaxNumFiles);
      if (NbFiles != num) {
          cerr << "weighted majority voting needs" << NbFiles << "comma separated entries: w1,w2,...,w" << NbFiles << endl;
          return EXIT_FAILURE ;
      }
      for(unsigned int i = 0 ; i < NbFiles ; i++)
    InputFiles.push_back(files[i]);
    }

  bool DanDistanceMapOn = ipExistsArgument(argv, "-danDistanceMap");
  char *mosaic = ipGetStringArgument(argv, "-mosaic", NULL);
  int mosaicStepSize = ipGetIntArgument(argv,"-mosaicStep",10);

  double locationTolerance = ipGetDoubleArgument(argv,"-setLocationTolerance",0.001);

  bool StdOn = ipExistsArgument(argv, "-std");
  if (StdOn)
  {
      NbFiles = ipGetStringMultipArgument(argv, "-std", files, MaxNumFiles);
      for(unsigned int i = 0 ; i < NbFiles ; i++)
    InputFiles.push_back(files[i]);
  }
  bool flip   = ipExistsArgument(argv, "-flip");
  tmp_str      = ipGetStringArgument(argv, "-flip", NULL);
  itk::Matrix< double , 3 , 3 > flipMatrix ;
  flipMatrix.SetIdentity() ;
  if( tmp_str )
  {
    bool x = false ;
    bool y = false ;
    bool z = false ;
    char * flipAxes[3];
    int num = ipExtractStringTokens(flipAxes, tmp_str, 3);
    for( int i = 0 ; i < num ; i++ )
    {
      if( !strcmp( flipAxes[ i ] , "x" ) )
      {
        if( x == true )//'x' appears multiple times
        {
          std::cout<<"x appears multiple times"<<std::endl;
          return EXIT_FAILURE ;
        }
        flipMatrix[ 0 ][ 0 ] = -1 ;
        x = true ;
      }
      else if( !strcmp( flipAxes[ i ] , "y" ) )
      {
        if( y == true )//'y' appears multiple times
        {
          std::cout<<"y appears multiple times"<<std::endl;
          return EXIT_FAILURE ;
        }
        flipMatrix[ 1 ][ 1 ] = -1 ;
        y = true ;
      }
      else if( !strcmp( flipAxes[ i ] , "z" ) )
      {
        if( z == true )//'z' appears multiple times
        {
          std::cout<<"z appears multiple times"<<std::endl;
          return EXIT_FAILURE ;
        }
        flipMatrix[ 2 ][ 2 ] = -1 ;
        z = true ;
      }
      else
      {
        std::cout<<"Error: can only flip along x, y, z or any combination of them"<<std::endl;
        return EXIT_FAILURE ;
      }
    }
  }
  bool nan   = ipExistsArgument(argv, "-NaNCor");

  bool rescalingOn = ipExistsArgument( argv , "-rescale" ) ;
  tmp_str = ipGetStringArgument( argv , "-rescale" , NULL ) ;
  float rescaling[2] ;
  if (tmp_str)
  {
    int num = ipExtractFloatTokens(rescaling, tmp_str, 2);
    if( 2 != num)
    {
      cerr << "Rescale needs 2 comma separated entries: min,max" << endl ;
      return EXIT_FAILURE ;
    }
  }

  bool rescalingPercOn = ipExistsArgument( argv , "-rescalePerc" ) ;
  tmp_str = ipGetStringArgument( argv , "-rescalePerc" , NULL ) ;
  float rescalingPercPara[4] ;
  if (tmp_str)
  {
    int num = ipExtractFloatTokens(rescalingPercPara, tmp_str, 4);
    if( 4 != num)
    {
      cerr << "Rescale needs 4 comma separated entries: lowPerc,lowVal,upPerc,upVal" << endl ;
      return EXIT_FAILURE ;
    }
  }
  bool rescalingMaskOn = ipExistsArgument( argv , "-rescaleMask" ) ;
  char *rescalingMask  = ipGetStringArgument(argv, "-rescaleMask", NULL);

  // **********************************************
  // **********************************************
  // **********************************************
  // Cmd line parsing done
  // **********************************************
  // **********************************************
  // **********************************************
  ImagePointer inputImage ;
  typedef itk::ImageBase< 3 > ImageBaseType ;
  ImageBaseType::Pointer inputBaseImage ;
  // Check the type of image that is loaded
  itk::ImageIOBase::IOPixelType pixelType ;
  itk::ImageIOBase::IOComponentType componentType ;
  GetImageType( inputFileName , pixelType , componentType ) ;
  bool diffusionImage = false ;
  //If tensor image
  typedef itk::Image< itk::DiffusionTensor3D< PixelType > , 3 > DiffusionImageType ;
  if (debug) cout << "Loading file " << inputFileName << endl;
  if( pixelType == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR
   || pixelType == itk::ImageIOBase::DIFFUSIONTENSOR3D
   || pixelType == itk::ImageIOBase::VECTOR
    )
  {
    if( !changeSpOn && !changeOrigOn && !editPixdimsOn && !cropOn && !center  && !maskOn )
    {
      std::cerr << "The only operations supported on Diffusion Tensor Images are: -editPixdims, -changeSp, -changeOrig, -crop , -mask and -center"<< std::endl ;
      return EXIT_FAILURE;
    }
    typedef itk::ImageFileReader< DiffusionImageType > DiffusionReaderType ;
    DiffusionReaderType::Pointer diffusionReader = DiffusionReaderType::New() ;
    diffusionReader->SetFileName( inputFileName ) ;
    try
    {
      diffusionReader->Update();
    }
    catch (ExceptionObject & err)
    {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputBaseImage = diffusionReader->GetOutput() ;
    diffusionImage = true ;
  }
  else
  {
 // load image
    if (debug) cout << "Loading file " << inputFileName << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(inputFileName) ;
    try
    {
      imageReader->Update();
    }
    catch (ExceptionObject & err)
    {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = imageReader->GetOutput();
    inputBaseImage = inputImage ;
  }



  // do something to InputImage
  if (combineFile) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_comb");

    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << combineFile << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(combineFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage2 = imageReader->GetOutput();

    IteratorType iterImage1 (inputImage, inputImage->GetBufferedRegion());
    IteratorType iterImage2 (inputImage2, inputImage->GetBufferedRegion());

    if (relabelOn) {
      if (debug) cout << "Relabeling image2 " << endl;
      // relabel image2 to contain consecutive labels starting at the max label in image1
      // get max in image1
      PixelType max = iterImage1.Get();
      while ( !iterImage1.IsAtEnd() )  {
      PixelType value =  iterImage1.Get();
      if (max < value) max = value;
      ++iterImage1;
      }
      iterImage1.GoToBegin();
      // analyze labels in image2
      typedef set<PixelType> labelSetType;
      labelSetType labelSet; // sets are automatically sorted and contain only unique entries
      labelSetType::iterator pos;
      while ( !iterImage2.IsAtEnd() )  {
        PixelType value =  iterImage2.Get();
        if (value) {
          labelSet.insert(value);
        }
        ++iterImage2;
      }
      iterImage2.GoToBegin();

      // relabel labels in image2
      if (debug) {
    std::cout << " Relabeling map : " << std::endl;
        pos = labelSet.begin();
        int intpos = 1;
        while (pos != labelSet.end()){
      std::cout << *pos << " ---> " << max + intpos << std::endl;
          intpos++;
          pos++;
        }

      }
      while ( !iterImage2.IsAtEnd() )  {
      PixelType value =  iterImage2.Get();
      if (value) {
        pos = labelSet.begin();
        int intpos = 1;
        while (pos != labelSet.end() && *pos != value){
          intpos++;
          pos++;
        }
        iterImage2.Set(max + intpos);

      }
      ++iterImage2;
      }
      iterImage2.GoToBegin();
    }

    // combine them now
    if (debug) cout << "Combining images " << endl;
    while ( !iterImage1.IsAtEnd() )  {
      PixelType value1 =  iterImage1.Get();
      PixelType value2 =  iterImage2.Get();

      if (!value1 && value2) {
       iterImage1.Set(value2);
      }
      ++iterImage1;
      ++iterImage2;
    }
  }
  else if (RegionAvgOn)
    {
      VolumeReaderType::Pointer maskReader = VolumeReaderType::New();
      maskReader->SetFileName(RegionMaskFile) ;
      castShortFilterType::Pointer castMask = castShortFilterType::New() ;
      castMask->SetInput( maskReader->GetOutput() ) ;
      castMask->Update() ;
      double avg_value = ImagesAverage(inputImage,castMask->GetOutput());
      cout<<"Average Intensity of the masked region is: "<<avg_value<<endl;
      return EXIT_SUCCESS ;
    }
  else if (RegionStdOn)
    {
      VolumeReaderType::Pointer maskReader = VolumeReaderType::New();
      maskReader->SetFileName(RegionMaskFile) ;
      castShortFilterType::Pointer castMask = castShortFilterType::New() ;
      castMask->SetInput( maskReader->GetOutput() ) ;
      castMask->Update() ;
      double std_value = ImagesStd(inputImage,castMask->GetOutput());
      cout<<"Standard deviation of the iamge intensity in the masked region is: "<<std_value<<endl;
      return EXIT_SUCCESS ;
    }
  else if( correlOn )
  {
    VolumeReaderType::Pointer image2Reader = VolumeReaderType::New();
    image2Reader->SetFileName( correlFiles[ 0 ] ) ;
    VolumeReaderType::Pointer maskReader = VolumeReaderType::New();
    maskReader->SetFileName(correlFiles[ 1 ] ) ;
    try
    {
      image2Reader->Update() ;
      maskReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    castShortFilterType::Pointer castMask = castShortFilterType::New() ;
    castMask->SetInput( maskReader->GetOutput() ) ;
    castMask->Update() ;
    double corr_value = ImagesCorrelation( inputImage , image2Reader->GetOutput() , castMask->GetOutput() ) ;
    cout << "Images correlation: " << corr_value << endl ;
    return EXIT_SUCCESS ;
  } else if (extractLabelOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_label");
    if (debug) cout << "extracting object " << extractLabel << endl;

    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(inputImage);
    threshFilter->SetLowerThreshold(extractLabel);
    threshFilter->SetUpperThreshold(extractLabel);
    threshFilter->SetOutsideValue (BGVAL);
    threshFilter->SetInsideValue (FGVAL);
    try {
      threshFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = threshFilter->GetOutput();

  } else if (OtsuOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_Otsu");

    if (debug) cout << "Otsu thresholding ..." << endl;

    OtsuFilterType::Pointer otsuFilter = OtsuFilterType::New();
    otsuFilter->SetInput(inputImage);

    otsuFilter->SetOutsideValue (FGVAL);
    otsuFilter->SetInsideValue (BGVAL);
    try {
      otsuFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    int threshold = otsuFilter->GetThreshold();
    if (debug) cout << "Otsu threshold " << threshold << endl;
    inputImage = otsuFilter->GetOutput();

  }
  //Added by Jean-Yves Yang
  //Apply an n classes Otsu thresholding
  else if(OtsuMultipleThresholdsOn) {
      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_OtsuMultipleThresholds");

      if (debug) cout << "OtsuMultipleThresholds thresholding ..." << endl;

      OtsuMultipleThresholdsFilterType::Pointer otsuMultipleThresholdsFilter = OtsuMultipleThresholdsFilterType::New();
      otsuMultipleThresholdsFilter->SetInput(inputImage);
      otsuMultipleThresholdsFilter->SetNumberOfHistogramBins(nbHistogramBins);
      otsuMultipleThresholdsFilter->SetNumberOfThresholds(nbThresholds);
      otsuMultipleThresholdsFilter->SetLabelOffset( labelOffset );
      try {
        otsuMultipleThresholdsFilter->Update();
      }
      catch (ExceptionObject & err) {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
      }
     typedef itk::OtsuMultipleThresholdsImageFilter< ImageType, ImageType > FilterType;
     FilterType::Pointer filter = FilterType::New();
     FilterType::ThresholdVectorType thresholds = filter->GetThresholds();
     if (debug)  {
         cout << "Otsu thresholds " << endl;
         for( size_t i = 0; i < thresholds.size(); i++ )
         {
             std:: cout << thresholds[i] << std::endl;
         }
     }
     /*typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleType;
     RescaleType::Pointer rescaler = RescaleType::New();
     rescaler->SetInput( otsuMultipleThresholdsFilter->GetOutput() );
     rescaler->SetOutputMinimum( 0 );
     rescaler->SetOutputMaximum( 255 );*/
     inputImage = otsuMultipleThresholdsFilter->GetOutput();
  }
  else if (thresholdOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_thresh");
    if (debug) cout << "thresholding image  " << tmin << "," << tmax << endl;

    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(inputImage);
    threshFilter->SetLowerThreshold(tmin);
    threshFilter->SetUpperThreshold(tmax);
    threshFilter->SetOutsideValue (BGVAL);
    threshFilter->SetInsideValue (FGVAL);
    try {
      threshFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = threshFilter->GetOutput();
  }  else if (threshMaskOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_threshMask");
    if (debug) cout << "threshold/mask image  " << tmaskmin << "," << tmaskmax << endl;

    maskThreshFilterType::Pointer threshFilter = maskThreshFilterType::New();
    threshFilter->SetInput(inputImage);
    threshFilter->SetOutsideValue (BGVAL);
    threshFilter->ThresholdOutside(tmaskmin,tmaskmax);
    try {
      threshFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = threshFilter->GetOutput();
  } else if (maskOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_mask");

    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << maskFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(maskFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "masking images  " << endl;

    if( diffusionImage )
    {
      DiffusionImageType::Pointer inputDiffusionImage =
          dynamic_cast< DiffusionImageType* >( inputBaseImage.GetPointer() ) ;
      ImageType::SizeType maskSize ;
      maskSize = inputImage2->GetLargestPossibleRegion().GetSize() ;
      ImageType::SizeType size ;
      size = inputDiffusionImage->GetLargestPossibleRegion().GetSize() ;
      //Check that diffusion input volume and mask volume have the same size
      for( int i = 0 ; i < 3 ; i++ )
      {
        if( size[ i ] != maskSize[ i ] )
        {
          std::cout<<"Mask and input diffusion volume must have the same size"<<std::endl ;
          return EXIT_FAILURE ;
        }
      }
      //Create output volume and fill it with null tensors
      DiffusionImageType::Pointer outputVolume = DiffusionImageType::New() ;
      outputVolume->CopyInformation( inputDiffusionImage ) ;
      outputVolume->SetRegions( size ) ;
      outputVolume->Allocate() ;
      itk::DiffusionTensor3D< PixelType > pixel ;
      pixel.Fill( (PixelType)0.0 ) ;
      outputVolume->FillBuffer( pixel ) ;
      //Create iterators
      typedef itk::ImageRegionIterator< DiffusionImageType > DiffusionIterator ;
      typedef itk::ImageRegionIterator< ImageType > MaskIterator ;
      DiffusionIterator it( inputDiffusionImage , inputDiffusionImage->GetLargestPossibleRegion() ) ;
      MaskIterator maskIt(  inputImage2 ,  inputImage2->GetLargestPossibleRegion() ) ;
      DiffusionIterator out( outputVolume , outputVolume->GetLargestPossibleRegion() ) ;
      //Copy tensors where mask is not null
      for( it.GoToBegin() , maskIt.GoToBegin() , out.GoToBegin() ;
           !it.IsAtEnd() ; ++it , ++maskIt , ++out )
      {
        if( maskIt.Get() )
        {
          out.Set( it.Get() ) ;
        }
      }
      outputVolume->SetMetaDataDictionary( inputBaseImage->GetMetaDataDictionary() ) ;
      inputBaseImage = outputVolume ;
    }
    else
    {
      maskFilterType::Pointer maskFilter = maskFilterType::New() ;
      maskFilter->SetInput1( inputImage ) ;
      maskFilter->SetInput2( inputImage2 ) ;
      maskFilter->SetCoordinateTolerance( locationTolerance ) ;
      maskFilter->SetDirectionTolerance( locationTolerance ) ;
      try
      {
        maskFilter->Update() ;
      }
      catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl ;
        cerr << err << endl ;
        return EXIT_FAILURE ;
      }
      inputImage = maskFilter->GetOutput();
    }
  } else if (addOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_add");

    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << addFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(addFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "adding images  " << endl;

    addFilterType::Pointer addFilter = addFilterType::New();
    addFilter->SetCoordinateTolerance( locationTolerance ) ;
    addFilter->SetDirectionTolerance( locationTolerance ) ;
    addFilter->SetInput1(inputImage);
    addFilter->SetInput2(inputImage2);
    try {
      addFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = addFilter->GetOutput();

  } else if (subOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_sub");

    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << subFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(subFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "subing images  " << endl;

    subFilterType::Pointer subFilter = subFilterType::New();
    subFilter->SetInput1(inputImage);
    subFilter->SetInput2(inputImage2);
    subFilter->SetCoordinateTolerance( locationTolerance ) ;
    subFilter->SetDirectionTolerance( locationTolerance ) ;
    try {
      subFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = subFilter->GetOutput();

  } else if (mulOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_mul");

    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << mulFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(mulFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "muling images  " << endl;

    mulFilterType::Pointer mulFilter = mulFilterType::New();
    mulFilter->SetInput1(inputImage);
    mulFilter->SetInput2(inputImage2);
    mulFilter->SetCoordinateTolerance( locationTolerance ) ;
    mulFilter->SetDirectionTolerance( locationTolerance ) ;
    try {
       mulFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = mulFilter->GetOutput();

  } else if (divOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_div");

    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << divFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(divFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "diving images  " << endl;

    divFilterType::Pointer divFilter = divFilterType::New();
    divFilter->SetInput1(inputImage);
    divFilter->SetInput2(inputImage2);
    divFilter->SetCoordinateTolerance( locationTolerance ) ;
    divFilter->SetDirectionTolerance( locationTolerance ) ;
    try {
      divFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = divFilter->GetOutput();

  } else if (matchHistogramOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_matchHisto");

    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << matchHistogramFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(matchHistogramFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "matching images  " << endl;

    matchHistogramFilterType::Pointer matchHistogramFilter = matchHistogramFilterType::New();
    matchHistogramFilter->SetSourceImage(inputImage);
    matchHistogramFilter->SetReferenceImage(inputImage2);
    matchHistogramFilter->SetThresholdAtMeanIntensity(matchHistoThresh);
    matchHistogramFilter->SetNumberOfHistogramLevels(matchHistoNumBins); // number of bins
    matchHistogramFilter->SetNumberOfMatchPoints(matchHistoNumPoints); // number of equally distributed match points
    try {
      matchHistogramFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = matchHistogramFilter->GetOutput();

  } else if (constOperOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_oper");
    if (debug) cout << "Operation  " << operID << "," << operVal << endl;

    IteratorType iterImage1 (inputImage, inputImage->GetBufferedRegion());
    while ( !iterImage1.IsAtEnd() )  {
      PixelType value1 =  iterImage1.Get();
      PixelType value2 ;
      if (operID == 0)
      {
        value2 = value1 + operVal;
      }
      else if (operID == 1)
      {
        value2 = value1 - operVal;
      }
      else if (operID == 2)
      {
        value2 = value1 * operVal;
      }
      else if (operID == 3)
      {
        value2 = value1 / operVal;
      }
      else
      {
        std::cerr << "opID for -constOper has to be between 0 and 3" << std::endl ;
        return EXIT_FAILURE ;
      }

      iterImage1.Set(value2);
      ++iterImage1;
    }
  } else if (connectiveCompOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_conComp");

    if (debug) cout << "Get all the  " <<Lbl  << " bigger elements " << endl;

    ConnectiveFilterType::Pointer Connective = ConnectiveFilterType::New();
    RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
    //Get the connectivity map of the image
    Connective->SetInput(inputImage);
    try {
      Connective->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }

    //Sort the labels according to their size, each labeled object has a different value
    relabelFilter->SetInput(Connective->GetOutput());
    try {
      relabelFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }

    if(Lbl == 0)
      inputImage = relabelFilter->GetOutput();
    else
      {
    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(relabelFilter->GetOutput());
    threshFilter->SetLowerThreshold(0.1);
    threshFilter->SetUpperThreshold(Lbl);
    threshFilter->SetOutsideValue(0);
    threshFilter->SetInsideValue(1);
    try {
      threshFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = threshFilter->GetOutput();
      }
  } else if (cleanCompOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_cleanComp");
    if (debug) cout << "cleaning component " << cleanLabel << " regions of smaller than " << cleanPercent << " Percent" << endl;

    inputImage = RunCleanComponent(inputImage, cleanLabel, cleanPercent);

  } else if (dilateOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_dilate");

    if (debug) cout << "dilate ball radius  " <<dilateRadius  << ", val " << dilateVal << endl;
    StructuringElementType structuringElement;
    structuringElement.SetRadius( dilateRadius );  // 3x3x3 structuring element
    structuringElement.CreateStructuringElement( );

    dilateFilterType::Pointer dilateFilter = dilateFilterType::New();
    dilateFilter->SetInput(inputImage);
    dilateFilter->SetDilateValue (dilateVal);
    dilateFilter->SetKernel( structuringElement );
    try {
      dilateFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }

    inputImage = dilateFilter->GetOutput();
  } else if (changeOrigOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_newOrig");
//    inputImage->SetOrigin(origCoor);
    inputBaseImage->SetOrigin(origCoor);
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if (changeSpOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_newSpacing");
    std::cout << "SPACING: " << spacingval[0] << " | " << spacingval[1] << " | " << spacingval[2]  << std::endl;
//    inputImage->SetSpacing(spacingval);
    inputBaseImage->SetSpacing(spacingval);
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if (center) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_centered");
    ImageType::SizeType size ;
    size = inputBaseImage->GetLargestPossibleRegion().GetSize() ;
    ImageType::PointType origin ;
    origin = inputBaseImage->GetOrigin() ;
    itk::Index< 3 > index ;
    for( int i = 0 ; i < 3 ; i++ )
    {
      index[ i ] = size[ i ] - 1 ;
    }
    itk::Point< double , 3 > corner ;
    inputBaseImage->TransformIndexToPhysicalPoint( index , corner ) ;
    itk::Point< double, 3 > newOrigin ;
    for( int i = 0 ; i < 3 ; i++ )
    {
      newOrigin[ i ] = ( origin[ i ] - corner[ i ] ) / 2.0 ;
    }
    inputBaseImage->SetOrigin( newOrigin ) ;
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if (erodeOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_erode");

    if (debug) cout << "dilate ball radius  " << erodeRadius << ", val " << erodeVal << endl;
    StructuringElementType structuringElement;
    structuringElement.SetRadius( erodeRadius );  // 3x3x3 structuring element
    structuringElement.CreateStructuringElement( );

    erodeFilterType::Pointer erodeFilter = erodeFilterType::New();
    erodeFilter->SetInput(inputImage);
    erodeFilter->SetErodeValue (erodeVal);
    erodeFilter->SetKernel( structuringElement );
    try {
      erodeFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }

    inputImage = erodeFilter->GetOutput();

  } else if(cropOn)
    {
    ImageRegionType extractRegion;
    extractRegion.SetIndex(0,cropParam[0]);
    extractRegion.SetIndex(1,cropParam[1]);
    extractRegion.SetIndex(2,cropParam[2]);
    extractRegion.SetSize(0,cropParam[3]);
    extractRegion.SetSize(1,cropParam[4]);
    extractRegion.SetSize(2,cropParam[5]);

    int dim[3];
    ImageRegionType imageRegion = inputBaseImage->GetLargestPossibleRegion();
    dim[0] = imageRegion.GetSize(0);
    dim[1] = imageRegion.GetSize(1);
    dim[2] = imageRegion.GetSize(2);
    if (debug) cout << "size of the original image " << dim[0] << "," << dim[1] << "," << dim[2]<< endl;
    if (debug) cout << "cropping (x,y,z,w,h,d) " << extractRegion.GetIndex(0) << ","
            << extractRegion.GetIndex(1) << "," << extractRegion.GetIndex(2) << ","
            << extractRegion.GetSize(0) << "," << extractRegion.GetSize(1) << ","
            << extractRegion.GetSize(2) <<  endl;
    if( diffusionImage )
    {
      typedef ExtractImageFilter< DiffusionImageType , DiffusionImageType > DiffusionCropFilterType ;
      DiffusionCropFilterType::Pointer cropFilter = DiffusionCropFilterType::New() ;
      cropFilter->SetInput( dynamic_cast< DiffusionImageType* >( inputBaseImage.GetPointer() ) ) ;
      cropFilter->SetExtractionRegion( extractRegion ) ;
      try
      {
        cropFilter->Update() ;
      }
      catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
      }
      cropFilter->GetOutput()->SetMetaDataDictionary( inputBaseImage->GetMetaDataDictionary() ) ;
      inputBaseImage = cropFilter->GetOutput() ;

    }
    else
    {
      CropFilterType::Pointer cropFilter = CropFilterType::New() ;
      cropFilter->SetInput( inputImage ) ;
      cropFilter->SetExtractionRegion( extractRegion ) ;
      try
      {
        cropFilter->Update() ;
      }
      catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
      }
      inputImage = cropFilter->GetOutput();
    }
  } else if (smoothOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_smooth");

    if (curveEvolOn) {
      if (smoothSize == -1) smoothSize = 0.125;
      cout << "smoothing  (curveEvolution): size " << smoothSize << ", iterations " << numIter << endl;
      curvFilterType::Pointer smoothFilter = curvFilterType::New();

      smoothFilter->SetInput(inputImage);
      smoothFilter->SetNumberOfIterations(numIter);
      smoothFilter->SetTimeStep(smoothSize);
      try {
    smoothFilter->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      inputImage = smoothFilter->GetOutput();
    } else if (gaussianOn) {
      gaussFilterType::Pointer smoothFilter = gaussFilterType::New();

      if (smoothSize == -1) smoothSize = 1.0;
      if (smoothSizes[0] != -1) {
    smoothFilter->SetVariance(smoothSizes);
    cout << "smoothing  (Gaussian): size " << smoothSizes[0] <<"," << smoothSizes[1] <<","
         << smoothSizes[2] <<", iterations " << numIter << endl;
      } else {
    smoothFilter->SetVariance(smoothSize);
    cout << "smoothing  (Gaussian): size " << smoothSize << ", iterations " << numIter << endl;
      }
      for (unsigned int i=0; i < numIter; i++) {
    smoothFilter->SetInput(inputImage);
    try {
      smoothFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = smoothFilter->GetOutput();
      }
    } else if (grayDilateOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayDilate): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      dilategrayFilterType::Pointer dilateFilter = dilategrayFilterType::New();

      for (unsigned int i=0; i < numIter; i++) {
        dilateFilter->SetInput(inputImage);
        structuringElement.SetRadius( (int) smoothSize );
        structuringElement.CreateStructuringElement( );
        dilateFilter->SetKernel( structuringElement );
        try {
          dilateFilter->Update();
        }
        catch (ExceptionObject & err) {
          cerr << "ExceptionObject caught!" << endl;
          cerr << err << endl;
          return EXIT_FAILURE;
        }
        inputImage = dilateFilter->GetOutput();
      }
    } else if (grayErodeOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayErode): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      erodegrayFilterType::Pointer erodeFilter = erodegrayFilterType::New();

      for (unsigned int i=0; i < numIter; i++) {
        erodeFilter->SetInput(inputImage);
        structuringElement.SetRadius( (int) smoothSize );
        structuringElement.CreateStructuringElement( );
        erodeFilter->SetKernel( structuringElement );
        try {
          erodeFilter->Update();
        }
        catch (ExceptionObject & err) {
          cerr << "ExceptionObject caught!" << endl;
          cerr << err << endl;
          return EXIT_FAILURE;
        }
        inputImage = erodeFilter->GetOutput();
      }
    } else if (grayCloseOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayClose): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      dilategrayFilterType::Pointer dilateFilter = dilategrayFilterType::New();
      erodegrayFilterType::Pointer erodeFilter = erodegrayFilterType::New();

      for (unsigned int i=0; i < numIter; i++) {
    dilateFilter->SetInput(inputImage);
    erodeFilter->SetInput(dilateFilter->GetOutput());
    structuringElement.SetRadius( (int) smoothSize );
    structuringElement.CreateStructuringElement( );
    dilateFilter->SetKernel( structuringElement );
    erodeFilter->SetKernel( structuringElement );
    try {
      erodeFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = erodeFilter->GetOutput();
      }
    } else if (grayOpenOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayOpen): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      dilategrayFilterType::Pointer dilateFilter = dilategrayFilterType::New();
      erodegrayFilterType::Pointer erodeFilter = erodegrayFilterType::New();

      for (unsigned int i=0; i < numIter; i++) {
    erodeFilter->SetInput(inputImage);
    dilateFilter->SetInput(erodeFilter->GetOutput());
    structuringElement.SetRadius( (int) smoothSize);
    dilateFilter->SetKernel( structuringElement );
    erodeFilter->SetKernel( structuringElement );
    try {
      dilateFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = dilateFilter->GetOutput();
      }
    } else if (grayFillHoleOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayFillHole): size " << smoothSize << ", iterations " << numIter << endl;

      fillholegrayFilterType::Pointer fillholeFilter = fillholegrayFilterType::New();
      for (unsigned int i=0; i < numIter; i++) {
        fillholeFilter->SetInput(inputImage);
        try {
          fillholeFilter->Update();
        }
        catch (ExceptionObject & err) {
          cerr << "ExceptionObject caught!" << endl;
          cerr << err << endl;
          return EXIT_FAILURE;
        }
        inputImage = fillholeFilter->GetOutput();
      }
    } else if (meanFilterOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (meanFilter): size " << smoothSize << ", iterations " << numIter << endl;
      meanFilterType::Pointer smoothFilter = meanFilterType::New();

      for (unsigned int i=0; i < numIter; i++) {
    smoothFilter->SetInput(inputImage);
    ImageType::SizeType size;
    size[0] = size[1] = size[2] = (int) smoothSize;
    smoothFilter->SetRadius( size );
    try {
      smoothFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = smoothFilter->GetOutput();
      }
    } else if (anisoDiffOn) {
      if (smoothSize == -1) smoothSize = 0.05;
      cout << "smoothing  (gradient anisotropic Diffusion ): size  " << smoothSize << ", iterations " << numIter << endl;
      anisoDiffFilterType::Pointer smoothFilter = anisoDiffFilterType::New();

      smoothFilter->SetInput(inputImage);
      smoothFilter->SetNumberOfIterations(numIter);
      smoothFilter->SetConductanceParameter(0.75);
      smoothFilter->SetTimeStep(smoothSize);
      try {
    smoothFilter->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      inputImage = smoothFilter->GetOutput();
    }
  } else if (normalizeEMSOn) {

    vector<ImagePointer> EMSImages;
    ImagePointer restImage;
    {
      for (int i=0; i<EMScount; i++){
    if (debug) cout << "Loading image " << EMSFiles[i]  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(EMSFiles[i].c_str()) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    EMSImages.push_back(imageReader->GetOutput());
      }
    }

    if (debug) cout << "normalizing images  " << endl;
    if(EMScount > 1){
      //add files together
      addFilterType::Pointer addFilter = addFilterType::New();
      addFilter->SetInput1(EMSImages[0]);
      addFilter->SetInput2(EMSImages[1]);
      try {
    addFilter->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      if(debug) cout << "adding first two images done " << endl;
      inputImage = addFilter->GetOutput();
      restImage = addFilter->GetOutput();

      for (int i = 2; i < EMScount; i++){
    addFilterType::Pointer add2Filter = addFilterType::New();
    add2Filter->SetInput1(inputImage);
    add2Filter->SetInput2(EMSImages[i]);
    try {
      add2Filter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = add2Filter->GetOutput();
    if(debug) cout << "adding on the "<<i<<"th image "<< EMSFiles[i]<<" done " << endl;
      }
    }
    else{
      inputImage = EMSImages[0];
    }
    vector<IteratorType> EMSIters;

    IteratorType restIter (restImage, restImage->GetBufferedRegion());
    IteratorType normIter (inputImage, inputImage->GetBufferedRegion());

    for (int i = 0; i < EMScount; i++){

      if(debug) cout << "EMSImage_" <<i << endl;
      IteratorType EMSIter (EMSImages[i], inputImage->GetBufferedRegion());
      EMSIters.push_back(EMSIter);
    }
    while ( !normIter.IsAtEnd() )  {
      PixelType normVal =  normIter.Get();
      if (normVal < 255) {
    restIter.Set(255 - normVal);
    for (int i = 0; i < EMScount; i++){
      ++EMSIters[i];
    }
      }else{
    double factor = 255/normVal;
    for (int i = 0; i < EMScount; i++){
      restIter.Set(0); //used to be 1, why? ?????????????
      EMSIters[i].Set(EMSIters[i].Get() * factor);
      ++EMSIters[i];
    }
      }
      ++restIter;
      ++normIter;
    }

    for (int i = 0; i < EMScount; i++){
      std::stringstream ss;
      ss << (i+1);
      outFileName.erase();
      if(debug) cout<<"base string is "<<base_string<<endl;
      outFileName.append(base_string);
      outFileName.append(ss.str());
      outFileName.append("_normalized");
      outFileName.append(format);
      castShortFilterType::Pointer castFilter = castShortFilterType::New();
      castFilter->SetInput(EMSImages[i]);
      try {
    castFilter->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
      if(!nocompOn) writer->UseCompressionOn();
      writer->SetFileName(outFileName.c_str());
      writer->SetInput(castFilter->GetOutput());
      writer->Write();
    }
    {
      std::stringstream ss;
      ss << (EMScount+1);
      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append(ss.str());
      outFileName.append("_normalized");
      outFileName.append(format);
      castShortFilterType::Pointer castFilter = castShortFilterType::New();
      castFilter->SetInput(restImage);
      try {
    castFilter->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
      if(!nocompOn) writer->UseCompressionOn();
      writer->SetFileName(outFileName.c_str());
      writer->SetInput(castFilter->GetOutput());
      writer->Write();
    }
    return EXIT_SUCCESS ;
  } else if (NormalizeOn) {

    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(inputFileName) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = imageReader->GetOutput();

    vector<ImagePointer> NormalizeImages;
    for (unsigned int i = 0; i< NbFiles; i++){
      if (debug) cout << "Loading image " << NormalizeFiles[i]  << endl;
      imageReader = VolumeReaderType::New();
      imageReader->SetFileName(NormalizeFiles[i].c_str()) ;
      try {
    imageReader->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      NormalizeImages.push_back(imageReader->GetOutput());
    }

    if (debug) cout << "normalizing images  " << endl;
    ImagePointer addImage;
    if(NbFiles > 1){
      //add files together
      addFilterType::Pointer addFilter = addFilterType::New();
      addFilter->SetInput1(NormalizeImages[0]);
      addFilter->SetInput2(NormalizeImages[1]);
      try {
    addFilter->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      if(debug) cout << "adding first two images done " << endl;
      addImage = addFilter->GetOutput();

      for (unsigned int i = 2; i < NbFiles; i++){
    addFilterType::Pointer add2Filter = addFilterType::New();
    add2Filter->SetInput1(addImage);
    add2Filter->SetInput2(NormalizeImages[i]);
    try {
      add2Filter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    addImage = add2Filter->GetOutput();
    if(debug) cout << "adding on the "<<i<<"th image "<< NormalizeFiles[i]<<" done " << endl;
      }
    }
    else{
      addImage = NormalizeImages[0];}

    vector<IteratorType> NormalizeIters;
    IteratorType addIter (addImage, addImage->GetBufferedRegion());
    IteratorType inputIter (inputImage, inputImage->GetBufferedRegion());

    for (unsigned int i = 0; i < NbFiles; i++){
      if(debug) cout << "NormalizeImage_" <<i << endl;
      IteratorType NormalizeIter (NormalizeImages[i], inputImage->GetBufferedRegion());
      NormalizeIters.push_back(NormalizeIter);}

    if(debug) cout << "NormValue: "<<NormValue<<endl;

    while ( !inputIter.IsAtEnd() )  {
      PixelType inputVal = inputIter.Get();
      if (inputVal != 0)
    {
      PixelType addVal = addIter.Get();
      double factor = NormValue/addVal;
      for (unsigned int i = 0; i < NbFiles; i++)
        NormalizeIters[i].Set(NormalizeIters[i].Get() * factor);
    }
      for (unsigned int i = 0; i < NbFiles; i++)
    ++NormalizeIters[i];
      ++addIter;
      ++inputIter;
    }

    for (unsigned int i = 0; i < NbFiles; i++){
      std::stringstream ss;
      ss << (i+1);
      outFileName.erase();
      if(debug) cout<<"base string is "<<base_string<<endl;
      outFileName.append(base_string);
      outFileName.append(ss.str());
      outFileName.append("_normalized");
      outFileName.append(format);
      castShortFilterType::Pointer castFilter = castShortFilterType::New();
      castFilter->SetInput(NormalizeImages[i]);
      try {
    castFilter->Update();
      }
      catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;
      }
      ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
      if(!nocompOn) writer->UseCompressionOn();
      writer->SetFileName(outFileName.c_str());
      writer->SetInput(castFilter->GetOutput());
      writer->Write();
    }
    return EXIT_SUCCESS ;
  } else if (editPixdimsOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_pixDims");

//    inputImage->SetSpacing(pixdims);
    inputBaseImage->SetSpacing(pixdims);
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if(imageCreationOn) {

    ImageType::Pointer image = ImageType::New();

    ImageType::SizeType ImDim;
    ImDim[0] = static_cast<int>(Dims[0]);
    ImDim[1] = static_cast<int>(Dims[1]);
    ImDim[2] = static_cast<int>(Dims[2]);
    //    inputImage->SetSize(ImDim);
    ImageType::SpacingType spacing;
    spacing[0] = Dims[3];
    spacing[1] = Dims[4];
    spacing[2] = Dims[5];
    image->SetSpacing(spacing);
    image->SetRegions(ImDim);
    image->Allocate();

    inputImage = image;

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_createIm");
  } else if (MaxOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_max");

    if (debug) cout << "Loading file2 " << MaxFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(MaxFile) ;
    if (debug) cout << "computing maximum  " << endl;

    MaximumImageFilterType::Pointer MaximumFilter = MaximumImageFilterType::New();
    MaximumFilter->SetInput1(inputImage);
    MaximumFilter->SetInput2(imageReader->GetOutput());
    try {
      MaximumFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = MaximumFilter->GetOutput();

  } else if (MinOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_min");

    if (debug) cout << "Loading file2 " << MinFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(MinFile) ;
    if (debug) cout << "computing maximum  " << endl;

    MinimumImageFilterType::Pointer MinimumFilter = MinimumImageFilterType::New();
    MinimumFilter->SetInput1(inputImage);
    MinimumFilter->SetInput2(imageReader->GetOutput());
    try {
      MinimumFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }
    inputImage = MinimumFilter->GetOutput();

  } else if(pwrOn) {

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_pwr");

    if(debug) std::cout << "Applying power: " << pwrval << std::endl;

    //Creating the output image
    ImageType::Pointer imagepwr = ImageType::New();
    ImageType::SizeType ImDimpwr;
    ImDimpwr[0] = (inputImage->GetLargestPossibleRegion()).GetSize()[0];
    ImDimpwr[1] = (inputImage->GetLargestPossibleRegion()).GetSize()[1];
    ImDimpwr[2] = (inputImage->GetLargestPossibleRegion()).GetSize()[2];
    ImageType::SpacingType spacingpwr;
    spacingpwr[0] = inputImage->GetSpacing()[0];
    spacingpwr[1] = inputImage->GetSpacing()[0];
    spacingpwr[2] = inputImage->GetSpacing()[0];
    imagepwr->SetSpacing(spacingpwr);
    imagepwr->SetRegions(ImDimpwr);
    imagepwr->Allocate();

    IteratorType iterImage1 (inputImage, inputImage->GetLargestPossibleRegion());
    IteratorType iterImage2 (imagepwr,imagepwr->GetLargestPossibleRegion());
    iterImage2.GoToBegin();
    iterImage1.GoToBegin();
    float outval;

    while ( !iterImage1.IsAtEnd() )  {
      PixelType value =  iterImage1.Get();
      outval = 0;
      if (value == 0) {
    outval = 0;
      } else {
    outval = pow(static_cast<float>(value),static_cast<float>(pwrval));
      }
      iterImage2.Set(outval);
      ++iterImage2;
      ++iterImage1;
    }
    inputImage = imagepwr;

  } else if(AvgOn) {

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_avg");

    if(debug) cout << "Computing average image " << endl;

    VolumeReaderType::Pointer ImageReader = VolumeReaderType::New();
    addFilterType::Pointer addFilter = addFilterType::New();
    for (unsigned int FileNumber = 0; FileNumber < NbFiles; FileNumber++)
      {
    // Reading image
    ImageReader->SetFileName(InputFiles[FileNumber].c_str());

    // Adding image
    addFilter->SetInput1(inputImage);
    addFilter->SetInput2(ImageReader->GetOutput());
  addFilter->SetCoordinateTolerance( locationTolerance ) ;
  addFilter->SetDirectionTolerance( locationTolerance ) ;
    try
      {
        addFilter->Update();
      }
    catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
      }
    inputImage = addFilter->GetOutput();
      }
    IteratorType iterImage1 (inputImage, inputImage->GetBufferedRegion());
    while ( !iterImage1.IsAtEnd() )
      {
    PixelType NewValue = iterImage1.Get() / (NbFiles+1);
    iterImage1.Set(NewValue);
    ++iterImage1;
      }

  } else if(StdOn) {

     outFileName.erase();
     outFileName.append(base_string);
     outFileName.append("_std");

     if(debug) cout << "Computing standard deviation image" << endl;
     ImagePointer squareSumImage ;

     squareFilterType::Pointer squareFilter = squareFilterType::New();
     squareFilter->SetInput(inputImage);
     try
     {
        squareFilter->Update();
     }
     catch (ExceptionObject & err)
     {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
     }
     squareSumImage = squareFilter->GetOutput();

     VolumeReaderType::Pointer ImageReader = VolumeReaderType::New();
     addFilterType::Pointer addFilter = addFilterType::New();
     addFilterType::Pointer addFilterSquare = addFilterType::New();

     for (unsigned int FileNumber = 0; FileNumber < NbFiles; FileNumber++)
     {
    // Reading image
        ImageReader->SetFileName(InputFiles[FileNumber].c_str());

    // Adding image
        addFilter->SetInput1(inputImage);
        addFilter->SetInput2(ImageReader->GetOutput());
        addFilter->SetCoordinateTolerance( locationTolerance ) ;
        addFilter->SetDirectionTolerance( locationTolerance ) ;
        try
        {
           addFilter->Update();
        }
        catch (ExceptionObject & err)
        {
           cerr << "ExceptionObject caught!" << endl;
           cerr << err << endl;
           return EXIT_FAILURE;
        }
        inputImage = addFilter->GetOutput();

    // sum of square computation
        squareFilter = squareFilterType::New();
        squareFilter->SetInput(ImageReader->GetOutput());
        try
        {
           squareFilter->Update();
        }
        catch (ExceptionObject & err)
        {
           cerr << "ExceptionObject caught!" << endl;
           cerr << err << endl;
           return EXIT_FAILURE;
        }

        addFilterSquare->SetInput1(squareSumImage);
        addFilterSquare->SetInput2(squareFilter->GetOutput());
        addFilterSquare->SetCoordinateTolerance( locationTolerance ) ;
        addFilterSquare->SetDirectionTolerance( locationTolerance ) ;
        try
        {
           addFilterSquare->Update();
        }
        catch (ExceptionObject & err)
        {
           cerr << "ExceptionObject caught!" << endl;
           cerr << err << endl;
           return EXIT_FAILURE;
        }
        squareSumImage = addFilterSquare->GetOutput();

     }

     IteratorType iterImage1 (inputImage, inputImage->GetBufferedRegion());
     while ( !iterImage1.IsAtEnd() )
     {
        PixelType NewValue = iterImage1.Get() / (NbFiles+1);
        iterImage1.Set(NewValue);
        ++iterImage1;
     }

     IteratorType iterImage2 (squareSumImage, squareSumImage->GetBufferedRegion());
     while ( !iterImage2.IsAtEnd() )
     {
        PixelType NewValue2 = iterImage2.Get() / (NbFiles+1);
        iterImage2.Set(NewValue2);
        ++iterImage2;
     }

     squareFilterType::Pointer squareFilterMean = squareFilterType::New();
     squareFilterMean->SetInput(inputImage);
     try
     {
        squareFilterMean->Update();
     }
     catch (ExceptionObject & err)
     {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
     }

     ImagePointer MeanSquareImage ;
     MeanSquareImage = squareFilterMean->GetOutput();

     IteratorType iterImageMean (MeanSquareImage, MeanSquareImage->GetBufferedRegion());
     IteratorType iterImageSquare (squareSumImage, squareSumImage->GetBufferedRegion());

     while ( !iterImageMean.IsAtEnd() && !iterImageSquare.IsAtEnd())
     {
        PixelType NewValue = iterImageSquare.Get() - iterImageMean.Get() ;
        iterImageMean.Set(NewValue);
        ++iterImageMean;
        ++iterImageSquare;
     }

     IteratorType iterImageSqrt (MeanSquareImage, MeanSquareImage->GetBufferedRegion());
     while ( !iterImageSqrt.IsAtEnd() )
     {
        PixelType NewValue = sqrt(iterImageSqrt.Get());
        iterImageSqrt.Set(NewValue);
        ++iterImageSqrt;
     }

     inputImage = MeanSquareImage ;

  }  else if(MajorityVotingOn) {

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_majVoting");

    if(debug) cout << "Majority voting " << endl;

    // Reading label images
    vector<ImagePointer> vLabelImages;
    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      {
    VolumeReaderType::Pointer LabelImageReader = VolumeReaderType::New();
    if (debug) cout << "Loading file " << InputFiles[LabelFileNumber] << endl;
    LabelImageReader->SetFileName(InputFiles[LabelFileNumber].c_str());
    try
      {
        LabelImageReader->Update();
      }
    catch (ExceptionObject & err)
      {
        cerr<<"ExceptionObject caught!"<<endl;
        cerr<<err<<endl;
        return EXIT_FAILURE;
      }
    vLabelImages.push_back(LabelImageReader->GetOutput());
      }
    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(vLabelImages[0]);
    duplicator->Update();

    // Creating output image
    ImagePointer MajVotingImage = const_cast<ImageType*>(duplicator->GetOutput());

    //  Iterators initialization
    vector<ConstIteratorType> vConstLabelIterator;
    ConstIteratorType ConstInputIterator(inputImage,inputImage->GetRequestedRegion());
    IteratorType OutputIterator(MajVotingImage, inputImage->GetRequestedRegion());

    ConstInputIterator.GoToBegin();
    OutputIterator.GoToBegin();
    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      {
    ConstIteratorType ConstLabelIterator(vLabelImages[LabelFileNumber],inputImage->GetRequestedRegion());
    vConstLabelIterator.push_back(ConstLabelIterator);
    vConstLabelIterator[LabelFileNumber].GoToBegin();
      }

    //  Compute the maximum value of intensity of the first parcellation image
    ShortPixelType MaxLabel;
    MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetImage(vLabelImages[0]);
    maxFilter->ComputeMaximum();
    MaxLabel = (ShortPixelType) maxFilter->GetMaximum();
    MaxLabel++;

    //  Filling output image
    while (!ConstInputIterator.IsAtEnd())
      {
    PixelType MaxVoxelValue = 0, MaxLabelValue = 0;
    PixelType *LabelArray;
    LabelArray = new PixelType[MaxLabel];

    for (int Label = 0; Label < MaxLabel; Label++)
      LabelArray[Label] = 0;

    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
          {
        LabelArray[(ShortPixelType)vConstLabelIterator[LabelFileNumber].Get()]++;
          }
    for (int Label = 0; Label < MaxLabel; Label++)
      {
        if (LabelArray[Label] > MaxVoxelValue)
          {
        MaxVoxelValue = LabelArray[Label];
        MaxLabelValue = Label;
          }
      }

        if(MaxLabelValue >= 0.5)
        OutputIterator.Set(MaxLabelValue);
        else
        OutputIterator.Set(0);

          OutputIterator.Set(MaxLabelValue);

    delete[] LabelArray;

    ++ConstInputIterator;
    ++OutputIterator;
    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      ++vConstLabelIterator[LabelFileNumber];
      }

    //  Relabeling: considering neighborhood
    NeighborhoodIteratorType::RadiusType Radius;
    Radius.Fill(1);
    NeighborhoodIteratorType NeighborhoodOutputIterator(Radius,MajVotingImage,MajVotingImage->GetRequestedRegion());
    NeighborhoodIteratorType::OffsetType offset1 = {{-1,0,0}};
    NeighborhoodIteratorType::OffsetType offset2 = {{1,0,0}};
    NeighborhoodIteratorType::OffsetType offset3 = {{0,-1,0}};
    NeighborhoodIteratorType::OffsetType offset4 = {{0,1,0 }};
    NeighborhoodIteratorType::OffsetType offset5 = {{0,0,-1}};
    NeighborhoodIteratorType::OffsetType offset6 = {{0,0,1}};

    for (OutputIterator.GoToBegin(), NeighborhoodOutputIterator.GoToBegin(); !OutputIterator.IsAtEnd(); ++OutputIterator, ++NeighborhoodOutputIterator)
    {
      PixelType MaxVoxelValue = 0, MaxLabelValue = 0;
      PixelType *LabelArray;
      LabelArray = new PixelType[MaxLabel];
      for (int Label = 0; Label < MaxLabel; Label++)
        LabelArray[Label] = 0;
      LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset1)]++;
      LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset2)]++;
      LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset3)]++;
      LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset4)]++;
      LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset5)]++;
      LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset6)]++;

      for (int Label = 0; Label < MaxLabel; Label++)
        {
          if (LabelArray[Label] > MaxVoxelValue)
        {
          MaxVoxelValue = LabelArray[Label];
          MaxLabelValue = Label;
        }
        }

      if ( (MaxVoxelValue >= 4) && (MaxLabelValue != 0) && (OutputIterator.Get() != MaxLabelValue) )
        OutputIterator.Set(MaxLabelValue);

      delete[] LabelArray;
    }

    inputImage = MajVotingImage;
  } else if(WeightedMajorityVotingOn) {

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_majVoting");

    if(debug) cout << "Majority voting " << endl;

    // Reading label images
    vector<ImagePointer> vLabelImages;
    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      {
    VolumeReaderType::Pointer LabelImageReader = VolumeReaderType::New();
    if (debug) cout << "Loading file " << InputFiles[LabelFileNumber] << endl;
    LabelImageReader->SetFileName(InputFiles[LabelFileNumber].c_str());
    try
      {
        LabelImageReader->Update();
      }
    catch (ExceptionObject & err)
      {
        cerr<<"ExceptionObject caught!"<<endl;
        cerr<<err<<endl;
        return EXIT_FAILURE;
      }
    vLabelImages.push_back(LabelImageReader->GetOutput());

      }

    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(vLabelImages[0]);
    duplicator->Update();

    // Creating output image
    ImagePointer MajVotingImage =  const_cast<ImageType*>(duplicator->GetOutput());

    //  Iterators initialization
    vector<ConstIteratorType> vConstLabelIterator;
    ConstIteratorType ConstInputIterator(inputImage,inputImage->GetRequestedRegion());
    IteratorType OutputIterator(MajVotingImage, inputImage->GetRequestedRegion());

    ConstInputIterator.GoToBegin();
    OutputIterator.GoToBegin();
    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      {
    ConstIteratorType ConstLabelIterator(vLabelImages[LabelFileNumber],inputImage->GetRequestedRegion());
    vConstLabelIterator.push_back(ConstLabelIterator);
    vConstLabelIterator[LabelFileNumber].GoToBegin();
      }

    //  Compute the maximum value of intensity of the first parcellation image
    ShortPixelType MaxLabel;
    MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetImage(vLabelImages[0]);
    maxFilter->ComputeMaximum();
    MaxLabel = (ShortPixelType) maxFilter->GetMaximum();
    MaxLabel++;

    //  Filling output image
    while (!ConstInputIterator.IsAtEnd())
      {
    PixelType MaxVoxelValue = 0, MaxLabelValue = 0;
    PixelType *LabelArray;
    LabelArray = new PixelType[MaxLabel];

    for (int Label = 0; Label < MaxLabel; Label++)
      LabelArray[Label] = 0;

    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
          {
      LabelArray[(ShortPixelType)vConstLabelIterator[LabelFileNumber].Get()] += wparameters[LabelFileNumber];
          }

    for (int Label = 0; Label < MaxLabel; Label++)
      {
        if (LabelArray[Label] > MaxVoxelValue)
          {
        MaxVoxelValue = LabelArray[Label];
        MaxLabelValue = Label;
          }
      }

        if(MaxLabelValue >= 0.5)
        OutputIterator.Set(MaxLabelValue);
        else
        OutputIterator.Set(0);
    delete[] LabelArray;

    ++ConstInputIterator;
    ++OutputIterator;
    for (unsigned int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      ++vConstLabelIterator[LabelFileNumber];
      }

    //  Relabeling: considering neighborhood
    NeighborhoodIteratorType::RadiusType Radius;
      Radius.Fill(1);
      NeighborhoodIteratorType NeighborhoodOutputIterator(Radius,MajVotingImage,MajVotingImage->GetRequestedRegion());
      NeighborhoodIteratorType::OffsetType offset1 = {{-1,0,0}};
      NeighborhoodIteratorType::OffsetType offset2 = {{1,0,0}};
      NeighborhoodIteratorType::OffsetType offset3 = {{0,-1,0}};
      NeighborhoodIteratorType::OffsetType offset4 = {{0,1,0 }};
      NeighborhoodIteratorType::OffsetType offset5 = {{0,0,-1}};
      NeighborhoodIteratorType::OffsetType offset6 = {{0,0,1}};

        for (OutputIterator.GoToBegin(), NeighborhoodOutputIterator.GoToBegin(); !OutputIterator.IsAtEnd(); ++OutputIterator, ++NeighborhoodOutputIterator){
        PixelType MaxVoxelValue = 0, MaxLabelValue = 0;
        PixelType *LabelArray;
        LabelArray = new PixelType[MaxLabel];
        for (int Label = 0; Label < MaxLabel; Label++)
            LabelArray[Label] = 0;
        LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset1)]++;
        LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset2)]++;
        LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset3)]++;
        LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset4)]++;
        LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset5)]++;
        LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset6)]++;

        for (int Label = 0; Label < MaxLabel; Label++){
            if (LabelArray[Label] > MaxVoxelValue){
            MaxVoxelValue = LabelArray[Label];
            MaxLabelValue = Label;
                }
        }

            if ( (MaxVoxelValue >= 4) && (MaxLabelValue != 0) && (OutputIterator.Get() != MaxLabelValue) )
      //      if ( (MaxLabelValue != 0) && (OutputIterator.Get() != MaxLabelValue) )
            OutputIterator.Set(MaxLabelValue);
           // else
        //    OutputIterator.Set(MaxLabelValue);

        delete[] LabelArray;
    }
        inputImage = MajVotingImage;
    }

 else if(DanDistanceMapOn) {
 //   typedef unsigned char InputPixelType;
  //  typedef unsigned short OutputPixelType;
   // typedef itk::Image<InputPixelType, 3> InputImageType;
//    typedef itk::Image<OutputPixelType, 3> OutputImageType;
    typedef itk::DanielssonDistanceMapImageFilter<ImageType, ImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();

    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
    RescalerType::Pointer scaler = RescalerType::New();
/*
    typedef otb::ImageFileReader<InputImageType> ReaderType;
    typedef otb::ImageFileWriter<OutputImageType> WriterType;

    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();

    reader->SetFileName(argv[1]);
    writer->SetFileName(argv[2]);
*/
    filter->SetInput(inputImage);
   // scaler->SetInput(filter->GetOutput());
   // scaler->SetInput(filter->GetVoronoiMap());
   // writer->SetInput(scaler->GetOutput());

    //scaler->SetOutputMaximum(65535L);
   // scaler->SetOutputMinimum(0L);

    //filter->InputIsBinaryOn();


    try {
        filter->Update();
    }
    catch (ExceptionObject & err) {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
    }
    inputImage = filter->GetVoronoiMap();
/*
    writer->Update();
    const char * voronoiMapFileName = argv[3];

    scaler->SetInput(filter->GetVoronoiMap());
    writer->SetFileName(voronoiMapFileName);
    writer->Update();
*/
 }
 else if( flip )
  {
    ImageType::SizeType size ;
    size = inputImage->GetLargestPossibleRegion().GetSize() ;
    ImageType::PointType origin ;
    origin = inputImage->GetOrigin() ;
    itk::Index< 3 > index ;
    for( int i = 0 ; i < 3 ; i++ )
    {
      index[ i ] = size[ i ] - 1 ;
    }
    itk::Point< double , 3 > corner ;
    inputImage->TransformIndexToPhysicalPoint( index , corner ) ;
    itk::Point< double , 3 > newOrigin ;
    for( int i = 0 ; i < 3 ; i++ )
    {
     if( (int)flipMatrix[ i ][ i ] == -1 )
     {
       newOrigin[ i ] = -corner[ i ] ;
     }
     else
     {
       newOrigin[ i ] = origin[ i ] ;
     }
    }
    typedef itk::AffineTransform< double , 3 > AffineTransformType ;
    AffineTransformType::Pointer flipTransform = AffineTransformType::New() ;
    flipTransform->SetMatrix( flipMatrix ) ;
    //Resample Image
    typedef itk::NearestNeighborInterpolateImageFunction< ImageType , double > NearestNeighborInterpolateType ;
    NearestNeighborInterpolateType::Pointer interpolator = NearestNeighborInterpolateType::New() ;
    itk::ResampleImageFilter< ImageType , ImageType >::Pointer resampler ;
    resampler = itk::ResampleImageFilter< ImageType , ImageType >::New() ;
    resampler->SetOutputParametersFromImage( inputImage ) ;
    resampler->SetOutputOrigin( newOrigin ) ;
    resampler->SetInput( inputImage ) ;
    resampler->SetInterpolator( interpolator ) ;
    resampler->SetTransform( flipTransform ) ;
    resampler->Update() ;
    inputImage = resampler->GetOutput() ;
    inputImage->SetOrigin( origin ) ;
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_flip");
  }else if( nan )
  {
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType it( radius , inputImage , inputImage->GetLargestPossibleRegion() ) ;
     long counter = 0 ;
     for( it.GoToBegin() ; !it.IsAtEnd() ; ++it )
     {
        if( isnan( it.GetCenterPixel() ) )
        {
           counter++ ;
           float avg = 0.0 ;
           int counter2 = 0 ;
           for( SizeValueType c = 0 ; c < it.Size() ; c++ )
           {
             if( !isnan( it.GetPixel( c ) ) )
             {
               avg += it.GetPixel( c ) ;
               counter2++ ;
             }
           }
           if( counter2 > 0)
           {
             it.SetCenterPixel( static_cast<PixelType>( avg/float(counter2) ) ) ;
           }
        }
     }
     std::cout<< "Number Of NaN found: " << counter << std::endl ;
     outFileName.erase();
     outFileName.append(base_string);
     outFileName.append("_NaN");
  }else if( rescalingOn )
  {
     typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescalingType ;
     RescalingType::Pointer rescalingFilter=RescalingType::New() ;
     rescalingFilter->SetInput( inputImage ) ;
     rescalingFilter->SetOutputMinimum( rescaling[ 0 ] ) ;
     rescalingFilter->SetOutputMaximum( rescaling[ 1 ] ) ;
     rescalingFilter->Update() ;
     inputImage = rescalingFilter->GetOutput() ;
     outFileName.erase();
     outFileName.append(base_string);
     outFileName.append("_rescale");
  }else if( rescalingPercOn )
   {
     double lowPerc = rescalingPercPara[0];
     double highPerc = rescalingPercPara[2];
    double lowPercInt = rescalingPercPara[1];
    double highPercInt = rescalingPercPara[3];
    if (debug) std::cout << "lower Percentile: " <<  lowPerc << ", upper Percentile: " << highPerc << std::endl;

    maskFilterType::Pointer maskFilter = maskFilterType::New() ;

    if (rescalingMaskOn) {
      VolumeReaderType::Pointer rescalemaskReader = VolumeReaderType::New();
      rescalemaskReader->SetFileName(rescalingMask) ;
      rescalemaskReader->Update();
      if (debug) std::cout << "Reading " << rescalingMask  << std::endl;

      maskFilter->SetInput1( inputImage ) ;
      maskFilter->SetInput2( rescalemaskReader->GetOutput() ) ;
      maskFilter->SetCoordinateTolerance( 0.01 ) ;
      maskFilter->SetDirectionTolerance( 0.01 ) ;
      try
      {
        maskFilter->Update() ;
      }
      catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl ;
        cerr << err << endl ;
        return EXIT_FAILURE ;
      }
    }
    const int numberOfBins = 1000;
    HistogramFilterType::HistogramType::SizeType size(1);
    size.Fill(numberOfBins);

    HistogramFilterType::Pointer histoFilter = HistogramFilterType::New();
    if (rescalingMaskOn) {
      histoFilter->SetInput(maskFilter->GetOutput());
    } else {
      histoFilter->SetInput(inputImage);
    }

    histoFilter->SetAutoMinimumMaximum(true);
    histoFilter->SetHistogramSize( size );
    try
      {
    histoFilter->Update();
      }
    catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl ;
        cerr << err << endl ;
        return EXIT_FAILURE ;
      }

    HistogramFilterType::HistogramType::Pointer histogram = histoFilter->GetOutput();

    HistogramFilterType::HistogramType::ConstIterator histogramIterator = histogram->Begin();
    HistogramFilterType::HistogramType::MeasurementVectorType mv;
    HistogramFilterType::HistogramType::AbsoluteFrequencyType f, total, cur_total;

    HistogramFilterType::HistogramType::MeasurementVectorType lowPercMV, highPercMV, MaxMV;

    ++histogramIterator ;
    // Skip first bin, not interested in background...

    total = 0;
    while( histogramIterator  != histogram->End() )
      {
    mv =  histogramIterator.GetMeasurementVector();
    f = histogramIterator.GetFrequency();

    //if (debug) std::cout << mv  << "," << f <<  std::endl;
    total = total + f;
    ++histogramIterator ;
      }
    if (debug) std::cout << "Total " << histogram->GetTotalFrequency () << ", using only " << 100.0 * total / histogram->GetTotalFrequency () << "% due to background & mask removal" << std::endl;

    // Get quantiles
    cur_total = 0;
    MaxMV = mv; // mv of last bin
    lowPercMV = MaxMV; // mv of last bin
    highPercMV = MaxMV; // mv of last bin

    histogramIterator = histogram->Begin();
    if (lowPerc == 0) {
      lowPercMV = histogramIterator.GetMeasurementVector();
    }
    ++histogramIterator ;

    while( histogramIterator  != histogram->End() )
      {
    mv =  histogramIterator.GetMeasurementVector();
    f = histogramIterator.GetFrequency();
    cur_total = cur_total + f;
    double percentile = 100.0 * cur_total / total;
    if ( (lowPercMV == MaxMV) && (percentile > lowPerc) ) {
      lowPercMV = mv;
    }
    if ( (highPercMV == MaxMV) && (percentile > highPerc) ) {
      highPercMV = mv;
    }

    ++histogramIterator ;
      }

    if (debug) std::cout << "Image lowPerc " << lowPercMV[0] << ", highPerc " << highPercMV[0] << std::endl;

    double scale = (highPercInt - lowPercInt) / (highPercMV[0] - lowPercMV[0]);
    double shift = highPercInt/scale - highPercMV[0];

    if (debug) std::cout << "shift: " << shift << ", scale " << scale << " - "
             <<  scale * (shift + lowPercMV[0]) << ", "
             <<  scale * (shift + highPercMV[0]) << std::endl;
    ShiftScaleFilterType::Pointer shiftScaleFilter = ShiftScaleFilterType::New();
    shiftScaleFilter->SetInput(inputImage);
    shiftScaleFilter->SetShift(shift);
    shiftScaleFilter->SetScale(scale);
    try
      {
    shiftScaleFilter->Update();
      }
    catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl ;
        cerr << err << endl ;
        return EXIT_FAILURE ;
      }

    inputImage = shiftScaleFilter->GetOutput();

    // Get Quantiles
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_rescalePerc");

  } else if (pixelLookupOn )
    {
      // IO 2012
      ImageType::IndexType idx ;
      idx[0] = pixelLookupIndex[0] ;
      idx[1] = pixelLookupIndex[1] ;
      idx[2] = pixelLookupIndex[2] ;
      std::cout << "Image value at pixel " << idx[0] << " " << idx[1] << " " << idx[2] << " is: " << std::endl ;
      std::cout <<  inputImage->GetPixel ( idx ) << std::endl ;
      return EXIT_SUCCESS ;

  } else if( mosaic )
  {
    if (debug) std::cout << "Loading: " << mosaic  << std::endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(mosaic) ;
    imageReader->Update();
    //Resample images to same size
    if (debug) std::cout << "Resampling: " << mosaic << std::endl;
    typedef itk::LinearInterpolateImageFunction< ImageType , double > LinearInterpolateType ;
    LinearInterpolateType::Pointer interpolator = LinearInterpolateType::New() ;
    itk::ResampleImageFilter< ImageType , ImageType >::Pointer resampler ;
    resampler = itk::ResampleImageFilter< ImageType , ImageType >::New() ;
    resampler->SetOutputParametersFromImage( inputImage ) ;
    resampler->SetInput( imageReader->GetOutput() ) ;
    resampler->SetInterpolator( interpolator ) ;
    resampler->Update() ;
    //iterates through both inputs and output image
    if (debug) std::cout << "Creating mosaic. Step size: " << mosaicStepSize << std::endl ;
    IteratorWithIndexType it1( inputImage , inputImage->GetLargestPossibleRegion() ) ;
    IteratorWithIndexType it2( resampler->GetOutput() , resampler->GetOutput()->GetLargestPossibleRegion() ) ;
    ImageType::Pointer outputImage = ImageType::New() ;
    outputImage->CopyInformation( inputImage ) ;
    outputImage->SetRegions( inputImage->GetLargestPossibleRegion() ) ;
    outputImage->Allocate() ;
    IteratorWithIndexType out( outputImage , outputImage->GetLargestPossibleRegion() ) ;
    ImageType::IndexType index ;
    for( it1.GoToBegin() , it2.GoToBegin() , out.GoToBegin() ; !out.IsAtEnd() ; ++it1, ++it2, ++out )
    {
      index = out.GetIndex() ;
      if( ( ( (index[ 0 ] / mosaicStepSize)%2 + (index[ 1 ] / mosaicStepSize) )%2 + index[ 2 ]/mosaicStepSize )%2 )
      {
        out.Set( it1.Get() ) ;
      }
      else
      {
        out.Set( it2.Get() ) ;
      }
    }
    inputImage = outputImage ;
  } else if( skeletonOn )
  {
    // inputImage is 
    SkeletonFilterType::Pointer skeletonFilter = SkeletonFilterType::New();
    skeletonFilter->SetInput (inputImage);
    skeletonFilter->Update();
    inputImage = skeletonFilter->GetOutput() ;
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_skeleton");
  }
  else {
    cout << "NOTHING TO DO, no operation selected..." << endl;
    return EXIT_FAILURE ;
  }

  // add the extension to the outputfile name
  outFileName.append(format);

  // when outputFileName is set
  if (outputFileName) {
     outFileName = string(outputFileName);
  }

  // write image
  if (debug) cout << "writing output data " << outFileName << endl;

  if( diffusionImage )
    {
      typedef itk::ImageFileWriter< DiffusionImageType > DiffusionImageWriter ;
      DiffusionImageType::Pointer outputDiffusionImage ;
      outputDiffusionImage = dynamic_cast< DiffusionImageType* >( inputBaseImage.GetPointer() ) ;
      if( !outputDiffusionImage )
    {
      std::cerr << "Error saving output diffusion image" << std::endl ;
      return EXIT_FAILURE ;
    }
      DiffusionImageWriter::Pointer writer = DiffusionImageWriter::New() ;
      if( !nocompOn )
    {
      writer->UseCompressionOn() ;
    }
      writer->SetFileName( outFileName.c_str() );
      writer->SetInput( outputDiffusionImage ) ;
      try
    {
      writer->Write() ;
    }
      catch( ExceptionObject & err )
    {
      cerr << "ExceptionObject caught!" << endl ;
      cerr << err << endl ;
      return EXIT_FAILURE ;
    }

    }
  else
    {
      if (writeByte){
	castBinaryFilterType::Pointer castFilter = castBinaryFilterType::New();
	castFilter->SetInput(inputImage);
    try {
      castFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    }

    BinaryVolumeWriterType::Pointer writer = BinaryVolumeWriterType::New();
    if(!nocompOn)
      {
        writer->UseCompressionOn();
      }
    writer->SetFileName(outFileName.c_str());
    writer->SetInput(castFilter->GetOutput());
    writer->Write();
      } else if (writeFloat) {
    VolumeWriterType::Pointer writer = VolumeWriterType::New();
    if(!nocompOn)
      {
        writer->UseCompressionOn();
      }
    writer->SetFileName(outFileName.c_str());
    writer->SetInput(inputImage);
    writer->Write();
      } else {
    castShortFilterType::Pointer castFilter = castShortFilterType::New();
    castFilter->SetInput(inputImage);
    try
      {
        castFilter->Update();
      }
    catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;
      }

    ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
    writer->SetFileName(outFileName.c_str());
    if(!nocompOn)
      {
        writer->UseCompressionOn();
      }
    writer->SetInput(castFilter->GetOutput());
    writer->Write();
      }
    }
  return EXIT_SUCCESS;
}
