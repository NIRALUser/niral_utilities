/* 
 * compute image math and combinations
 *
 * author:  Martin Styner 
 * 
 * changes:
 *
 */

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4503 )
#endif

#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>

//#include "Disclaimer.h" 

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

#include <itkThresholdImageFilter.h> 
#include <itkBinaryThresholdImageFilter.h> 
#include <itkImageRegionIterator.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkOrientImageFilter.h>
#include <itkCastImageFilter.h>


#include "argio.h"
#include "MouseBrainOrient.h" 
#include "ImageMomentsCalculator.h"

#define DEFAULT_SAMP 2
// number of samples by default

using namespace std;
using namespace itk;

typedef unsigned short PixelType;
typedef float FloatType;
typedef short ShortPixelType;
typedef unsigned char BinaryPixelType;
enum { ImageDimension = 3 };
typedef Image<PixelType,ImageDimension>       ImageType;
typedef Image<PixelType,ImageDimension>       ShortImageType;
typedef Image<FloatType,ImageDimension>  	  FloatImageType;
typedef Image<BinaryPixelType,ImageDimension> BinaryImageType;
typedef ImageType::RegionType                 ImageRegionType;
typedef ImageRegionIterator< ImageType >      Iterator;
typedef ImageType::Pointer                    ImagePointer;
typedef ShortImageType::Pointer			      ShortImagePointer;

typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileWriter< ImageType >          VolumeWriterType;
typedef ImageFileWriter< BinaryImageType >    BinaryVolumeWriterType;
typedef ImageFileReader< ShortImageType >     ShortVolumeReaderType;
typedef ImageFileWriter< ShortImageType >     ShortVolumeWriterType;

typedef ExtractImageFilter<ImageType, ImageType> CropFilterType;
typedef ThresholdImageFilter< ImageType > maskThreshFilterType;

typedef BinaryThresholdImageFilter<ImageType, ImageType> ThreshFilterType;
typedef BinaryThresholdImageFilter<FloatImageType, ImageType> FloatThreshFilterType;
typedef ConnectedComponentImageFilter<FloatImageType, FloatImageType> ConnectiveFilterType;
typedef RelabelComponentImageFilter<FloatImageType, FloatImageType> RelabelFilterType;

typedef ImageMomentsCalculator<ImageType> ImageMomentsCalculatorType;
typedef LabelStatisticsImageFilter<ImageType, ImageType> LabelStatisticsImageFilterType;
typedef OrientImageFilter<ImageType, ImageType> OrientImageFilterType;
typedef CastImageFilter<ImageType, FloatImageType> castFloatImageType;
typedef CastImageFilter<FloatImageType, ImageType> castShortImageType;

static int debug = 0;

int main(const int argc, const char **argv)
{
	if (argc <= 1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help")) {
		cout << "MouseBrainOrient 1.1 version (Oct 2008)" << endl;
		cout << " computes base statistics of an image" << endl;
		cout << "usage: MouseBrainOrient infile [-showinfo] -tmin tmin -tmax tmax -idwi idwifile -outfile outfile" << endl;
		cout << "infile                 input dataset" << endl;;
		cout << "showinfo				show information about direction cosine, principal axes, and center of gravity" << endl;
	    cout << "idwi					reference IDWI file" << endl;
		cout << "tmin                   lower threshold to find mouse brain" << endl;
		cout << "tmax                   upper threshold to find mouse brain" << endl;
		cout << "outfile                output dataset" << endl;;
		cout << endl << endl;
		exit(0);
	}

	const int BGVAL = 0;
	const int FGVAL = 1;

	char *inputFileName = strdup(argv[1]);
	char *idwiFileName = ipGetStringArgument(argv, "-idwi", NULL);  
	char *outputFileName = ipGetStringArgument(argv, "-outfile", NULL);  
	bool bShowInfo = ipExistsArgument(argv, "-showinfo");
	char *base_string;

	int tmin = ipGetIntArgument(argv, "-tmin", 9000);
	int tmax = ipGetIntArgument(argv, "-tmax", 13000);

	ImagePointer inputImage, idwiImage, thresholdImage, conCompImage;

	// load image
	if (debug) {
		cout << "Loading file " << inputFileName << endl;
	}

	VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
	imageReader->SetFileName(inputFileName) ;
	imageReader->Update();
	inputImage = imageReader->GetOutput();

	VolumeReaderType::Pointer idwiReader = VolumeReaderType::New();
	idwiReader->SetFileName(idwiFileName);
	idwiReader->Update();
	idwiImage = idwiReader->GetOutput();


	// do threshold image
	ThreshFilterType::Pointer threshFilter = ThreshFilterType::New();
	threshFilter->SetInput(idwiImage);
	threshFilter->SetLowerThreshold(tmin);
	threshFilter->SetUpperThreshold(tmax);
	threshFilter->SetOutsideValue(BGVAL);
	threshFilter->SetInsideValue(FGVAL);
	threshFilter->Update();
	thresholdImage = threshFilter->GetOutput();

	ConnectiveFilterType::Pointer connective = ConnectiveFilterType::New();
	RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	FloatThreshFilterType::Pointer thresFilter2 = FloatThreshFilterType::New();

	castFloatImageType::Pointer floatType = castFloatImageType::New();
	floatType->SetInput(thresholdImage);
	floatType->Update();

	connective->SetInput(floatType->GetOutput());
	relabelFilter->SetInput(connective->GetOutput());
	thresFilter2->SetInput(relabelFilter->GetOutput());


	thresFilter2->SetLowerThreshold(0.1);
	thresFilter2->SetUpperThreshold(1);
	thresFilter2->SetOutsideValue(0.0);
	thresFilter2->SetInsideValue(1.0);
	thresFilter2->Update();
	conCompImage = thresFilter2->GetOutput();

	VolumeWriterType::Pointer writer = VolumeWriterType::New();
	writer->SetFileName("conComp.gipl.gz");
	writer->SetInput(conCompImage);
	writer->Write();

	//writer->SetFileName("connective.gipl.gz");
	//writer->SetInput(connective->GetOutput());
	//writer->Write();

	ImageMomentsCalculatorType::Pointer momentCalculator = ImageMomentsCalculatorType::New();
	momentCalculator->SetImage(thresholdImage);
	momentCalculator->Compute();

	ImageMomentsCalculatorType::MatrixType principalAxes = momentCalculator->GetPrincipalAxes();

	int dim[3];
	ImageRegionType imageRegion = inputImage->GetLargestPossibleRegion();
	dim[0] = imageRegion.GetSize(0);
	dim[1] = imageRegion.GetSize(1);
	dim[2] = imageRegion.GetSize(2);

	int cornerFrom[3] = { 0, 0, 0 };
	int cornerTo[3] = { dim[0], dim[1], dim[2] } ;
	LabelStatisticsImageFilterType::Pointer labelFilter = LabelStatisticsImageFilterType::New();
	labelFilter->SetInput(conCompImage);
	labelFilter->SetLabelInput(conCompImage);
	labelFilter->Update();
	LabelStatisticsImageFilterType::BoundingBoxType boundingBox = labelFilter->GetBoundingBox(1);

	cornerFrom[0] = boundingBox[0];
	cornerFrom[1] = boundingBox[2];
	cornerFrom[2] = boundingBox[4];
	cornerTo[0]   = boundingBox[1];
	cornerTo[1]   = boundingBox[3];
	cornerTo[2]   = boundingBox[5];

	ImageMomentsCalculatorType::VectorType firstMoment = momentCalculator->GetFirstMoments();
	for (int i = 0; i < 3; i++) {
		firstMoment[i] -= ((cornerFrom[i] + cornerTo[i]) / 2.0);
	}

	int columnOrder[3] = { 1, 0, 2};
	char orientation[4] = { 0, 0, 0, 0 };
	char orientationReference[3] = { 'P', 'R', 'S' };
	char orientationReference2[3] = { 'A', 'L', 'I' };

	ImageType::DirectionType direction;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (abs(principalAxes[j][i]) > 0.9) {
				if (firstMoment[i] > 0) {
					orientation[i] = orientationReference[j];
				} else {
					orientation[i] = orientationReference2[j];
				}
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		if (orientation[i] == 'R' || orientation[i] == 'L') {
			if (orientation[i] == 'L') {
				direction(i, 0) = 1;
			} else {
				direction(i, 0) = 1;
			}
			direction(i, 1) = 0;	
			direction(i, 2) = 0;	
		} else if (orientation[i] == 'P' || orientation[i] == 'A') {
			direction(i, 0) = 0; 
			if (orientation[i] == 'A') {
				direction(i, 1) = 1;	
			} else {
				direction(i, 1) = -1;	
			}
			direction(i, 2) = 0;	
		} else if (orientation[i] == 'S' || orientation[i] == 'I') {
			direction(i, 0) = 0; 
			direction(i, 1) = 0;	
			if (orientation[i] == 'I') {
				direction(i, 2) = 1;	
			} else {
				direction(i, 2) = -1;	
			}
		}
	}

	direction = direction.GetTranspose();
	cout << "Image Orientation: " << orientation << endl;

	if (bShowInfo) {
		cout << "Direction Cosine: " << endl << direction << endl;
		cout << "Principal Axes: " << endl << principalAxes << endl;
		cout << "Center of Gravity: " << endl << momentCalculator->GetCenterOfGravity() << endl;	
		return 0;
	}


	ImagePointer orientedImage = NULL;

	OrientImageFilterType::Pointer orienter = OrientImageFilterType::New();
    orienter->SetInput(inputImage);
	orienter->UseImageDirectionOff();
	orienter->SetGivenCoordinateDirection(direction);
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
    orienter->Update();
	orientedImage = orienter->GetOutput();

	writer->SetFileName(outputFileName);
	writer->SetInput(orientedImage);
	writer->Write();

	return 0;
}

