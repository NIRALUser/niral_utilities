/* 
 * compute image statistics 
 *
 * author:  Martin Styner 
 * 
 * changes:  this tool is programmed very UGLY, sorry in advance to anyone having to edit this thing
 *  - we should probably nicely re-write this whole thing, with classes, using more ITK, etc 
 *   AND make it GenerateCLP compatible...
 */

#include <iostream>
#include <fstream>

using namespace std;

#include <string.h>
#include <sys/types.h>
#include <stdlib.h>    // for exit, system
#include <math.h>
#include <errno.h>
#include <vector>
#include <algorithm>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkExtractImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageIOBase.h>

#include "itkVector.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkMembershipSample.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "itkSubsample.h"
#include "itkStatisticsAlgorithm.h"

#include "itkListSample.h"

#include "argio.h"
#include "ImageStat.h" 

#include "itkHausdorffDistanceImageFilter.h"

#define DEFAULT_SAMP 2
// number of samples by default

using namespace std;
using namespace itk;

typedef float PixelType;
typedef float ExtractPixelType;
typedef unsigned char Write2DPixelType;
typedef unsigned char BinaryPixelType;
enum { ImageDimension3 = 3, ImageDimension2 = 2 };
typedef Image<PixelType,ImageDimension3>       ImageType;
typedef Image<BinaryPixelType,ImageDimension3> BinaryImageType;
typedef ImageType::RegionType                 ImageRegionType;
typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileReader< BinaryImageType >    BinaryVolumeReaderType;
typedef ImageRegionIterator< ImageType >      Iterator;
typedef ImageRegionIterator< BinaryImageType> BinaryIterator;
typedef ImageType::Pointer                    ImagePointer;
typedef BinaryImageType::Pointer              BinaryImagePointer;

typedef Image<ExtractPixelType,ImageDimension2>      ExtractImageType;
typedef ExtractImageFilter< ImageType, ExtractImageType > ExtractFilterType ;
typedef Image<Write2DPixelType,ImageDimension2>      Write2DImageType;
typedef RescaleIntensityImageFilter< ExtractImageType , Write2DImageType > RescaleIntImageType;
typedef ImageFileWriter< Write2DImageType >   ExtractImageWriterType;

typedef itk::Statistics::ImageToListSampleAdaptor<ImageType>	ImageSampleType;
typedef ImageSampleType::Iterator				ImageSampleIterator;

typedef itk::Vector<PixelType,2> MeasurementVectorType;
typedef itk::Statistics::ListSample <MeasurementVectorType> SampleType;
typedef itk::Statistics::MembershipSample< SampleType >	MembershipSampleType;

typedef itk::Vector<PixelType, 3 > VectorType;
typedef itk::MinimumMaximumImageCalculator<ImageType > MinMaxCalcType;
typedef itk::Statistics::Subsample<MembershipSampleType::ClassSampleType> SubsampleType;
		
static void histo_vol(ImagePointer inputImage, BinaryImagePointer mask,char * histofile,
		      double min, double max,double  smin,double  smax, int samp);

static void histo_probvol(ImagePointer inputImage, ImagePointer mask, short probFactor, char * histofile,
			  double min, double max,double  smin,double  smax, int samp);
			  
bool searchList(int label, vector<int> list)
{
	int size = list.size();
	bool found = 0;
	for (int i=0; (i<size) && (found == 0) ; i++) if (list[i]==label) found = 1;
	return found;
}

static int debug = 0;

int main(const int argc, const char **argv)
{
	if (argc <=1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help")) {
		std::cout << "ImageStat 1.6 version (4/17/14)" << std::endl;
		std::cout << " computes base statistics of an image" << std::endl;
		std::cout << "usage: ImageStat infile [-outbase outbase] [-mask maskfile]  [-v]"<< std::endl;
		std::cout << "       [-histo [-min xmin] [-max xmax] [-samp s]] " << std::endl;
		std::cout << std::endl;
		std::cout << "infile         input dataset" << std::endl;;
		std::cout << "-outbase outbase     base-outputfilename, if omitted then the same as base of input" << std::endl; 
		std::cout << "-mask file     Mask [bytedata, 0/255 value] to mask out image data, where mask is set" << std::endl;
		std::cout << "-probmask file Mask [shortdata, 0/62000 value] to mask out probabilitically image data, where mask is > 0" << std::endl;
		std::cout << "-probfactor f  Maximum probability value/Probability normalization factor, needs to be set if using -probmask" << std::endl;
		std::cout << "               The proper probability of all values in probmask is computed by division with this value" << std::endl;
		std::cout << "-verbose       verbose mode " << std::endl << std::endl; 
		std::cout << "-histo         Histogram corrected and uncorrected for pixeldimension" << std::endl;
		std::cout << "-min xmin      Minimal value for histogram/volume/labels [DEFAULT masked_minval]" << std::endl; 
		std::cout << "-max xmax      Maximal value for histogram/volume/labels [DEFAULT masked_maxval]" << std::endl;
		std::cout << "-samp s        Number of Samples for histogram/volume [DEFAULT xmax-xmin]" << std::endl << std::endl; 
		std::cout << "-threeSlice[C x,y,z]   create three orthogonal views through the image center or the optionally supplied coordinate" << std::endl;
		std::cout << "-info          print image info on standard out" << std::endl;
		std::cout << "-label labelfile    give volume, mean intensity, standard deviation, min, max, quantiles for each label" << std::endl;
		std::cout << "      -display   display the results on the standard output (if -volumeSummary and/or -intensitySummary used, display theresults of these options)" << std::endl;
		std::cout << "      -quantile [list of numbers separate by a coma]   quantile value [DEFAULT: 1 5 33 50 66 95 99]" << std::endl;
		std::cout << "      -maskForLabel maskForLabelfile   binary mask" << std::endl;
		std::cout << "      -probabilityMap probabilitymap   use probability map to compute volume, mean intensity etc." << std::endl;
		std::cout << "      -volumeSummary   give volume for each label" << std::endl;
		std::cout << "      -fillEmpty       provide feedback for all possible labels from min to max (see -min, -max)"<< std::endl;
		std::cout << "      -intensitySummary   give mean intensity, standard deviation, min, max, quantiles for each label" << std::endl;
		std::cout << "      -volOverlap labelFile  compute Tanimoto (intersection over union), Dice coefficient and Hausdorff distance of input label image with additional image" << std::endl;
		std::cout << "      Caution: option min/max currently only work for -label" << std::endl;

		std::cout << std::endl << std::endl;
		exit(0);
	}
	
	debug                  = ipExistsArgument(argv, "-verbose");
	bool histOn            = ipExistsArgument(argv, "-histo");
	bool threeSliceOn      = ipExistsArgument(argv, "-threeSlice");
	bool infoOn            = ipExistsArgument(argv, "-info");
	bool labelOn           = ipExistsArgument(argv, "-label");
	bool displayOn         = ipExistsArgument(argv, "-display");
	bool volOverlapOn      = ipExistsArgument(argv, "-volOverlap");
	bool quantileOn        = ipExistsArgument(argv, "-quantile");
	bool maskForLabelOn    = ipExistsArgument(argv, "-maskForLabel");
	bool probabilityMapOn  = ipExistsArgument(argv, "-probabilityMap");
	bool volumeSummaryOn   = ipExistsArgument(argv, "-volumeSummary");
	bool fillEmptyOn       = ipExistsArgument(argv, "-fillEmpty");
	bool minOn             = ipExistsArgument(argv, "-min");
	bool maxOn             = ipExistsArgument(argv, "-max");
	bool intensitySummaryOn  = ipExistsArgument(argv, "-intensitySummary");

	char *inputFileName = strdup(argv[1]);
	char *outbase    = ipGetStringArgument(argv, "-outbase", NULL);  
	char *maskfile   = ipGetStringArgument(argv, "-mask", NULL);
	char *probmaskfile   = ipGetStringArgument(argv, "-probmask", NULL);
	char *label = ipGetStringArgument(argv, "-label", NULL);

	double xmin = ipGetDoubleArgument(argv, "-min", 32000);
	double xmax = ipGetDoubleArgument(argv, "-max", 0);

	double *quant_tab = NULL; 
	int tabSize = 0;
	if (quantileOn) {
		char * tmp_str      = ipGetStringArgument(argv, "-quantile", NULL);
		if (tmp_str) {
			quant_tab = new double[100];
			tabSize = ipExtractDoubleTokens(quant_tab, tmp_str, 100);
		}
	}
	else
	{
		tabSize = 7;
		quant_tab = new double[tabSize];
		quant_tab[0] = 1;
		quant_tab[1] = 5;
		quant_tab[2] = 33;
		quant_tab[3] = 50;
		quant_tab[4] = 66;
		quant_tab[5] = 95;
		quant_tab[6] = 99;
	}
	
	char *maskForLabel = ipGetStringArgument(argv, "-maskForLabel", NULL);
	
	char *probabilityMap = ipGetStringArgument(argv, "-probabilityMap", NULL);
	
	short probfactor = ipGetIntArgument(argv, "-probfactor", 0);
	if (probmaskfile && probfactor == 0) {
		std::cout << " when using -probmask for a probability mask, then -probfactor has to be supplied for normalization of the probability values" << std::endl;
		exit(-1);
	}
	char *volOverlapFile = ipGetStringArgument(argv, "-volOverlap", NULL);
	
	char * base_string;
	if (!outbase) {
		base_string = strdup(ipGetBaseName(inputFileName));
	} else {
		base_string = outbase;
	}

	ImagePointer inputImage ;
	BinaryImagePointer maskImage;
	ImagePointer  probmaskImage;
	ImagePointer  labelImage;
	ImagePointer  volOverlapImage;
	ImagePointer maskForLabelImage;
	ImagePointer probabilityMapImage;
	ImageSampleType::Pointer imageSample;
	
	vector<int> labelList;
  
	int sliceX, sliceY, sliceZ;
	sliceX = sliceY = sliceZ = -1;
	if (threeSliceOn) {
		// are there any arguments to threeSlice
		char * tmp_str      = ipGetStringArgument(argv, "-threeSliceC", NULL);
		int textend[3];
		if (tmp_str) {
			int num = ipExtractIntTokens(textend, tmp_str, 3);
			if (3 == num) {
				sliceX = textend[0];
				sliceY = textend[1];
				sliceZ = textend[2];
			}
		}
	}

	try
	{
		// load image
		if (debug) std::cout << "Loading file " << inputFileName << std::endl;
		VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
		imageReader->SetFileName(inputFileName) ;
		imageReader->Update();
		if (debug) std::cout << "Loading file done " << std::endl;
		inputImage = imageReader->GetOutput();
		//read mask
		if (maskfile) {
			if (debug) std::cout << "Loading mask " << maskfile  << std::endl;
			BinaryVolumeReaderType::Pointer maskReader = BinaryVolumeReaderType::New();
			maskReader->SetFileName(maskfile) ;
			maskReader->Update();
			maskImage = maskReader->GetOutput();
		}

		//read prob mask
		if (probmaskfile) {
			if (debug) std::cout << "Loading probability mask " << maskfile  << std::endl;
			VolumeReaderType::Pointer probaReader = VolumeReaderType::New();
			probaReader->SetFileName(probmaskfile) ;
			probaReader->Update();
			probmaskImage = probaReader->GetOutput();
		}
		//read label
		if (labelOn) {
			if (debug) std::cout << "Loading Label File " << label << std::endl;
			VolumeReaderType::Pointer imageReader2 = VolumeReaderType::New();
			imageReader2->SetFileName(label) ;
			try
			{
				imageReader2->Update();
			}
			catch (ExceptionObject & err) 
			{
				cerr<<"ExceptionObject caught!"<<endl;
				cerr<<err<<endl;
				return EXIT_FAILURE;	
			}
			labelImage = imageReader2->GetOutput();
		}
		if (maskForLabelOn) {
			if (debug) std::cout << "Loading Mask File " << maskForLabel << std::endl;
			VolumeReaderType::Pointer imageReader3 = VolumeReaderType::New();
			imageReader3->SetFileName(maskForLabel) ;
			try 
			{
				imageReader3->Update();
			}
			catch (ExceptionObject & err) 
			{
				cerr<<"ExceptionObject caught!"<<endl;
				cerr<<err<<endl;
				return EXIT_FAILURE;	
			}
			maskForLabelImage = imageReader3->GetOutput();
		}
		if (probabilityMapOn) {
			if (debug) std::cout << "Loading Probability Map" << probabilityMap << std::endl;
			
			VolumeReaderType::Pointer imageReader4 = VolumeReaderType::New();
			imageReader4->SetFileName(probabilityMap) ;
			try 
			{
				imageReader4->Update();
			}
			catch (ExceptionObject & err) 
			{
				cerr<<"ExceptionObject caught!"<<endl;
				cerr<<err<<endl;
				return EXIT_FAILURE;	
			}
			probabilityMapImage = imageReader4->GetOutput();
			
		}
		if (volOverlapFile) {
			if (debug) std::cout << "Loading Volume Overlap Label File " << volOverlapFile  << std::endl;
			VolumeReaderType::Pointer overlapReader = VolumeReaderType::New();
			overlapReader->SetFileName(volOverlapFile) ;
			overlapReader->Update();
			volOverlapImage = overlapReader->GetOutput();
		}
		if (histOn) {
			//Get extrema over volume/slice/mask
		        // TODO: RESPECT min and max provided on command line
			double min = +1000000, max = -1000000;
			double smin = min;
			double smax = max;
    
			Iterator iterImage (inputImage, inputImage->GetBufferedRegion());
			//initialize
			iterImage.GoToBegin();
    
			if (maskImage.IsNotNull()) {
				BinaryIterator iterMask (maskImage, maskImage->GetBufferedRegion());
				while ( !iterImage.IsAtEnd() )  {
					PixelType value =  iterImage.Get();
					BinaryPixelType Maskvalue =  iterMask.Get();
					if (min > value) min = value;
					if (max < value) max = value;
					if (Maskvalue) {
						if (smin > value) smin = value;
						if (smax < value) smax = value;
					}
					++iterMask;
					++iterImage;
				}
			} else {
				// for both with or without probmaskImage this needs to be done
				while ( !iterImage.IsAtEnd() )  {
					PixelType value =  iterImage.Get();
					if (min > value) min = value;
					if (max < value) max = value;
					++iterImage;
				}
				if (min == max) max = min + 1;
				smin = min;
				smax = max;
			}
			smin = ipGetDoubleArgument(argv, "-min", smin);
			smax = ipGetDoubleArgument(argv, "-max", smax);
    

			if (smin > smax)  {double val = smax; smax = smin; smin = val; }
			if (smin < min)   smin = min;
			if (smax > max)   smax = max;
			if (debug) std::cout << "Min: " << smin << ", Max: " << smax << std::endl;

			char histofile [5024];
			sprintf(histofile,"%s_vol.txt", base_string);
    
			int samp        = ipGetIntArgument(argv, "-samp", (int) (smax-smin+1));
			if (samp <= 1)  samp = DEFAULT_SAMP;
    
    
			if (probmaskImage.IsNotNull()) {
				histo_probvol(inputImage,probmaskImage,probfactor,histofile,min,max,smin,smax,samp);
			} else {
				histo_vol(inputImage,maskImage,histofile,min,max,smin,smax,samp);
			}
		}
		if (threeSliceOn) {
			ImageRegionType extractRegion;
			ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
			extractFilter->SetInput(inputImage);
			RescaleIntImageType::Pointer rescaleFilter = RescaleIntImageType::New();
			rescaleFilter->SetInput(extractFilter->GetOutput());
			rescaleFilter->SetOutputMinimum(0);
			rescaleFilter->SetOutputMaximum(255);
			ExtractImageWriterType::Pointer writer = ExtractImageWriterType::New();
			writer->SetInput(rescaleFilter->GetOutput());

			string outfileBase(base_string);
			string outfileName;
			int sliceDir, numSlices; 

			sliceDir = 0;
			extractRegion = (inputImage->GetLargestPossibleRegion());
			numSlices = extractRegion.GetSize(sliceDir);
			extractRegion.SetIndex(0,0);
			extractRegion.SetIndex(1,0);
			extractRegion.SetIndex(2,0);
			if (sliceX >= 0) {
				extractRegion.SetIndex(sliceDir, sliceX);
			} else {
				extractRegion.SetIndex(sliceDir, numSlices / 2);
			}
			extractRegion.SetSize(sliceDir, 0);
			extractFilter->SetExtractionRegion(extractRegion);

			outfileName = outfileBase + "_sliceX.tiff";
			writer->SetFileName(outfileName.c_str()); 
			writer->Write();
    
			sliceDir = 1;
			extractRegion = (inputImage->GetLargestPossibleRegion());
			numSlices = extractRegion.GetSize(sliceDir);
			extractRegion.SetIndex(0,0);
			extractRegion.SetIndex(1,0);
			extractRegion.SetIndex(2,0);
			if (sliceY >= 0) {
				extractRegion.SetIndex(sliceDir, sliceY);
			} else {
				extractRegion.SetIndex(sliceDir, numSlices / 2);
			}
			extractRegion.SetSize(sliceDir, 0);
			extractFilter->SetExtractionRegion(extractRegion);
			outfileName = outfileBase + "_sliceY.tiff";
			writer->SetFileName(outfileName.c_str()); 
			writer->Write();

			sliceDir = 2;
			extractRegion = (inputImage->GetLargestPossibleRegion());
			numSlices = extractRegion.GetSize(sliceDir);
			extractRegion.SetIndex(0,0);
			extractRegion.SetIndex(1,0);
			extractRegion.SetIndex(2,0);
			if (sliceZ >= 0) {
				extractRegion.SetIndex(sliceDir, sliceZ);
			} else {
				extractRegion.SetIndex(sliceDir, numSlices / 2);
			}
			extractRegion.SetSize(sliceDir, 0);
			extractFilter->SetExtractionRegion(extractRegion);
			outfileName = outfileBase + "_sliceZ.tiff";
			writer->SetFileName(outfileName.c_str()); 
			writer->Write();
		}
		if (infoOn) {
      
			MinMaxCalcType::Pointer minmaxCalc = MinMaxCalcType::New();
			minmaxCalc->SetImage(inputImage);
			minmaxCalc->Compute();
			double min = minmaxCalc->GetMinimum();
			double max = minmaxCalc->GetMaximum();

      			//Get the voxel type
			itk::ImageIOBase::IOComponentType input_component_type = imageReader->GetImageIO()->GetComponentType();
			std::string input_component_string = imageReader->GetImageIO()->GetComponentTypeAsString(input_component_type);
      
			ImageRegionType region;
			region = (inputImage->GetLargestPossibleRegion());
			std::cout << "Filename: " << inputFileName << std::endl;
			std::cout << "Dims: " << region.GetSize(0) << " "
					<< region.GetSize(1) << " "
					<< region.GetSize(2) << " " << std::endl;
			std::cout << "Pixdims: " << (inputImage->GetSpacing())[0] << " "
					<< (inputImage->GetSpacing())[1] << " "
					<< (inputImage->GetSpacing())[2] << std::endl;
			std::cout << "Origin: " << (inputImage->GetOrigin())[0] << " "
					<< (inputImage->GetOrigin())[1] << " "
					<< (inputImage->GetOrigin())[2] << std::endl;
			std::cout << "Minimum .. Maximum : " << min << " .. "  << max << std::endl;
			std::cout << "Voxel type: " << input_component_string << std::endl;
		}

		if (labelOn){

		  	MinMaxCalcType::Pointer minmaxCalc = MinMaxCalcType::New();
			minmaxCalc->SetImage(labelImage);
			minmaxCalc->Compute();
			int maxLabel = minmaxCalc->GetMaximum();
			int minLabel = minmaxCalc->GetMinimum();

			if (minOn) minLabel = xmin;
			if (maxOn) maxLabel = xmax;

			ImageSampleType::Pointer labelSample = ImageSampleType::New();
			labelSample->SetImage(labelImage); //initializing to avoid compilation warning
			PixelType labelValue = 0;
			ImageSampleIterator iterLabel = labelSample->Begin() ;
			if(maskForLabelOn)
			{
				Iterator itMaskForLabel( maskForLabelImage, maskForLabelImage->GetRequestedRegion());
				while( iterLabel != labelSample->End() )
				{				
					labelValue = iterLabel.GetMeasurementVector()[0];
					if ((labelValue != 0) && (!searchList(labelValue,labelList))
					    && (itMaskForLabel.Get() != 0)) {
					  labelList.push_back( labelValue );
					}
					++iterLabel;
					++itMaskForLabel;
				}
			}
			else
			{
				while(iterLabel != labelSample->End())
				{					
				  labelValue = iterLabel.GetMeasurementVector()[0];
				  if ((labelValue != 0) && (!searchList(labelValue,labelList))) {
				    labelList.push_back(labelValue);
				  }
				  ++iterLabel;
				}
			}
			sort(labelList.begin(), labelList.end());
			int nbLabels = labelList.size();
			SampleType::Pointer sample = SampleType::New();
			MeasurementVectorType mv;
			Iterator itIm( inputImage, inputImage->GetRequestedRegion() );
			itIm.GoToBegin();
			if(maskForLabelOn)
			{
				Iterator itMask( maskForLabelImage, maskForLabelImage->GetRequestedRegion());
				itMask.GoToBegin();
				while((!itIm.IsAtEnd())&&(!itMask.IsAtEnd()))
				{
					mv[0] = itIm.Get();
					mv[1] = itMask.Get();
					sample -> PushBack(mv);
					++itIm;
					++itMask;
				}
			}
			else
			{
				if(probabilityMapOn)
				{
					Iterator itProbMap( probabilityMapImage, probabilityMapImage->GetRequestedRegion());
					itProbMap.GoToBegin();
					while((!itIm.IsAtEnd())&&(!itProbMap.IsAtEnd()))
					{
						mv[0] = itIm.Get();
						mv[1] = itProbMap.Get();
						sample -> PushBack(mv);
						++itIm;
						++itProbMap;
					}
				}
				else
				{
					while(!itIm.IsAtEnd())
					{
						mv[0] = itIm.Get();
						mv[1] = 1;
						sample -> PushBack(mv);
						++itIm;
					}
				}
			}
			MembershipSampleType::Pointer membershipSample = MembershipSampleType::New();
			membershipSample->SetSample(sample);
			membershipSample->SetNumberOfClasses(nbLabels);

			SampleType::Iterator iter = sample->Begin();
			iterLabel = labelSample->Begin();
			while( iter != sample->End() )
			{				
			  labelValue = iterLabel.GetMeasurementVector()[0];
			  if ((labelValue != 0 ) && ((iter.GetMeasurementVector()[1] != 0) ))
			    {
			      membershipSample->AddInstance(labelValue,iter.GetInstanceIdentifier());
			    }
			  ++iter ;
			  ++iterLabel;
			}

			//Write stat infos in the output file
			int nbPix;
			int l;
			VectorType spacing = inputImage->GetSpacing ();
			float pixVol = spacing[0] * spacing[1] * spacing[2];

			double tab [5][labelList.size()];
			double quantiles_tab[labelList.size()][tabSize];

			for (int j=0; j<(int)labelList.size() ; j++ )
			{
				l=labelList[j];
									      
				nbPix = membershipSample->GetClassSample(l)->Size();
				MembershipSampleType::ClassSampleType::ConstPointer classSample;
				classSample = membershipSample->GetClassSample(l);
		
				SubsampleType::Pointer Subsample = SubsampleType::New();
	
				Subsample->SetSample(classSample);
				Subsample->InitializeWithAllInstances();
		
				SubsampleType::Iterator SubsampleIter = Subsample->Begin();
	
				itk::Statistics::Algorithm::HeapSort<SubsampleType>(Subsample,0,0,Subsample->Size());
				SubsampleIter = Subsample->Begin();

				// mean computation
				double totalProb = 0;
				double meanSum = 0;
				for (unsigned int k=0; k<(Subsample->Size()); k++)
				{
					totalProb = totalProb + (Subsample->GetMeasurementVectorByIndex(k)[1]);
					meanSum = meanSum + ((Subsample->GetMeasurementVectorByIndex(k)[0])*(Subsample->GetMeasurementVectorByIndex(k)[1]));
				}
				double mean = (1/totalProb)*meanSum ;
			
				// standart deviation computation
		
				double sdSum = 0;
				for (unsigned int k=0; k<(Subsample->Size()); k++)
				{
					sdSum = sdSum + ((Subsample->GetMeasurementVectorByIndex(k)[0]) - mean)*((Subsample->GetMeasurementVectorByIndex(k)[0]) - mean)*(Subsample->GetMeasurementVectorByIndex(k)[1]);
				}
				double standardDeviation = sqrt((1/totalProb)*sdSum) ;
			
				tab[0][j]=totalProb*pixVol;
				tab[1][j]=mean;
				tab[2][j]=standardDeviation;
				tab[3][j]=Subsample->GetMeasurementVectorByIndex(0)[0];
				tab[4][j]=Subsample->GetMeasurementVectorByIndex(Subsample->Size()-1)[0];

				//quantiles
				float size = Subsample->Size();
				for(int i=0; i<tabSize; i++)
				{
					float quantValue = size*(quant_tab[i])/100;
					if(floor(quantValue) == 0)
					{
						quantiles_tab[j][i]=Subsample->GetMeasurementVectorByIndex((int)(floor(quantValue)))[0];
					}
					else
					{
						if (((quantValue - floor(quantValue)) == 0)||((quantValue - floor(quantValue)) < 0.5))
						{
							quantiles_tab[j][i]=Subsample->GetMeasurementVectorByIndex((int)((floor(quantValue))-1))[0];
						}
						else
						{
							if((quantValue - floor(quantValue)) > 0.5)
							{
								quantiles_tab[j][i]=Subsample->GetMeasurementVectorByIndex((int)((floor(quantValue))))[0];
							}
							else
							{
								quantiles_tab[j][i]=((Subsample->GetMeasurementVectorByIndex((int)(floor(quantValue)) -1 )[0]) + (Subsample->GetMeasurementVectorByIndex((int)(floor(quantValue)))[0]))/2;
							}
						}
					}
				}
			}

			std::string input= inputFileName;
			size_t found;
			found=input.find_last_of("/\\");
			std::string inputFileNameTail =input.substr(found+1);
			if (!(volumeSummaryOn || intensitySummaryOn))
			{
			  if (debug) std::cout << "Computing Label stats" << label  << std::endl;
			  char statfile [5024];
			  sprintf(statfile,"%s_stat.txt", base_string);
			  ofstream efile(statfile, ios::out);  
			  if (!efile) {
			    cerr << "Error: open of file \"" << statfile << "\" failed." << std::endl;
			    exit(-1);
			  }
			  efile.precision(6);
			  efile<<"##########################################################################"<<endl;
			  efile<<"# File format :"<<endl;
			  efile<<"# LABEL \t VOLUME \t MEAN \t STD \t MIN \t MAX \t QUANTILES"<<endl;
			  efile<<"# Fields :"<<endl;
			  efile<<"#\t LABEL         Label number"<<endl;
			  efile<<"#\t VOLUME        Volume of voxels that have that label in cubic mm"<<endl;
			  efile<<"#\t MEAN          Mean intensity of those voxels"<<endl;
			  efile<<"#\t STD            Standard deviation of those voxels"<<endl;
			  efile<<"#\t MIN           Min intensity of those voxels"<<endl;
			  efile<<"#\t MAX           Max intensity of those voxels"<<endl;
			  efile<<"#\t QUANTILES     Quantile values [DEFAULT: 1 5 33 50 66 95 99]"<<endl;
			  efile<<"##########################################################################"<<endl<<endl; 
			  for (int j=0; j<(int)labelList.size() ; j++ )
			    {
			      l=labelList[j];
			      if ( (!minOn || (minOn && l >= minLabel)) && 
				   (!maxOn || (maxOn && l <= maxLabel)))
				{
				  efile << l;
				  efile << "\t \t" << tab[0][j];
				  efile << "\t \t" << tab[1][j];
				  efile << "\t \t" << tab[2][j];
				  efile << "\t \t" << tab[3][j];
				  efile << "\t \t" << tab[4][j];
				  if(displayOn)
				    {
				      std::cout << l;
				      std::cout << "\t \t" << tab[0][j];
				      std::cout << "\t \t" << tab[1][j];
				      std::cout << "\t \t" << tab[2][j];
				      std::cout << "\t \t" << tab[3][j];
				      std::cout << "\t \t" << tab[4][j];
				    }
				  for(int i=0; i<tabSize; i++)
				    {
				      efile << "\t \t" << quantiles_tab[j][i];
				      if(displayOn)
					efile << "\t \t" << quantiles_tab[j][i];	
				    }
				  efile << std::endl;
				  if(displayOn)
				    std::cout<<endl;
				}
			      efile.close();
			    }
			}
			if (volumeSummaryOn){
			  char statfile [5024];
			  sprintf(statfile,"%s_volumeSummary.csv", base_string);
			  ofstream efile(statfile, ios::out);  
			  if (!efile) {
			    cerr << "Error: open of file \"" << statfile << "\" failed." << std::endl;
			    exit(-1);
			  }
			  efile.precision(6);
			  efile << inputFileNameTail << ",VOLUME";
			  if(displayOn)
			    std::cout << inputFileNameTail << ",VOLUME";
			  // display results for all labels, even if some of them are empty
			  int i = 0;
			  for (int Label = 1; Label <= maxLabel; Label++)
			    {
			      if (labelList[i] == Label)
				{
				  efile <<","<< tab[0][i];
				  if(displayOn)
				    std::cout <<","<< tab[0][i];
				  i++;
				}
			      else
				{
				  efile <<",0";
				  if(displayOn)
				    std::cout <<",0";
				}
			    }
			  efile << std::endl;
			  if(displayOn)
			    std::cout<<endl;
			  efile.close();
			}
			if(intensitySummaryOn){
			  char statfile [5024];
			  sprintf(statfile,"%s_intensitySummary.csv", base_string);
			  ofstream efile(statfile, ios::out);  
			  if (!efile) {
			    cerr << "Error: open of file \"" << statfile << "\" failed." << std::endl;
			    exit(-1);
			  }
			  efile.precision(6);
			  std::string rowname [11]={"MEAN","STD","MIN","MAX","QUANTILES 1%","QUANTILES 5%","QUANTILES 33%","QUANTILES 50%","QUANTILES 66%","QUANTILES 95%","QUANTILES 99%"};
			  for (int i=0;i<4;i++)
			    {
			      efile << inputFileNameTail;
			      efile<<","<<rowname [i];
			      if(displayOn)
				{
				  std::cout << inputFileNameTail;
				  std::cout<<","<<rowname [i];
				}
			      // display results for all labels, even if some of them are empty
			      int j = 0;
			      for (int Label = 1; Label <= maxLabel; Label++ )
				{
				  if (labelList[j] == Label)
				    {
				      efile<<","<<tab[i+1][j];
				      if(displayOn)
					std::cout<<","<<tab[i+1][j];
				      j++;
				    }
				  else
				    {
				      efile<<",0";
				      if(displayOn)
					std::cout<<",0";
				    }
				}
			      efile<<endl;
			      if(displayOn)
				std::cout<<endl;
			    }
			  for(int k=0; k<tabSize; k++)
			    {
			      efile << inputFileNameTail;
			      efile<<","<<rowname [k+4];
			      if(displayOn)
				{
				  std::cout << inputFileNameTail;
				  std::cout<<","<<rowname [k+4];
				}
			      // display results for all labels, even if some of them are empty
			      int j = 0;
			      for (int Label = 1; Label <= maxLabel; Label++ )
				{
				  if (labelList[j] == Label)
				    {
				      efile<<","<<quantiles_tab[j][k];
				      if(displayOn)
					std::cout<<","<<quantiles_tab[j][k];
				      j++;
				    }
				  else
				    {
				      efile<<",0";
				      if(displayOn)
					std::cout<<",0";
				    }
				}
			      efile<<endl;
			      if(displayOn)
				std::cout<<endl;
			    }
			  efile.close();
			}
			delete[] quant_tab;
		}
		if (volOverlapOn){
		  if (debug) std::cout << "Computing Overlap stats" << label  << std::endl;
		  char statfile [5024];
		  sprintf(statfile,"%s_volOverlap.txt", base_string);
		  ofstream efile(statfile, ios::out);
		  if (!efile) {
		    cerr << "Error: open of file \"" << statfile << "\" failed." << std::endl;
		    exit(-1);
		  }
		  efile.precision(6);
		  efile<<"##########################################################################"<<endl;
		  efile<<"# file: "<<statfile << std::endl;
		  efile<<"# contains volume overlap between " << inputFileName << std::endl;
		  efile<<"                              and " << volOverlapFile<< std::endl;
		  efile<<"##########################################################################"<<endl;
		  
		  int volIntersection, volUnion, volA, volB;
      
			volIntersection = volUnion = volA = volB = 0;
      
			Iterator iterImageA (inputImage, inputImage->GetBufferedRegion());
			Iterator iterImageB (volOverlapImage, volOverlapImage->GetBufferedRegion());
    
			while ( !iterImageA.IsAtEnd() )  {
				PixelType valueA =  iterImageA.Get();
				PixelType valueB =  iterImageB.Get();
	
				if (valueA != 0)  volA++;
				if (valueB != 0)  volB++;
				if (valueA != 0 && valueB != 0) volIntersection++;
				if ((valueA != 0 && valueB == 0) ||
								 (valueA == 0 && valueB != 0) ) volUnion++;

				++iterImageA;
				++iterImageB;
			}

			typedef itk::HausdorffDistanceImageFilter<ImageType,ImageType> HausdorffFilterType;
			HausdorffFilterType::Pointer hausdorff = HausdorffFilterType::New();

			hausdorff->SetInput1( inputImage );
			hausdorff->SetInput2( volOverlapImage );
			hausdorff->Update();
			double hausdorffDistance = hausdorff->GetHausdorffDistance() ;
			

			volUnion = volUnion + volIntersection;
			efile << "Union = " << volUnion << std::endl;
			efile << "Intersection = " << volIntersection << std::endl;
			efile << "volA = " << volA << std::endl;
			efile << "volB = " << volB << std::endl;
			efile<< "DiceCoef = " << 2.0 * (double) volIntersection / (volA + volB) << std::endl;
			efile<< "TanimotoCoef = " << (double) volIntersection / volUnion << std::endl;
			efile << "HausdorffDist = " << hausdorffDistance << std::endl ;
			efile.close();

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

//**********************************************************************
//**********************************************************************
//**********************************************************************

static void histo_vol(ImagePointer inputImage, BinaryImagePointer maskImage,char * histofile,
		      double min, double max,double  smin,double  smax, int samp)
{
	double avg = 0 ;
	double savg = 0, snum = 0, num = 0;
	double var = 0, svar = 0;
  
	if (debug) std::cout << "computing histogram" << std::endl;

	Iterator iterImage (inputImage, inputImage->GetBufferedRegion());
	//Initialize
	iterImage.GoToBegin();

  // compute average
	if (maskImage.IsNotNull()) {
		BinaryIterator iterMask (maskImage, maskImage->GetBufferedRegion());	
		//Initialize
		iterMask.GoToBegin();
		while ( !iterImage.IsAtEnd() )  {
			PixelType value =  iterImage.Get();
			BinaryPixelType Maskvalue =  iterMask.Get();
			avg += value;
			num++;
			if (Maskvalue && value >= smin && value <= smax) {
				savg += value;
				snum++;
			}
			++iterMask;
			++iterImage;
		}
	} else {
		while ( !iterImage.IsAtEnd() )  {
			PixelType value =  iterImage.Get();
			avg += value;
			num++;
			if (value >= smin && value <= smax) {
				savg += value;
				snum++;
			}
			++iterImage;
		}
	}
	
	avg /= num;
	if (snum) savg /= snum;
    
  // compute variance and histogram
	iterImage.GoToBegin();
	float histogram[40000];
	double inc;
	inc = (smax - smin) / (double) (samp-1);
	int index;
    
	if (maskImage.IsNotNull()) {
		BinaryIterator iterMask (maskImage, maskImage->GetBufferedRegion());
		while ( !iterImage.IsAtEnd() )  {
			PixelType value =  iterImage.Get();
			BinaryPixelType Maskvalue =  iterMask.Get();
			var += (avg - value) * (avg - value);
			if (Maskvalue && value >= smin && value <= smax) {
				svar += (savg - value) * (savg - value);
				index = (int) ((value - smin) / inc+0.5);
				if (index < 0 || index > (samp - 1)) {
					cerr << "index !!!! " << index << " " << inc;
					cerr << " " << value << " " << smax << " " << smin << std::endl;
				} else {
					histogram[index] += 1.0;
				}
			}
			++iterMask;
			++iterImage;
		}
	} else {
		while ( !iterImage.IsAtEnd() )  {
			PixelType value =  iterImage.Get();
			var += (avg - value) * (avg - value);
			if (value >= smin && value <= smax) {
				svar += (savg - value) * (savg - value);
				index = (int) ((value - smin) / inc+0.5);
				if (index < 0 || index > (samp - 1)) {
					cerr << "index !!!! " << index << " " << inc;
					cerr << " " << value << " " << smax << " " << smin << std::endl;
				} else {
					histogram[index] += 1.0;
				}
			}
			++iterImage;
		}
	}
	var /= (num-1); 
	if (snum) {
		svar /= (snum-1);
	}
    
  // write statistics
	ofstream efile(histofile, ios::out);  
	if (!efile) {
		cerr << "Error: open of file \"" << histofile << "\" failed." << std::endl;
		exit(-1);
	}

	double scaling_factor = (inputImage->GetSpacing())[0] * (inputImage->GetSpacing())[1] * 
				(inputImage->GetSpacing())[2];

	efile.precision(10);
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "volumes = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," 
				<< (double) histogram[index] * scaling_factor ;
	}
	efile << "}};" << std::endl << std::endl;
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "histo = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," << histogram[index] ;
	}
	efile << "}};" << std::endl << std::endl;

	double sumHisto = 0.0;
	for (index = 0; index < samp; index ++) {
		sumHisto = sumHisto + histogram[index];
	}
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "Percentile = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," << (double) histogram[index]/sumHisto ;
	}
	efile << "}};" << std::endl << std::endl;

	double sumPerc = 0.0;
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "CumulativePercentile = {{";
		else {
			efile << "},{";
		}
		sumPerc = sumPerc +(double) histogram[index]/sumHisto;
		efile << index*inc + smin << "," << sumPerc;
	}
	efile << "}};" << std::endl << std::endl;

	double quantileIncrement = (double) 1.0/samp;
	for (double quantile = 0.0; quantile <= 1.0; quantile +=  quantileIncrement) {

		sumPerc = 0.0;
		double quantileVal = 0;
		for (index = 0; index < samp; index ++) {
			sumPerc = sumPerc +(double) histogram[index]/sumHisto;
			if (quantile > sumPerc) {
				quantileVal = index*inc + smin;
			} 
		}
		if (quantile==0)
			efile << "Quantile = {{";
		else {
			efile << "},{";
		}
		efile <<  quantile << "," << quantileVal;
	}
	efile << "}};" << std::endl << std::endl;

	efile << "Fullmin = " << min << ";" << std::endl;
	efile << "Fullmax = " << max << ";" << std::endl;
	efile << "Fullavg = " << avg << ";" << std::endl;
	efile << "Fullvar = " << var << ";" << std::endl;
	efile << "Fullstdev = " << sqrt(var) << ";" << std::endl << std::endl;
	efile << "Selmin = " << smin << ";" << std::endl;
	efile << "Selmax = " << smax << ";" << std::endl;
	efile << "Selavg = " << savg << ";" << std::endl;
	efile << "Selvar = " << svar << ";" << std::endl;
	efile << "Selstdev = " << sqrt(svar) << ";" << std::endl << std::endl;
	efile << "NofFullVoxels = " << num << ";" << std::endl;
	efile << "NofSelVoxels = " << snum << ";"<< std::endl;
	efile << "Nofsamples = " << samp << ";" << std::endl;
	efile << "SumFullVoxels = " << avg * num << ";" << std::endl;
	efile << "SumSelVoxels = " << savg * snum << ";"<< std::endl;
	efile << "VolumeSumFullVoxels = " << avg * num * scaling_factor << ";" << std::endl;
	efile << "VolumeSumSelVoxels = " << savg * snum * scaling_factor << ";"<< std::endl;
	efile << "SumFullVoxelsRatio = " << avg * num / max << ";" << std::endl;
	efile << "SumSelVoxelsRatio = " << savg * snum / smax << ";"<< std::endl;
	efile << "VolumeFullRatio = " << avg * num / max * scaling_factor << ";" << std::endl;
	std::cout<<avg * num / max * scaling_factor<<endl;
	efile << "VolumeSelRatio = " << savg * snum / smax * scaling_factor << ";"<< std::endl;
    
	efile.close();
}

static void histo_probvol(ImagePointer inputImage, ImagePointer maskImage, short probFactor, char * histofile,
			  double min, double max,double  smin,double  smax, int samp)
{
	double avg = 0 ;
	double savg = 0, snum = 0, num = 0;
	double var = 0, svar = 0;
  
	Iterator iterImage (inputImage, inputImage->GetBufferedRegion());
	Iterator iterMask (maskImage, maskImage->GetBufferedRegion());
	//Initialize
	iterImage.GoToBegin();
	iterMask.GoToBegin();
	if (debug) std::cout << "computing histogram" << std::endl;

	while ( !iterImage.IsAtEnd() )  {
		PixelType value =  iterImage.Get();
		double Maskvalue =  (double) iterMask.Get() / probFactor;
		avg += (double) value * Maskvalue;
		num += (double) Maskvalue;
		if (value >= smin && value <= smax) {
			savg += (double)  value * Maskvalue;
			snum += (double) Maskvalue;;
		}
		++iterMask;
		++iterImage;
	}

	avg /= num;
	if (snum) savg /= snum;

  // compute variance and histogram
	iterImage.GoToBegin();
	iterMask.GoToBegin();
	float histogram[40000];
	double inc;
	inc = (smax - smin) / (double) (samp-1);
	int index;

	while ( !iterImage.IsAtEnd() )  {
		PixelType value =  iterImage.Get();
		double Maskvalue =  (double) iterMask.Get()/ probFactor;;
		var += (avg - value) * Maskvalue * (avg - value) * Maskvalue;
		if ( value >= smin && value <= smax) {
			svar += (savg - value) * Maskvalue * (savg - value) * Maskvalue;
			index = (int) ((value - smin) / inc+0.5);
			if (index < 0 || index > (samp - 1)) {
				cerr << "index !!!! " << index << " " << inc;
				cerr << " " << value << " " << smax << " " << smin << std::endl;
			} 
			else {
				histogram[index] += (double) Maskvalue;
			}
		}
		++iterMask;
		++iterImage;
	}
	var /= (num-1); 
	if (snum) {
		svar /= (snum-1);
	}

	// write statistics
	ofstream efile(histofile, ios::out);  
	if (!efile) {
		cerr << "Error: open of file \"" << histofile << "\" failed." << std::endl;
		exit(-1);
	}

	double scaling_factor = (inputImage->GetSpacing())[0] * (inputImage->GetSpacing())[1] * 
				(inputImage->GetSpacing())[2];

	efile.precision(10);
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "volumes = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," 
				<< (double) histogram[index] * scaling_factor ;
	}
	efile << "}};" << std::endl << std::endl;
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "histo = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," << histogram[index] ;
	}
	efile << "}};" << std::endl << std::endl;

	double sumHisto = 0.0;
	for (index = 0; index < samp; index ++) {
		sumHisto = sumHisto + histogram[index];
	}
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "Percentile = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," << (double) histogram[index]/sumHisto ;
	}
	efile << "}};" << std::endl << std::endl;

	double sumPerc = 0.0;
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "CumulativePercentile = {{";
		else {
			efile << "},{";
		}
		sumPerc = sumPerc +(double) histogram[index]/sumHisto;
		efile << index*inc + smin << "," << sumPerc;
	}
	efile << "}};" << std::endl << std::endl;

	double quantileIncrement = (double) 1.0/ (samp - 1);
	for (double quantile = 0.0; quantile <= 1.0; quantile +=  quantileIncrement) {

		sumPerc = 0.0;
		double quantileVal = 0;
		for (index = 0; index < samp; index ++) {
			sumPerc = sumPerc +(double) histogram[index]/sumHisto;
			if (quantile > sumPerc) {
				quantileVal = index*inc + smin;
			} 
		}
		if (quantile==0)
			efile << "Quantile = {{";
		else {
			efile << "},{";
		}
		efile <<  quantile << "," << quantileVal;
	}
	efile << "}};" << std::endl << std::endl;

	efile << "Fullmin = " << min << ";" << std::endl;
	efile << "Fullmax = " << max << ";" << std::endl;
	efile << "Fullavg = " << avg << ";" << std::endl;
	efile << "Fullvar = " << var << ";" << std::endl;
	efile << "Fullstdev = " << sqrt(var) << ";" << std::endl << std::endl;
	efile << "Selmin = " << smin << ";" << std::endl;
	efile << "Selmax = " << smax << ";" << std::endl;
	efile << "Selavg = " << savg << ";" << std::endl;
	efile << "Selvar = " << svar << ";" << std::endl;
	efile << "Selstdev = " << sqrt(svar) << ";" << std::endl << std::endl;
	efile << "NofFullVoxels = " << num << ";" << std::endl;
	efile << "NofSelVoxels = " << snum << ";"<< std::endl;
	efile << "Nofsamples = " << samp << ";" << std::endl;
	efile << "SumFullVoxels = " << avg * num << ";" << std::endl;
	efile << "SumSelVoxels = " << savg * snum << ";"<< std::endl;
	efile << "VolumeSumFullVoxels = " << avg * num * scaling_factor << ";" << std::endl;
	efile << "VolumeSumSelVoxels = " << savg * snum * scaling_factor << ";"<< std::endl;
	efile << "SumFullVoxelsRatio = " << avg * num / max << ";" << std::endl;
	efile << "SumSelVoxelsRatio = " << savg * snum / smax << ";"<< std::endl;
	efile << "VolumeFullRatio = " << avg * num / max * scaling_factor << ";" << std::endl;
	std::cout<<avg * num / max * scaling_factor<<endl;
	efile << "VolumeSelRatio = " << savg * snum / smax * scaling_factor << ";"<< std::endl;
    
	efile.close();
}



