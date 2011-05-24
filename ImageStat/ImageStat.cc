/* 
 * compute image statistics 
 *
 * author:  Martin Styner 
 * 
 * changes:
 *
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

#define DEFAULT_SAMP 2
// number of samples by default

using namespace std;
using namespace itk;

typedef short PixelType;
typedef short ExtractPixelType;
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

//typedef itk::Statistics::ScalarImageToListAdaptor<ImageType>	ImageSampleType;
typedef itk::Statistics::ImageToListSampleAdaptor<ImageType>	ImageSampleType;
typedef ImageSampleType::Iterator				ImageSampleIterator;

typedef itk::Vector<short,2> MeasurementVectorType;
typedef itk::Statistics::ListSample <MeasurementVectorType> SampleType;
typedef itk::Statistics::MembershipSample< SampleType >	MembershipSampleType;

typedef itk::Vector<float, 3 > VectorType;
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
		cout << "ImageStat 1.5 version (10.May.99)" << endl;
		cout << " computes base statistics of an image" << endl;
		cout << "usage: ImageStat infile [-outbase outbase] [-mask maskfile]  [-v]"<< endl;
		cout << "       [-histo [-min xmin] [-max xmax] [-samp s]] " << endl;
		cout << endl;
		cout << "infile         input dataset" << endl;;
		cout << "-outbase outbase     base-outputfilename, if omitted then the same as base of input" << endl; 
		cout << "-mask file     Mask [bytedata, 0/255 value] to mask out image data, where mask is set" << endl;
		cout << "-probmask file Mask [shortdata, 0/62000 value] to mask out probabilitically image data, where mask is > 0" << endl;
		cout << "-probfactor f  Maximum probability value/Probability normalization factor, needs to be set if using -probmask" << endl;
		cout << "               The proper probability of all values in probmask is computed by division with this value" << endl;
		cout << "-verbose       verbose mode " << endl << endl; 
		cout << "-histo         Histogram corrected and uncorrected for pixeldimension" << endl;
		cout << "-min xmin      Minimal x value for histogram/volume [DEFAULT masked_minval]" << endl; 
		cout << "-max xmax      Maximal x value for histogram/volume [DEFAULT masked_maxval]" << endl;
		cout << "-samp s        Number of Samples for histogram/volume [DEFAULT xmax-xmin]" << endl << endl; 
		cout << "-threeSlice[C x,y,z]   create three orthogonal views through the image center or the optionally supplied coordinate" << endl;
		cout << "-info          print image info on standard out" << endl;
		cout << "-label labelfile    give volume, mean intensity, standard deviation, min, max, quantiles for each label" << endl;
		cout << "      -display   display the results on the standard output (if -volumeSummary and/or -intensitySummary used, display theresults of these options)" << endl;
		cout << "      -quantile [list of numbers separate by a coma]   quantile value [DEFAULT: 1 5 33 50 66 95 99]" << endl;
		cout << "      -maskForLabel maskForLabelfile   binary mask" << endl;
		cout << "      -probabilityMap probabilitymap   use probability map to compute volume, mean intensity etc." << endl;
		cout << "      -volumeSummary   give volume for each label" << endl;
		cout << "      -intensitySummary   give mean intensity, standard deviation, min, max, quantiles for each label" << endl;
		cout << "-volOverlap labelFile  compute both Tanimoto (intersection over union) and Dice coefficient of input label image with additional image" << endl;
		cout << endl << endl;
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
	bool volumeSummaryOn  = ipExistsArgument(argv, "-volumeSummary");
	bool intensitySummaryOn  = ipExistsArgument(argv, "-intensitySummary");
	
	char *inputFileName = strdup(argv[1]);
	char *outbase    = ipGetStringArgument(argv, "-outbase", NULL);  
	char *maskfile   = ipGetStringArgument(argv, "-mask", NULL);
	char *probmaskfile   = ipGetStringArgument(argv, "-probmask", NULL);
	char *label = ipGetStringArgument(argv, "-label", NULL);
	
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
	/*char *probabilityMap[1000];
	vector<string> probMap;
	int probabilityMapNb = 0;
	if (probabilityMapOn)
	{
		probabilityMapNb = ipGetStringMultipArgument(argv, "-probabilityMap",probabilityMap,1000);
		for(int i = 0 ; i < probabilityMapNb ; i++)
			probMap.push_back(probabilityMap[i]);
	}*/
	
	char *probabilityMap = ipGetStringArgument(argv, "-probabilityMap", NULL);
	
	short probfactor = ipGetIntArgument(argv, "-probfactor", 0);
	if (probmaskfile && probfactor == 0) {
		cout << " when using -probmask for a probability mask, then -probfactor has to be supplied for normalization of the probability values" << endl;
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
	//vector<ImagePointer> vprobaMaps;
	ImagePointer probabilityMapImage;
	ImageSampleType::Pointer imageSample;
	ImageSampleType::Pointer labelSample;
	
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
		if (debug) cout << "Loading file " << inputFileName << endl;
		VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
		imageReader->SetFileName(inputFileName) ;
		imageReader->Update();
		if (debug) cout << "Loading file done " << endl;
		inputImage = imageReader->GetOutput();
		//read mask
		if (maskfile) {
			if (debug) cout << "Loading mask " << maskfile  << endl;
			BinaryVolumeReaderType::Pointer imageReader = BinaryVolumeReaderType::New();
			imageReader->SetFileName(maskfile) ;
			imageReader->Update();
			maskImage = imageReader->GetOutput();
		}

		//read mask
		if (probmaskfile) {
			if (debug) cout << "Loading probability mask " << maskfile  << endl;
			VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
			imageReader->SetFileName(probmaskfile) ;
			imageReader->Update();
			probmaskImage = imageReader->GetOutput();
		}
		//read label
		if (labelOn) {
			if (debug) cout << "Loading Label File " << label << endl;
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
			if (debug) cout << "Loading Mask File " << maskForLabel << endl;
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
			if (debug) cout << "Loading Probability Map" << probabilityMap << endl;
			
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
			
		/*	for (int probaMapNumber = 0; probaMapNumber < probabilityMapsNb; probaMapNumber++)
			{
				VolumeReaderType::Pointer probaMapReader = VolumeReaderType::New();
				if (debug) cout << "Loading file " << probMaps[probaMapNumber] << endl;
				probaMapReader->SetFileName(probMaps[probaMapNumber].c_str());
				try 
				{
					probaMapReader->Update();
				}
				catch (ExceptionObject & err) 
				{
					cerr<<"ExceptionObject caught!"<<endl;
					cerr<<err<<endl;
					return EXIT_FAILURE;	
				}
				vprobaMaps.push_back(probaMapReader->GetOutput());
			}*/
		}
		if (volOverlapFile) {
			if (debug) cout << "Loading Volume Overlap Label File " << volOverlapFile  << endl;
			VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
			imageReader->SetFileName(volOverlapFile) ;
			imageReader->Update();
			volOverlapImage = imageReader->GetOutput();
		}
		if (histOn) {
			//Get extrema over volume/slice/mask
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
			if (debug) cout << "Min: " << smin << ", Max: " << smax << endl;

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

			ImageSampleType::Pointer labelSample = ImageSampleType::New();
			labelSample->SetImage(labelImage);
			PixelType label;
			ImageSampleIterator iterLabel = labelSample->Begin() ;
			if(maskForLabelOn)
			{
				Iterator itMaskForLabel( maskForLabelImage, maskForLabelImage->GetRequestedRegion());
				while( iterLabel != labelSample->End() )
				{
					label = iterLabel.GetMeasurementVector()[0];
					if ((label != 0) && (!searchList(label,labelList)) && (itMaskForLabel.Get() != 0)) labelList.push_back( label );
					++iterLabel;
					++itMaskForLabel;
				}
			}
			else
			{
				while(iterLabel != labelSample->End())
				{
					label = iterLabel.GetMeasurementVector()[0];
					if ((label != 0) && (!searchList(label,labelList))) labelList.push_back(label);
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
				label = iterLabel.GetMeasurementVector()[0];
				if ((label != 0 ) && ((iter.GetMeasurementVector()[1] != 0) ))
				{
					membershipSample->AddInstance(label,iter.GetInstanceIdentifier());
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
				nbPix =membershipSample->GetClassSample(l)->Size(); 
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
				if (debug) cout << "Computing Label stats" << label  << endl;
				char statfile [5024];
				sprintf(statfile,"%s_stat.txt", base_string);
				ofstream efile(statfile, ios::out);  
				if (!efile) {
					cerr << "Error: open of file \"" << statfile << "\" failed." << endl;
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
					efile << l;
					efile << "\t \t" << tab[0][j];
					efile << "\t \t" << tab[1][j];
					efile << "\t \t" << tab[2][j];
					efile << "\t \t" << tab[3][j];
					efile << "\t \t" << tab[4][j];
					if(displayOn)
					{
						cout << l;
						cout << "\t \t" << tab[0][j];
						cout << "\t \t" << tab[1][j];
						cout << "\t \t" << tab[2][j];
						cout << "\t \t" << tab[3][j];
						cout << "\t \t" << tab[4][j];
					}
					for(int i=0; i<tabSize; i++)
					{
						efile << "\t \t" << quantiles_tab[j][i];
						if(displayOn)
							efile << "\t \t" << quantiles_tab[j][i];	
					}
					efile << endl;
					if(displayOn)
						cout<<endl;
				}
				efile.close();
			}
			if (volumeSummaryOn){
				char statfile [5024];
				sprintf(statfile,"%s_volumeSummary.csv", base_string);
				ofstream efile(statfile, ios::out);  
				if (!efile) {
					cerr << "Error: open of file \"" << statfile << "\" failed." << endl;
					exit(-1);
				}
				efile.precision(6);
				efile << inputFileNameTail << ",VOLUME";
				if(displayOn)
				  cout << inputFileNameTail << ",VOLUME";
				// display results for all labels, even if some of them are empty
				int i = 0;
				for (int Label = 1; Label <= maxLabel; Label++)
				{
				  if (labelList[i] == Label)
				    {
				      efile <<","<< tab[0][i];
				      if(displayOn)
					cout <<","<< tab[0][i];
				      i++;
				    }
				  else
				    {
				      efile <<",0";
				      if(displayOn)
					cout <<",0";
				    }
				}
				efile << endl;
				if(displayOn)
				  cout<<endl;
				efile.close();
			}
			if(intensitySummaryOn){
				char statfile [5024];
				sprintf(statfile,"%s_intensitySummary.csv", base_string);
				ofstream efile(statfile, ios::out);  
				if (!efile) {
					cerr << "Error: open of file \"" << statfile << "\" failed." << endl;
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
				      cout << inputFileNameTail;
				      cout<<","<<rowname [i];
				    }
				  // display results for all labels, even if some of them are empty
				  int j = 0;
				  for (int Label = 1; Label <= maxLabel; Label++ )
				    {
				      if (labelList[j] == Label)
					{
					  efile<<","<<tab[i+1][j];
					  if(displayOn)
					    cout<<","<<tab[i+1][j];
					  j++;
					}
				      else
					{
					  efile<<",0";
					  if(displayOn)
					    cout<<",0";
					}
				    }
				  efile<<endl;
				  if(displayOn)
				    cout<<endl;
				}
				for(int k=0; k<tabSize; k++)
				{
				  efile << inputFileNameTail;
				  efile<<","<<rowname [k+4];
				  if(displayOn)
				    {
				      cout << inputFileNameTail;
				      cout<<","<<rowname [k+4];
				    }
				  // display results for all labels, even if some of them are empty
				  int j = 0;
				  for (int Label = 1; Label <= maxLabel; Label++ )
				    {
				      if (labelList[j] == Label)
					{
					  efile<<","<<quantiles_tab[j][k];
					  if(displayOn)
					    cout<<","<<quantiles_tab[j][k];
					  j++;
					}
				      else
					{
					  efile<<",0";
					  if(displayOn)
					    cout<<",0";
					}
				    }
				  efile<<endl;
				  if(displayOn)
				    cout<<endl;
				}
				efile.close();
			}
		delete[] quant_tab;
	}
		if (volOverlapOn){
			if (debug) cout << "Computing Overlap stats" << label  << endl;
			char statfile [5024];
			sprintf(statfile,"%s_volOverlap.txt", base_string);
			ofstream efile(statfile, ios::out);
			if (!efile) {
				cerr << "Error: open of file \"" << statfile << "\" failed." << endl;
				exit(-1);
			}
			efile.precision(6);
			efile<<"##########################################################################"<<endl;
			efile<<"# file: "<<statfile << endl;
			efile<<"# contains volume overlap between " << inputFileName << endl;
			efile<<"                              and " << volOverlapFile<< endl;
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
			volUnion = volUnion + volIntersection;
			efile << "Union = " << volUnion << endl;
			efile << "Intersection = " << volIntersection << endl;
			efile << "volA = " << volA << endl;
			efile << "volB = " << volB << endl;
			efile<< "DiceCoef = " << 2.0 * (double) volIntersection / (volA + volB) << endl;
			efile<< "TanimotoCoef = " << (double) volIntersection / volUnion << endl;

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
  
	if (debug) cout << "computing histogram" << endl;

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
					cerr << " " << value << " " << smax << " " << smin << endl;
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
					cerr << " " << value << " " << smax << " " << smin << endl;
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
		cerr << "Error: open of file \"" << histofile << "\" failed." << endl;
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
	efile << "}};" << endl << endl;
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "histo = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," << histogram[index] ;
	}
	efile << "}};" << endl << endl;

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
	efile << "}};" << endl << endl;

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
	efile << "}};" << endl << endl;

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
	efile << "}};" << endl << endl;

	efile << "Fullmin = " << min << ";" << endl;
	efile << "Fullmax = " << max << ";" << endl;
	efile << "Fullavg = " << avg << ";" << endl;
	efile << "Fullvar = " << var << ";" << endl;
	efile << "Fullstdev = " << sqrt(var) << ";" << endl << endl;
	efile << "Selmin = " << smin << ";" << endl;
	efile << "Selmax = " << smax << ";" << endl;
	efile << "Selavg = " << savg << ";" << endl;
	efile << "Selvar = " << svar << ";" << endl;
	efile << "Selstdev = " << sqrt(svar) << ";" << endl << endl;
	efile << "NofFullVoxels = " << num << ";" << endl;
	efile << "NofSelVoxels = " << snum << ";"<< endl;
	efile << "Nofsamples = " << samp << ";" << endl;
	efile << "SumFullVoxels = " << avg * num << ";" << endl;
	efile << "SumSelVoxels = " << savg * snum << ";"<< endl;
	efile << "VolumeSumFullVoxels = " << avg * num * scaling_factor << ";" << endl;
	efile << "VolumeSumSelVoxels = " << savg * snum * scaling_factor << ";"<< endl;
	efile << "SumFullVoxelsRatio = " << avg * num / max << ";" << endl;
	efile << "SumSelVoxelsRatio = " << savg * snum / smax << ";"<< endl;
	efile << "VolumeFullRatio = " << avg * num / max * scaling_factor << ";" << endl;
	cout<<avg * num / max * scaling_factor<<endl;
	efile << "VolumeSelRatio = " << savg * snum / smax * scaling_factor << ";"<< endl;
    
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
	if (debug) cout << "computing histogram" << endl;

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
				cerr << " " << value << " " << smax << " " << smin << endl;
			} else {
			  
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
		cerr << "Error: open of file \"" << histofile << "\" failed." << endl;
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
	efile << "}};" << endl << endl;
	for (index = 0; index < samp; index ++) {
		if (index==0)
			efile << "histo = {{";
		else {
			efile << "},{";
		}
		efile << index*inc + smin << "," << histogram[index] ;
	}
	efile << "}};" << endl << endl;

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
	efile << "}};" << endl << endl;

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
	efile << "}};" << endl << endl;

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
	efile << "}};" << endl << endl;

	efile << "Fullmin = " << min << ";" << endl;
	efile << "Fullmax = " << max << ";" << endl;
	efile << "Fullavg = " << avg << ";" << endl;
	efile << "Fullvar = " << var << ";" << endl;
	efile << "Fullstdev = " << sqrt(var) << ";" << endl << endl;
	efile << "Selmin = " << smin << ";" << endl;
	efile << "Selmax = " << smax << ";" << endl;
	efile << "Selavg = " << savg << ";" << endl;
	efile << "Selvar = " << svar << ";" << endl;
	efile << "Selstdev = " << sqrt(svar) << ";" << endl << endl;
	efile << "NofFullVoxels = " << num << ";" << endl;
	efile << "NofSelVoxels = " << snum << ";"<< endl;
	efile << "Nofsamples = " << samp << ";" << endl;
	efile << "SumFullVoxels = " << avg * num << ";" << endl;
	efile << "SumSelVoxels = " << savg * snum << ";"<< endl;
	efile << "VolumeSumFullVoxels = " << avg * num * scaling_factor << ";" << endl;
	efile << "VolumeSumSelVoxels = " << savg * snum * scaling_factor << ";"<< endl;
	efile << "SumFullVoxelsRatio = " << avg * num / max << ";" << endl;
	efile << "SumSelVoxelsRatio = " << savg * snum / smax << ";"<< endl;
	efile << "VolumeFullRatio = " << avg * num / max * scaling_factor << ";" << endl;
	cout<<avg * num / max * scaling_factor<<endl;
	efile << "VolumeSelRatio = " << savg * snum / smax * scaling_factor << ";"<< endl;
    
	efile.close();
}



