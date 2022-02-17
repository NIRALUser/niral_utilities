//#if defined(_MSC_VER)
//#pragma warning ( disable : 4786 )
//#endif

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkCastImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>

#include <itkImageToImageFilter.h>
#include <itkDiffusionTensor3D.h>
#include <itkImageDuplicator.h>

#include <itkExtractImageFilter.h>
#include <itkConstantPadImageFilter.h>

#include "argio.h"
#include <itkConnectedComponentImageFilter.h>
#include <itkDiffusionTensor3DReconstructionImageFilter.h>

#include <itkHistogramMatchingImageFilter.h>

using namespace itk;
using namespace std;

int help_function()
{
  std::cout<<"Usage: " << std::endl;
  std::cout << "DTITensorNormalization <Input_DTI_file> -atlasdti <Atlas_DTI_file> [options]\n" << std::endl;
  std::cout << "Input_DTI_file, Atlas_DTI_file : .nhdr or .nrrd file\n" << std::endl;
  std::cout << "Output_DTI_file: (Input_DTI_file)_Norm_DTI.nhdr file\n" << std::endl;
  
  std::cout << "With the options: " << std::endl;
  std::cout << " -matchHistoPara bins,points,thresh  optional parameters for matchHistogram \n (bins= number of histogram bins; points = number of control points; thresh = boolean for mean intensity threshold [0/1]). \n Default: 256,50,1 \n" << std::endl;
  std::cout << " --help: \tdisplay this message \n" << std::endl;
  std::cout << "This tool performs FA and cdf based tensor normalization of the input DTI image to the atlas DTI image. \n" << std::endl;
  std::cout << "Authors: Aditya Gupta, Martin Styner.\n" << std::endl;
  return 0;
//The normalization is based on mapping the individual tensor eigen values of the input DTI to particular cumulative distribution function plane in the atlas space and matching the atlas FA values on each plane.

}

int main(int argc, const char* argv[])
{
	if (argc<2)
    		{
      		cout<<"Please provide the input DTI file name (.nhdr or .nrrd format) and the atlas file name (.nhdr or .nrrd format) OR use --help for options\nExit without generating results!\n";
      		return -1;
    		}

  	if( (strcmp(argv[1],"--help")==0) || (strcmp(argv[1],"--HELP")==0) || (strcmp(argv[1],"--Help")==0) || (argc == 1) )
    		{
      		help_function();
      		return -1;
    		}

	string fileName(argv[1]);
	string atlasDTIFileName(ipGetStringArgument(argv, "-atlasdti", ""));
	//cout << "test1010"<< endl;

	char *inputFileName = strdup(argv[1]);
	char *base_string, *base_string_atlas;
	base_string = strdup(ipGetBaseName(inputFileName));
	base_string_atlas = strdup(ipGetBaseName(atlasDTIFileName.c_str()));
	string inputFAFileName ("temp");
	string inputnormFAFileName ("temp");
	string atlasFAFileName ("temp");
	string outFileName ("temp");
	string outDTIFileName ("temp");
	
	char * formatChar = ipGetStringArgument(argv, "-extension", ".nhdr");
  	string format;
  	if (! strchr(formatChar, '.')) 
		format = string(".") + string(formatChar);
	else 
		format = string(formatChar);

	// Definitions
	typedef short int ShortPixelType;
	typedef double DoublePixelType;
	typedef float PixelType;

	const int ImageDimension = 3;
	typedef Image<DoublePixelType,ImageDimension>      ImageType;
	typedef Image<ShortPixelType, ImageDimension> 	   ImageType_short;
	typedef Image<DoublePixelType, ImageDimension> 	   ImageType_d;
	typedef ImageType::Pointer                    	   ImagePointer;
	typedef ImageSliceIteratorWithIndex< ImageType >      	   IteratorType;
	typedef ImageSliceIteratorWithIndex< ImageType_short >     IteratorType_short;
	typedef itk::Point< short int, ImageDimension >       PointType;

	//Histogram definitions
	typedef HistogramMatchingImageFilter< ImageType , ImageType > matchHistogramFilterType;
	int matchHistoNumBins, matchHistoNumPoints, matchHistoThresh;
	bool Histoparameters = ipExistsArgument(argv, "-matchHistoPara");
	if (Histoparameters)
		{
		matchHistoNumBins = ipGetIntArgument(argv,"-matchHistoPara",1);
		matchHistoNumPoints = ipGetIntArgument(argv,"-matchHistoPara",2);
		matchHistoThresh = ipGetIntArgument(argv,"-matchHistoPara",3);
		}
	else
		{
		matchHistoNumBins = 256;
        	matchHistoNumPoints = 50;
        	matchHistoThresh = 1;
		}

	typedef ConnectedComponentImageFilter< ImageType, ImageType_short > ConnectedComponentImageFilterType;

	//DTI definitions
	typedef itk::ImageBase< 3 > 			ImageBaseType ;
	ImageBaseType::Pointer inputBaseImage ;

	typedef DiffusionTensor3D< DoublePixelType > 	InputTensorDataType;

	typedef itk::Image<InputTensorDataType, 3>	DiffusionImageType;
	DiffusionImageType::Pointer inputDTImage;

	typedef ImageFileReader< ImageType >          	VolumeReaderType;
	typedef ImageFileWriter< ImageType >          	VolumeWriterType;

	typedef ImageFileReader< ImageType_short >          	VolumeReaderType_short;
	typedef ImageFileWriter< ImageType_short >          	VolumeWriterType_short;

	typedef itk::ImageFileReader< DiffusionImageType > DiffusionReaderType ;
	typedef itk::ImageFileWriter< DiffusionImageType > DiffusionWriterType ;
  	
	DiffusionTensor3D< double >::EigenValuesArrayType eigenValues ;
    	DiffusionTensor3D< double >::EigenVectorsMatrixType eigenVectors ;
	DiffusionTensor3D< double >::EigenValuesArrayType aeigenValues ;
    	DiffusionTensor3D< double >::EigenVectorsMatrixType aeigenVectors ;
	DiffusionTensor3D< double >::EigenValuesArrayType feigenValues ;
    	DiffusionTensor3D< double >::EigenVectorsMatrixType feigenVectors ;
	DiffusionTensor3D< double >::EigenValuesArrayType testeigenValues ;
    	DiffusionTensor3D< double >::EigenVectorsMatrixType testeigenVectors ;
    	
	//Read input DTI image
	DiffusionReaderType::Pointer diffusionReader1 = DiffusionReaderType::New() ;
	DiffusionReaderType::Pointer diffusionReader2 = DiffusionReaderType::New() ;
	try
	{
		diffusionReader1->SetFileName(fileName.c_str()) ;
		diffusionReader1->Update() ;
	}
	catch (ExceptionObject & err)
	{
		cerr << " Error reading Input DTI File... ExceptionObject caught!" << err << endl;
	  	return EXIT_FAILURE;
	}
	inputDTImage = diffusionReader1->GetOutput() ;	

	ImageType::SizeType DTIImage_size;
	DTIImage_size = inputDTImage->GetLargestPossibleRegion().GetSize();
	int DimensionX, DimensionY, DimensionZ;
	DimensionX = DTIImage_size[0];
	DimensionY = DTIImage_size[1];
	DimensionZ = DTIImage_size[2];
	cout << "Image Dimension X : " << DimensionX << "  Dimension Y : " << DimensionY << "  Dimension Z : " << DimensionZ << endl;

	//Read atlas DTI image
	DiffusionImageType::Pointer atlasDTImage;
	try
	{
		diffusionReader2->SetFileName(atlasDTIFileName.c_str()) ;
		diffusionReader2->Update() ;
	}
	catch (ExceptionObject & err)
	{
		cerr << " Error reading Atlas Image... ExceptionObject caught!" << err << endl;
	  	return EXIT_FAILURE;
	}
	atlasDTImage = diffusionReader2->GetOutput() ;

	ImageType::SizeType atlasDTImage_size;
	atlasDTImage_size = atlasDTImage->GetLargestPossibleRegion().GetSize();
	int AtlasDimensionX, AtlasDimensionY, AtlasDimensionZ;
	AtlasDimensionX = atlasDTImage_size[0];
	AtlasDimensionY = atlasDTImage_size[1];
	AtlasDimensionZ = atlasDTImage_size[2];
	cout << "Atlas Dimension X : " << AtlasDimensionX << "  Dimension Y : " << AtlasDimensionY << "  Dimension Z : " << AtlasDimensionZ << endl;

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////GENERATING FA IMAGES FOR INPUT AND ATLAS AND THE INTENSITY NORM FA/////////
////////////////////////////////////////////////////////////////////////////////////////////
	//generate FA image for INPUT DTI
	inputFAFileName.erase();
    	inputFAFileName.append(base_string);
    	inputFAFileName.append("_FA");
	inputFAFileName.append(format);

	string str;
	str += "dtiprocess --dti_image  ";
	str += inputFileName;
  	//str += " --RD_output ";
	str += " -f ";
  	str += inputFAFileName;
  	
//	std::cout << str << std::endl;

	system(str.c_str());

	//generate FA image for ATLAS DTI
	atlasFAFileName.erase();
    	atlasFAFileName.append(base_string_atlas);
    	atlasFAFileName.append("_FA");
	atlasFAFileName.append(format);

	string str1;
	str1 += "dtiprocess --dti_image  ";
	str1 += atlasDTIFileName;
  //	str1 += " --RD_output ";
	str1 += " -f ";
  	str1 += atlasFAFileName;
  	
//	std::cout << str1 << std::endl;
	system(str1.c_str());

	//read input FA and atlas FA files
	ImagePointer FAorigImage;
	VolumeReaderType::Pointer FAorigImageReader = VolumeReaderType::New();
	FAorigImageReader->SetFileName(inputFAFileName.c_str());
	try
	{
		FAorigImageReader->Update();
	}
	catch (ExceptionObject & err)
	{
	  	cerr << " Error reading image FA... ExceptionObject caught!" << err << endl;
	  	return EXIT_FAILURE;
	}
	FAorigImage = FAorigImageReader->GetOutput();

	ImagePointer FAatlasImage;
	VolumeReaderType::Pointer FAatlasImageReader = VolumeReaderType::New();
	FAatlasImageReader->SetFileName(atlasFAFileName.c_str());
	try
	{
		FAatlasImageReader->Update();
	}
	catch (ExceptionObject & err)
	{
	  	cerr << " Error reading atlas... ExceptionObject caught!" << err << endl;
	  	return EXIT_FAILURE;
	}
	FAatlasImage = FAatlasImageReader->GetOutput();

	IteratorType iterFAorigImage (FAorigImage, FAorigImage->GetLargestPossibleRegion());
	IteratorType iterFAatlasImage (FAatlasImage, FAatlasImage->GetLargestPossibleRegion());
	
	//Generate intensity normalized FA image
	typedef itk::ImageDuplicator< ImageType > DuplicatorType;
     	DuplicatorType::Pointer duplicator_FAnorm = DuplicatorType::New();
     	duplicator_FAnorm->SetInputImage(FAorigImage);
     	duplicator_FAnorm->Update();

	ImageType::Pointer ImageNormFA = duplicator_FAnorm->GetOutput();
	
	matchHistogramFilterType::Pointer matchHistogramFilter = matchHistogramFilterType::New();
    	matchHistogramFilter->SetSourceImage(FAorigImage);
    	matchHistogramFilter->SetReferenceImage(FAatlasImage);
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
    	ImageNormFA = matchHistogramFilter->GetOutput();

	IteratorType iterImageNormFA (ImageNormFA, ImageNormFA->GetLargestPossibleRegion());
		
	//Write the intensity normalized FA image
	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_normFA");
	outFileName.append(format);
	VolumeWriterType::Pointer writer_normFA = VolumeWriterType::New();
	writer_normFA->SetFileName(outFileName.c_str()); 
    	writer_normFA->SetInput(ImageNormFA);
    	writer_normFA->Write();

///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////CREATING LAMBDA IMAGES////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	//Create Images
	//copy and create image using itk Duplicator
	DuplicatorType::Pointer duplicator_lambda1 = DuplicatorType::New();
     	duplicator_lambda1->SetInputImage(FAorigImage);
     	duplicator_lambda1->Update();
	DuplicatorType::Pointer duplicator_lambda2 = DuplicatorType::New();
     	duplicator_lambda2->SetInputImage(FAorigImage);
     	duplicator_lambda2->Update();
	DuplicatorType::Pointer duplicator_lambda3 = DuplicatorType::New();
     	duplicator_lambda3->SetInputImage(FAorigImage);
     	duplicator_lambda3->Update();

     	ImageType::Pointer imagelambda1 = duplicator_lambda1->GetOutput();
	ImageType::Pointer imagelambda2 = duplicator_lambda2->GetOutput();
	ImageType::Pointer imagelambda3 = duplicator_lambda3->GetOutput();
	IteratorType iterImagelambda1 (imagelambda1, imagelambda1->GetLargestPossibleRegion());	
	IteratorType iterImagelambda2 (imagelambda2, imagelambda2->GetLargestPossibleRegion());
	IteratorType iterImagelambda3 (imagelambda3, imagelambda3->GetLargestPossibleRegion());	

    	//computing eigen values
	typedef itk::ImageSliceIteratorWithIndex<DiffusionImageType> IterDTI;
	IterDTI iterDTImage1 (inputDTImage, inputDTImage->GetLargestPossibleRegion());
	InputTensorDataType inputTensor1;
	
	float max0 = 0.0;
	float max1 = 0.0;
	float max2 = 0.0;
	float min0 = 0.0;
	float min1 = 0.0;
	float min2 = 0.0;
	float scaling_fact = 100000;
	//obtaining the eigen values and assigning to different images and then scaling the eigen values
	int count = 0;

	iterDTImage1.GoToBegin();
	iterImagelambda1.GoToBegin();
	iterImagelambda2.GoToBegin();
	iterImagelambda3.GoToBegin();

	iterDTImage1.SetFirstDirection(0);
	iterImagelambda1.SetFirstDirection(0);
	iterImagelambda2.SetFirstDirection(0);
	iterImagelambda3.SetFirstDirection(0);
	
	iterDTImage1.SetSecondDirection(1);
	iterImagelambda1.SetSecondDirection(1);
	iterImagelambda2.SetSecondDirection(1);
	iterImagelambda3.SetSecondDirection(1);

	while (!iterDTImage1.IsAtEnd())
		{
		while (!iterDTImage1.IsAtEndOfSlice())
			{
			while (!iterDTImage1.IsAtEndOfLine())
				{
				inputTensor1 = iterDTImage1.Get();
				inputTensor1.ComputeEigenAnalysis(eigenValues, eigenVectors);
				if (count>1)
				{
					if (eigenValues[0] > max0)
						max0 = eigenValues[0];
					if (eigenValues[1] > max1)
						max1 = eigenValues[1];
					if (eigenValues[2] > max2)
						max2 = eigenValues[2];
					if (min0 > eigenValues[0])
						min0 = eigenValues[0];
					if (min1 > eigenValues[1])
						min1 = eigenValues[1];
					if (min2 > eigenValues[2])
						min2 = eigenValues[2];	
				}
				else
				{
					max0 = eigenValues[0];
					max1 = eigenValues[1];
					max2 = eigenValues[2];
					min0 = eigenValues[0];
					min1 = eigenValues[1];
					min2 = eigenValues[2];
				}
				iterImagelambda3.Set(eigenValues[0]*scaling_fact);
				iterImagelambda2.Set(eigenValues[1]*scaling_fact);
				iterImagelambda1.Set(eigenValues[2]*scaling_fact);
				++iterDTImage1;
				++iterImagelambda1;
				++iterImagelambda2;
				++iterImagelambda3;
				count++;
				}
			iterDTImage1.NextLine();
			iterImagelambda1.NextLine();
			iterImagelambda2.NextLine();
			iterImagelambda3.NextLine();
			}	
		iterDTImage1.NextSlice();
		iterImagelambda1.NextSlice();
		iterImagelambda2.NextSlice();
		iterImagelambda3.NextSlice();
		}
	
/*	cout << " Max eigen value 0 is : " << max0*scaling_fact << endl;
	cout << " Max eigen value 1 is : " << max1*scaling_fact << endl;
	cout << " Max eigen value 2 is : " << max2*scaling_fact << endl;
	cout << " For loop count is : " << count << endl;*/

///////////////////////////////////////detecting background////////////////////////////////////// 


	ConnectedComponentImageFilterType::Pointer connectedcomponentFilter1 = ConnectedComponentImageFilterType::New(); 
	ConnectedComponentImageFilterType::Pointer connectedcomponentFilter2 = ConnectedComponentImageFilterType::New();	

	ImageType_short::Pointer foregroundImage = ImageType_short::New();
	ImageType_short::SizeType ImDim;
	ImDim[0] = DimensionX;
	ImDim[1] = DimensionY;
	ImDim[2] = DimensionZ;
//	cout << "ImageDimensions 0: " << ImDim[0] << "  1: " << ImDim[1] << "  2  " << ImDim[2] << endl; 
	ImageType_short::SpacingType spacinginfo;
	spacinginfo[0] = FAorigImage->GetSpacing()[0];
	spacinginfo[1] = FAorigImage->GetSpacing()[0];
	spacinginfo[2] = FAorigImage->GetSpacing()[0];
//	cout << " ImageSpacing 0: " << spacinginfo[0] << "  1:  " << spacinginfo[1] << "  2  " << spacinginfo[2] << endl;
	foregroundImage->SetSpacing(spacinginfo);
	foregroundImage->SetRegions(ImDim);
	foregroundImage->Allocate();

/*	foregroundImagetemp->SetSpacing(spacinginfo);
	foregroundImagetemp->SetRegions(ImDim);
	foregroundImagetemp->Allocate();*/

	connectedcomponentFilter1->SetInput(FAorigImage);
	connectedcomponentFilter1->Update();
	foregroundImage = connectedcomponentFilter1->GetOutput();

	//Write the foreground image
	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_foreground");
	outFileName.append(format);
	VolumeWriterType_short::Pointer writer_foregroundFA = VolumeWriterType_short::New();
	writer_foregroundFA->SetFileName(outFileName.c_str()); 
    	writer_foregroundFA->SetInput(foregroundImage);
    	writer_foregroundFA->Write();

	IteratorType_short iterforegroundImage (foregroundImage, foregroundImage->GetLargestPossibleRegion());
	
	vector< vector < vector<int> > > foreground_matrix (DimensionX, vector< vector<int> > (DimensionY, vector<int>(DimensionZ,0)));
	int ImageMaskcount = 0;
	iterforegroundImage.GoToBegin();
	iterforegroundImage.SetFirstDirection(0);
	iterforegroundImage.SetSecondDirection(1);

	while (!iterforegroundImage.IsAtEnd())
		{
		while (!iterforegroundImage.IsAtEndOfSlice())
			{
			while (!iterforegroundImage.IsAtEndOfLine())
				{
				int temp = iterforegroundImage.Get();
				//foreground_matrix[ifore][jfore][kfore] = temp;
				if (temp == 1)
					{ImageMaskcount++;}
				++iterforegroundImage;
				}
			iterforegroundImage.NextLine();
			}	
		iterforegroundImage.NextSlice();
		}
	
	
	//cout << " ImageMaskCount is : " << ImageMaskcount << endl;
	
////////////////////ATLAS FOREGROUND///////////////////
	ImageType_short::Pointer foregroundAtlas = ImageType_short::New();
	ImageType_short::SizeType ImDimA;
	ImDimA[0] = AtlasDimensionX;
	ImDimA[1] = AtlasDimensionY;
	ImDimA[2] = AtlasDimensionZ;
	//cout << "AtlasDimensions 0: " << ImDimA[0] << "  1: " << ImDimA[1] << "  2  " << ImDimA[2] << endl; 
	ImageType_short::SpacingType spacinginfoA;
	spacinginfoA[0] = FAatlasImage->GetSpacing()[0];
	spacinginfoA[1] = FAatlasImage->GetSpacing()[0];
	spacinginfoA[2] = FAatlasImage->GetSpacing()[0];
	//cout << " AtlasSpacing 0: " << spacinginfoA[0] << "  1:  " << spacinginfoA[1] << "  2  " << spacinginfoA[2] << endl;
	foregroundAtlas->SetSpacing(spacinginfoA);
	foregroundAtlas->SetRegions(ImDimA);
	foregroundAtlas->Allocate();

	connectedcomponentFilter2->SetInput(FAatlasImage);
	connectedcomponentFilter2->Update();
	foregroundAtlas = connectedcomponentFilter2->GetOutput();

	//Write the foreground image
	outFileName.erase();
    	outFileName.append(base_string_atlas);
    	outFileName.append("_foreground");
	outFileName.append(format);
	VolumeWriterType_short::Pointer writer_foregroundAtlasFA = VolumeWriterType_short::New();
	writer_foregroundAtlasFA->SetFileName(outFileName.c_str()); 
    	writer_foregroundAtlasFA->SetInput(foregroundAtlas);
    	writer_foregroundAtlasFA->Write();

	IteratorType_short iterforegroundAtlas (foregroundAtlas, foregroundAtlas->GetLargestPossibleRegion());
	
	iterforegroundAtlas.GoToBegin();
	iterforegroundAtlas.SetFirstDirection(0);
	iterforegroundAtlas.SetSecondDirection(1);

	/////////////////////////assigning lambda values to matrices///////////////////////////////
	//Assigning Lambda Values to matrices
	vector< vector < vector<float> > > lambda0 (DimensionX, vector< vector<float> > (DimensionY, vector<float>(DimensionZ,0)));
	vector< vector < vector<float> > > lambda1 (DimensionX, vector< vector<float> > (DimensionY, vector<float>(DimensionZ,0)));
	vector< vector < vector<float> > > lambda2 (DimensionX, vector< vector<float> > (DimensionY, vector<float>(DimensionZ,0)));
	int casecount = 0;
	int atlascount = 0;
	iterImagelambda1.GoToBegin();
	iterImagelambda2.GoToBegin();
	iterImagelambda3.GoToBegin();
	iterforegroundImage.GoToBegin();

	iterforegroundImage.SetFirstDirection(0);
	iterforegroundImage.SetSecondDirection(1);

	for ( int i = 0 ; i< DimensionX ; i++)
		{
		for (int j = 0 ; j < DimensionY ; j++)
			{
			for (int k = 0 ; k < DimensionZ ; k++)
				{
				int background_value = iterforegroundImage.Get();
				//cout << "  " << background_value;
				if (background_value != 0)
					{
					lambda0[i][j][k] = iterImagelambda3.Get();//*scaling_fact;
					lambda1[i][j][k] = iterImagelambda2.Get();//*scaling_fact;
					lambda2[i][j][k] = iterImagelambda1.Get();//*scaling_fact;
					casecount++;
					}
				else
					{
					lambda0[i][j][k] = 0;
					lambda1[i][j][k] = 0;
					lambda2[i][j][k] = 0;
					}
					++iterImagelambda1;
					++iterImagelambda2;
					++iterImagelambda3;
					++iterforegroundImage;
					
				}
			}
		}
			
	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_orig_lambda1");
	outFileName.append(format);
	//cout << "Output file Original Lambda1 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writer1 = VolumeWriterType::New();
	writer1->SetFileName(outFileName.c_str()); 
    	writer1->SetInput(imagelambda1);
    	writer1->Write();

	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_orig_lambda2");
	outFileName.append(format);
	//cout << "Output file Original Lambda2 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writer2 = VolumeWriterType::New();
	writer2->SetFileName(outFileName.c_str()); 
    	writer2->SetInput(imagelambda2);
    	writer2->Write();

	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_orig_lambda3");
	outFileName.append(format);
	//cout << "Output file Original Lambda3 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writer3 = VolumeWriterType::New();
	writer3->SetFileName(outFileName.c_str()); 
    	writer3->SetInput(imagelambda3);
    	writer3->Write();

	////////////////////////////////////////////////ATLAS//////////////////////////////////////////
	//Create Images
	//copy and create image using itk Duplicator
	DuplicatorType::Pointer duplicator_atlas1 = DuplicatorType::New();
     	duplicator_atlas1->SetInputImage(FAatlasImage);
     	duplicator_atlas1->Update();
	DuplicatorType::Pointer duplicator_atlas2 = DuplicatorType::New();
     	duplicator_atlas2->SetInputImage(FAatlasImage);
     	duplicator_atlas2->Update();
	DuplicatorType::Pointer duplicator_atlas3 = DuplicatorType::New();
     	duplicator_atlas3->SetInputImage(FAatlasImage);
     	duplicator_atlas3->Update();
 
    	ImageType::Pointer atlaslambda1 = duplicator_atlas1->GetOutput();
	ImageType::Pointer atlaslambda2 = duplicator_atlas2->GetOutput();
	ImageType::Pointer atlaslambda3 = duplicator_atlas3->GetOutput();
	
	//computing eigen images for atlas
	IterDTI iteratlasDTImage1 (atlasDTImage, atlasDTImage->GetLargestPossibleRegion());
	InputTensorDataType inputTensor2;
	IteratorType iteratlaslambda1 (atlaslambda1, atlaslambda1->GetLargestPossibleRegion());	
	IteratorType iteratlaslambda2 (atlaslambda2, atlaslambda2->GetLargestPossibleRegion());
	IteratorType iteratlaslambda3 (atlaslambda3, atlaslambda3->GetLargestPossibleRegion());

	count = 0;
	float amax0 = 0.0;
	float amax1 = 0.0;
	float amax2 = 0.0;
	float amin0 = 0.0;
	float amin1 = 0.0;
	float amin2 = 0.0;
	
	iteratlasDTImage1.GoToBegin();
	iteratlaslambda1.GoToBegin();
	iteratlaslambda2.GoToBegin();
	iteratlaslambda3.GoToBegin();
	
	iteratlasDTImage1.SetFirstDirection(0);
	iteratlaslambda1.SetFirstDirection(0);
	iteratlaslambda2.SetFirstDirection(0);
	iteratlaslambda3.SetFirstDirection(0);
	
	iteratlasDTImage1.SetSecondDirection(1);
	iteratlaslambda1.SetSecondDirection(1);
	iteratlaslambda2.SetSecondDirection(1);
	iteratlaslambda3.SetSecondDirection(1);

	while (!iteratlasDTImage1.IsAtEnd())
		{
		while (!iteratlasDTImage1.IsAtEndOfSlice())
			{
			while (!iteratlasDTImage1.IsAtEndOfLine())
				{
				count++;
				inputTensor2 = iteratlasDTImage1.Get();
				inputTensor2.ComputeEigenAnalysis(aeigenValues, aeigenVectors);
				if (count>1)
					{
					if (aeigenValues[0] > amax0)
						amax0 = aeigenValues[0];
					if (aeigenValues[1] > amax1)
						amax1 = aeigenValues[1];
					if (aeigenValues[2] > amax2)
						amax2 = aeigenValues[2];
					if (amin0 > aeigenValues[0])
						amin0 = aeigenValues[0];
					if (amin1 > aeigenValues[1])
						amin1 = aeigenValues[1];
					if (amin2 > aeigenValues[2])
						amin2 = aeigenValues[2];
					}
				else
					{
					amax0 = aeigenValues[0];
					amax1 = aeigenValues[1];
					amax2 = aeigenValues[2];
					amin0 = aeigenValues[0];
					amin1 = aeigenValues[1];
					amin2 = aeigenValues[2];	
					//cout << "Assigning atlas first max values " << endl; 
					}
				
				iteratlaslambda3.Set(aeigenValues[0]*scaling_fact);
				iteratlaslambda2.Set(aeigenValues[1]*scaling_fact);
				iteratlaslambda1.Set(aeigenValues[2]*scaling_fact);
    				++iteratlasDTImage1;
				++iteratlaslambda1;
				++iteratlaslambda2;
				++iteratlaslambda3;
				}
			iteratlasDTImage1.NextLine();
			iteratlaslambda1.NextLine();
			iteratlaslambda2.NextLine();
			iteratlaslambda3.NextLine();
			}	
		iteratlasDTImage1.NextSlice();
		iteratlaslambda1.NextSlice();
		iteratlaslambda2.NextSlice();
		iteratlaslambda3.NextSlice();
		}


	/*cout << " Max Atlas eigen value 0 is : " << amax0*scaling_fact << endl;
	cout << " Max Atlas eigen value 1 is : " << amax1*scaling_fact << endl;
	cout << " Max Atlas eigen value 2 is : " << amax2*scaling_fact << endl;
	cout << " For loop count is : " << count << endl;*/

	vector< vector < vector<float> > > alambda0 (AtlasDimensionX, vector< vector<float> > (AtlasDimensionY, vector<float>(AtlasDimensionZ,0)));
	vector< vector < vector<float> > > alambda1 (AtlasDimensionX, vector< vector<float> > (AtlasDimensionY, vector<float>(AtlasDimensionZ,0)));
	vector< vector < vector<float> > > alambda2 (AtlasDimensionX, vector< vector<float> > (AtlasDimensionY, vector<float>(AtlasDimensionZ,0)));
	
	iteratlaslambda1.GoToBegin();
	iteratlaslambda2.GoToBegin();
	iteratlaslambda3.GoToBegin();
	iterforegroundAtlas.GoToBegin();

	iterforegroundAtlas.SetFirstDirection(0);
	iterforegroundAtlas.SetSecondDirection(1);

	for ( int i = 0 ; i < AtlasDimensionX ; i++)
		{
		for (int j =0; j < AtlasDimensionY ;j++)
			{
			for (int k = 0; k < AtlasDimensionZ; k++)
				{
				PixelType not_background_value_Atlas = iterforegroundAtlas.Get();
				if (not_background_value_Atlas == 1)
					{
					alambda0[i][j][k] = iteratlaslambda3.Get(); // * scaling_fact;
					alambda1[i][j][k] = iteratlaslambda2.Get(); // * scaling_fact;
					alambda2[i][j][k] = iteratlaslambda1.Get(); // * scaling_fact;
					atlascount++;
					}
				else
					{
					alambda0[i][j][k] = 0;
					alambda1[i][j][k] = 0;
					alambda2[i][j][k] = 0;
					}
				++iteratlaslambda1;
				++iteratlaslambda2;
				++iteratlaslambda3;
				++iterforegroundAtlas;
				}
			}
		}
	
	outFileName.erase();
    	outFileName.append("atlas_lambda1");
	outFileName.append(format);
	//cout << "Output file Atlas Lambda1 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writer4 = VolumeWriterType::New();
	writer4->SetFileName(outFileName.c_str()); 
    	writer4->SetInput(atlaslambda1);
    	writer4->Write();

	outFileName.erase();
    	outFileName.append("atlas_lambda2");
	outFileName.append(format);
	//cout << "Output file Atlas Lambda2 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writer5 = VolumeWriterType::New();
	writer5->SetFileName(outFileName.c_str()); 
    	writer5->SetInput(atlaslambda2);
    	writer5->Write();

	outFileName.erase();
    	outFileName.append("atlas_lambda3");
	outFileName.append(format);
	//cout << "Output file Atlas Lambda3 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writer6 = VolumeWriterType::New();
	writer6->SetFileName(outFileName.c_str()); 
    	writer6->SetInput(atlaslambda3);
    	writer6->Write();

	//////////////forming the 3D bins

	float MAXIMUM_lambda = 0.0;
	if (max0 > max1 && max0 > max2)
		MAXIMUM_lambda = max0;
	else if (max1 > max0 && max1 > max2)
		MAXIMUM_lambda = max1;
	else 
		MAXIMUM_lambda = max2;
	
	int bin_size = 10;
	int number_bins = round((MAXIMUM_lambda*scaling_fact*2)/bin_size) + 1 ; //replace 500 by max eigen value
	
	//vector< vector < vector<int> > > hist_image (number_bins, vector< vector<int> > (number_bins, vector<int>(number_bins,0)));
	//vector< vector < vector<int> > > hist_atlas (number_bins, vector< vector<int> > (number_bins, vector<int>(number_bins,0)));
	
	//normalization

	//equation of plane
	//plane_intersection_point represented by variable plane
	//three points defining the plane are ppt1, ppt2, ppt3
	int plane;
	float t, atlast;
	int pdf[number_bins]; //number_bins
	int atlaspdf[number_bins];
	for (int i =0; i<number_bins ; i++)
		{
		pdf[i] = 0;
		atlaspdf[i] = 0;
		}
	

	for ( int i = 0 ; i< DimensionX ; i++)
		{
		for (int j = 0 ; j < DimensionY ; j++)
			{
			for (int k = 0 ; k < DimensionZ ; k++)
				{
				if (lambda0[i][j][k] != 0 && lambda1[i][j][k] != 0 && lambda2[i][j][k] != 0)
					{
					for (int l = 0 ; l < (number_bins) ; l++) //50 planes 500 max each plane at 10
						{
						plane = (l+1)*bin_size;
						t = plane / (lambda0[i][j][k] + lambda1[i][j][k] + lambda2[i][j][k]);
						//cout << "t: " << t << "		";
						if (t > 1) 
							{
							pdf[l]++;
							break;
							}
						}
					}
	 			}
			}
		}
		
	for ( int i = 0 ; i< AtlasDimensionX ; i++)
		{
		for (int j = 0 ; j < AtlasDimensionY ; j++)
			{
			for (int k = 0 ; k < AtlasDimensionZ ; k++)
				{
				if (alambda0[i][j][k] != 0 && alambda1[i][j][k] != 0 && alambda2[i][j][k] != 0)
					{
					for (int l = 0 ; l < (number_bins) ; l++) //50 planes 500 max each plane at 10
						{
						plane = (l+1)*bin_size;
						atlast = plane / (alambda0[i][j][k] + alambda1[i][j][k]	+ alambda2[i][j][k]);
						//cout << "atlast: " << atlast << "	";
						if (atlast > 1) 
							{
							atlaspdf[l]++;
							break;
							}
						}
					}
				}
			}
		}
	
	double probability_density_fn[number_bins];
	double cdf[number_bins];
	double atlas_probability_density_fn[number_bins];
	double atlas_cdf[number_bins];
	for (int l = 0 ; l < number_bins ; l++)
		{
		probability_density_fn[l] = (double)pdf[l]/casecount;
		atlas_probability_density_fn[l] = (double)atlaspdf[l]/atlascount;
		if (l == 0)
			{
			cdf[l] = probability_density_fn[l];
			atlas_cdf[l] = atlas_probability_density_fn[l];	
			//cout << "initialize cdf assignment " <<endl;
			}
		else	
			{
			cdf[l] = probability_density_fn[l] + cdf[l-1];
			atlas_cdf[l] = atlas_probability_density_fn[l] + atlas_cdf[l-1];
			}
		}

	//closest case plane to points
	vector< vector < vector<float> > > min_dist (DimensionX, vector< vector<float> > (DimensionY, vector<float>(DimensionZ,0)));
	vector< vector < vector<int> > > closest_plane (DimensionX, vector< vector<int> > (DimensionY, vector<int>(DimensionZ,0)));
	float dist;
	for ( int i = 0 ; i< DimensionX ; i++)
		{
		for (int j = 0 ; j < DimensionY ; j++)
			{
			for (int k = 0 ; k < DimensionZ ; k++)
				{
				if (lambda0[i][j][k] != 0 && lambda1[i][j][k] != 0 && lambda2[i][j][k] != 0)
					{
					for (int l = 0; l < (number_bins) ; l++)
						{
						plane = (l + 1)*bin_size;
						dist = abs(plane - (lambda0[i][j][k] + lambda1[i][j][k] + lambda2[i][j][k]));
						if (l>0)	
							{
							if (dist < min_dist[i][j][k])
								{
								min_dist[i][j][k] = dist;
								closest_plane[i][j][k] = l;
								}
							}
						else
							{
							min_dist[i][j][k] = dist;
							closest_plane[i][j][k] = 0;
							}
						}
					//cout << "Closest plane of point " << lambda0[i][j][k] <<"  " << lambda1[i][j][k] <<  "  " << lambda2[i][j][k] << " is " << closest_plane[i][j][k] << endl;
					//cout << " check closest planes of points and enter to continue" <<" i: " << i << " j: " <<  j << " k : " << k << endl;
					//cout << " enter to continue: " << endl;
					//getchar(); 
					}
				}
			}
		}
	
	//create FA values look up table
	//vector< vector < vector<float> > > FAlut (500, vector< vector<float> > (500, vector<float>(500,0.0)));
	//std::vector<planemap> cdfmap;
	float cdf_dif, tA, closest_cdf_ptX, closest_cdf_ptY, closest_cdf_ptZ, BestFA;
	float min_cdf = 0.0;
	vector< vector < vector<int> > > closest_atlas_plane (DimensionX, vector< vector<int> > (DimensionY, vector<int>(DimensionZ,0)));
	

	int temp1;
	int *selected_Atlas_planes = new int [number_bins * number_bins * number_bins * 1000];
	//int *selected_Atlas_planes;
	int count_Aselected_planes = 0;
	for ( int i = 0 ; i< DimensionX ; i++)
		{
		for (int j = 0 ; j < DimensionY ; j++)
			{
			for (int k = 0 ; k < DimensionZ ; k++)
				{
				
				temp1 = closest_plane[i][j][k];
				for (int l = 0; l < number_bins ; l++)
					{
					cdf_dif = abs(cdf[temp1] - atlas_cdf[l]);
					if (l>0)
						{
						if (cdf_dif < min_cdf)
							{
							min_cdf = cdf_dif;  
							closest_atlas_plane[i][j][k] = l;
							}
						}
					else
						{
						min_cdf = cdf_dif;
						closest_atlas_plane[i][j][k] = 0;
						}
					}
				*(selected_Atlas_planes+count_Aselected_planes) = closest_atlas_plane[i][j][k];
			//	cout << *(selected_Atlas_planes+count_Aselected_planes) << "  ";
				count_Aselected_planes++;
				}
			}
		}  

	/*cout << "Selected Atlas Planes : " << endl;
 
	for (int i =0; i<count_Aselected_planes; i++)
		{
		cout << "    " << *(selected_Atlas_planes+i);
		} 
	getchar();*/

	int *selected_temp;
	selected_temp = selected_Atlas_planes;
	//count_Aselected_planes = 7;
	//remove duplicated values from *selected_Atlas_planes
	int *current;
	int *end;
	end = selected_Atlas_planes + count_Aselected_planes -1;

	for (current = selected_Atlas_planes+1; selected_Atlas_planes < end; selected_Atlas_planes++)
		{
		current = selected_Atlas_planes+1;
		while (current < end)
			{
	//		cout << *selected_Atlas_planes << "  " << *current << "  " << *end << endl;
			if (*current == *selected_Atlas_planes)
				{
				*current = *end--;
				count_Aselected_planes--;
//				for (int j = 0 ; j< count_Aselected_planes; j++)
//					cout << " Values: " << *(selected_temp + j);
				}
			else
				{
				current++;
				}
			}
		}

	/*cout << "Selected Atlas Planes after removing duplicates: " << endl;
	for (int j = 0 ; j< count_Aselected_planes; j++)
	cout << "  " << *(selected_temp + j);
	
	cout << "test : " << *(selected_temp + 1) << endl;
	getchar();*/
	
	
	int *cdfnumpoints = new int [number_bins];
	float *cdflambda1 = new float [number_bins * number_bins * number_bins * 1000];
	float *cdflambda2 = new float [number_bins *number_bins * number_bins * 1000];
	float *cdflambda3 = new float [number_bins *number_bins * number_bins * 1000];
	float *cdfFAvalue = new float [number_bins *number_bins * number_bins * 1000];

	float FAlut;
	/*ofstream FAvaluesOutput;
	FAvaluesOutput.open("Favaluesfile.txt");*/
	//int countFAlut = 0;
	int count_eachplane = 0;
	cout << "Computing FA on CDF planes.........."<<endl;
	//for (int planeaxis = 10; planeaxis < (number_bins*bin_size) ; planeaxis = planeaxis+ 10)
	//	{
	/*	for (int i =0 ; i < 1000; i++)
		{
		if (selected_Atlas_planes*/	
		//count_eachplane = 0;
	//int count_eachplane_actual[count_Aselected_planes] = 0;
	int planeaxis;
	for (int i = 0 ; i < count_Aselected_planes ; i++)
		{
		planeaxis = (*(selected_temp+i)) * 10;
		for (float t1 = 0; t1 < (number_bins*bin_size); t1=t1+1)
			{
			for (float t2 = 0; t2 < (number_bins*bin_size); t2=t2+1)
				{
				for (float t3 = 0; t3 < (number_bins*bin_size); t3=t3+1)
					{
						if (t1 + t2 +t3 - planeaxis == 0)
						{
						FAlut = sqrt(0.5) * sqrt((pow((t1-t2),2) + pow((t1-t3),2) + pow((t2-t3),2)) / (pow(t1,2) + pow(t2,2) + pow(t3,2)));
						//FAlut = (t2 + t3)/3;
						*(cdfFAvalue+count_eachplane) = FAlut;
						*(cdflambda1+count_eachplane) = t1;
						*(cdflambda2+count_eachplane) = t2;
						*(cdflambda3+count_eachplane) = t3;
						
	//					cout << "Lambda Values: " << *(cdflambda1+count_eachplane);
	//					cout << " " << *(cdflambda2+count_eachplane);
	//					cout << " " << *(cdflambda3+count_eachplane);
	//					cout << "  FA:" << *(cdfFAvalue+count_eachplane) << endl;
						count_eachplane++;
						}
					}
				}
			}
		//cout << "Counts in each plane:" << count_eachplane << endl;
		*(cdfnumpoints + i) = count_eachplane;
		//*(cdfnumpoints + (planeaxis/10)) = count_eachplane;
		//cdfnumpoints++;
		}

	//cout << "enter to continue: " <<  endl;
	//getchar();
	//closest atlas plane to case plane and closest point on plane to each point
	//substituting the cdf values
	int temp, FAmatchfoundcount = 0;
	//float adjusted_value_temp;
//	int count_negative_lambda1 = 0;
//	int count_negative_lambda2 = 0;
//	int count_negative_lambda3 = 0;
	int number_voxels_not_norm = 0;
	float maxtA = 0.0;
	float dist_FA;
	float mindist = 0.0;
	float bestlambda1 = 0.0;
	float bestlambda2 = 0.0;
	float bestlambda3 = 0.0;
	bool FAmatchfound;
	vector< vector < vector<float> > > newlambda0 (DimensionX, vector< vector<float> > (DimensionY, vector<float>(DimensionZ,0)));
	vector< vector < vector<float> > > newlambda1 (DimensionX, vector< vector<float> > (DimensionY , vector<float>(DimensionZ,0)));
	vector< vector < vector<float> > > newlambda2 (DimensionX, vector< vector<float> > (DimensionY, vector<float>(DimensionZ,0)));
	iterImageNormFA.GoToBegin();
	iterImageNormFA.SetFirstDirection(0);
	iterImageNormFA.SetSecondDirection(1);
	int testcount = 0;
	int tempi, tempj, tempk;

	
for ( int i = 0 ; i< DimensionX ; i++)
		{
		for (int j = 0 ; j < DimensionY ; j++)
			{
			for (int k = 0 ; k < DimensionZ ; k++)
				{
				FAmatchfound = false;
				if (lambda0[i][j][k] != 0 && lambda1[i][j][k] != 0 && lambda2[i][j][k] != 0)
					{
					testcount++;
					temp = closest_plane[i][j][k];
					for (int l = 0; l < number_bins ; l++)
						{
						cdf_dif = abs(cdf[temp] - atlas_cdf[l]);
						if (l>0)
							{
							if (cdf_dif < min_cdf)
								{
								min_cdf = cdf_dif;  
								closest_atlas_plane[i][j][k] = l;
								}
							}
						else
							{
							min_cdf = cdf_dif;
							closest_atlas_plane[i][j][k] = 0;
							}
						}
					//cout << " cdf of plane closest to point: " << cdf[temp] << " cdf of atlas plane " << atlas_cdf[closest_atlas_plane[i][j][k]] << endl;
					int closest_atlas_plane_value = closest_atlas_plane[i][j][k];
					//cout << "Closest atlas plane: " << closest_atlas_plane_value << endl;
					//cout << " minimum cdf : " << min_cdf << endl;
				//intersection point of closest atlas plane to point
				//equation of plane x + y + z = plane
				//parametric equation of line normal to plane x = x(lambda) + t; y = y(lambda) + t; z = z(lambda) + t;
					int G = (closest_atlas_plane_value+1)*bin_size;
					tA = (G - (lambda0[i][j][k] + lambda1[i][j][k] + lambda2[i][j][k]))/(3);
					
					closest_cdf_ptX = lambda0[i][j][k]+tA;
					closest_cdf_ptY = lambda1[i][j][k]+tA;
					closest_cdf_ptZ = lambda2[i][j][k]+tA;
					
					BestFA = (double)(iterImageNormFA.Get() / 10000);
					
					int count_simFA = 0;
					float FAsim[10000] = {0.0};
					int seq_start, seq_end, num_points_onplane;
					//get number of points on the plane and then check for similar FA
					int pos_in_Array = 0;
					int y;
					for (y = 0; y < count_Aselected_planes; y++)
						{
						if (*(selected_temp + y) == closest_atlas_plane_value)
							{
							pos_in_Array = y;
							break;
							}
						}

					if (pos_in_Array == 0)
						{
						num_points_onplane = *(cdfnumpoints);
						seq_start = 0;
						seq_end = *(cdfnumpoints);
						}
					else	
						{
						num_points_onplane = *(cdfnumpoints+pos_in_Array) - *(cdfnumpoints+pos_in_Array-1); 
					//	cout << " Number of points on plane: " << num_points_onplane;
						seq_start = *(cdfnumpoints+y-1) + 1;
						seq_end =  *(cdfnumpoints+y);
						}
					
					float difference_FA = 0.0;
					float min_dif = 0.0;
					for (int l = seq_start ; l < seq_end; l++)
						{
						float MatchFAval = *(cdfFAvalue+l);
						difference_FA = abs (MatchFAval - BestFA);
						if (l == seq_start)
							{
							min_dif = difference_FA;
							}
						else if (difference_FA < min_dif)
							{
							min_dif = difference_FA;
							}
						}
					for (int l = seq_start ; l < seq_end; l++)
						{
						float MatchFAval = *(cdfFAvalue+l);
						float MatchLambda1 = *(cdflambda1+l);
						float MatchLambda2 = *(cdflambda2+l);
						float MatchLambda3 = *(cdflambda3+l);
						
						if ((abs(MatchFAval -BestFA) == min_dif) && (MatchLambda1 > 0 && MatchLambda2 > 0  && MatchLambda3 > 0))
						//if ((abs(MatchFAval -BestFA) < 0.01) && (MatchLambda1 > 0 && MatchLambda2 > 0  && MatchLambda3 > 0))
							{
					//		cout << " BestFA: " << BestFA << " FA match found: " << MatchFAval << endl;
					//		cout << " lambda1: " << MatchLambda1 << " lambda2: " << MatchLambda2 << " lambda3: " << MatchLambda3 << endl;
							FAmatchfound = true;
							
							count_simFA++;
							FAsim[count_simFA] = MatchFAval;
							dist_FA = pow((closest_cdf_ptX - MatchLambda1),2) + pow((closest_cdf_ptY - MatchLambda2),2) + pow((closest_cdf_ptZ - MatchLambda3),2);
							if (count_simFA == 1)
								{
								mindist = dist_FA;
								bestlambda1 = MatchLambda1;
								bestlambda2 = MatchLambda2;
								bestlambda3 = MatchLambda3; 
								}
								else
								{
								if (dist_FA < mindist)
									{
									mindist = dist_FA;
									bestlambda1 = MatchLambda1;
									bestlambda2 = MatchLambda2;
									bestlambda3 = MatchLambda3;
									}
								}
							}
						}
						
					if (FAmatchfound == true)
					{
				//	cout << " Close FA value found and lambda values substituted" << endl;
					newlambda0[i][j][k] = bestlambda1;
					newlambda1[i][j][k] = bestlambda2;
					newlambda2[i][j][k] = bestlambda3;
					FAmatchfoundcount++;
					}
					else
					{
						{
						newlambda0[i][j][k] = lambda0[i][j][k];
						newlambda1[i][j][k] = lambda1[i][j][k];
						newlambda2[i][j][k] = lambda2[i][j][k];
						number_voxels_not_norm++;
						}
					}
					
					}
				++iterImageNormFA;
				}
			}
		}
		

/*	cout << "FA match found for " << FAmatchfoundcount << "number of voxels." << endl;
	cout << "Maximum adjustment to lambda values is : " << maxtA << endl;
	//cout << "Values mapped to negative lambda1: " << count_negative_lambda1 << "  lambda2: " << count_negative_lambda2 << " lambda3: " << count_negative_lambda3 << endl;
	cout << "Number voxels not normalized are: " << number_voxels_not_norm << endl;*/

	DuplicatorType::Pointer duplicatorN1 = DuplicatorType::New();
     	duplicatorN1->SetInputImage(FAorigImage);
     	duplicatorN1->Update();
	DuplicatorType::Pointer duplicatorN2 = DuplicatorType::New();
     	duplicatorN2->SetInputImage(FAorigImage);
     	duplicatorN2->Update();
	DuplicatorType::Pointer duplicatorN3 = DuplicatorType::New();
     	duplicatorN3->SetInputImage(FAorigImage);
     	duplicatorN3->Update();

     	ImageType::Pointer imagelambdanew1 = duplicatorN1->GetOutput();
	ImageType::Pointer imagelambdanew2 = duplicatorN2->GetOutput();
	ImageType::Pointer imagelambdanew3 = duplicatorN3->GetOutput();
	
    	//computing eigen values
	IteratorType iterImagelambdaNew1 (imagelambdanew1, imagelambdanew1->GetLargestPossibleRegion());	
	IteratorType iterImagelambdaNew2 (imagelambdanew2, imagelambdanew2->GetLargestPossibleRegion());
	IteratorType iterImagelambdaNew3 (imagelambdanew3, imagelambdanew3->GetLargestPossibleRegion());

	IteratorType iterImagelambdar1 (imagelambda1, imagelambda1->GetLargestPossibleRegion());	
	IteratorType iterImagelambdar2 (imagelambda2, imagelambda2->GetLargestPossibleRegion());
	IteratorType iterImagelambdar3 (imagelambda3, imagelambda3->GetLargestPossibleRegion());

	iterImagelambdar1.GoToBegin();
	iterImagelambdar2.GoToBegin();
	iterImagelambdar3.GoToBegin();
	
	iterImagelambdaNew1.GoToBegin();
	iterImagelambdaNew2.GoToBegin();
	iterImagelambdaNew3.GoToBegin();

	iterImagelambdar1.SetFirstDirection(0);
	iterImagelambdar2.SetFirstDirection(0);
	iterImagelambdar3.SetFirstDirection(0);
	iterImagelambdaNew1.SetFirstDirection(0);
	iterImagelambdaNew2.SetFirstDirection(0);
	iterImagelambdaNew3.SetFirstDirection(0);

	iterImagelambdar1.SetSecondDirection(1);
	iterImagelambdar2.SetSecondDirection(1);
	iterImagelambdar3.SetSecondDirection(1);
	iterImagelambdaNew1.SetSecondDirection(1);
	iterImagelambdaNew2.SetSecondDirection(1);
	iterImagelambdaNew3.SetSecondDirection(1);

	for ( int i = 0 ; i< DimensionX ; i++)
		{
		for (int j = 0 ; j < DimensionY ; j++)
			{
			for (int k = 0 ; k < DimensionZ ; k++)
				{
				if (lambda0[i][j][k] != 0 && lambda1[i][j][k] != 0 && lambda2[i][j][k] != 0)
					{
					iterImagelambdaNew3.Set(newlambda0[i][j][k]); //scaling_fact);
					iterImagelambdaNew2.Set(newlambda1[i][j][k]); //scaling_fact);
					iterImagelambdaNew1.Set(newlambda2[i][j][k]); //scaling_fact);
					}
				else
					{
					iterImagelambdaNew1.Set(iterImagelambdar1.Get());
					iterImagelambdaNew2.Set(iterImagelambdar2.Get());
					iterImagelambdaNew3.Set(iterImagelambdar3.Get());
					}
				++iterImagelambdar1;
				++iterImagelambdar2;
				++iterImagelambdar3;

				++iterImagelambdaNew1;
				++iterImagelambdaNew2;
				++iterImagelambdaNew3;
				}
			}
		}

	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_norm_lambda3");
	outFileName.append(format);
	//cout << "Output file Normalized Lambda3 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writerN1 = VolumeWriterType::New();
	writerN1->SetFileName(outFileName.c_str()); 
    	writerN1->SetInput(imagelambdanew1);
    	writerN1->Write();

	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_norm_lambda2");
	outFileName.append(format);
	//cout << "Output file Normalized Lambda2 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writerN2 = VolumeWriterType::New();
	writerN2->SetFileName(outFileName.c_str()); 
    	writerN2->SetInput(imagelambdanew2);
    	writerN2->Write();

	outFileName.erase();
    	outFileName.append(base_string);
    	outFileName.append("_norm_lambda1");
	outFileName.append(format);
	//cout << "Output file Normalized Lambda1 : " << outFileName.c_str() <<endl;
	VolumeWriterType::Pointer writerN3 = VolumeWriterType::New();
	writerN3->SetFileName(outFileName.c_str()); 
    	writerN3->SetInput(imagelambdanew3);
    	writerN3->Write();


	Matrix< double , 3 , 3 > mat; 
    	Matrix< double , 3 , 3 > mat1;
	Matrix< double,  3 , 3 > mat2;
	mat.Fill(0);

	typedef itk::ImageDuplicator< DiffusionImageType > DTIDuplicatorType;
	DTIDuplicatorType::Pointer duplicatorDTI1 = DTIDuplicatorType::New();
     	duplicatorDTI1->SetInputImage(inputDTImage);
     	duplicatorDTI1->Update();

	DiffusionImageType::Pointer outputDTImage = duplicatorDTI1->GetOutput();
	IterDTI iterDTImage2 (outputDTImage, outputDTImage->GetLargestPossibleRegion());
	InputTensorDataType inputTensorf;
	InputTensorDataType outputTensorf;

	IteratorType iterImagelambdaNew11 (imagelambdanew1, imagelambdanew1->GetLargestPossibleRegion());	
	IteratorType iterImagelambdaNew21 (imagelambdanew2, imagelambdanew2->GetLargestPossibleRegion());
	IteratorType iterImagelambdaNew31 (imagelambdanew3, imagelambdanew3->GetLargestPossibleRegion());

	iterDTImage2.SetFirstDirection(0);
	iterImagelambdaNew11.SetFirstDirection(0);
	iterImagelambdaNew21.SetFirstDirection(0);
	iterImagelambdaNew31.SetFirstDirection(0);

	iterDTImage2.SetSecondDirection(1);
	iterImagelambdaNew11.SetSecondDirection(1);
	iterImagelambdaNew21.SetSecondDirection(1);
	iterImagelambdaNew31.SetSecondDirection(1);

	iterDTImage2.GoToBegin();
	iterImagelambdaNew11.GoToBegin();
	iterImagelambdaNew21.GoToBegin();
	iterImagelambdaNew31.GoToBegin();
	
	while (!iterDTImage2.IsAtEnd())
		{
		while (!iterDTImage2.IsAtEndOfSlice())
			{
			while (!iterDTImage2.IsAtEndOfLine())
				{
				inputTensorf = iterDTImage2.Get();
				inputTensorf.ComputeEigenAnalysis(feigenValues, feigenVectors);
				mat[0][0] = (iterImagelambdaNew31.Get()/scaling_fact);
				mat[1][1] = (iterImagelambdaNew21.Get()/scaling_fact);
				mat[2][2] = (iterImagelambdaNew11.Get()/scaling_fact);

				/*mat[0][0] = eigenValues[0];
				mat[1][1] = eigenValues[1];
				mat[2][2] = eigenValues[2];*/

				mat2 = feigenVectors.GetTranspose();
				mat1 = mat2 * mat * feigenVectors;
				for (int i=0; i<3; i++)
					{
					for (int j=i; j<3; j++)
						{
						outputTensorf(i,j) = mat1[i][j];
 						}
					}
				iterDTImage2.Set(outputTensorf);
				++iterDTImage2;
				++iterImagelambdaNew11;
				++iterImagelambdaNew21;
				++iterImagelambdaNew31;
				}
			iterDTImage2.NextLine();
			iterImagelambdaNew11.NextLine();
			iterImagelambdaNew21.NextLine();
			iterImagelambdaNew31.NextLine();
			}	
		iterDTImage2.NextSlice();
		iterImagelambdaNew11.NextSlice();
		iterImagelambdaNew21.NextSlice();
		iterImagelambdaNew31.NextSlice();
		}
			

	DiffusionWriterType::Pointer diffusionWriter = DiffusionWriterType::New() ;
	
	outDTIFileName.erase();
    	outDTIFileName.append(base_string);
    	outDTIFileName.append("_Norm_DTI");
	outDTIFileName.append(format);

	cout << "Output NORMALIZED DTI file name : " << outDTIFileName.c_str() <<endl;

	//VolumeWriterType::Pointer writer = VolumeWriterType::New();
	diffusionWriter->SetFileName(outDTIFileName.c_str()); 
    	diffusionWriter->SetInput(outputDTImage);
    	try
		{
		diffusionWriter->Write();
		}
	catch (ExceptionObject & err) 
		{
      		cerr << "Error in writing DTImage... ExceptionObject caught!" << endl;
      		cerr << err << endl;
      		return EXIT_FAILURE;	
    		}



}
