/*=========================================================================

Program:   DTI crop 
Module:    $RCSfile: CropDTI.cxx,v $
Language:  C++
Date:      $Date: 2009/04/14 11:53:47 $
Version:   $Revision: 1.2 $

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
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
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>

#include <itkExtractImageFilter.h>
#include <itkConstantPadImageFilter.h>

#include <itkDiffusionTensor3DReconstructionImageFilter.h>

#include "argio.h"

using namespace itk;
using namespace std;

static int debug = 1;

int main(int argc, const char* argv[])
{


	if (argc <= 1 || ipExistsArgument(argv,"-h ") || ipExistsArgument(argv,"-help") || 
		ipExistsArgument(argv,"-usage"))
	{
		cout << "usage: " << argv[0] <<" infile [-o outfile] [-region  px,py,pz,w,h,d | -size w,h,d] [-v] " << endl
			 << "  -region px,py,pz,w,h,d crop/uncrop(embed/pad with 0-matrix) DTI image " << endl
			 << "          including starting point index px, py, pz" << endl
			 << "          with dimensions width(w),height(h) and depth(d)" << endl
			 << "  -size   w,h,d crop/uncrop(embed/pad with 0-matrix) DTI image with" << endl
			 << "          the input data put in center and with sizes of width(w),height(h) and depth(d)" << endl
			 << "  -v      verbose mode" << endl
			 << endl;
		exit(0) ;
	}

	// argument processing
	string fileName(argv[1]);

	debug =  ipExistsArgument(argv,"-v"); 

	string outfileName(ipGetStringArgument(argv, "-o", ""));

	if (outfileName.empty()) {
		outfileName = "output.mha";
		cout << "no outputname specified using " << outfileName << endl;
	} 

	if (debug) cout << "Output filename: " << outfileName << endl;

	bool regionOn = ipExistsArgument(argv,"-region"); 
	const int numCropParam = 6;
	int cropParam[numCropParam];

	if (regionOn)
	{  
		char *tmp_str    = ipGetStringArgument(argv, "-region", NULL);
		int numDim       = ipExtractIntTokens(cropParam, tmp_str, numCropParam);
		if (numDim != numCropParam) {              
			cerr << argv[0] << ": region needs "<< numCropParam << " parameters.\n";
			exit(-1);
		}
		else
			if (debug) cout << "region parameters: "<< cropParam[0]<< " " << cropParam[1]<< " "<< cropParam[2]
			<< " "<< cropParam[3]<< " "<< cropParam[4]<< " "<< cropParam[5]<< endl;
		free(tmp_str);
	}

	bool sizeOn = ipExistsArgument(argv,"-size"); 
	const int numSizeParam = 3;
	int sizeParam[numSizeParam];

	if (sizeOn)
	{  
		char *tmp_str    = ipGetStringArgument(argv, "-size", NULL);
		int numDim       = ipExtractIntTokens(sizeParam, tmp_str, numSizeParam);
		if (numDim != numSizeParam) {              
			cerr << argv[0] << ": size needs "<< numSizeParam << " parameters.\n";
			exit(-1);
		}
		else
			if (debug) cout << "region parameters: "<< sizeParam[0]<< " " << sizeParam[1]<< " "<< sizeParam[2]	<< endl;
			free(tmp_str);
	}

	if (regionOn && sizeOn) {
		cerr << "region and size cannot be used at the same time" << endl;
		exit(-1);
	}  

	if( !regionOn && !sizeOn)
	{
		cerr << "-region or -size, one of them must be set" << endl;
		exit(-1);
	}

	// argument processing done, do real work

	// Definitions
	typedef short ShortPixelType;
	typedef double DoublePixelType;

	typedef DiffusionTensor3DReconstructionImageFilter< ShortPixelType, ShortPixelType, DoublePixelType >  TensorReconstructionImageFilterType;
	typedef TensorReconstructionImageFilterType::GradientDirectionContainerType    GradientDirectionContainerType;
	typedef TensorReconstructionImageFilterType::TensorPixelType    TensorPixelType;

	const int Dimension = 3;
	typedef itk::Point< double, Dimension >			PointType;
	typedef itk::Image< TensorPixelType, Dimension> DtiImageType;
	typedef DtiImageType::RegionType				DtiImageRegionType;
	typedef TensorPixelType::RealValueType          RealValueType;

	typedef itk::ImageFileReader<DtiImageType>  DTIReaderType;
	typedef itk::ImageFileWriter<DtiImageType>  DTIWriterType;

	// read file
	DTIReaderType::Pointer imageReader = DTIReaderType::New();
	if (debug) cout << "Loading file " << endl;
	try
	{
		imageReader->SetFileName(fileName.c_str()) ;
		imageReader->Update() ;
	}
	catch (ExceptionObject e)
	{
		e.Print(cout) ;
		exit(0) ;
	}

	DtiImageType::RegionType imageOriginalRegion = imageReader->GetOutput()->GetLargestPossibleRegion();
	DtiImageType::SizeType imageOriginalSize = imageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

	DtiImageType::RegionType	outputRegion;
	DtiImageType::IndexType		outputIndex;
	DtiImageType::SizeType		outputSize;

	if(sizeOn)
	{
		cropParam[0] = (imageOriginalSize[0] - sizeParam[0])/2; 
		cropParam[1] = (imageOriginalSize[1] - sizeParam[1])/2; ; 
		cropParam[2] = (imageOriginalSize[2] - sizeParam[2])/2; ; 

		cropParam[3] = sizeParam[0]; 
		cropParam[4] = sizeParam[1]; 
		cropParam[5] = sizeParam[2]; 
	}

	outputIndex[0] = 0;//cropParam[0]; 
	outputIndex[1] = 0;//cropParam[1]; 
	outputIndex[2] = 0;//ropParam[2]; 

	outputSize[0] = cropParam[3]; 
	outputSize[1] = cropParam[4]; 
	outputSize[2] = cropParam[5]; 

	outputRegion.SetIndex(outputIndex);
	outputRegion.SetSize( outputSize);

	DtiImageType::Pointer	PaddedDTI = DtiImageType::New();
	PaddedDTI->SetLargestPossibleRegion(outputRegion);
	PaddedDTI->SetBufferedRegion(outputRegion);
	PaddedDTI->SetSpacing( imageReader->GetOutput()->GetSpacing() );
	PaddedDTI->SetOrigin(imageReader->GetOutput()->GetOrigin() );
	PaddedDTI->Allocate();

	double tensorPixel0[6]={0.0,0.0,0.0,0.0,0.0,0.0};
	//double tensorPixel1[6]={1.0,0.0,0.0,1.0,0.0,1.0};
	PaddedDTI->FillBuffer(tensorPixel0);

	typedef ImageRegionConstIterator< DtiImageType > constDTIIterator;
	typedef ImageRegionIterator< DtiImageType > DTIIterator;

	DTIIterator iOutput( PaddedDTI,outputRegion);
	iOutput.GoToBegin();

	constDTIIterator iInput(  imageReader->GetOutput(), imageReader->GetOutput()->GetLargestPossibleRegion());
	iInput.GoToBegin();

	DtiImageType::IndexType		index;
	TensorPixelType pixel;
	while( !iOutput.IsAtEnd() )
	{
		index = iOutput.GetIndex();
		index[0]=index[0]+cropParam[0]; 
		index[1]=index[1]+cropParam[1]; 
		index[2]=index[2]+cropParam[2]; 
		//std::cout<< "index: "<< index[0] <<" "<< index[1] <<" " << index[2] <<" "  <<std::endl;
		if(imageOriginalRegion.IsInside(index))
		{
			pixel = imageReader->GetOutput()->GetPixel(index);
			iOutput.Set( pixel  );
			//std::cout<< "in index: "<< index[0] <<" "<< index[1] <<" " << index[2] <<" "  <<std::endl;
			//std::cout<< pixel[0] <<" "<< pixel[1] <<" " << pixel[2] <<std::endl;
			//std::cout<< pixel[1] <<" "<< pixel[3] <<" " << pixel[4] <<std::endl;
			//std::cout<< pixel[2] <<" "<< pixel[4] <<" " << pixel[5] <<std::endl;
		}
		else
		{	
			iOutput.Set( tensorPixel0 );
			//std::cout<< "out index: "<< index[0] <<" "<< index[1] <<" " << index[2] <<" "  <<std::endl;
		}
		++iOutput;	
	}

	if (debug) cout << "writing output data " << outfileName << endl;

	DTIWriterType::Pointer DTIWriter = DTIWriterType::New();
	DTIWriter->SetFileName(outfileName.c_str());
	DTIWriter->SetInput(PaddedDTI);
	DTIWriter->Update();

	return 0 ;
}

