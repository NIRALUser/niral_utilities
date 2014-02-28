/*=========================================================================
 * * Program :   $Insight Segmentation & Registration Toolkit $ * Module  :   $DMDBiomarkerTool: DMDBase.h $ * Purpose :   $The base class of Duchenne Muscle Dystrophy biomarker tools $
 * Language:   $C++ $
 * Date    :   $Date: 2010-06-04 12:36:34 $
 * Version :   $Revision: 0.10 $
 * Authors :   $Jiahui Wang, Martin Styner $
 * Update  :   1. Created the data class (06-04-10)
 * Copyright (c) Nero Image Research and Analysis Lab.@UNC All rights reserved.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even 
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.
 *
=========================================================================*/

//#include "itkOrientedImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
//#include "itkImageDuplicator.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"

#include "itkExtractImageFilter.h"
#include "itkSliceBySliceImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "itkOtsuThresholdImageFilter.h"

#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"

#include <math.h>
#include <iostream>
#include <fstream>

//using namespace std;

#ifndef DMDDATA    // prevent for redefine
#define DMDDATA    
#define SEARCHED_PIXEL_VALUE -1
#define TEMP_PIXEL_VALUE -2


class DMDData
{
    private:
        
    public:
        typedef float                                                      PixelType;
        typedef long			                                   uncharPixelType;
	typedef std::string                                                StrITKType;
        typedef int                                                        IntITKType;
        typedef float                                                      FITKType;
        typedef double                                                     LfITKType;
	DMDData() {

	}
        DMDData( IntITKType d, StrITKType inputDirec ) {
 	    inputDirectory = inputDirec;
	}
        // declare the pixel type and image dimension
        typedef itk::Image< PixelType, 2 >                         OrientedImage2DType;
        typedef itk::Image< PixelType, 3 >                         OrientedImageType;
        typedef itk::Image< uncharPixelType, 3 >                   uncharOrientedImageType;
	// create image file reader and writer
        typedef itk::ImageSeriesReader< OrientedImageType >                OrientedSeriesReaderType;
        typedef itk::ImageFileWriter< OrientedImageType >                  OrientedWriterType;
	typedef itk::Image< PixelType, 3 >                                 ImageReaderInputType;
	typedef itk::ImageFileReader< ImageReaderInputType >               segMuscleReaderType;  
	//	typedef itk::Image< PixelType, 3 >                                 ImageReaderInputType;
        // dicomIO was created to connect GDCMImage (provide services for reading and writing DICOM files) and reader
        typedef itk::GDCMImageIO                                           ImageIOType;
        // to identify from a given directory the set of filenames that belong together to the same volumetric image
        typedef itk::RecursiveGaussianImageFilter<OrientedImageType, OrientedImageType>    FilterType;
	typedef itk::GDCMSeriesFileNames                                   NamesGeneratorType;    
        typedef itk::MetaDataDictionary                                    DictionaryType;            
        typedef std::vector< StrITKType >                                  FileNamesContainer;        
	typedef itk::ImageFileReader< OrientedImageType >                  OrientedImageReaderType;
        typedef itk::ImageRegionConstIterator< OrientedImageType >         ConstIteratorType;
        typedef itk::ImageRegionIterator< OrientedImageType >              IteratorType;
	typedef itk::ConstShapedNeighborhoodIterator< ImageReaderInputType > ShapedNeighborhoodIteratorType;	

        // connected component filter
	typedef itk::ConnectedComponentImageFilter<uncharOrientedImageType, uncharOrientedImageType> itkConnectedComponentFilterType;
        typedef itk::RelabelComponentImageFilter<uncharOrientedImageType, uncharOrientedImageType>    itkRelabelComponentFilterType;
 


	FileNamesContainer                                                 fileNames;
	StrITKType                                                         inputDirectory;
        itkConnectedComponentFilterType::Pointer labeller; //= itkConnectedComponentFilterType::New();
        itkRelabelComponentFilterType::Pointer relabeller;// = itkRelabelComponentFilterType::New();
	
        // tells the GDCMSereiesFileNames object to use additional DICOM information to distinguish unique volumes within the directory.
        // note that SetUseSeriesDetails(true) must be called before SetDirectory()	          
	void dataReader ( StrITKType fileName, OrientedImageReaderType::Pointer & reader );
	void seriesDataReader ( StrITKType seriesIdentifie, OrientedSeriesReaderType::Pointer & reader, ImageIOType::Pointer & dicomIO );
	void dataWriter ( OrientedWriterType::Pointer writer, StrITKType outputFilename );
	void readDicomTag ( ImageIOType::Pointer & dicomIO, std::vector<StrITKType> & tagvalue, std::vector<StrITKType> & tagkey );
	void imageErod ( );
	void smoothGradAnisoDiff( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, int smoothItr, float smoothTimeStep, float smoothConductance );
	void smoothGaussian( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, float sigma );
	//void imageInitialize (OrientedSeriesReaderType::Pointer input, OrientedImageType::Pointer &data, float &voxelVol );
        //void imageInitialize ( OrientedSeriesReaderType::Pointer input, OrientedImageType::Pointer &data );
	void imageInitialize (OrientedImageType::Pointer input, OrientedImageType::Pointer &data, float &voxelVol );
        void imageInitialize ( OrientedImageType::Pointer input, OrientedImageType::Pointer &data );
        float muscleMean ( OrientedSeriesReaderType::Pointer input, OrientedImageReaderType::Pointer mask, int &volume ) ;     
        void morphErod ( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, IntITKType size ) ;
        void morphDilat ( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, IntITKType size ) ;
        void morphMultiGrayErod2DIn3D ( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, IntITKType size ) ;
        void outHistogram ( OrientedImageType::Pointer data, OrientedImageType::Pointer roi, OrientedImageType::Pointer mask, StrITKType Identify );
        //void interpolate3D (OrientedImageType::Pointer inputImage, OrientedImageType::Pointer outputImage);
        void interpolate3D (const OrientedImageType::Pointer inputImage, StrITKType);
        float overlap (OrientedImageType::Pointer, OrientedImageType::Pointer, int, bool, float );
        void connectedLabeling2D (const OrientedImage2DType::Pointer input );
        void connectedLabeling3D (const OrientedImageType::Pointer input );
        void waterTubeSegmentation (const OrientedImageType::Pointer input );
        void fatSuppressBiasIdentify (const OrientedImageType::Pointer inputT2, OrientedImageType::Pointer inputT2FS);
        void histogramMatchinMask (const OrientedImageType::Pointer input, OrientedImageType::Pointer inputRef, OrientedImageType::Pointer output, int bin, int point);
        int connectedComponentRegionGrowing (float lower, float upper, OrientedImageType::Pointer & input, OrientedImageType::IndexType index);
        int connectedComponentLabeling (OrientedImageType::Pointer & input);


        
};
////////////////////////////////////l//////////////////////////////////////////////////////////////////////////////
void DMDData::seriesDataReader( StrITKType seriesIdentifier, OrientedSeriesReaderType::Pointer & reader, ImageIOType::Pointer & dicomIO )
{

    NamesGeneratorType::Pointer        nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->SetDirectory(inputDirectory.c_str());
    reader = OrientedSeriesReaderType::New();
    dicomIO = ImageIOType::New();
    //   std::cout << "Now reading series: " << seriesIdentifier.c_str() << std::endl; 
    fileNames = nameGenerator->GetFileNames( seriesIdentifier );
    
    std::cout << "Now reading series: " << seriesIdentifier << std::endl; 
    reader->SetImageIO( dicomIO );
    std::cout << "Now reading series: " << fileNames[0] << std::endl; 
    reader->SetFileNames( fileNames );
    
    
    try {
        reader->Update();
       
    }
    catch (itk::ExceptionObject &ex) {
       std::cout << ex << std::endl;
       //  continue;
       exit(0);
    }            
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::dataWriter(OrientedWriterType::Pointer writer, StrITKType outputFilename )
{
    writer->SetFileName(outputFilename);   
    //    std::cout  << "Writing the image as " << outputFilename << std::endl << std::endl;
    try {
        writer->Update();
    }
    catch (itk::ExceptionObject &ex){
  	std::cout << "write:" << ex << std::endl;
	//return EXIT_FAILURE;
	exit(0);
    }
}
////////////////////////////////////l//////////////////////////////////////////////////////////////////////////////
void DMDData::readDicomTag( ImageIOType::Pointer & dicomIO, std::vector<StrITKType> & tagvalue, std::vector<StrITKType> & tagkey )
{
    const                                           DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
    DictionaryType::ConstIterator                   itr = dictionary.Begin();
    DictionaryType::ConstIterator                   end = dictionary.End();
    //    std::string                                     outputFilename, fileExtension, lastCh;
    while( itr != end ) {
        typedef itk::MetaDataObject<StrITKType> MetaDataStringType;
        itk::MetaDataObjectBase::Pointer  entry = itr->second;
        MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;
        // take DICOM tag and pass it to the GetLabelFromTag() method of the GDCMImageIO class. 
        // This method checks the DICOM dictionary and returns the string label associated to the tag that we are providing in the tagkey variable. 
        // The method itself return false if the tagkey is not found in the dictionary.
        if( entryvalue ) {
	    tagkey.push_back(itr->first);
            StrITKType labelId;
            //  bool found =  itk::GDCMImageIO::GetLabelFromTag( tagkey, labelId ); 
            // The actual value of the dictionary entry is obtained as a string with the GetMetaDataObjectValue() method.
	    tagvalue.push_back(entryvalue->GetMetaDataObjectValue()); 
            // print out an entry by concatenating the DICOM Name or label, the numeric tag and its actual value.
		    /*  if( found ) {
                std::cout << "(" << tagkey << ") " << labelId;
                std::cout << " = " << tagvalue.c_str() << std::endl;
            }
	            else {
              std::cout << "(" << tagkey <<  ") " << "Unknown";
              std::cout << " = " << tagvalue.c_str() << std::endl;
		      }*/                
	    
        }
        // Finally we just close the loop that will walk through all the Dictionary entries.
        ++itr;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::imageErod ( )
{
    unsigned int element_radius = 1;
    ShapedNeighborhoodIteratorType::RadiusType radius;
    radius.Fill(element_radius);

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::smoothGradAnisoDiff( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, int smoothItr, float smoothTimeStep, float smoothConductance )
{
    typedef itk::GradientAnisotropicDiffusionImageFilter< OrientedImageType, OrientedImageType >  FilterType;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput( input );
    filter->SetNumberOfIterations( smoothItr ); // 1
    filter->SetTimeStep( smoothTimeStep );                 // 0.05
    filter->SetConductanceParameter( smoothConductance );     // 0.75

    filter->Update();

    output = filter->GetOutput();
    

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::smoothGaussian( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, float sigma )
{
   // typedef itk::GradientAnisotropicDiffusionImageFilter< OrientedImageType, OrientedImageType >  FilterType;
    typedef itk::RecursiveGaussianImageFilter< OrientedImageType, OrientedImageType > FilterType;

    FilterType::Pointer filterX = FilterType::New();
    FilterType::Pointer filterY = FilterType::New();
    FilterType::Pointer filterZ = FilterType::New();

    filterX->SetDirection( 0 );
    filterY->SetDirection( 1 );
    filterZ->SetDirection( 2 );

    filterX->SetOrder( FilterType::ZeroOrder );
    filterY->SetOrder( FilterType::ZeroOrder );
    filterZ->SetOrder( FilterType::ZeroOrder );

    filterX->SetNormalizeAcrossScale( false );
    filterY->SetNormalizeAcrossScale( false );
    filterZ->SetNormalizeAcrossScale( false );


    filterX->SetInput( input );
    filterY->SetInput( filterX->GetOutput() );
    filterZ->SetInput( filterY->GetOutput() );

    filterX->SetSigma( sigma );
    filterY->SetSigma( sigma );
    filterZ->SetSigma( sigma );

    filterZ->Update();

    output = filterZ->GetOutput();
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void DMDData::imageInitialize ( OrientedSeriesReaderType::Pointer input, OrientedImageType::Pointer &data, float &voxelVol )
void DMDData::imageInitialize ( OrientedImageType::Pointer input, OrientedImageType::Pointer &data, float &voxelVol )
{
    OrientedImageType::RegionType                  region = input->GetLargestPossibleRegion() ;
    OrientedImageType::SpacingType                 spacing = input->GetSpacing() ;
    OrientedImageType::DirectionType               direction = input->GetDirection() ;
    OrientedImageType::PointType                   origin = input->GetOrigin() ;
   // std::cout << "resolution of the image:  " << spacing[0] << "  " << spacing[1] << "  " << spacing[2] << std::endl;
    voxelVol = spacing[0] * spacing[1] * spacing[2] ;
    data->SetRegions( region );
    data->SetSpacing( spacing );
    data->SetDirection( direction );
    data->SetOrigin( origin );                    
    data->Allocate();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void DMDData::imageInitialize ( OrientedSeriesReaderType::Pointer input, OrientedImageType::Pointer &data )
void DMDData::imageInitialize ( OrientedImageType::Pointer input, OrientedImageType::Pointer &data )
{
    OrientedImageType::RegionType                  region = input->GetLargestPossibleRegion() ;
    OrientedImageType::SpacingType                 spacing = input->GetSpacing() ;
    OrientedImageType::DirectionType               direction = input->GetDirection() ;
    OrientedImageType::PointType                   origin = input->GetOrigin() ;
   // std::cout << "resolution of the image:  " << spacing[0] << "  " << spacing[1] << "  " << spacing[2] << std::endl;
    data->SetRegions( region );
    data->SetSpacing( spacing );
    data->SetDirection( direction );
    data->SetOrigin( origin );                    
    data->Allocate();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::dataReader ( StrITKType fileName, OrientedImageReaderType::Pointer & reader )
{
    reader = OrientedImageReaderType::New();
    reader->SetFileName( fileName );
    reader->Update();

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float DMDData::muscleMean ( OrientedSeriesReaderType::Pointer input, OrientedImageReaderType::Pointer mask, int &volume )
{
    float mean = 0;
    
    DMDData::OrientedImageType::Pointer origImage = input->GetOutput();
    ConstIteratorType constOrigIterator( origImage, origImage->GetRequestedRegion() );
    OrientedImageType::Pointer maskImage = mask->GetOutput();
    ConstIteratorType constMaskIterator( maskImage, maskImage->GetRequestedRegion() );
    if ( origImage->GetLargestPossibleRegion().GetSize()[0] == maskImage->GetLargestPossibleRegion().GetSize()[0] && origImage->GetLargestPossibleRegion().GetSize()[1] == maskImage->GetLargestPossibleRegion().GetSize()[1] && origImage->GetLargestPossibleRegion().GetSize()[2] == maskImage->GetLargestPossibleRegion().GetSize()[2] ) {
        for ( constOrigIterator.GoToBegin(), constMaskIterator.GoToBegin(); !constOrigIterator.IsAtEnd(); ++constOrigIterator, ++constMaskIterator ) {
            if ( constMaskIterator.Get() > 0 ) {
                mean += constOrigIterator.Get();
                volume++;
            }
        }
    }
    
    return mean;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::morphErod ( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, IntITKType size ) 
{
    typedef unsigned char   InputPixelType;
    typedef itk::Image< PixelType, 3 >   InputImageType;
    typedef itk::BinaryBallStructuringElement< PixelType, 3 >                             StructuringElementType;        
    typedef itk::GrayscaleErodeImageFilter< InputImageType, InputImageType, StructuringElementType >  ErodeFilterType;
    typedef itk::GrayscaleDilateImageFilter< InputImageType, InputImageType, StructuringElementType > DilateFilterType;
    
    ErodeFilterType::Pointer  grayscaleErode  = ErodeFilterType::New();
    DilateFilterType::Pointer grayscaleDilate = DilateFilterType::New();
   
    StructuringElementType  structuringElement;

    structuringElement.SetRadius( size );  // ( size x 2 + 1 ) x ( size x 2 + 1 ) structuring element

    structuringElement.CreateStructuringElement();

    grayscaleErode->SetKernel( structuringElement );
    //grayscaleDilate->SetKernel( structuringElement );
    grayscaleErode->SetInput( input );
    
    typedef itk::ImageFileWriter< InputImageType >  WriterType;
    WriterType::Pointer writerErosion  = WriterType::New();
    writerErosion->SetFileName( "erosion.nrrd" );
    writerErosion->SetInput( grayscaleErode->GetOutput() );
    writerErosion->Update();
     
    output = grayscaleErode->GetOutput();
    output->Update(); 
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::morphDilat ( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, IntITKType size ) 
{
    typedef unsigned char   InputPixelType;
    typedef itk::Image< PixelType, 3 >   InputImageType;
    typedef itk::BinaryBallStructuringElement< PixelType, 3 >                             StructuringElementType;        
    typedef itk::GrayscaleDilateImageFilter< InputImageType, InputImageType, StructuringElementType > DilateFilterType;
    
    DilateFilterType::Pointer grayscaleDilate = DilateFilterType::New();
   
    StructuringElementType  structuringElement;

    structuringElement.SetRadius( size );  // ( size x 2 + 1 ) x ( size x 2 + 1 ) structuring element

    structuringElement.CreateStructuringElement();

    grayscaleDilate->SetKernel( structuringElement );
    //grayscaleDilate->SetKernel( structuringElement );
    grayscaleDilate->SetInput( input );
    
    typedef itk::ImageFileWriter< InputImageType >  WriterType;
    WriterType::Pointer writerErosion  = WriterType::New();
    writerErosion->SetFileName( "dilation.nrrd" );
    writerErosion->SetInput( grayscaleDilate->GetOutput() );
    writerErosion->Update();
     
    output = grayscaleDilate->GetOutput();
    output->Update(); 
} 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void DMDData::morphMultiGrayErod2DIn3D ( OrientedImageType::Pointer input, OrientedImageType::Pointer &output, IntITKType size ) 
{
    // apply 2D morphological erosion to each slice of 3D volumetric data
    typedef unsigned char   InputPixelType;
    typedef itk::Image< PixelType, 3 >   InputImageType;
    typedef itk::Image< PixelType, 2 >   OutputImageType;
    typedef itk::BinaryBallStructuringElement< PixelType, 2 >                             StructuringElementType;        
//    typedef itk::GrayscaleErodeImageFilter< OutputImageType, OutputImageType, StructuringElementType >  ErodeFilterType;
//    typedef itk::GrayscaleDilateImageFilter< OutputImageType, OutputImageType, StructuringElementType > DilateFilterType;

    typedef itk::ExtractImageFilter< InputImageType, OutputImageType > FilterType;
   // int numberofslice = 0;

    typedef itk::SliceBySliceImageFilter< InputImageType, InputImageType > slicebysliceFilterType;
    slicebysliceFilterType::Pointer slicebyslicefilter = slicebysliceFilterType::New();
    slicebyslicefilter->SetInput( input );

    typedef itk::GrayscaleErodeImageFilter< slicebysliceFilterType::InternalInputImageType, slicebysliceFilterType::InternalOutputImageType, StructuringElementType >  ErodeFilterType;
    typedef itk::GrayscaleDilateImageFilter< slicebysliceFilterType::InternalInputImageType, slicebysliceFilterType::InternalOutputImageType, StructuringElementType > DilateFilterType;
    ErodeFilterType::Pointer  grayscaleErode  = ErodeFilterType::New();
    DilateFilterType::Pointer grayscaleDilate = DilateFilterType::New();
   
    StructuringElementType  structuringElement;

    structuringElement.SetRadius( size );  // ( size x 2 + 1 ) x ( size x 2 + 1 ) structuring element

    structuringElement.CreateStructuringElement();

    grayscaleErode->SetKernel( structuringElement );
    //grayscaleDilate->SetKernel( structuringElement );
   // grayscaleErode->SetInput( filter->GetOutput() );
    slicebyslicefilter->SetFilter(grayscaleErode);

    /*
    typedef itk::ImageFileWriter< InputImageType >  WriterType;
    WriterType::Pointer writerErosion  = WriterType::New();
    writerErosion->SetFileName( "erosion.nrrd" );
    writerErosion->SetInput( grayscaleErode->GetOutput() );
    writerErosion->Update();*/
    
    //output = grayscaleErode->GetOutput();
    output = slicebyslicefilter->GetOutput();
    output->Update(); 
 
    typedef itk::ImageFileWriter< InputImageType >  WriterType;
    WriterType::Pointer writerErosion  = WriterType::New();
    writerErosion->SetFileName( "../data/erosion.nrrd" );
    writerErosion->SetInput( output );
    writerErosion->Update();
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::outHistogram ( OrientedImageType::Pointer data, OrientedImageType::Pointer roi, OrientedImageType::Pointer mask, StrITKType Identify )
{
    ConstIteratorType constIterator( data, data->GetRequestedRegion() ) ;
    ConstIteratorType constMaskIterator( mask, mask->GetRequestedRegion() ) ;
    IteratorType roiIterator( roi, roi->GetRequestedRegion() ) ;
    StrITKType histofile, musclename[3] = {"semit", "rect", "cran"};
    float min = 0, max = 0;    
    int minBin = 0, maxBin = 0 ;
    int label[3] = {SEMIT, REC_FEM, CRAN_SART} ;
 
    for ( int i = 0; i < 3; i++ ) {
        histofile = Identify + musclename[i] + "_histogram";
        std::ofstream efile( histofile.c_str(), std::ios::out );
        min = 1000;
        max = 0;
        for ( constIterator.GoToBegin(), constMaskIterator.GoToBegin(), roiIterator.GoToBegin() ; !constIterator.IsAtEnd(); ++constIterator, ++constMaskIterator, ++roiIterator ) {
            if ( constMaskIterator.Get() == label[i] ) {
                if (constIterator.Get() < min)
                    min = constIterator.Get(); 
                if (constIterator.Get() > max)
                    max = constIterator.Get(); 
                roiIterator.Set( constIterator.Get() );
            }
            else {
                roiIterator.Set( 0 );
            }
        }
        
  //      cout << Identify << "  " << musclename[i] << ":  " << min << "  " << max << endl; 
        minBin = (int)min;
        maxBin = (int)max;
        if(min < 0) {
            for ( constIterator.GoToBegin(), constMaskIterator.GoToBegin(), roiIterator.GoToBegin() ; !constIterator.IsAtEnd(); ++constIterator, ++constMaskIterator, ++roiIterator ) {
                if ( constMaskIterator.Get() == label[i] ) {
                    roiIterator.Set( constIterator.Get() - min + 1);
                   // if(constIterator.Get() == min)
                     //   std::cout <<  roiIterator.Get() << endl;
                }
                else {
                    roiIterator.Set( 0 );
                }   
            }
            max = max - min + 1;
            min = min - min + 1;
        }
   //     cout << Identify << "  " << musclename[i] << ":  " << min << "  " << max << "  " << minBin << "  " << maxBin  << endl; 
         
        typedef itk::Statistics::ScalarImageToHistogramGenerator< OrientedImageType >   HistogramGeneratorType;
        HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
        histogramGenerator->SetInput( roi );
        histogramGenerator->SetNumberOfBins( (int)max );
        histogramGenerator->SetMarginalScale( 1.0 );
        histogramGenerator->Compute();
        typedef HistogramGeneratorType::HistogramType  HistogramType;
        const HistogramType *histogram = histogramGenerator->GetOutput();
        const unsigned int histogramSize = histogram->Size();

        for( unsigned int bin = 1; bin < histogramSize; bin++ ) {
            if (minBin < 0)
                efile <<  minBin + (int) bin << "  " << histogram->GetFrequency( bin, 0 ) << "\n";
            else
                efile << bin << "  " << histogram->GetFrequency( bin, 0 ) << "\n";
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::interpolate3D ( const OrientedImageType::Pointer input, StrITKType caseID)
{
    OrientedImageType::Pointer output = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceSemit = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceRecFem = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceCranSart = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceVolume = OrientedImageType::New() ;
    OrientedImageType::Pointer segSemit = OrientedImageType::New() ;
    OrientedImageType::Pointer segRecFem = OrientedImageType::New() ;
    OrientedImageType::Pointer segCranSart = OrientedImageType::New() ;
    OrientedImageType::Pointer segVolume = OrientedImageType::New() ;
    DMDData::OrientedWriterType::Pointer                writer = DMDData::OrientedWriterType::New();
    writer->UseCompressionOn();
    OrientedImageType::RegionType                  region = input->GetLargestPossibleRegion() ;
    OrientedImageType::SpacingType                 spacing = input->GetSpacing() ;
    OrientedImageType::DirectionType               direction = input->GetDirection() ;
    OrientedImageType::PointType                   origin = input->GetOrigin() ;
    OrientedImageType::SizeType                    size = input->GetLargestPossibleRegion().GetSize();
    OrientedImageType::IndexType                   indexInput, indexOutput;
    float interpolatedVal = 0, upperVal = 0, bottomVal = 0, pixelVal = 0; 
    float distance = 5.0, upperIndex = 0, bottomIndex = 0;
    float upperFactor = 0, bottomFactor = 0, interpolatThresh = 0.5, inplane = 0, upperSum = 0, bottomSum = 0;
    float currentLabel = 0; 
    std::string outputfilename;

    std::cout << "size of the image:  " <<  size[0] << "  " << size[1] << "  " << size[2] << std::endl;

    output->SetRegions( region );
    output->SetSpacing( spacing );
    output->SetDirection( direction );
    output->SetOrigin( origin );                    
    output->Allocate();
    referenceVolume->SetRegions( region );
    referenceVolume->SetSpacing( spacing );
    referenceVolume->SetDirection( direction );
    referenceVolume->SetOrigin( origin );                    
    referenceVolume->Allocate();

    referenceSemit->SetRegions( region );
    referenceSemit->SetSpacing( spacing );
    referenceSemit->SetDirection( direction );
    referenceSemit->SetOrigin( origin );                    
    referenceSemit->Allocate();
    referenceRecFem->SetRegions( region );
    referenceRecFem->SetSpacing( spacing );
    referenceRecFem->SetDirection( direction );
    referenceRecFem->SetOrigin( origin );                    
    referenceRecFem->Allocate();
    referenceCranSart->SetRegions( region );
    referenceCranSart->SetSpacing( spacing );
    referenceCranSart->SetDirection( direction );
    referenceCranSart->SetOrigin( origin );                    
    referenceCranSart->Allocate();
    segVolume->SetRegions( region );
    segVolume->SetSpacing( spacing );
    segVolume->SetDirection( direction );
    segVolume->SetOrigin( origin );                    
    segVolume->Allocate();
    segSemit->SetRegions( region );
    segSemit->SetSpacing( spacing );
    segSemit->SetDirection( direction );
    segSemit->SetOrigin( origin );                    
    segSemit->Allocate();
    segRecFem->SetRegions( region );
    segRecFem->SetSpacing( spacing );
    segRecFem->SetDirection( direction );
    segRecFem->SetOrigin( origin );                    
    segRecFem->Allocate();
    segCranSart->SetRegions( region );
    segCranSart->SetSpacing( spacing );
    segCranSart->SetDirection( direction );
    segCranSart->SetOrigin( origin );                    
    segCranSart->Allocate();

    for(int l = CRAN_SART; l < CAUD_SART; l++){
        currentLabel = l;
        //search the range of slices of a muscle
        int maxSlice = 0, minSlice = size[2], roundRange = 0;
        bool jumpout = 0;
        for(indexInput[2] = 0; indexInput[2] < (int)size[2]; indexInput[2]++) {
            for(indexInput[1] = 0; indexInput[1] < (int)size[1]; indexInput[1]++) {
                for(indexInput[0] = 0; indexInput[0] < (int)size[0]; indexInput[0]++) {
                    pixelVal = input->GetPixel(indexInput);
                    if(pixelVal == currentLabel){
                        if(indexInput[2] < minSlice)
                            minSlice = indexInput[2];
                        if(indexInput[2] > maxSlice)
                            maxSlice = indexInput[2];
                        jumpout = 1;
                        break;
                    }
                }
                if(jumpout)
                    break;
            }
            if(jumpout)
                jumpout = 0;
        }
        roundRange = (maxSlice - minSlice + 1);
        roundRange = (int)(floor(roundRange / distance) * distance); //rounded range of slices
        std::cout << "the slice range of muscle " << currentLabel << " is " << minSlice << " - " << maxSlice << "  " << roundRange << std::endl;
        for(indexOutput[2] = 0; indexOutput[2] < roundRange; indexOutput[2]++) { // start from the 3rd slice
            for(indexOutput[1] = 1; indexOutput[1] < (int)size[1] - 1; indexOutput[1]++) {
                for(indexOutput[0] = 1; indexOutput[0] < (int)size[0] - 1; indexOutput[0]++) {
                    int k = 3;
                    inplane = 0; upperSum = 0; bottomSum = 0;
                    for (int i = -k / 2 ; i <= k / 2; i++){
                        for (int j = -k / 2 ; j <= k / 2; j++){
                            indexInput[0] = indexOutput[0] + i; 
                            indexInput[1] = indexOutput[1] + j; 
                            upperIndex = (indexOutput[2] / (int)distance + 1) * distance + minSlice;
                            indexInput[2] = (int)upperIndex;
                            upperVal = input->GetPixel(indexInput);
                            if(upperVal == currentLabel)
                                upperVal = input->GetPixel(indexInput);
                            else
                                upperVal = 0;
                            bottomIndex = indexOutput[2] / (int)distance * distance + minSlice;
                            indexInput[2] = (int)bottomIndex;
                            bottomVal = input->GetPixel(indexInput);
                            if(bottomVal == currentLabel){
                                bottomVal = input->GetPixel(indexInput);
                            }
                            else
                                bottomVal = 0;
                            upperFactor = ((float)(indexOutput[2] + minSlice) - bottomIndex) / distance;
                            bottomFactor =  (upperIndex - (float)(indexOutput[2] + minSlice)) / distance;
                            inplane += (upperVal * upperFactor + bottomVal * bottomFactor);
                        }
                    }
                    interpolatedVal = inplane /(k * k);
                    if(interpolatedVal >= (interpolatThresh * currentLabel))
                        interpolatedVal = currentLabel;
                    else
                        interpolatedVal = 0;
                    indexOutput[2] += minSlice;
                    if(output->GetPixel(indexOutput) == 0){
                        output->SetPixel(indexOutput, interpolatedVal);
                    }
                    indexOutput[2] -= minSlice;
                }
            }
        }
    }
    //output interpolated muscle volume
    writer->SetInput( output );  
    std::stringstream fs;
    fs << distance;
    //outputfilename = "../data/" + caseID + "_intervolfull_left_femur_mask" + ".nrrd";
    outputfilename = "../data/" + caseID + "_intervolfull_every" + fs.str() + ".nrrd";
    dataWriter(writer, outputfilename);

    float overlapRate = 0;
    //calculate overlap of semitendinaus
    overlapRate = overlap(input, output, SEMIT, 1, distance);    
    std::cout << "the overlap of semit is:  " << overlapRate << std::endl;
    //calculate overlap of rec fem
    overlapRate = overlap(input, output, REC_FEM, 1, distance);    
    std::cout << "the overlap of rec fem is:  " << overlapRate << std::endl;
    //calculate overlap of cran sart
    overlapRate = overlap(input, output, CRAN_SART, 1, distance);    
    std::cout << "the overlap of cran sart is:  " << overlapRate << std::endl;

    for(indexOutput[2] = 1; indexOutput[2] < (int)size[2] - 1; indexOutput[2]++) {
        for(indexOutput[1] = 1; indexOutput[1] < (int)size[1] - 1; indexOutput[1]++) {
            for(indexOutput[0] = 1; indexOutput[0] < (int)size[0] - 1; indexOutput[0]++) {
                if(indexOutput[2] % (int)distance != 0){
                    segVolume->SetPixel(indexOutput, output->GetPixel(indexOutput));
                    referenceVolume->SetPixel(indexOutput, input->GetPixel(indexOutput));
                    if(output->GetPixel(indexOutput) == SEMIT)
                        segSemit->SetPixel(indexOutput, SEMIT);    
                    if(input->GetPixel(indexOutput) == SEMIT)
                        referenceSemit->SetPixel(indexOutput, SEMIT);
                    if(output->GetPixel(indexOutput) == REC_FEM)
                        segRecFem->SetPixel(indexOutput, REC_FEM);    
                    if(input->GetPixel(indexOutput) == REC_FEM)
                        referenceRecFem->SetPixel(indexOutput, REC_FEM);
                    if(output->GetPixel(indexOutput) == CRAN_SART)
                        segCranSart->SetPixel(indexOutput, CRAN_SART);    
                    if(input->GetPixel(indexOutput) == CRAN_SART)
                        referenceCranSart->SetPixel(indexOutput, CRAN_SART);
                }
                else { // remove the manually segmented slices
                   // reference->SetPixel(indexOutput, input->GetPixel(indexOutput));
                    segVolume->SetPixel(indexOutput, output->GetPixel(indexOutput));
                    referenceVolume->SetPixel(indexOutput, input->GetPixel(indexOutput));
                }
            }
        }
    }         
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void DMDData::interpolate3D ( const OrientedImageType::Pointer input, StrITKType caseID)
{
    OrientedImageType::Pointer output = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceSemit = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceRecFem = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceCranSart = OrientedImageType::New() ;
    OrientedImageType::Pointer referenceVolume = OrientedImageType::New() ;
    OrientedImageType::Pointer segSemit = OrientedImageType::New() ;
    OrientedImageType::Pointer segRecFem = OrientedImageType::New() ;
    OrientedImageType::Pointer segCranSart = OrientedImageType::New() ;
    OrientedImageType::Pointer segVolume = OrientedImageType::New() ;
    DMDData::OrientedWriterType::Pointer                writer = DMDData::OrientedWriterType::New();
    writer->UseCompressionOn();
    OrientedImageType::RegionType                  region = input->GetLargestPossibleRegion() ;
    OrientedImageType::SpacingType                 spacing = input->GetSpacing() ;
    OrientedImageType::DirectionType               direction = input->GetDirection() ;
    OrientedImageType::PointType                   origin = input->GetOrigin() ;
    OrientedImageType::SizeType                    size = input->GetLargestPossibleRegion().GetSize();
    OrientedImageType::IndexType                   indexInput, indexOutput;
    float interpolatedVal = 0, upperVal = 0, bottomVal = 0, pixelVal = 0; 
    float distance = 5.0, upperIndex = 0, bottomIndex = 0;
    float upperFactor = 0, bottomFactor = 0, interpolatThresh = 0.5, inplane = 0, upperSum = 0, bottomSum = 0;
    float currentLabel = 0; 
    std::string outputfilename;

    std::cout << "size of the image:  " <<  size[0] << "  " << size[1] << "  " << size[2] << std::endl;

    output->SetRegions( region );
    output->SetSpacing( spacing );
    output->SetDirection( direction );
    output->SetOrigin( origin );                    
    output->Allocate();
    referenceVolume->SetRegions( region );
    referenceVolume->SetSpacing( spacing );
    referenceVolume->SetDirection( direction );
    referenceVolume->SetOrigin( origin );                    
    referenceVolume->Allocate();

    referenceSemit->SetRegions( region );
    referenceSemit->SetSpacing( spacing );
    referenceSemit->SetDirection( direction );
    referenceSemit->SetOrigin( origin );                    
    referenceSemit->Allocate();
    referenceRecFem->SetRegions( region );
    referenceRecFem->SetSpacing( spacing );
    referenceRecFem->SetDirection( direction );
    referenceRecFem->SetOrigin( origin );                    
    referenceRecFem->Allocate();
    referenceCranSart->SetRegions( region );
    referenceCranSart->SetSpacing( spacing );
    referenceCranSart->SetDirection( direction );
    referenceCranSart->SetOrigin( origin );                    
    referenceCranSart->Allocate();
    segVolume->SetRegions( region );
    segVolume->SetSpacing( spacing );
    segVolume->SetDirection( direction );
    segVolume->SetOrigin( origin );                    
    segVolume->Allocate();
    segSemit->SetRegions( region );
    segSemit->SetSpacing( spacing );
    segSemit->SetDirection( direction );
    segSemit->SetOrigin( origin );                    
    segSemit->Allocate();
    segRecFem->SetRegions( region );
    segRecFem->SetSpacing( spacing );
    segRecFem->SetDirection( direction );
    segRecFem->SetOrigin( origin );                    
    segRecFem->Allocate();
    segCranSart->SetRegions( region );
    segCranSart->SetSpacing( spacing );
    segCranSart->SetDirection( direction );
    segCranSart->SetOrigin( origin );                    
    segCranSart->Allocate();

    for(int l = CRAN_SART; l < CAUD_SART; l++){
        currentLabel = l;
        //search the range of slices of a muscle
        int maxSlice = 0, minSlice = size[1], roundRange = 0;
        bool jumpout = 0;
        for(indexInput[1] = 0; indexInput[1] < (int)size[1]; indexInput[1]++) {
            for(indexInput[2] = 0; indexInput[2] < (int)size[2]; indexInput[2]++) {
                for(indexInput[0] = 0; indexInput[0] < (int)size[0]; indexInput[0]++) {
                    pixelVal = input->GetPixel(indexInput);
                    if(pixelVal == currentLabel){
                        if(indexInput[1] < minSlice)
                            minSlice = indexInput[1];
                        if(indexInput[1] > maxSlice)
                            maxSlice = indexInput[1];
                        jumpout = 1;
                        break;
                    }
                }
                if(jumpout)
                    break;
            }
            if(jumpout)
                jumpout = 0;
        }
        roundRange = (maxSlice - minSlice + 1);
        roundRange = (int)(floor(roundRange / distance) * distance); //rounded range of slices
        std::cout << "the slice range of muscle " << currentLabel << " is " << minSlice << " - " << maxSlice << "  " << roundRange << std::endl;
        for(indexOutput[1] = 0; indexOutput[1] < roundRange; indexOutput[1]++) { // start from the 3rd slice
            for(indexOutput[2] = 1; indexOutput[2] < (int)size[2] - 1; indexOutput[2]++) {
                for(indexOutput[0] = 1; indexOutput[0] < (int)size[0] - 1; indexOutput[0]++) {
                    int k = 3;
                    inplane = 0; upperSum = 0; bottomSum = 0;
                    for (int i = -k / 2 ; i <= k / 2; i++){
                        for (int j = -k / 2 ; j <= k / 2; j++){
                            indexInput[0] = indexOutput[0] + i; 
                            indexInput[2] = indexOutput[2] + j; 
                            upperIndex = (indexOutput[1] / (int)distance + 1) * distance + minSlice;
                            indexInput[1] = (int)upperIndex;
                            upperVal = input->GetPixel(indexInput);
                            if(upperVal == currentLabel)
                                upperVal = input->GetPixel(indexInput);
                            else
                                upperVal = 0;
                            bottomIndex = indexOutput[1] / (int)distance * distance + minSlice;
                            indexInput[1] = (int)bottomIndex;
                            bottomVal = input->GetPixel(indexInput);
                            if(bottomVal == currentLabel){
                                bottomVal = input->GetPixel(indexInput);
                            }
                            else
                                bottomVal = 0;
                            upperFactor = ((float)(indexOutput[1] + minSlice) - bottomIndex) / distance;
                            bottomFactor =  (upperIndex - (float)(indexOutput[1] + minSlice)) / distance;
                            inplane += (upperVal * upperFactor + bottomVal * bottomFactor);
                        }
                    }
                    interpolatedVal = inplane /(k * k);
                    if(interpolatedVal >= (interpolatThresh * currentLabel))
                        interpolatedVal = currentLabel;
                    else
                        interpolatedVal = 0;
                    indexOutput[1] += minSlice;
                    if(output->GetPixel(indexOutput) == 0){
                        output->SetPixel(indexOutput, interpolatedVal);
                    }
                    indexOutput[1] -= minSlice;
                }
            }
        }
    }
    //output interpolated muscle volume
    writer->SetInput( output );  
    std::stringstream fs;
    fs << distance;
    //outputfilename = "../data/" + caseID + "_intervolfull_left_femur_mask" + ".nrrd";
    outputfilename = "../data/" + caseID + "_intervolfull_every" + fs.str() + ".nrrd";
    dataWriter(writer, outputfilename);

    float overlapRate = 0;
    //calculate overlap of semitendinaus
    overlapRate = overlap(input, output, SEMIT, 1, distance);    
    std::cout << "the overlap of semit is:  " << overlapRate << std::endl;
    //calculate overlap of rec fem
    overlapRate = overlap(input, output, REC_FEM, 1, distance);    
    std::cout << "the overlap of rec fem is:  " << overlapRate << std::endl;
    //calculate overlap of cran sart
    overlapRate = overlap(input, output, CRAN_SART, 1, distance);    
    std::cout << "the overlap of cran sart is:  " << overlapRate << std::endl;

    for(indexOutput[1] = 1; indexOutput[1] < (int)size[1] - 1; indexOutput[1]++) {
        for(indexOutput[2] = 1; indexOutput[2] < (int)size[2] - 1; indexOutput[2]++) {
            for(indexOutput[0] = 1; indexOutput[0] < (int)size[0] - 1; indexOutput[0]++) {
                if(indexOutput[1] % (int)distance != 0){
                    segVolume->SetPixel(indexOutput, output->GetPixel(indexOutput));
                    referenceVolume->SetPixel(indexOutput, input->GetPixel(indexOutput));
                    if(output->GetPixel(indexOutput) == SEMIT)
                        segSemit->SetPixel(indexOutput, SEMIT);    
                    if(input->GetPixel(indexOutput) == SEMIT)
                        referenceSemit->SetPixel(indexOutput, SEMIT);
                    if(output->GetPixel(indexOutput) == REC_FEM)
                        segRecFem->SetPixel(indexOutput, REC_FEM);    
                    if(input->GetPixel(indexOutput) == REC_FEM)
                        referenceRecFem->SetPixel(indexOutput, REC_FEM);
                    if(output->GetPixel(indexOutput) == CRAN_SART)
                        segCranSart->SetPixel(indexOutput, CRAN_SART);    
                    if(input->GetPixel(indexOutput) == CRAN_SART)
                        referenceCranSart->SetPixel(indexOutput, CRAN_SART);
                }
                else { // remove the manually segmented slices
                   // reference->SetPixel(indexOutput, input->GetPixel(indexOutput));
                    segVolume->SetPixel(indexOutput, output->GetPixel(indexOutput));
                    referenceVolume->SetPixel(indexOutput, input->GetPixel(indexOutput));
                }
            }
        }
    }         
}
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float DMDData::overlap ( OrientedImageType::Pointer input, OrientedImageType::Pointer output, int label, bool exception, float interval)
{
    // exception: TRUE: the overlap rate may except some slices; FALSE: otherwise
    // interval: indicate the slices should be excepted. for example, exclude every third slices, then interval = 3
    OrientedImageType::SizeType                    size = input->GetLargestPossibleRegion().GetSize();
    OrientedImageType::IndexType                   index;
    float referVal = 0, segVal = 0, interVol = 0, unionVol = 0, referVol = 0, segVol = 0, numberofoutside = 0;

    for(index[2] = 1; index[2] < (int)size[2] - 1; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                if(exception){
                    if( index[2] % (int)interval != 0){
                        referVal = input->GetPixel(index);
                        segVal = output->GetPixel(index);
                        if((referVal == label) && (segVal == label))
                            interVol++;
                        if((referVal == label) || (segVal == label))
                            unionVol++;
                        if((referVal != label) && (segVal == label))
                            numberofoutside++;
                        if((referVal == label))
                            referVol++;
                        if((segVal == label))
                            segVol++;
                    }                
                }
                else{
                    referVal = input->GetPixel(index);
                    segVal = output->GetPixel(index);
                    if((referVal == label) && (segVal == label))
                        interVol++;
                    if((referVal == label) || (segVal == label))
                        unionVol++;
                    if((referVal == label))
                        referVol++;
                    if((segVal == label))
                        segVol++;
                }
            }
        }
    }
   
    std::cout << "number of pixels outside reference standard:  " << numberofoutside << "  " << 100 * numberofoutside / segVol << std::endl;
    return interVol / unionVol * 100;  // percentage
   
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::connectedLabeling2D ( const OrientedImage2DType::Pointer input)
{
   
 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::connectedLabeling3D ( const OrientedImageType::Pointer input)
{
    typedef itk::ImageFileWriter< uncharOrientedImageType >                  WriterType;

    // Labeller
    WriterType::Pointer   writer = WriterType::New();
    OrientedImageType::RegionType                  region = input->GetLargestPossibleRegion() ;
    OrientedImageType::SpacingType                 spacing = input->GetSpacing() ;
    OrientedImageType::DirectionType               direction = input->GetDirection() ;
    OrientedImageType::PointType                   origin = input->GetOrigin() ;
    OrientedImageType::SizeType                    size = input->GetLargestPossibleRegion().GetSize();
    OrientedImageType::IndexType                   index;
    OrientedImageType::Pointer                     morphOutput = OrientedImageType::New() ;
    uncharOrientedImageType::Pointer               output = uncharOrientedImageType::New() ;
    
    morphOutput->SetRegions( region );
    morphOutput->SetSpacing( spacing );
    morphOutput->SetDirection( direction );
    morphOutput->SetOrigin( origin );                    
    morphOutput->Allocate();

    output->SetRegions( region );
    output->SetSpacing( spacing );
    output->SetDirection( direction );
    output->SetOrigin( origin );                    
    output->Allocate();

    morphErod ( input, morphOutput, 2 ) ;
    morphDilat ( morphOutput, morphOutput, 2 ) ;

    // initialized the volume data
    for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                output->SetPixel(index, (unsigned char)morphOutput->GetPixel(index));
            }
        }
    }     

//    itkMetaImageIOType::Pointer		labellerMetaWriterIO  = itkMetaImageIOType::New();

    labeller = itkConnectedComponentFilterType::New();
    labeller->SetInput( output );
    labeller->Update();
   
    relabeller = itkRelabelComponentFilterType::New();
    relabeller->SetInput( labeller->GetOutput() );
    relabeller->Update();

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::waterTubeSegmentation (const OrientedImageType::Pointer input )
{
    typedef itk::ImageFileWriter< uncharOrientedImageType >                  WriterType;
    WriterType::Pointer   uncharwriter = WriterType::New();

    OrientedImageType::RegionType                  region = input->GetLargestPossibleRegion() ;
    typedef itk::OtsuThresholdImageFilter<OrientedImageType, OrientedImageType >    FilterType;
    OrientedWriterType::Pointer                                                     writer = OrientedWriterType::New();
    FilterType::Pointer                                                             filter = FilterType::New();
    OrientedImageType::IndexType                   index;
    OrientedImageType::SizeType                    size = input->GetLargestPossibleRegion().GetSize();
    float actualVolWaterTube = 50000, diffVolActualLabel = 50000;
    int selectedLabel = 0; 
    
    filter->SetInput( input );
    filter->SetOutsideValue( 1 );
    filter->SetInsideValue( 0 );
    filter->SetNumberOfHistogramBins( 1024 );
    filter->Update();
    writer->SetInput(filter->GetOutput());
    dataWriter(writer, "../data/otsuThresholdOutput.nrrd");
    
    connectedLabeling3D ( filter-> GetOutput());
 
    unsigned long nObjects = relabeller->GetNumberOfObjects();
    const std::vector<unsigned long> sizeOfObjects = relabeller->GetSizeOfObjectsInPixels();
    std::cout << "Number of objects : " << nObjects << std::endl;
    std::vector<unsigned long>::const_iterator it;
    //remove small objectes
    int smallObjThresh;
    int i;
    for( i = 0, it = sizeOfObjects.begin(); it != sizeOfObjects.end(); ++i, ++it){
        std::cout << "Size of object " << i + 1 << ": " << (*it) << " pixels" << std::endl;
        if(*it < 1000){
            smallObjThresh = i + 1;
            break;
        }
    }

    std::cout << "the number of objects remained:  " << smallObjThresh << std::endl;
    float cylinderVolume[(const int)nObjects], cylinderVolumeTheory[(const int)nObjects], cylinderRadius[(const int)nObjects], cylinderLength[(const int)nObjects];
    for (int i = 0; i < nObjects; i++) {
        cylinderVolume[i] = 0;
        cylinderVolumeTheory[i] = 0;
        cylinderRadius[i] = 0;
        cylinderLength[i] = 0;
    }
     
    for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        float radius[(const int)nObjects]; 
        for (int i = 0; i < nObjects; i++) {
            radius[i] = 0; }
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                int pixelLabel = relabeller->GetOutput()->GetPixel(index);
                if(pixelLabel <= smallObjThresh){ 
                    labeller->GetOutput()->SetPixel(index, pixelLabel); 
                }
                else{
                    labeller->GetOutput()->SetPixel(index, 0); 
                }
            }
        }
    }
   

    std::cout << "small object threshold: " <<  smallObjThresh << std::endl;
    typedef itk::Image< uncharPixelType, 3 >   LabelImageType;
    typedef itk::LabelGeometryImageFilter< LabelImageType > LabelGeometryType;
    LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();
    labelGeometryFilter->SetInput( labeller->GetOutput() );
    typedef LabelGeometryType::VectorType   VectorType;
  // These generate optional outputs.
    labelGeometryFilter->CalculatePixelIndicesOn();
    labelGeometryFilter->CalculateOrientedBoundingBoxOn();;
    labelGeometryFilter->CalculateOrientedLabelRegionsOn();
    labelGeometryFilter->Update();
    std::cout << "\n\nRUNNING THE FILTER WITHOUT AN INTENSITY IMAGE..." << std::endl;
    labelGeometryFilter->Print(std::cout);

    for( int i = 1; i <= smallObjThresh; i++) { 
        float axelength[3], volume = labelGeometryFilter->GetVolume(i), rater = 0 ; 
        const VectorType EigenValue = labelGeometryFilter->GetEigenvalues( i );
        VectorType::const_iterator itr = EigenValue.begin();

        int j = 0;
        while( itr != EigenValue.end() ){
            std::cout << "Index = " << *itr << "  " << 4 * sqrt(*itr) << std::endl;
            axelength[j] = 4 * sqrt(*itr);
            ++itr;
            j++;
        }
       // rater = axelength[0] / 2 * axelength[1] / 2 * 3.1415926 * axelength[2] / volume;
        rater = volume / actualVolWaterTube;
        std::cout << "score for lable " << i << " is " << rater << "  volume:  " << volume << std::endl;
        if (diffVolActualLabel > fabs(volume - actualVolWaterTube)) {
            diffVolActualLabel = fabs(volume - actualVolWaterTube);
            selectedLabel = i;
        }
    }
    std::cout << "selected label is " << selectedLabel << "   difference in volume is " << diffVolActualLabel << std::endl;
    std::cout << "labeling finished! " << std::endl;
    uncharwriter->SetInput( labeller->GetOutput());
    uncharwriter->SetFileName("../data/labeledImage.nrrd");
    try {
        uncharwriter->Update();
    }
    catch (itk::ExceptionObject &ex){
        std::cout << "write:" << ex << std::endl;
        //return EXIT_FAILURE;
        exit(0);
    }
    std::cout << "output labeld image" << std::endl;

    for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                int pixelLabel = labeller->GetOutput()->GetPixel(index);
                if(pixelLabel == selectedLabel){ 
                    labeller->GetOutput()->SetPixel(index, 1); 
                }
                else{
                    labeller->GetOutput()->SetPixel(index, 0); 
                }
            }
        }
    }
    uncharwriter->SetInput( labeller->GetOutput());
    uncharwriter->SetFileName("../data/selectedlabeledImage.nrrd");
    try {
        uncharwriter->Update();
    }
    catch (itk::ExceptionObject &ex){
        std::cout << "write:" << ex << std::endl;
        //return EXIT_FAILURE;
        exit(0);
    }   
    std::cout << "output selected label" << std::endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::fatSuppressBiasIdentify (const OrientedImageType::Pointer inputT2, OrientedImageType::Pointer inputT2FS)
{
    float                                          meanIntT2 = 0, meanIntT2FS = 0, volume = 0;
    OrientedImageType::IndexType                   index;
    OrientedImageType::SizeType                    size = inputT2->GetLargestPossibleRegion().GetSize();
    OrientedWriterType::Pointer                    writer = OrientedWriterType::New();
   
    for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                int pixelLabel = labeller->GetOutput()->GetPixel(index);
                if(pixelLabel == 1){ 
                    meanIntT2 += inputT2->GetPixel(index);
                    meanIntT2FS += inputT2FS->GetPixel(index);
                    volume++;
                }
            }
        }
    }
   
    meanIntT2 = meanIntT2 / volume;
    meanIntT2FS = meanIntT2FS / volume;
    std::cout << "mean intensity is  " << meanIntT2 << "  " << meanIntT2FS << std::endl;
    if( meanIntT2FS < meanIntT2) {
        std::cout << "t2FS correction is needed" << std::endl;
        for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
            for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
                for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                    float correctPixel = inputT2FS->GetPixel(index) * meanIntT2 / meanIntT2FS ; 
                    inputT2FS->SetPixel(index, 1); 
                }
            }
        }
        writer->SetInput(inputT2FS);
        dataWriter(writer, "../data/correctedT2FS.nrrd");
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DMDData::histogramMatchinMask (const OrientedImageType::Pointer inputMove, OrientedImageType::Pointer inputRef, OrientedImageType::Pointer output, int bin, int point)
{
    OrientedImageType::IndexType                   index;
    OrientedImageType::SizeType                    size = inputMove->GetLargestPossibleRegion().GetSize();
    
    
    // build histogram
    for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
           
            }
        }
    }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int  DMDData::connectedComponentRegionGrowing (float lower, float upper, OrientedImageType::Pointer & input, OrientedImageType::IndexType index)
{
    int componentVolume = 0;

    StrITKType outputFilename = "../data/regionGrowing.nrrd"    ;
    OrientedWriterType::Pointer bwriter = OrientedWriterType::New();

    bwriter->SetFileName(outputFilename);
    typedef itk::ConnectedThresholdImageFilter< OrientedImageType, OrientedImageType > ConnectedFilterType;
    ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
    ConnectedFilterType::ConnectivityEnumType connectivity = ConnectedFilterType::FaceConnectivity;
    connectedThreshold->SetInput( input );
    bwriter->SetInput( connectedThreshold->GetOutput() );
    connectedThreshold->SetLower( lower );
    connectedThreshold->SetUpper( upper );
    /*connectivity = ConnectedFilterType::FaceConnectivity;
    connectedThreshold->SetConnectivity( connectivity );
    std::cout << connectedThreshold->GetConnectivity() << std::endl;
*/
    connectivity = ConnectedFilterType::FullConnectivity; 
    connectedThreshold->SetConnectivity( connectivity );
 //   std::cout << connectedThreshold->GetConnectivity() << std::endl;
    connectedThreshold->SetReplaceValue( TEMP_PIXEL_VALUE );
    connectedThreshold->SetSeed( index );
    try
    {
        bwriter->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
    }

    OrientedImageType::SizeType   size = connectedThreshold->GetOutput()->GetLargestPossibleRegion().GetSize();
    OrientedImageType::IndexType  tmpindex;
    //calculate the volume of the component
    for(tmpindex[2] = 0 ; tmpindex[2] < size[2]; tmpindex[2]++) {
        for(tmpindex[1] = 0; tmpindex[1] < size[1]; tmpindex[1]++) {
            for(tmpindex[0] = 0; tmpindex[0] < size[0]; tmpindex[0]++) {
                if (connectedThreshold->GetOutput()->GetPixel(tmpindex) == TEMP_PIXEL_VALUE){
                    componentVolume++;
                    input->SetPixel(tmpindex, SEARCHED_PIXEL_VALUE);
                }
            }
        }
    }

    return componentVolume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int  DMDData::connectedComponentLabeling (OrientedImageType::Pointer & input)
{
    //typedef itk::ImageFileWriter< uncharOrientedImageType >                  WriterType;
    typedef itk::ImageFileWriter< OrientedImageType >                  WriterType;
    // Labeller
    WriterType::Pointer   writer = WriterType::New();
    OrientedImageType::RegionType                  region = input->GetLargestPossibleRegion() ;
    OrientedImageType::SpacingType                 spacing = input->GetSpacing() ;
    OrientedImageType::DirectionType               direction = input->GetDirection() ;
    OrientedImageType::PointType                   origin = input->GetOrigin() ;
    OrientedImageType::SizeType                    size = input->GetLargestPossibleRegion().GetSize();
    OrientedImageType::IndexType                   index;
    OrientedImageType::Pointer                     morphOutput = OrientedImageType::New() ;
    uncharOrientedImageType::Pointer               output = uncharOrientedImageType::New() ;

    output->SetRegions( region );
    output->SetSpacing( spacing );
    output->SetDirection( direction );
    output->SetOrigin( origin );                    
    output->Allocate();

    // initialized the volume data
    for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                output->SetPixel(index, (unsigned char)input->GetPixel(index));
            }
        }
    }   
    //    itkMetaImageIOType::Pointer		labellerMetaWriterIO  = itkMetaImageIOType::New();
    labeller = itkConnectedComponentFilterType::New();
    labeller->SetInput( output );
    labeller->Update();

    for(index[2] = 0; index[2] < (int)size[2]; index[2]++) {
        for(index[1] = 0; index[1] < (int)size[1]; index[1]++) {
            for(index[0] = 0; index[0] < (int)size[0]; index[0]++) {
                input->SetPixel(index, labeller->GetOutput()->GetPixel(index));
            }
        }
    }
//    rescaleFilter->SetOutputMinimum(0);
//    rescaleFilter->SetOutputMaximum(255);
//    rescaleFilter->SetInput(labelFilter->GetOutput());

//    relabeller = itkRelabelComponentFilterType::New();
 //   relabeller->SetInput( labeller->GetOutput() );
 //   relabeller->SetOutputMinimum(0);
  //  relabeller->SetOutputMaximum(1024);
  //  relabeller->Update();
//    std::cout << labeller->GetObjectCount() << std::endl;
    writer->SetInput(input);
    writer->SetFileName("../data/labelling_o.nrrd");
    writer->Update();

    return labeller->GetObjectCount();
}
#endif
