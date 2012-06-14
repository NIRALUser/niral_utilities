/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/ResampleDTI.cxx $
  Language:  C++
  Date:      $Date: 2010/08/06 16:23:18 $
  Version:   $Revision: 1.12 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/



#include "itkDiffusionTensor3DResample.h"
#include "itkDiffusionTensor3DRigidTransform.h"
#include "itkDiffusionTensor3DFSAffineTransform.h"
#include "itkDiffusionTensor3DPPDAffineTransform.h"
#include "itkDiffusionTensor3DNonRigidTransform.h"
#include "itkDiffusionTensor3DNearestNeighborInterpolateFunction.h"
#include "itkDiffusionTensor3DLinearInterpolateFunction.h"
#include "itkDiffusionTensor3DWindowedSincInterpolateImageFunction.h"
#include "itkDiffusionTensor3DBSplineInterpolateImageFunction.h"
//#include "itkDiffusionTensor3DTransform.h"
#include <itkPeriodicBoundaryCondition.h>
#include <itkMatrix.h>
#include <itkPoint.h>
#include <itkImageFileReader.h>
#include <itkImageIOBase.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkTransformFileReader.h>
#include "concatenateHFieldCLP.h"
#include <itkDiffusionTensor3D.h>
#include "itkDiffusionTensor3DLogImageFilter.h"
#include "itkDiffusionTensor3DExpImageFilter.h"
#include "itkDiffusionTensor3DZeroCorrection.h"
#include "itkDiffusionTensor3DAbsCorrection.h"
#include "itkDiffusionTensor3DNearestCorrection.h"
#include "itkDiffusionTensor3DRead.h"
#include "itkDiffusionTensor3DWrite.h"
#include "dtiprocessFiles/deformationfieldio.h"
#include "itkWarpTransform3D.h"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <list>
#include <itkVectorLinearInterpolateImageFunction.h>
#include "itkTransformDeformationFieldFilter.h"
#include "vnl/vnl_math.h"
#include <itkVectorResampleImageFilter.h>
//#include <itkBSplineDeformableTransform.h>
#include "dtiprocessFiles/itkDeformationFieldToHFieldImageFilter.h"
#include <vector>
#include <string>
// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace {

#define RADIUS 3



//Structure of the command lign parameters
struct parameters
{
  std::vector<std::string> deffieldtoaverage;
  std::vector<std::string> typeOfFieldtoaverage;
  std::string referenceVolume ;
  std::string deffield1 ;
  std::string deffield2 ;
  std::string typeOfField1 ;
  std::string typeOfField2 ;
  std::string outputHfield ;
  std::string outputHfieldAvg;
};



//resamples field to reference image size; local filter so that the memory is freed once it has run
void ResampleDeformationField( DeformationImageType::Pointer &field ,
                               const itk::Point< double , 3 > &origin ,
                               const itk::Vector< double , 3 > &spacing ,
                               const itk::Size< 3 > &size ,
                               const itk::Matrix< double , 3 , 3 > &direction
                             )
{
  //Check if the field does not already have the same properties as the output image:
  //It would save some time if we did not have to resample the field
  itk::Point< double , 3 > fieldOrigin ;
  itk::Vector< double , 3 > fieldSpacing ;
  itk::Size< 3 > fieldSize ;
  itk::Matrix< double , 3 , 3 > fieldDirection ;
  fieldOrigin = field->GetOrigin() ;
  fieldSpacing = field->GetSpacing() ;
  fieldSize = field->GetLargestPossibleRegion().GetSize() ;
  fieldDirection = field->GetDirection() ;
  if(  fieldSize == size
    && fieldSpacing == spacing
    && fieldDirection == direction
    && fieldOrigin == origin
    )
  {
    return ;
  }
  typedef itk::VectorLinearInterpolateImageFunction< DeformationImageType > VectorInterpolatorType ;
  VectorInterpolatorType::Pointer linearVectorInterpolator = VectorInterpolatorType::New() ;
  typedef itk::VectorResampleImageFilter< DeformationImageType ,
                                          DeformationImageType ,
                                          double
                                        > ResampleImageFilter ;
  ResampleImageFilter::Pointer resampleFieldFilter = ResampleImageFilter::New() ;
  DeformationPixelType defaultPixel ;
  defaultPixel.Fill( 0.0 ) ;
  resampleFieldFilter->SetDefaultPixelValue( defaultPixel ) ;
  resampleFieldFilter->SetInput( field ) ;
  resampleFieldFilter->SetInterpolator( linearVectorInterpolator ) ;
  resampleFieldFilter->SetOutputDirection( direction ) ;
  resampleFieldFilter->SetSize( size ) ;
  resampleFieldFilter->SetOutputSpacing( spacing ) ;
  resampleFieldFilter->SetOutputOrigin( origin ) ;
  resampleFieldFilter->Update() ;
  field = resampleFieldFilter->GetOutput() ;
}


template< class PixelType >
int Do( parameters list )
{
  
    //read the reference volume 
    typedef itk::Image< unsigned char , 3 > ImageType ;
    typedef itk::ImageFileReader< ImageType > ReaderType ;
    typedef typename ReaderType::Pointer ReaderTypePointer ;
    ReaderTypePointer readerReference ;
    readerReference = ReaderType::New() ;
    readerReference->SetFileName( list.referenceVolume.c_str() ) ;
    readerReference->UpdateOutputInformation() ;

      //set the size, origin, spacing,direction of the reference volume
      itk::Point< double , 3 > m_Origin ;
      itk::Vector< double , 3 > m_Spacing ;
      itk::Size< 3 > m_Size ;
      itk::Matrix< double , 3 , 3 > m_Direction;
       m_Spacing = readerReference->GetOutput()->GetSpacing() ;
       m_Size= readerReference->GetOutput()->GetLargestPossibleRegion().GetSize() ;
       m_Origin = readerReference->GetOutput()->GetOrigin() ;
       m_Direction = readerReference->GetOutput()->GetDirection() ;

    typename DeformationImageType::Pointer fieldPointer1 ;
    typename DeformationImageType::Pointer fieldPointer2 ;

       //set if the field is a displacement or a H- field
       DeformationFieldType dftype1 = HField ;
       DeformationFieldType dftype2 = HField ;
       if( !list.typeOfField1.compare( "displacement" ) )
       {
         dftype1 = Displacement ;
       }
	if( !list.typeOfField2.compare( "displacement" ) )
       {
         dftype2 = Displacement ;
       }
       //reads deformation field and if it is a h-field, it transforms it to a displacement field
       fieldPointer1 = readDeformationField( list.deffield1 , dftype1 ) ;
       fieldPointer2 = readDeformationField( list.deffield2 , dftype2 ) ;

      //Create warp transform
      typedef itk::WarpTransform3D< double > WarpTransformType ;
      
      typename WarpTransformType::Pointer warpTransform2 = WarpTransformType::New() ;

        //Resample the deformation field so that it has the same properties as the reference image 
      ResampleDeformationField( fieldPointer1 ,m_Origin , m_Spacing ,m_Size ,m_Direction) ;
        //Resample the deformation field so that it has the same properties as the reference image 
      ResampleDeformationField( fieldPointer2 ,m_Origin ,m_Spacing , m_Size ,m_Direction) ;

      //Compute the transformation field adding all the deformation field together
	warpTransform2->SetDeformationField( fieldPointer2 ) ;
        typedef itk::TransformDeformationFieldFilter< double , double , 3 > itkTransformDeformationFieldFilterType ;
        typename itkTransformDeformationFieldFilterType::Pointer transformDeformationFieldFilter = itkTransformDeformationFieldFilterType::New() ;

        transformDeformationFieldFilter->SetInput( fieldPointer1 ) ;
        transformDeformationFieldFilter->SetTransform(  warpTransform2.GetPointer()  ) ;
        transformDeformationFieldFilter->Update() ;
        fieldPointer1 = transformDeformationFieldFilter->GetOutput() ;


      //Save HField
        typedef itk::ImageFileWriter< DeformationImageType > HFieldWriterType ;
        typedef typename HFieldWriterType::Pointer HFieldWriterTypePointer ;
        HFieldWriterTypePointer hfwriter = HFieldWriterType::New() ; 

        hfwriter->SetInput( fieldPointer1 ) ;
        hfwriter->SetFileName( list.outputHfield ) ;
        try
        {
          hfwriter->Update() ;
        }
        catch( itk::ExceptionObject & Except )
        {
          std::cerr << "Writing output HField: Exception caught!"
                    << std::endl ;
          std::cerr << Except << std::endl ;
          return EXIT_FAILURE ;
        }

    
    return EXIT_SUCCESS ;
}

} // end of anonymous namespace

template< class PixelType >
int averageDeformationField(parameters list )
{
//read the reference volume 
    typedef itk::Image< unsigned char , 3 > ImageType ;
    typedef itk::ImageFileReader< ImageType > ReaderType ;
    typedef typename ReaderType::Pointer ReaderTypePointer ;
    ReaderTypePointer readerReference ;
    readerReference = ReaderType::New() ;
    readerReference->SetFileName( list.referenceVolume.c_str() ) ;
    readerReference->UpdateOutputInformation() ;

      //set the size, origin, spacing,direction of the reference volume
      itk::Point< double , 3 > m_Origin ;
      itk::Vector< double , 3 > m_Spacing ;
      itk::Size< 3 > m_Size ;
      itk::Matrix< double , 3 , 3 > m_Direction;
       m_Spacing = readerReference->GetOutput()->GetSpacing() ;
       m_Size= readerReference->GetOutput()->GetLargestPossibleRegion().GetSize() ;
       m_Origin = readerReference->GetOutput()->GetOrigin() ;
       m_Direction = readerReference->GetOutput()->GetDirection() ;

       int n=list.deffieldtoaverage.size();
	typename DeformationImageType::Pointer fieldPointer1 ;
	DeformationFieldType dftype1 = HField ;
	if( !list.typeOfField1.compare( "displacement" ) )
		{
			dftype1 = Displacement ;
		}
	fieldPointer1 = readDeformationField( list.deffield1 , dftype1 ) ;
	//Resample the deformation field so that it has the same properties as the reference image 
      	ResampleDeformationField( fieldPointer1 ,m_Origin , m_Spacing ,m_Size ,m_Direction) ;
	typedef itk::ImageRegionIterator< DeformationImageType > DeformationConstIteratorType;
	DeformationConstIteratorType it1 (fieldPointer1,fieldPointer1->GetLargestPossibleRegion());
	for(int i=0;i<n;i++)
	{
		
		typename DeformationImageType::Pointer fieldPointer2 ;
		//set if the field is a displacement or a H- field
		DeformationFieldType dftype2 = HField ;
		std::cout<<"verification of the "<<i<<"eme type of DFiled : "<<list.typeOfFieldtoaverage[i]<<std::endl;
		if( !list.typeOfFieldtoaverage[i].compare( "displacement" ) )
		{
			dftype2 = Displacement ;
		}
		//reads deformation field and if it is a h-field, it transforms it to a displacement field
		fieldPointer2 = readDeformationField( list.deffieldtoaverage[i] , dftype2 ) ;
	 	
        	//Resample the deformation field so that it has the same properties as the reference image 
      		ResampleDeformationField( fieldPointer2 ,m_Origin ,m_Spacing , m_Size ,m_Direction) ;
		
		it1.GoToBegin();
		DeformationConstIteratorType it2 (fieldPointer2,fieldPointer2->GetLargestPossibleRegion());
		it2.GoToBegin();
		itk::Vector<double,3> sum;sum[0]=0;sum[1]=0;sum[2]=0;
		while(!it1.IsAtEnd() && !it2.IsAtEnd())
		{
			sum+=(it1.Get())+(it2.Get());	
			it1.Set(sum);
			sum[0]=0;sum[1]=0;sum[2]=0;
			++it1;++it2;
		}

	}
        //DeformationConstIteratorType it1 (fieldPointer1,fieldPointer1->GetLargestPossibleRegion());
	it1.GoToBegin();
	itk::Vector<double,3> temp;
	while(!it1.IsAtEnd())
        {
		temp=(it1.Get())/(n+1);
		it1.Set(temp);
		++it1;
	}
	//Save HField
        typedef itk::ImageFileWriter< DeformationImageType > HFieldWriterType ;
        typedef typename HFieldWriterType::Pointer HFieldWriterTypePointer ;
        HFieldWriterTypePointer hfwriterAvg = HFieldWriterType::New() ; 

        hfwriterAvg->SetInput( fieldPointer1 ) ;
        hfwriterAvg->SetFileName( list.outputHfieldAvg ) ;
        try
        {
          hfwriterAvg->Update() ;
        }
        catch( itk::ExceptionObject & Except )
        {
          std::cerr << "Writing output HField: Exception caught!"
                    << std::endl ;
          std::cerr << Except << std::endl ;
          return EXIT_FAILURE ;
        }

return EXIT_SUCCESS ;
}


int main( int argc , char * argv[] )
{
  PARSE_ARGS ;
  parameters list ;
  list.deffieldtoaverage=deffieldtoaverage;
  list.typeOfFieldtoaverage=typeOfFieldtoaverage;
  list.referenceVolume = referenceVolume ;
  list.deffield1 = deffield1 ;
  list.deffield2 = deffield2 ;
  list.typeOfField1 = typeOfField1 ;
  list.typeOfField2 = typeOfField2 ;
  list.outputHfield = outputHfield ;
  list.outputHfieldAvg = outputHfieldAvg;
 // Do< float > ( list ) ;

  if (list.typeOfFieldtoaverage[0].compare( "" ) && list.deffieldtoaverage[0].compare( "" ) && list.outputHfieldAvg.compare( "" ))
  {
	averageDeformationField< float >(list ) ;
  }
	

  return EXIT_SUCCESS ;
}
