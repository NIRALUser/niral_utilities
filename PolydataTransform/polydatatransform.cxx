
#include <string>
#include <iostream>
#include <fstream>

// VTK includes
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkLine.h>
#include <vtkVersion.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

// ITK includes
#include <itkDiffusionTensor3D.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkVersion.h>

#include "fiberfileIO.hxx"
// #include "deformationfieldoperations.h"
#include "deformationfieldio.h"
#include "dtitypes.h"
#include "polydatatransformCLP.h"


int main(int argc, char* argv[])
{
  PARSE_ARGS;

  if( fiberFile.empty() )
  {
    std::cerr << "An input fiber file has to be specified" << std::endl;
    return EXIT_FAILURE;
  }
  if( fiberOutput.empty() )
  {
    std::cerr << "An output fiber file has to be specified" << std::endl;
    return EXIT_FAILURE;
  }

  if (hField.empty() && displacementField.empty()){
    /* Linear transforms */

    vtkSmartPointer<vtkPolyData> polydata;
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    int noTransformCount = 0;

    if(fiberFile.rfind(".vtk")!= std::string::npos || fiberFile.rfind(".vtp")!= std::string::npos)
    {
      polydata = readVTKFile(fiberFile.c_str());
    }

    else if(fiberFile.rfind(".fcsv")!= std::string::npos)
    {
      polydata = readFCSVFile(fiberFile.c_str());
    }
    else
    {
      throw itk::ExceptionObject("Unknown file format for input polydata must be .vtp, .vtk or .fcsv");
    }

    if ((translationVector[0] != 0) || (translationVector[1] != 0) || (translationVector[2] != 0)){
      transform->Translate(translationVector[0],translationVector[1],translationVector[2]);
    }
    else{noTransformCount++;}
    if((rotationVector[0] != 0) && ((rotationVector[1] != 0) || (rotationVector[2] != 0) || (rotationVector[3] != 0))){
      transform->RotateWXYZ(rotationVector[0],rotationVector[1],rotationVector[2],rotationVector[3]);
    }
    else{noTransformCount++;}
    if(rotationX != 0){
      transform->RotateX(rotationX);
    }
    else{noTransformCount++;}
    if(rotationY != 0){
      transform->RotateY(rotationY);
    }
    else{noTransformCount++;}
    if(rotationZ != 0){
      transform->RotateZ(rotationZ);
    }
    else{noTransformCount++;}
    if((scaling[0] != 1) || (scaling[1] != 1) || (scaling[2] != 1)){
      transform->Scale(scaling[0], scaling[1], scaling[2]);
    }
    else{noTransformCount++;}

    if (noTransformCount == 6){
      std::cout << "A transform has to be specified" << std::endl;
      return EXIT_FAILURE;
    }
    else{
      vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New(); 
      transformFilter->SetInputData(polydata); 
      transformFilter->SetTransform(transform); 
      transformFilter->Update();
      polydata = transformFilter->GetOutput();
      if(fiberOutput.rfind(".vtk")!= std::string::npos || fiberFile.rfind(".vtp")!= std::string::npos){
        writeVTKFile(fiberOutput.c_str(), polydata);
      }

      else if(fiberOutput.rfind(".fcsv")!= std::string::npos)
      {
        writeFCSVFile(fiberOutput.c_str(), polydata);
      }
      else
      {
        throw itk::ExceptionObject("Unknown file format for output polydata must be .vtp, .vtk or .fcsv");
        return EXIT_FAILURE;
      }
      std::cout<<"Done!"<<std::endl;
      return EXIT_SUCCESS;
    }
    
  }


  // vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  // std::cout << "Reading " << fiberFile<< std::endl;
  // reader->SetFileName(fiberFile.c_str());
  // reader->Update();

  // // Extract the polydata
  // vtkSmartPointer<vtkPolyData> polydata =
  //   reader->GetOutput();
  
  DeformationImageType::Pointer deformationfield;
  if( !hField.empty() )
  {
    deformationfield = readDeformationField(hField, HField);
  }
  else if( !displacementField.empty() )
  {
    std::cout << "1"<< std::endl;
    deformationfield = readDeformationField(displacementField, Displacement);
    std::cout << "2"<< std::endl;
  }
  else
  {
    std::cout << "3"<< std::endl;
    deformationfield = NULL;
  }
  typedef itk::VectorLinearInterpolateImageFunction<DeformationImageType, double> DeformationInterpolateType;
  DeformationInterpolateType::Pointer definterp;
  if(deformationfield)
  {
    std::cout << "4"<< std::endl;
    definterp = DeformationInterpolateType::New();
    definterp->SetInputImage(deformationfield);
  }
  else
  {
    std::cout << "5"<< std::endl;
    std::cout << "A deformation field has to be specified" << std::endl;
    return EXIT_FAILURE;
  }
  vtkPoints* inPts;
  vtkIdType numPoints;
  vtkPoints    * points = vtkPoints::New();
  vtkSmartPointer<vtkPolyData> polydata;
  if(fiberFile.rfind(".vtk")!= std::string::npos || fiberFile.rfind(".vtp")!= std::string::npos)
  {
    polydata = readVTKFile(fiberFile.c_str());
  }

  else if(fiberFile.rfind(".fcsv")!= std::string::npos)
  {
    polydata = readFCSVFile(fiberFile.c_str());
  }
  else
  {
    throw itk::ExceptionObject("Unknown file format for input fibers must be .vtp, .vtk or .fcsv");
  }


    typedef DeformationInterpolateType::ContinuousIndexType ContinuousIndexType;
    ContinuousIndexType ci, origci;

    // For each point along the fiber
    numPoints= polydata->GetNumberOfPoints();
    inPts = polydata->GetPoints();

    typedef DTIPointType::PointType PointType;
    PointType fiberpoint;
    double fiberpointtemp[3];

    for(int i=0;i<numPoints;i++)
    {	
        inPts->GetPoint( i, fiberpointtemp );
	if (invx) {
	  fiberpoint[0] = - fiberpointtemp[0];
	} else {
	  fiberpoint[0] = + fiberpointtemp[0];
	}
	if (invy) {
	  fiberpoint[1] = - fiberpointtemp[1];
	} else {
	  fiberpoint[1] = + fiberpointtemp[1];
	}
	if (invz) {
	  fiberpoint[2] = - fiberpointtemp[2];
	} else {
	  fiberpoint[2] = + fiberpointtemp[2];	
	}
        deformationfield->TransformPhysicalPointToContinuousIndex(fiberpoint,ci);
	if( !deformationfield->GetLargestPossibleRegion().IsInside( ci ) )
        {
          std::cerr << "Fiber is outside deformation field image. Deformation field has to be in the fiber space. Warning: Original position will be used" << std::endl ;
        }
        DeformationPixelType warp(definterp->EvaluateAtContinuousIndex(ci).GetDataPointer());
        for(unsigned int j =0; j < 3; j++)
        {
          fiberpoint[j] +=warp[j];
        }	//convert LPS to RAS (vtk)
	if (invx)
        {
	  fiberpoint[0] = - fiberpoint[0];
	}
        else
        {
	  fiberpoint[0] = + fiberpoint[0];
	}
	if (invy)
        {
	  fiberpoint[1] = - fiberpoint[1];
	}
        else
        {
	  fiberpoint[1] = + fiberpoint[1];
	}
	if (invz)
        {
	  fiberpoint[2] = - fiberpoint[2];
	}
        else
        {
	  fiberpoint[2] = + fiberpoint[2];	
	}

        points->InsertPoint(i,fiberpoint[0],fiberpoint[1],fiberpoint[2]);
	
    }
    polydata->SetPoints(points);

  if(fiberOutput.rfind(".vtk")!= std::string::npos || fiberFile.rfind(".vtp")!= std::string::npos)
  {
    writeVTKFile(fiberOutput.c_str(), polydata);
  }

  else if(fiberOutput.rfind(".fcsv")!= std::string::npos)
  {
    writeFCSVFile(fiberOutput.c_str(), polydata);
  }
  else
  {
    throw itk::ExceptionObject("Unknown file format for output fibers must be .vtp, .vtk or .fcsv");
  }
  std::cout<<"Done!"<<std::endl;
  std::cout<<"Number of fibers saved: "<<polydata->GetNumberOfLines()<<std::endl;

  // std::cout<<std::endl;
  // std::cout<<"Saving fibers...."<<std::endl;
  // std::cout<<fiberOutput<<std::endl;
  // vtkPolyDataWriter * fiberwriter = vtkPolyDataWriter::New();
  // fiberwriter->SetFileName(fiberOutput.c_str());
  // #if (VTK_MAJOR_VERSION < 6)
  // fiberwriter->SetInput(polydata);
  // #else
  // fiberwriter->SetInputData(polydata);
  // #endif
  // fiberwriter->SetFileTypeToBinary();
  // fiberwriter->Update();
  // try
  //   {
  //   fiberwriter->Write();
  //   std::cout<<"Done!"<<std::endl;
  //   std::cout<<"Number of fibers saved: "<<polydata->GetNumberOfLines()<<std::endl;
  //   }
  // catch(...)
  //   {
  //   std::cout << "Error while saving fiber file." << std::endl;
  //   throw;
  //   }

  return EXIT_SUCCESS;
}

