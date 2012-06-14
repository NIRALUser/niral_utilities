
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

// ITK includes
#include <itkDiffusionTensor3D.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkVersion.h>

#include "deformationfieldoperations.h"
#include "dtitypes.h"
#include "polydatatransformCLP.h"


int main(int argc, char* argv[])
{
  PARSE_ARGS;

  if(fiberFile == "")
    {
      std::cerr << "A fiber file has to be specified" << std::endl;
      return EXIT_FAILURE;
    }
 
  
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  std::cout << "Reading " << fiberFile<< std::endl;
  reader->SetFileName(fiberFile.c_str());
  reader->Update();
 
  // Extract the polydata
  vtkSmartPointer<vtkPolyData> polydata =
    reader->GetOutput();
  
 
  DeformationImageType::Pointer deformationfield(NULL);
  if(hField != "")
    deformationfield = readDeformationField(hField, HField);
  else if(displacementField != "")
    deformationfield = readDeformationField(displacementField, Displacement);
  else
    deformationfield = NULL;
  
  typedef itk::VectorLinearInterpolateImageFunction<DeformationImageType, double> DeformationInterpolateType;
  DeformationInterpolateType::Pointer definterp(NULL);
  if(deformationfield)
  {
    definterp = DeformationInterpolateType::New();
    definterp->SetInputImage(deformationfield);
  } else {
    std::cout << "A deformation field has to be specified" << std::endl;
      return EXIT_FAILURE;
  }
  
    
    typedef DeformationInterpolateType::ContinuousIndexType ContinuousIndexType;
    ContinuousIndexType ci, origci;
    // For each point along the fiber
     vtkIdType numPoints= polydata->GetNumberOfPoints();
     vtkPoints* inPts = polydata->GetPoints();
     vtkPoints    * points = vtkPoints::New();
  

    typedef DTIPointType::PointType PointType;
    PointType fiberpoint;
    double fiberpointtemp[3];

    for(int i=0;i<numPoints;i++)
    {	
        inPts->GetPoint( i, fiberpointtemp );
	//convert RAS to LPS (vtk)
	fiberpoint[0] = - fiberpointtemp[0];
	fiberpoint[1] = - fiberpointtemp[1];
	fiberpoint[2] = + fiberpointtemp[2];	

        deformationfield->TransformPhysicalPointToContinuousIndex(fiberpoint,ci);
	if( !deformationfield->GetLargestPossibleRegion().IsInside( ci ) )
        {
          std::cerr << "Fiber is outside deformation field image. Deformation field has to be in the fiber space. Warning: Original position will be used" << std::endl ;
        }
        DeformationPixelType warp(definterp->EvaluateAtContinuousIndex(ci).GetDataPointer());
        for(unsigned int j =0; j < 3; j++)
          fiberpoint[j] +=warp[j]; 

	//convert LPS to RAS (vtk)
	fiberpoint[0] = - fiberpoint[0];
	fiberpoint[1] = - fiberpoint[1];
	fiberpoint[2] = + fiberpoint[2];

        points->InsertPoint(i,fiberpoint[0],fiberpoint[1],fiberpoint[2]);
	
    }
    polydata->SetPoints(points);
    

  
  if(fiberOutput != "")
  {
  std::cout<<std::endl;
  std::cout<<"Saving fibers...."<<std::endl;
  std::cout<<fiberOutput<<std::endl;
  vtkPolyDataWriter * fiberwriter = vtkPolyDataWriter::New();
  fiberwriter->SetFileName(fiberOutput.c_str());
  fiberwriter->SetInput(polydata);
  fiberwriter->SetFileTypeToBinary();
  fiberwriter->Update();
  try
    {
    fiberwriter->Write();
    std::cout<<"Done!"<<std::endl;
    std::cout<<"Number of fibers saved: "<<polydata->GetNumberOfLines()<<std::endl;
    }
  catch(...)
    {
    std::cout << "Error while saving fiber file." << std::endl;
    throw;
    }
	}

  return EXIT_SUCCESS;
}

