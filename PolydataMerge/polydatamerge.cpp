#include <string>
#include <iostream>
#include <fstream>

// VTK includes

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
#include <vtkPolyLine.h>
#include <vtkAppendPolyData.h>
#include <vtkVersion.h>

#include "polydatamergeCLP.h"

int main(int argc, char* argv[])
{
  PARSE_ARGS;

  if(fiberFiles[0] == "")
    {
      std::cerr << "An first input fiber file has to be specified" << std::endl;
      return EXIT_FAILURE;
    }

 if(fiberFiles[1] == "")
    {
	std::cerr << "An second input fiber file has to be specified" << std::endl;
      return EXIT_FAILURE;
    }

 if(fiberOutput== "")
   {
	std::cerr << "An output fiber file name has to be specified" << std::endl;
      return EXIT_FAILURE;
    }

  vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New();
  std::cout << "Reading " << fiberFiles[0]<< std::endl;
  reader1->SetFileName(fiberFiles[0].c_str());
  reader1->Update();
 
  // Extract the polydata 1
  vtkSmartPointer<vtkPolyData> polydata1 =
    reader1->GetOutput();

  vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
  std::cout << "Reading " << fiberFiles[1]<< std::endl;
  reader2->SetFileName(fiberFiles[1].c_str());
  reader2->Update();
 
  // Extract the polydata 2
  vtkSmartPointer<vtkPolyData> polydata2 =
    reader2->GetOutput();

 vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

 vtkAppendPolyData *apd = vtkAppendPolyData::New(); 
 #if (VTK_MAJOR_VERSION < 6)
 apd->AddInput (polydata2);
 apd->AddInput (polydata1);
 #else
 apd->AddInputData (polydata2);
 apd->AddInputData (polydata1);
 #endif
 polydata=apd->GetOutput ();

  if(fiberOutput != "")
  {
  std::cout<<std::endl;
  std::cout<<"Saving fibers...."<<std::endl;
  std::cout<<fiberOutput<<std::endl;
  vtkPolyDataWriter * fiberwriter = vtkPolyDataWriter::New();
  fiberwriter->SetFileName(fiberOutput.c_str());
  #if (VTK_MAJOR_VERSION < 6)
  fiberwriter->SetInput(polydata);
  #else
  fiberwriter->SetInputData(polydata);
  #endif
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


