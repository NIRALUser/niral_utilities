#include <itkMetaMeshConverter.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "itkCorrespondenceEvaluator.h"

//////////////////////////////////////////////////////////////////////////////////

typedef itk::DefaultStaticMeshTraits <double, 3, 3, double, double> traitsType ;
typedef itk::Mesh <double, 3, traitsType> meshType ;
typedef itk::MeshSpatialObject <meshType> meshSOType ;
typedef itk::MetaMeshConverter <3, double, traitsType> MeshConverterType ;
typedef itk::CorrespondenceEvaluator <meshType> evaluatorType ;

//////////////////////////////////////////////////////////////////////////////////

bool readShapes ( std::string listName, evaluatorType::Pointer evaluator, unsigned int &nSamples )
{
  std::ifstream shapeListFile;
  shapeListFile.open( listName.c_str(), std::ios_base::in );
  if (!shapeListFile)
  {
    std::cerr << "Unable to open shape list \"" << listName << "\"!" << std::endl;
    return false ;
  }

  // parse the shape list file
  int numSamples = 0 ;
  std::string currentFileName ;
  std::vector < std::string > fileNames ;

  while (!shapeListFile.eof())
  {
    std::getline( shapeListFile, currentFileName );
    if (currentFileName.empty() || currentFileName[0]=='#') { continue; }

    numSamples++;
    fileNames.resize ( numSamples ) ;
    fileNames[numSamples - 1] = currentFileName ;
  }

  shapeListFile.close();

  evaluator->SetNumberOfInputs ( numSamples ) ;

  // read the shapes
  MeshConverterType * itkConverter = new MeshConverterType() ;
  for ( int i = 0 ; i < numSamples ; i++ )
  {
    std::cout << "reading " << fileNames[i] << std::endl;
    meshSOType::Pointer meshSO = itkConverter->ReadMeta ( fileNames[i].c_str() ) ;
    meshType::Pointer mesh = meshSO->GetMesh() ;
    evaluator->SetInput ( i, mesh ) ;
  }
  delete (itkConverter);

  nSamples = numSamples ;

  return true ;
}

int main( int argc, char *argv[] )
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " MeshListFile OutputFile [-gen] [-spec] [-uniform]" << std::endl << std::endl ;
    std::cout << "MeshListFile: List of shapes in the population" << std::endl;
    std::cout << "OutputFile: Where to store the evaluation results" << std::endl;
    std::cout << "-gen: Compute generalization" << std::endl;
    std::cout << "-spec: Compute specificity" << std::endl;
    std::cout << "-uniform: Use uniform distribution for random numbers (default is Gaussian)" << std::endl ;
    std::cout << "-numSpec: number of samples created for specificity" << std::endl ;
    std::cout << std::endl ;
    return -1;
  }

  // parse command line arguments
  std::string meshListFileName = argv[1];
  std::string outputFileName = argv[2];

  bool computeGeneralization = false ;
  bool computeSpecificity = false ;
  bool gaussian = true ;

  int N = 1000 ;

  if ( argc > 3 )
  {
    for ( short int i = 3 ; i < argc ; i++ )
    {
      if ( (!computeGeneralization) && ( !strcmp ( "-gen", argv[i] ) ) )
      {
        computeGeneralization = true ;
      }
      if ( (!computeSpecificity) && ( !strcmp ( "-spec", argv[i] ) ) )
      {
        computeSpecificity = true ;
      }
      if ( (gaussian) && ( !strcmp ( "-uniform", argv[i] ) ) )
	{
	  gaussian = false ;
	}
      if (( !strcmp ( "-numSpec", argv[i] ) ) )
	{
	  N = atoi(argv[i+1]);
	}
    }
  }

  // read input meshes
  unsigned int numSamples ;
  evaluatorType::Pointer evaluator = evaluatorType::New() ;

  if ( ! readShapes ( meshListFileName, evaluator, numSamples ) )
  {
    return -1 ;
  }
    std::cout << "reading done " << std::endl;

  // evaluate correspondence
  std::vector < double > generalization, specificity ;
  std::vector < double > generalizationError, specificityError ;
  generalization.resize ( numSamples - 1 ) ;
  specificity.resize ( numSamples ) ;
  generalizationError.resize ( numSamples - 1 ) ;
  specificityError.resize ( numSamples ) ;

  if ( computeGeneralization )
  {
  for ( unsigned int i = 0  ; i < numSamples - 1 ; i++ )
  {
    std::cout << "Starting generalization computation with " << i << " shape parameters... " ;
    generalization[i] = evaluator->GetGeneralization ( i, generalizationError[i] ) ;
    std::cout << "Finished." << std::endl ;
  }
  }
  if ( computeSpecificity )
  {
  evaluator->SetDistribution ( gaussian ) ;
    for ( unsigned int i = 0 ; i < numSamples ; i++ )
  {
    std::cout << "Starting specificity computation with " << i << " shape parameters... " ;
    specificity[i] = evaluator->GetSpecificity ( N, i, specificityError[i] ) ;
    std::cout << "Finished." << std::endl ;
    }
  }

  // write out evaluation results
  std::ofstream outputFile ;
  outputFile.open ( outputFileName.c_str() ) ;
  outputFile << "Correspondence evaluation for: " << argv[1] << std::endl << std::endl ;
  if ( computeGeneralization )
  {
    outputFile << "Generalization: " << std::endl ;
    for ( unsigned int i = 0 ; i < numSamples - 1 ; i++ )
    {
      outputFile << i << " " << generalization[i] << " " << generalizationError[i] << std::endl ;
    }
  }
  if ( computeSpecificity )
  {
    outputFile << "Specificity: " << std::endl ;
    for ( unsigned int i = 0 ; i < numSamples ; i++ )
    {
      outputFile << i << " " << specificity[i] << " " << specificityError[i] << std::endl ;
    }
  }
  outputFile.close () ;

  return 0 ;
}

