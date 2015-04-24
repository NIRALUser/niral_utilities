#ifndef _itkCorrespondenceEvaluator_txx
#define _itkCorrespondenceEvaluator_txx

#include "itkCorrespondenceEvaluator.h"
#include <fstream>
#include <stdlib.h>
#include <time.h>

namespace itk
{

  template<class TOutputMeshType>
  CorrespondenceEvaluator<TOutputMeshType>
  ::CorrespondenceEvaluator()
  {
    this->m_NumberOfInputs = 0 ;
    srand ( time(0) ) ;
    this->gaussianDistribution = true ;
  }

  template<class TOutputMeshType>
  CorrespondenceEvaluator<TOutputMeshType>
  ::~CorrespondenceEvaluator()
  {
  }

  template<class TOutputMeshType>
  double CorrespondenceEvaluator<TOutputMeshType>
  ::GetGeneralization(unsigned int numberOfShapeParameters, double &stdError) 
  {
    double mean = 0, stddev = 0 ;
    std::vector < double > sample ;

    if ( numberOfShapeParameters > this->m_NumberOfInputs - 2 ) 
    {
      itkExceptionMacro ("Invalid number of shape parameters for generalization computation.") ;
      return -1 ;
    }

	sample.resize ( this->m_NumberOfInputs ) ;

    for ( int i = 0 ; i < this->m_NumberOfInputs ; i++ )
    {
      sample[i] = this->ComputeGeneralizationError ( i, numberOfShapeParameters ) ;
      mean += sample[i] ;
    }

	// compute sample mean
    mean /= this->m_NumberOfInputs ;

	// compute sample standard deviation
	for ( int i = 0 ; i < this->m_NumberOfInputs ; i++ )
	{
		sample[i] -= mean ;
		stddev += sample[i] * sample[i] ;
	}
	stddev = sqrt ( stddev / this->m_NumberOfInputs ) ;

	// compute standard error of generalization
	stdError = stddev / sqrt ( this->m_NumberOfInputs - 1 ) ;
	

    return mean ;
  }

  template<class TOutputMeshType>
  double CorrespondenceEvaluator<TOutputMeshType>
  ::GetSpecificity(unsigned int N, unsigned int numberOfShapeParameters, double &stdError) 
  {
    double specificity = 0 ;

    if ( numberOfShapeParameters > this->m_NumberOfInputs - 1 ) 
    {
      itkExceptionMacro ("Invalid number of shape parameters for generalization computation.") ;
      return -1 ;
    }

    // build model 
    MatrixType pcaMatrix, mean ;
    itk::Point <double, 3> sample ;
    int nPoints = this->m_Meshes[0]->GetNumberOfPoints() ;
    pcaMatrix.set_size( 3*nPoints, this->m_NumberOfInputs );
    mean.set_size ( 3, nPoints ) ;

    for (unsigned int row=0, pointId=0; pointId<nPoints; row+=3, pointId++)
    {
      // compute mean of the training set
      mean[0][pointId] = mean[1][pointId] = mean[2][pointId] = 0 ;

      for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ )
      {
        this->m_Meshes[sampleIdx]->GetPoint ( pointId, &sample ) ;

        mean[0][pointId] += sample[0] ;
        mean[1][pointId] += sample[1] ;
        mean[2][pointId] += sample[2] ;
      }
      mean[0][pointId] /= ( this->m_NumberOfInputs ) ;
      mean[1][pointId] /= ( this->m_NumberOfInputs ) ;
      mean[2][pointId] /= ( this->m_NumberOfInputs ) ;

      // build pca matrix, the matrix where we are going to perform pca to find principal modes of variation of the model
      for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ )
      {
        this->m_Meshes[sampleIdx]->GetPoint ( pointId, &sample ) ;
              
        (pcaMatrix)[row  ][sampleIdx] = sample[0] - mean[0][pointId];
        (pcaMatrix)[row+1][sampleIdx] = sample[1] - mean[1][pointId];
        (pcaMatrix)[row+2][sampleIdx] = sample[2] - mean[2][pointId];
      }
    }

    std::vector<double> SqrtEigenValues;
    MatrixType PrincipalDirections;

    SqrtEigenValues.resize ( this->m_NumberOfInputs ) ;

    vnl_svd<double> svd( pcaMatrix );
    for (unsigned int m=0; m<this->m_NumberOfInputs-1; m++) 
    {
      double sv = svd.W( m );
      SqrtEigenValues[m] = sqrt ( sv ) ;
    }
    PrincipalDirections = svd.U();

  	std::vector < double > sampleErrors ;
	  sampleErrors.resize ( N ) ;

    for ( int i = 0 ; i < N ; i++ )
    {
        // generate random shape parameters
        std::vector < double > b ;
        b.resize ( numberOfShapeParameters ) ;

        for ( unsigned int m = 0 ; m < numberOfShapeParameters ; m++ )
        {
            //b[m] = 3 * ( ( rand () / (double) RAND_MAX ) * 2 * SqrtEigenValues[m] - SqrtEigenValues[m] ) ;
            if ( this->gaussianDistribution )
            {
                b[m] = this->GenerateGaussianRandomNumber() * SqrtEigenValues[m] ;
            }
            else
            {
                b[m] = this->GenerateUniformRandomNumber() * 3 * SqrtEigenValues[m] ; 
            }
        }

        // construct random shape
        MatrixType reconstruction ;

        reconstruction.set_size ( 3, nPoints ) ;

        for ( unsigned int pt = 0 ; pt< nPoints ; pt++ )
        {
          reconstruction[0][pt] = mean[0][pt] ;
          reconstruction[1][pt] = mean[1][pt] ;
          reconstruction[2][pt] = mean[2][pt] ;

          for ( unsigned int m = 0 ; m < numberOfShapeParameters ; m++ )
          {
            // add the component associated with the mth mode of variation
            reconstruction[0][pt] += PrincipalDirections[3*pt][m] * b[m] ;
            reconstruction[1][pt] += PrincipalDirections[3*pt+1][m] * b[m] ;
            reconstruction[2][pt] += PrincipalDirections[3*pt+2][m] * b[m] ;
          }
        }

        // find nearest shape in the training set and compute error
        int minJ = -1 ; 
        double minDist = 999999999 ;
        double dist ;
		MatrixType currentMesh ;
		currentMesh.set_size ( 3, nPoints ) ;
        for ( int j = 0 ; j < this->m_NumberOfInputs ; j++ )
        {
			for ( unsigned int pt = 0 ; pt < nPoints ; pt++ )
            {
                this->m_Meshes[j]->GetPoint ( pt, &sample ) ;
				currentMesh[0][pt] = sample[0] ;
				currentMesh[1][pt] = sample[1] ;
				currentMesh[2][pt] = sample[2] ;
            }
            //if ( ( i == 0 ) && ( j == 0 ) )
       
			      dist = this->SurfaceDistance ( &currentMesh, &reconstruction, nPoints ) ;
            
            if ( dist < minDist )
            {
                minJ = j ;
                minDist = dist ;
            }
        }

        specificity += minDist ;
		sampleErrors[i] = minDist ;
    }

	// compute mean error
    specificity /= N ;

	// compute sample standard deviation
	double stdDev = 0 ;
	for ( int i = 0 ; i < N ; i++ )
	{
		sampleErrors[i] -= specificity ;
		stdDev += sampleErrors[i] * sampleErrors[i] ;
	}
	stdDev = sqrt ( stdDev / N ) ;

	// compute standard error of specificity
	stdError = stdDev / sqrt ( N ) ;

    return specificity ;
  }

  template<class TOutputMeshType>
  void CorrespondenceEvaluator<TOutputMeshType>
  ::SetNumberOfInputs( unsigned int number )
  {
    if ( number < 0 ) 
    {
      itkExceptionMacro( "Invalid number of inputs for correspondence evaluation." );
      return ;
    }
    this->m_NumberOfInputs = number ;
    this->m_Meshes.resize ( number ) ;
  }

  template<class TOutputMeshType>
  void CorrespondenceEvaluator<TOutputMeshType>
  ::SetInput( unsigned int idx, MeshPointer input )
  {
    if ( ( idx < 0 ) || ( idx >= this->m_NumberOfInputs ) ) 
    {
      itkExceptionMacro ("Invalid input id for correspondence evalution.") ;
      return ;
    }
    
    this->m_Meshes[idx] = input ;
    this->BringInputToOrigin ( idx ) ;
  }

  template<class TOutputMeshType>
  void CorrespondenceEvaluator<TOutputMeshType>
  ::BringInputToOrigin (unsigned int id)
  {
    int nPoints = this->m_Meshes[id]->GetNumberOfPoints() ;
    double center[3] ;
    itk::Point <double, 3> sample ;

    center[0] = center[1] = center[2] = 0 ;

    for (unsigned int pointId=0; pointId < nPoints; pointId++)
    {
      this->m_Meshes[id]->GetPoint ( pointId, &sample ) ;
      center[0] += sample[0] ;
      center[1] += sample[1] ;
      center[2] += sample[2] ;
    }  
    center[0] /= nPoints ;  
    center[1] /= nPoints ;
    center[2] /= nPoints ;
    std::cout << center[0] << " " << center[1] << " " << center[2] << std::endl ;

    for (unsigned int pointId=0; pointId < nPoints; pointId++)
    {
      this->m_Meshes[id]->GetPoint ( pointId, &sample ) ;
      sample[0] -= center[0] ;
      sample[1] -= center[1] ;
      sample[2] -= center[2] ;
      this->m_Meshes[id]->SetPoint ( pointId, sample ) ;
    }
  }

  template<class TOutputMeshType>
  double CorrespondenceEvaluator<TOutputMeshType>
  ::ComputeGeneralizationError ( unsigned int i, unsigned int M ) 
  {
    if ( ( i < 0 ) || ( i >= this->m_NumberOfInputs ) ) 
    {
      itkExceptionMacro ("Invalid input id for correspondence evalution.") ;
      return -1 ;
    }
   
    //////////////////////////////////////////////////////////////////////
    // 1. BUILD THE MODEL FROM THE TRAINING SET, WITH Xi REMOVED        //
    //////////////////////////////////////////////////////////////////////

    MatrixType pcaMatrix, mean ;
    itk::Point <double, 3> sample ;
    int nPoints = this->m_Meshes[0]->GetNumberOfPoints() ;
    pcaMatrix.set_size( 3*nPoints, this->m_NumberOfInputs-1 );
    mean.set_size ( 3, nPoints ) ;

    for (unsigned int row=0, pointId=0; pointId<nPoints; row+=3, pointId++)
    {
      // compute mean of the training set, with Xi removed
      mean[0][pointId] = mean[1][pointId] = mean[2][pointId] = 0 ;

      for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ )
      {
        if ( sampleIdx == i ) 
        {
          continue ;
        }
        this->m_Meshes[sampleIdx]->GetPoint ( pointId, &sample ) ;

        mean[0][pointId] += sample[0] ;
        mean[1][pointId] += sample[1] ;
        mean[2][pointId] += sample[2] ;
      }
      mean[0][pointId] /= ( this->m_NumberOfInputs - 1 ) ;
      mean[1][pointId] /= ( this->m_NumberOfInputs - 1 ) ;
      mean[2][pointId] /= ( this->m_NumberOfInputs - 1 ) ;

      // build pca matrix, the matrix where we are going to perform pca to find principal modes of variation of the model
      unsigned int skip = 0 ;
      for ( unsigned int sampleIdx = 0 ; sampleIdx < this->m_NumberOfInputs ; sampleIdx++ )
      {
        if ( sampleIdx == i ) 
        {
          skip = 1 ;
          continue ;
        }
        
        this->m_Meshes[sampleIdx]->GetPoint ( pointId, &sample ) ;
              
        (pcaMatrix)[row  ][sampleIdx-skip] = sample[0] - mean[0][pointId];
        (pcaMatrix)[row+1][sampleIdx-skip] = sample[1] - mean[1][pointId];
        (pcaMatrix)[row+2][sampleIdx-skip] = sample[2] - mean[2][pointId];
      }
    }

    std::vector<double> EigenValues;
    MatrixType PrincipalDirections;

    EigenValues.resize ( this->m_NumberOfInputs ) ;

    vnl_svd<double> svd( pcaMatrix );
    for (unsigned int m=0; m<this->m_NumberOfInputs-1; m++) 
    {
      double sv = svd.W( m );
      EigenValues[m] = sv ;
    }
    PrincipalDirections = svd.U();
    
    //////////////////////////////////////////////////////////////////////
    // 2. ESTIMATE THE MODEL PARAMETERS FOR Xi                          //
    //////////////////////////////////////////////////////////////////////

    std::vector < double > b ;
    b.resize ( M ) ;
    
    for ( unsigned int m = 0 ; m < M ; m++ )
    {
        b[m] = 0 ;
        for ( unsigned int pt = 0 ; pt < nPoints ; pt++ )
        {
            this->m_Meshes[i]->GetPoint ( pt, &sample ) ;
            for ( unsigned int internal = 0 ; internal < 3 ; internal++ )
            {
                b[m] += PrincipalDirections[pt*3+internal][m] * ( sample[internal] - mean[internal][pt] ) ;
            }
            // (3*642)x57
            // (3*642)
        }
    }
    
    //////////////////////////////////////////////////////////////////////
    // 3. RECONSTRUCT Xi USING M SHAPE PARAMETERS                       //
    //////////////////////////////////////////////////////////////////////

    MatrixType reconstruction ;

    reconstruction.set_size ( 3, nPoints ) ;

    for ( unsigned int pt = 0 ; pt< nPoints ; pt++ )
    {
      reconstruction[0][pt] = mean[0][pt] ;
      reconstruction[1][pt] = mean[1][pt] ;
      reconstruction[2][pt] = mean[2][pt] ;

      for ( unsigned int m = 0 ; m < M ; m++ )
      {
        // add the component associated with the mth mode of variation
        reconstruction[0][pt] += PrincipalDirections[pt*3][m] * b[m] ;
        reconstruction[1][pt] += PrincipalDirections[pt*3+1][m] * b[m] ;
        reconstruction[2][pt] += PrincipalDirections[pt*3+2][m] * b[m] ;
      }
    }
    
    //////////////////////////////////////////////////////////////////////
    // 4. CALCULATE THE APPROXIMATION ERROR                             //
    //////////////////////////////////////////////////////////////////////
    
	MatrixType currentMesh ;
	currentMesh.set_size ( 3, nPoints ) ;
	for ( unsigned int pt = 0 ; pt < nPoints ; pt++ )
	{
		this->m_Meshes[i]->GetPoint ( pt, &sample ) ;
		currentMesh[0][pt] = sample[0] ;
		currentMesh[1][pt] = sample[1] ;
		currentMesh[2][pt] = sample[2] ;
	}

	return this->SurfaceDistance ( &currentMesh, &reconstruction, nPoints ) ;

  }

  template<class TOutputMeshType>
  double CorrespondenceEvaluator<TOutputMeshType>
  ::SurfaceDistance (MatrixType *from, MatrixType *to, unsigned int nPoints) 
  {
	double distance[3] ;
	double error = 0 ;

	for ( unsigned int pt = 0 ; pt < nPoints ; pt++ )
	{
	  distance[0] = (*from)[0][pt] - (*to)[0][pt] ;
      distance[1] = (*from)[1][pt] - (*to)[1][pt] ;
      distance[2] = (*from)[2][pt] - (*to)[2][pt] ;
      error += sqrt ( distance[0] * distance[0] + distance[1] * distance[1] + distance[2] * distance[2] ) ;	
	}
	return error / nPoints ;
  }

  template<class TOutputMeshType>
  void CorrespondenceEvaluator<TOutputMeshType>
  ::WriteSurfaceForDebug (MatrixType *surface, unsigned int nPoints, std::string filename) 
  {
	  std::ofstream out ;
    out.open (filename.c_str()) ;
    for ( unsigned int pt = 0 ; pt < nPoints ; pt++ )
	  {
      out << (*surface)[0][pt] << " " ;
      out << (*surface)[1][pt] << " " ;
      out << (*surface)[2][pt] << std::endl  ;
	  }
	  out.close () ;
  }

  template<class TOutputMeshType>
  double CorrespondenceEvaluator<TOutputMeshType>
  ::GenerateUniformRandomNumber () 
  {
    // generate a random number uniformly distributed between [-1..1]
    double uniform = ( rand () / (double) RAND_MAX ) * 2 - 1.0 ;
    return uniform ;
  }

  template<class TOutputMeshType>
  double CorrespondenceEvaluator<TOutputMeshType>
  ::GenerateGaussianRandomNumber () 
  {
    // generate a random number normally distributed between [-1..1]

    // from http://www.taygeta.com/random/gaussian.html
    // Algorithm by Dr. Everett (Skip) Carter, Jr.

    float x1, x2, w, y1 ;

    do {
        x1 = this->GenerateUniformRandomNumber() ;
        x2 = this->GenerateUniformRandomNumber() ;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;

    return y1 ;
  } 

}

#endif
