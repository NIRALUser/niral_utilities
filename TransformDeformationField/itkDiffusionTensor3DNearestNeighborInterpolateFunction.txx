/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DNearestNeighborInterpolateFunction.txx $
  Language:  C++
  Date:      $Date: 2010-04-29 11:58:49 -0400 (Thu, 29 Apr 2010) $
  Version:   $Revision: 13073 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DNearestNeighborInterpolateFunction_txx
#define __itkDiffusionTensor3DNearestNeighborInterpolateFunction_txx

#include "itkDiffusionTensor3DNearestNeighborInterpolateFunction.h"

namespace itk
{


template< class TData , class TCoordRep >
typename DiffusionTensor3DNearestNeighborInterpolateFunction< TData , TCoordRep >
::TensorDataType
DiffusionTensor3DNearestNeighborInterpolateFunction< TData , TCoordRep >
::EvaluateAtContinuousIndex( const ContinuousIndexType & index ) const
//::Evaluate( const PointType &point )
{
  if( this->m_Image.IsNotNull() )
    {
      typename DiffusionImageType::IndexType pixelIndex ;
      this->ConvertContinuousIndexToNearestIndex( index , pixelIndex ) ;
      return this->m_Image->GetPixel( pixelIndex ) ;
    }
  else
    {
      itkExceptionMacro( << " No InputImage set" ) ;
    }
}

}//end namespace itk

#endif
