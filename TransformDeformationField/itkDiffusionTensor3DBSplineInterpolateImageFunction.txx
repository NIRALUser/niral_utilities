/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DBSplineInterpolateImageFunction.txx $
  Language:  C++
  Date:      $Date: 2010-04-29 11:58:49 -0400 (Thu, 29 Apr 2010) $
  Version:   $Revision: 13073 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DBSplineInterpolateImageFunction_txx
#define __itkDiffusionTensor3DBSplineInterpolateImageFunction_txx

#include "itkDiffusionTensor3DBSplineInterpolateImageFunction.h"

namespace itk
{
    
template< class TData , class TCoordRep >
DiffusionTensor3DBSplineInterpolateImageFunction< TData , TCoordRep >
::DiffusionTensor3DBSplineInterpolateImageFunction()
{
  m_SplineOrder = 1 ;
}    
    
template< class TData , class TCoordRep >
void
DiffusionTensor3DBSplineInterpolateImageFunction< TData , TCoordRep >
::AllocateInterpolator()
{
  for( int i = 0 ; i < 6 ; i++ )
    {
    bSplineInterpolateFunction[ i ] = BSplineInterpolateFunction::New() ;
    bSplineInterpolateFunction[ i ]->SetSplineOrder( m_SplineOrder ) ;
    this->m_Interpol[ i ] = bSplineInterpolateFunction[ i ] ;
    }
}


}//end itk namespace

#endif

