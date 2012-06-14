/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DInterpolateImageFunction.txx $
  Language:  C++
  Date:      $Date: 2010-04-29 11:58:49 -0400 (Thu, 29 Apr 2010) $
  Version:   $Revision: 13073 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DInterpolateImageFunction_txx
#define __itkDiffusionTensor3DInterpolateImageFunction_txx

#include "itkDiffusionTensor3DInterpolateImageFunction.h"

namespace itk
{

template< class TData , class TCoordRep >
DiffusionTensor3DInterpolateImageFunction< TData , TCoordRep >
::DiffusionTensor3DInterpolateImageFunction()
{
//  m_InputImage = 0 ;
  latestTime = 0 ;
//  SetDefaultPixelValue( ZERO ) ;
}





/*
template< class TData >
void
DiffusionTensor3DInterpolateImageFunction< TData >
::SetDefaultPixelValue( TensorRealType defaultPixelValue )
{
  m_DefaultPixelValue = defaultPixelValue ;
  m_DefaultPixel.SetIdentity() ;
  for( unsigned int i = 0 ; i < 3 ; i++ ) 
    {
    m_DefaultPixel( i , i ) *= static_cast< TData >( this->m_DefaultPixelValue ) ;
    }
}
*/


}//end namespace itk

#endif
