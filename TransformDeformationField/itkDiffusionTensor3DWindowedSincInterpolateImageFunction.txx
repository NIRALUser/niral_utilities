/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DWindowedSincInterpolateImageFunction.txx $
  Language:  C++
  Date:      $Date: 2010-04-29 11:58:49 -0400 (Thu, 29 Apr 2010) $
  Version:   $Revision: 13073 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DWindowedSincInterpolateFunction_txx
#define __itkDiffusionTensor3DWindowedSincInterpolateFunction_txx

#include "itkDiffusionTensor3DWindowedSincInterpolateImageFunction.h"

namespace itk
{
    
    
template< class TData ,
          unsigned int VRadius ,
          class TWindowFunction ,
          class TBoundaryCondition ,
          class TCoordRep
        >
void
DiffusionTensor3DWindowedSincInterpolateImageFunction< TData ,
                                                       VRadius ,
                                                       TWindowFunction ,
                                                       TBoundaryCondition ,
                                                       TCoordRep
                                                     >
::AllocateInterpolator()
{
  for( int i = 0 ; i < 6 ; i++ )
    {
    windowedSincInterpolator[ i ] = WindowedSincInterpolateImageFunctionType::New() ;
    this->m_Interpol[ i ] = windowedSincInterpolator[ i ] ;
    }
}



}//end itk namespace

#endif

