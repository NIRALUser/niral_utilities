/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DTransform.txx $
  Language:  C++
  Date:      $Date: 2010/03/09 16:34:14 $
  Version:   $Revision: 1.2 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DTransform_txx
#define __itkDiffusionTensor3DTransform_txx

#include "itkDiffusionTensor3DTransform.h"

namespace itk
{

template< class TData >
DiffusionTensor3DTransform< TData >
::DiffusionTensor3DTransform()
{
  //Initialize the Measurement Frame to Identity
  m_MeasurementFrame.SetIdentity() ;
}


}//end namespace itk
#endif
