/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDeformationFieldToHFieldImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2010/03/09 18:01:02 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDeformationFieldToHFieldImageFilter_h
#define __itkDeformationFieldToHFieldImageFilter_h

#include <itkImageToImageFilter.h>

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.

/** \class DeformationFieldToHFieldImageFilter
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template <typename TInputImage,
          typename TOutputImage=TInputImage>
class ITK_EXPORT DeformationFieldToHFieldImageFilter :
    public
ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DeformationFieldToHFieldImageFilter  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef typename TOutputImage::PixelType        OutputPixelType;
  typedef typename Superclass::InputImageType     InputImageType;
  typedef typename TInputImage::PixelType         InputPixelType;
  typedef typename InputPixelType::ValueType      InputValueType;

  typedef typename TInputImage::SpacingType SpacingType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }

  // need to override GenerateData (This should be threaded)
  void GenerateData();

  OutputPixelType ComputeH(typename InputImageType::ConstPointer input,
                           typename InputImageType::IndexType ind,
                           typename InputImageType::PixelType dvec
                          );

  
protected:
  DeformationFieldToHFieldImageFilter() {};
  virtual ~DeformationFieldToHFieldImageFilter() {};

private:
  DeformationFieldToHFieldImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformationFieldToHFieldImageFilter.txx"
#endif
  
  
#endif
