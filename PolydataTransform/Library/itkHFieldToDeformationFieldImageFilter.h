/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHFieldToDeformationFieldImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHFieldToDeformationFieldImageFilter_h
#define __itkHFieldToDeformationFieldImageFilter_h

#include <itkImageToImageFilter.h>

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.

/** \class HFieldToDeformationFieldImageFilter
 * \brief Computes the Mean Diffusivity for every pixel of a input tensor image.
 *
 * HFieldToDeformationFieldImageFilter applies pixel-wise the invokation for
 * computing the mean diffusivity of every pixel. The pixel type of the
 * input image is expected to implement a method GetTrace(), and
 * to specify its return type as RealValueType.
 *
 * \sa TensorRelativeAnisotropyImageFilter
 * \sa TensorFractionalAnisotropyImageFilter
 * \sa DiffusionTensor3D
 *
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template <typename TInputImage,
          typename TOutputImage=TInputImage>
class ITK_EXPORT HFieldToDeformationFieldImageFilter :
    public
ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef HFieldToDeformationFieldImageFilter  Self;
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
  void PrintSelf(std::ostream& os, Indent indent) //const override
    { this->Superclass::PrintSelf( os, indent ); }

  // need to override GenerateData (This should be threaded)
  void GenerateData();// override;

  OutputPixelType ComputeDisplacement(typename InputImageType::ConstPointer input,
                                     typename InputImageType::IndexType ind,
                                     typename InputImageType::PixelType hvec);


protected:
  HFieldToDeformationFieldImageFilter() {};
  virtual ~HFieldToDeformationFieldImageFilter() {};

private:
  HFieldToDeformationFieldImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHFieldToDeformationFieldImageFilter.txx"
#endif


#endif
