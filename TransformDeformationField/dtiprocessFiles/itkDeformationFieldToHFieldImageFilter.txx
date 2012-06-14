/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDeformationFieldToHFieldImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2010/03/09 18:01:02 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

namespace itk
{

//const double eps = 1e-16;
template <typename TInputImage, typename TOutputImage>
void DeformationFieldToHFieldImageFilter<TInputImage,TOutputImage>::GenerateData()
{
//  Superclass::GenerateInputRequestedRegion();
//  outputImage->SetRequrestedRegion(inputImage
  this->AllocateOutputs();

  const typename InputImageType::ConstPointer input(this->GetInput() );

  typename OutputImageType::Pointer output(this->GetOutput());
  
  typename InputImageType::RegionType inputRequestedRegion(input->GetRequestedRegion());

  typename OutputImageType::RegionType outputRequestedRegion(output->GetRequestedRegion());

  ImageRegionConstIteratorWithIndex<InputImageType> it = ImageRegionConstIteratorWithIndex<InputImageType>(
    input, inputRequestedRegion);
  
  ImageRegionIterator<OutputImageType> oit = ImageRegionIterator<OutputImageType>(
    output, outputRequestedRegion);
  
  //  typename InputImageType::SpacingType spacing = input->GetSpacing();
  for(it.GoToBegin(), oit.GoToBegin();!it.IsAtEnd(); ++it, ++oit)
    {
    InputPixelType dvec = it.Get();
    typename OutputImageType::IndexType index = it.GetIndex();

    oit.Set(this->ComputeH(input,index,dvec));
    }

}

template <typename TInputImage, typename TOutputImage>
typename TOutputImage::PixelType 
DeformationFieldToHFieldImageFilter<TInputImage,TOutputImage>::ComputeH(typename InputImageType::ConstPointer input,
                                                                        typename InputImageType::IndexType ind,
                                                                        typename InputImageType::PixelType dvec
                                                                       )
{
  typedef typename InputPixelType::ValueType CoordRepType;
  const unsigned int Dimension = InputImageType::ImageDimension;

  typedef itk::Point<CoordRepType, Dimension> PointType;
  PointType dpt;


  typedef ContinuousIndex<CoordRepType, Dimension> ContinuousIndexType;
  ContinuousIndexType dci;
  for(unsigned int i = 0; i < Dimension; ++i)
    {
    dpt[i] = dvec[i] ;
    }
//  input->TransformIndexToPhysicalPoint(ind,ipt);
  input->TransformPhysicalPointToContinuousIndex(dpt,dci);
  typename OutputImageType::PixelType pix ;
  for( unsigned int i = 0 ; i < Dimension ; ++i )
  {
    pix[ i ] = dci[ i ] + ind[ i ] ;
  }
  return pix ;
}

} // namespace itk
