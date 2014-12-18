#ifndef DTITYPES_H
#define DTITYPES_H

// ITK Data types
#include <itkImage.h>
#include <itkVectorContainer.h>
#include <itkVector.h>
#include <itkDiffusionTensor3D.h>
#include <itkAffineTransform.h>
#include <itkDTITubeSpatialObject.h>
#include <itkVectorImage.h>

// VNL Includes
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>

// Define necessary types for images
typedef double RealType;
typedef double TransformRealType;
typedef unsigned char LabelType;
const unsigned int DIM = 3;

typedef unsigned short ScalarPixelType;
typedef itk::DiffusionTensor3D<double> TensorPixelType;
typedef itk::Vector<double,3> DeformationPixelType;

typedef itk::VectorImage<ScalarPixelType, DIM> VectorImageType;
typedef itk::Image<TensorPixelType, DIM> TensorImageType;

typedef itk::Image<DeformationPixelType, DIM> DeformationImageType;

typedef itk::Image<RealType, DIM> RealImageType;
typedef itk::Image<ScalarPixelType, DIM> IntImageType;


typedef TensorImageType::SizeType ImageSizeType;
typedef TensorImageType::SpacingType ImageSpacingType;

typedef itk::AffineTransform<TransformRealType,3> AffineTransformType;

typedef vnl_vector_fixed<double, 3> GradientType;
typedef itk::VectorContainer<unsigned int, GradientType> GradientListType;

enum InterpolationType {NearestNeighbor, Linear, Cubic};
enum TensorReorientationType {FiniteStrain, PreservationPrincipalDirection};

enum EigenValueIndex {Lambda1 = 0, Lambda2, Lambda3};

typedef itk::DTITubeSpatialObject<3> DTITubeType;
typedef DTITubeType::TubePointType DTIPointType;
typedef DTITubeType::PointListType DTIPointListType;


#endif
