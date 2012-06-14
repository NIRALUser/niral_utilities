/*The Log-Euclidean framework is protected by a patent: for any commercial use, please contact {Vincent.Arsigny, Xavier.Pennec, Nicholas.Ayache}@Sophia.Inria.fr
*/
#ifndef __itkDiffusionTensor3DExpImageFilter_h
#define __itkDiffusionTensor3DExpImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"
#include <itkMatrix.h>
#include <itkDiffusionTensor3D.h>
#include "itkDiffusionTensor3DExtended.h"

namespace itk
{
  
/** \class DiffusionTensor3DExpImageFilter
 * \brief Computes pixel-wise the exponential of diffusion tensors.
 *
 * This filter is templated over the pixel type of the input image
 * and the pixel type of the output image. 
 *
 * The filter will walk over all the pixels in the input image, and for
 * each one of them it will do the following: 
 *
 * - cast the pixel value to \c DiffusionTensor3DExtended<double>, 
 * - Compute the eigenvalues and the eigenvectors
 * - apply the \c vcl_exp() function to the eigenvalues 
 *   of the output image 
 * - store the casted value into the output image.
 * 
 * The filter expect both images to have the same dimension (e.g. both 2D, 
 * or both 3D, or both ND).
 * The filter needs DiffusionTensor3D images to work
 *
 */
namespace Functor {  
  
template< class TInput , class TOutput >
class DiffusionTensor3DExp
{
public:
  DiffusionTensor3DExp() {} ;
  ~DiffusionTensor3DExp() {} ;
  bool operator!=( const DiffusionTensor3DExp & other ) const
    {
    return ( *this != other ) ;
    }
  bool operator==( const DiffusionTensor3DExp & other ) const
    {
    return !( *this != other ) ;
    }
  inline DiffusionTensor3D< TOutput > operator()( const DiffusionTensor3D< TInput > & A )
    {
    DiffusionTensor3D< TOutput > tensor ;
    Matrix< double , 3 , 3 > mat ;
    Matrix< double , 3 , 3 > matexp ;
    typename DiffusionTensor3DExtended< double >::EigenValuesArrayType eigenValues ;
    typename DiffusionTensor3DExtended< double >::EigenVectorsMatrixType eigenVectors ;
    DiffusionTensor3DExtended< double > tensorDouble ;
    tensorDouble=( DiffusionTensor3DExtended< TInput > )A ;
    tensorDouble.ComputeEigenAnalysis( eigenValues,eigenVectors ) ;
    //mat.Fill(0);
    for( int i = 0 ; i < 3 ; i++ )
      {
      mat[ i ][ i ] = vcl_exp( eigenValues[ i ] ) ;
      }
    eigenVectors = eigenVectors.GetTranspose();
    matexp = eigenVectors * mat * eigenVectors.GetInverse() ;
    tensorDouble.SetTensorFromMatrix( matexp ) ;
    for( int i = 0 ; i < 6 ; i++ )
      { tensor[ i ] = ( TOutput ) tensorDouble[ i ] ; }
    return tensor ;
    }
}; 
}//end of Functor namespace

template < class TInputImage , class TOutputImage >
class DiffusionTensor3DExpImageFilter :
    public
    UnaryFunctorImageFilter< TInputImage,TOutputImage , 
         Functor::DiffusionTensor3DExp< 
    typename TInputImage::PixelType::ComponentType, 
    typename TOutputImage::PixelType::ComponentType >   >
{
public:
  /** Standard class typedefs. */
  typedef DiffusionTensor3DExpImageFilter  Self ;
  typedef UnaryFunctorImageFilter< TInputImage , TOutputImage, 
       Functor::DiffusionTensor3DExp< typename TInputImage::PixelType , 
                                      typename TOutputImage::PixelType> >
                                                                 Superclass ;
  typedef SmartPointer< Self >   Pointer ;
  typedef SmartPointer< const Self >  ConstPointer ;

  /** Method for creation through the object factory. */
  itkNewMacro( Self ) ;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputTensorTypeCheck ,
    ( Concept::SameType< DiffusionTensor3D< typename TInputImage::PixelType::ComponentType > ,
                         typename TInputImage::PixelType > ) ) ;
  itkConceptMacro( OutputTensorTypeCheck ,
    ( Concept::SameType< DiffusionTensor3D< typename TOutputImage::PixelType::ComponentType > ,
                         typename TOutputImage::PixelType> ) ) ;

  /** End concept checking */
#endif

protected:
  DiffusionTensor3DExpImageFilter() {}
  virtual ~DiffusionTensor3DExpImageFilter() {}

private:
  DiffusionTensor3DExpImageFilter( const Self& ) ; //purposely not implemented
  void operator=( const Self & ) ; //purposely not implemented

};

} // end namespace itk


#endif
