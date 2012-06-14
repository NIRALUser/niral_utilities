/*The Log-Euclidean framework is protected by a patent: for any commercial use, please contact {Vincent.Arsigny, Xavier.Pennec, Nicholas.Ayache}@Sophia.Inria.fr
*/

#ifndef __itkDiffusionTensor3DLogImageFilter_h
#define __itkDiffusionTensor3DLogImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"
#include <itkMatrix.h>
#include "itkDiffusionTensor3DExtended.h"
#include "define.h"
//#define NEG (-25)

namespace itk
{
  
/** \class DiffusionTensor3DLogImageFilter
 * \brief Computes the log pixel-wise of diffusion tensors.
 *
 * This filter is templated over the pixel type of the input image
 * and the pixel type of the output image. 
 *
 * The filter will walk over all the pixels in the input image, and for
 * each one of them it will do the following: 
 *
 * - cast the pixel value to \c DiffusionTensor3DExtended<double>, 
 * - Compute the eigenvalues and the eigenvectors
 * - apply the \c vcl_log() function to the eigenvalues 
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
class DiffusionTensor3DLog
{
public:
  DiffusionTensor3DLog() {}
  ~DiffusionTensor3DLog() {}
  bool operator!=( const DiffusionTensor3DLog & other ) const
  {
    return ( *this != other ) ;
  }
  bool operator==( const DiffusionTensor3DLog & other ) const
  {
    return !( *this != other ) ;
  }
  inline DiffusionTensor3D< TOutput > operator()
       ( const DiffusionTensor3D< TInput > & A )
  {
    DiffusionTensor3D< TOutput > tensor;
    Matrix< double , 3 , 3 > mat ;
    Matrix< double , 3 , 3 > matlog ;
    typename DiffusionTensor3DExtended< double >::EigenValuesArrayType eigenValues ;
    typename DiffusionTensor3DExtended< double >::EigenVectorsMatrixType eigenVectors ;
    DiffusionTensor3DExtended< double > tensorDouble ;
    tensorDouble = ( DiffusionTensor3DExtended< TInput > ) A ;
    tensorDouble.ComputeEigenAnalysis( eigenValues , eigenVectors ) ;
    for( int i = 0 ; i < 3 ; i++ )
      {
//      if( eigenValues[ i ] > 0 )
      if( eigenValues[ i ] > ZERO )
        { mat[ i ][ i ] = vcl_log( eigenValues[ i ] ) ; }
      else
        { mat [i][i] = LOGNEG ; }
      }
    eigenVectors = eigenVectors.GetTranspose();
    matlog = eigenVectors * mat * eigenVectors.GetInverse() ;
    tensorDouble.SetTensorFromMatrix( matlog ) ;
    for( int i = 0 ; i < 6 ; i++ )
      { tensor[ i ] = ( TOutput ) tensorDouble[ i ] ; }
    return tensor ;
  }
}; 
}//end of Functor namespace

template < class TInputImage, class TOutputImage >
class DiffusionTensor3DLogImageFilter :
    public
  UnaryFunctorImageFilter< TInputImage , TOutputImage , 
                        Functor::DiffusionTensor3DLog< 
  typename TInputImage::PixelType::ComponentType , 
  typename TOutputImage::PixelType::ComponentType> >
{
public:
  /** Standard class typedefs. */
  typedef DiffusionTensor3DLogImageFilter Self ;
  typedef UnaryFunctorImageFilter< TInputImage , TOutputImage , 
           Functor::DiffusionTensor3DLog< typename TInputImage::PixelType , 
                                typename TOutputImage::PixelType> >  Superclass ;
  typedef SmartPointer< Self > Pointer ;
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
                                            typename TOutputImage::PixelType>) ) ;

  /** End concept checking */
#endif

protected:
  DiffusionTensor3DLogImageFilter() {}
  virtual ~DiffusionTensor3DLogImageFilter() {}

private:
  DiffusionTensor3DLogImageFilter( const Self& ) ; //purposely not implemented
  void operator=( const Self& ) ; //purposely not implemented

};

} // end namespace itk


#endif
