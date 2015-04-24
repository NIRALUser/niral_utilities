#ifndef _itkCorrespondenceEvaluator_h
#define _itkCorrespondenceEvaluator_h

#include "itkObject.h"
#include "vnl/vnl_matrix.h"

namespace itk
{
  template<class TOutputMeshType>
    class CorrespondenceEvaluator : public Object
  {

  public:

    /** Standard class typedefs. */
    typedef CorrespondenceEvaluator           Self;
    typedef Object                            Superclass;
    typedef SmartPointer<Self>                Pointer;
    typedef SmartPointer<const Self>          ConstPointer;

    /** Convenient typedefs. */
    typedef TOutputMeshType                   MeshType;
    typedef typename MeshType::Pointer        MeshPointer;
    typedef std::vector < MeshPointer >       MeshPointerArray;
    typedef vnl_matrix<double>                MatrixType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(CorrespondenceEvaluator, Object);

    double GetGeneralization(unsigned int numberOfShapeParameters, double &stdError) ;
    double GetSpecificity(unsigned int N, unsigned int numberOfShapeParameters, double &stdError) ;

    /** Sets the number of meshes to be analyzed. */
    void SetNumberOfInputs( unsigned int number );

    /** Returns the number of analyzed meshes. */
    itkGetConstMacro( NumberOfInputs, unsigned int );

    /** Sets the maximum number of allowed modes for the PCA analysis */
    itkGetConstMacro( MaxModes, unsigned int );

    /** Sets the specified parameterized input mesh. */
    void SetInput( unsigned int idx, MeshPointer input );



    /** Sets the random number distribution to be used: gaussian(true) or uniform(false) **/
    void SetDistribution ( bool gaussian )
    {
      gaussianDistribution = gaussian ;
    }

  protected:

    CorrespondenceEvaluator();
    ~CorrespondenceEvaluator();

    MeshPointerArray m_Meshes ;
    unsigned int m_NumberOfInputs ;
    unsigned int m_MaxModes ;
    bool gaussianDistribution ;

    double ComputeGeneralizationError ( unsigned int i, unsigned int M ) ;
    double GenerateUniformRandomNumber () ;
    double GenerateGaussianRandomNumber () ;

    double SurfaceDistance ( MatrixType *from, MatrixType *to, unsigned int nPoints ) ;
    void WriteSurfaceForDebug (MatrixType *surface, unsigned int nPoints, std::string filename) ;
    void BringInputToOrigin (unsigned int id) ;
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCorrespondenceEvaluator.txx"
#endif

#endif
