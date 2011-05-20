#include <string>
#include <fstream>
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkFSSurfaceReader.h"
#include "vtkFSSurfaceScalarReader.h"
#include "KWMeshVisuAttributeIO.h"
#include <itkDefaultDynamicMeshTraits.h>
#include <itkMetaMeshConverter.h>
#include "itkMeshTovtkPolyData.h"

class MeshIO
{
 public: 
  MeshIO () ;
  ~MeshIO () ;

  //void SetAttributes ( vtkFloatArray *attributes )
  //{
  //  this->m_attributes = attributes ;
  //}

  vtkFloatArray *GetAttributes ()
    {
      return this->m_attributes ;
    } 

  vtkPolyData *GetMesh () 
    {
      vtkPolyData *result = vtkPolyData::New () ;
      result->DeepCopy ( m_mesh ) ;
      return result ;
    }

  bool Read ( std::string fileName ) ;

  void WriteAttributes ( vtkFloatArray *attributes, std::string fileName ) ;

 protected:
  void ProcessFileName ( std::string fileName ) ;
  bool ReadAttributes () ;
  //void WriteAttributes () ;
  bool ReadMesh () ;

 private:
  std::string m_fileExtension ;
  std::string m_fileName ;
  vtkFloatArray *m_attributes ;
  vtkPolyData *m_mesh ;
  vtkPolyDataReader *m_vtkPolyDataReader ;
  vtkFSSurfaceReader *m_vtkFSSurfaceReader ;
};
