#include "MeshIO.h"

MeshIO::MeshIO ()
{
  this->m_fileName = "" ;
  this->m_fileExtension = "" ;
  this->m_attributes = NULL ;
  this->m_mesh = NULL ;
  this->m_vtkPolyDataReader = NULL ;
  this->m_vtkFSSurfaceReader = NULL ;
}

MeshIO::~MeshIO ()
{
  if ( this->m_vtkPolyDataReader ) 
    this->m_vtkPolyDataReader->Delete () ;
  if ( this->m_vtkFSSurfaceReader ) 
    this->m_vtkFSSurfaceReader->Delete () ;
  /*if ( this->m_attributes )
    m_attributes->Delete () ;
  if ( this->m_mesh ) 
    m_mesh->Delete () ;
  */
}

void MeshIO::ProcessFileName ( std::string fileName ) 
{
  this->m_fileName = fileName ;
  this->m_fileExtension = "" ;

  std::string infile_str ( fileName ) ;

  int lastPoint = infile_str.rfind( '.' );
  if ( lastPoint > 0 )
    {
      this->m_fileExtension = infile_str.substr( lastPoint );
    }
}

bool MeshIO::Read ( std::string fileName )
{
  this->ProcessFileName ( fileName ) ;

  if ( this->ReadMesh () )
    return true ;
  if ( this->ReadAttributes () )
    return true ;
  return false ;
}

bool MeshIO::ReadMesh ()
{
  if ( this->m_fileExtension.compare ( ".vtk" ) == 0 )
    {
      if ( this->m_vtkPolyDataReader ) 
	this->m_vtkPolyDataReader->Delete () ;
      this->m_vtkPolyDataReader = vtkPolyDataReader::New() ;
      this->m_vtkPolyDataReader->SetFileName ( this->m_fileName.c_str() ) ;
      this->m_vtkPolyDataReader->Update () ;
      this->m_mesh = this->m_vtkPolyDataReader->GetOutput () ;
      return true ;
    }
  else if ( this->m_fileExtension.compare ( ".meta" ) == 0 ) 
    {
      typedef itk::DefaultDynamicMeshTraits < double, 3, 3, double, double > MeshTraitsType ;
      typedef itk::Mesh < double, 3, MeshTraitsType > itkMeshType ;
      typedef itk::MeshSpatialObject < itkMeshType > itkMeshSOType ;
      typedef itk::MetaMeshConverter < 3, double, MeshTraitsType > MeshConverterType ;

      MeshConverterType *itkConverter = new MeshConverterType () ;
      itkMeshSOType::Pointer meshSO = itkConverter->ReadMeta ( this->m_fileName.c_str() ) ;
      itkMeshType::Pointer mesh = meshSO->GetMesh () ;
      delete ( itkConverter ) ;

      itkMeshTovtkPolyData * ITK_VTK_Converter = new itkMeshTovtkPolyData ;
      ITK_VTK_Converter->SetInput ( mesh ) ;
      this->m_mesh = ITK_VTK_Converter->GetOutput () ;
      mesh->DisconnectPipeline () ;
      delete ( ITK_VTK_Converter ) ;
      return true ;
    }

  else if ( ( this->m_fileExtension.compare ( ".white" ) == 0 )  
	    || ( this->m_fileExtension.compare ( ".inflated" ) == 0 )
	    || ( this->m_fileExtension.compare ( ".orig" ) == 0 )
	    || ( this->m_fileExtension.compare ( ".pial" ) == 0 )
	    || ( this->m_fileExtension.compare ( ".nofix" ) == 0 ) )

    {
      if ( this->m_vtkFSSurfaceReader )
	this->m_vtkFSSurfaceReader->Delete () ;
      this->m_vtkFSSurfaceReader = vtkFSSurfaceReader::New () ;
      this->m_vtkFSSurfaceReader->SetFileName ( this->m_fileName.c_str() ) ;
      this->m_vtkFSSurfaceReader->Update () ;
      this->m_mesh = this->m_vtkFSSurfaceReader->GetOutput () ;
      return true ;
    }
  return false ;
}

bool MeshIO::ReadAttributes () 
{
  if ( this->m_fileExtension.compare ( ".txt" ) == 0 )
    {
      // KWMeshVisu style scalars
      KWMeshVisuAttributeIO *attributeReader = new ( KWMeshVisuAttributeIO ) ;
      attributeReader->SetFileName ( this->m_fileName.c_str() ) ;
      attributeReader->ReadAttributes () ;
      this->m_attributes = attributeReader->GetAttributes () ;
      delete ( attributeReader ) ;
      return true ;
    }
  else if ( ( this->m_fileExtension.compare ( ".sulc" ) == 0 ) 
	    || ( this->m_fileExtension.compare ( ".curv" ) == 0 ) 
	    || ( this->m_fileExtension.compare ( ".area" ) == 0 ) 
	    || ( this->m_fileExtension.compare ( ".thickness" ) == 0 ) )
    {
      // FreeSurfer style scalars
      this->m_attributes = vtkFloatArray::New () ;
      vtkFSSurfaceScalarReader *fsScalarReader = vtkFSSurfaceScalarReader::New () ;
      fsScalarReader->SetFileName ( this->m_fileName.c_str() ) ;
      fsScalarReader->SetOutput ( this->m_attributes ) ;
      fsScalarReader->ReadFSScalars() ; 
      this->m_attributes = fsScalarReader->GetOutput () ;
      fsScalarReader->Delete () ;
      return true ;
    }

  return false ;
}

void MeshIO::WriteAttributes ( vtkFloatArray *attributes, std::string fileName ) 
{
  this->ProcessFileName ( fileName ) ;

  if ( this->m_fileExtension.compare ( ".txt" ) == 0 )
   {
     // KWMeshVisu style scalars
     KWMeshVisuAttributeIO *attributeWriter = new ( KWMeshVisuAttributeIO ) ;
     attributeWriter->SetFileName ( this->m_fileName.c_str() ) ;
     attributeWriter->SetAttributes ( attributes ) ;
     attributeWriter->WriteAttributes () ;
     delete ( attributeWriter ) ;
     return ;
   }
}
