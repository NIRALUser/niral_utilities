#include <qapplication.h>
#include "QtGlSliceView.h"
#include "QtSlicer.h"
#include <qfiledialog.h>
#include <qslider.h>

#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <qwindowsstyle.h>
#include <qplatinumstyle.h>
#include <qmotifstyle.h>
#include <qmotifplusstyle.h>
#include <qcdestyle.h>

int main( int argc, char* argv[] ) 
{

  QApplication myApp( argc, argv );

  QtSlicer m_GUI( 0, 0, TRUE );
  myApp.setMainWidget(&m_GUI);

  m_GUI.setCaption( "Insight Qt Slicer" );
  myApp.setStyle( new QPlatinumStyle );
  QPalette p( QColor( 239, 239, 239 ) );
  myApp.setPalette( p, TRUE );





  try
    {
    m_GUI.exec();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception during GUI execution" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
 
  return 0;

}

