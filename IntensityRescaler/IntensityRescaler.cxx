#include "ImageManager.h"
#include "QtImageLoader.h"
#include "IntensityRescalerGUIControls.h"
#include "IntensityRescalerCommandLine.h"
#include <iostream>
#include <qstring.h>
#include <qfiledialog.h>
#include <qapplication.h>
#include <qplatinumstyle.h>
  
int main( int argc, char* argv[] ) 
{
  if ((argc==2) && (!strcmp(argv[1],"-version")))
  {
    std::cout << "IntensityRescaler 1.0 - Compiled on: " << __DATE__ << " - " __TIME__  << std::endl;
    return 0;
  }

  if (argc < 2)
  {
    QApplication myApp( argc, argv );
    
    IntensityRescalerGUIControls m_GUI( 0, 0, TRUE );
      myApp.setMainWidget(&m_GUI);

    m_GUI.setCaption( "Intensity Rescaler 1.0 (Qt)" );
    myApp.setStyle( new QPlatinumStyle );
    QPalette p( QColor( 239, 239, 239 ) );
    myApp.setPalette( p, TRUE );

    m_GUI.show();
    myApp.exec();
    return 0;
  }
  else
  {
     IntensityRescalerCommandLine m_commandline;
     if ((argc>1) && ((!strcmp(argv[1],"-create")) || (!strcmp(argv[1],"-c"))))
     {
       if (argc>2)
       {
          m_commandline.Create(argv[2]);
       }
       else
       {
         m_commandline.Create("IRescalerExample.irs"); 
       }
     }
     else
     if ((argc>2) && ((!strcmp(argv[1],"-input")) || (!strcmp(argv[1],"-i"))))
     {
       m_commandline.Load(argv[2]);
       
       if ((argc==4) && (!strcmp(argv[3],"-v")))
       {
         m_commandline.DisplayOptions();
         m_commandline.Run(true);
       }
       else
        m_commandline.Run();
     }
     else
     {
       std::cout << "IntensityRescaler 1.0 options:" << std::endl;
       std::cout << "-create (-c) [filename]: Create a example script file" << std::endl;
       std::cout << "-input (-i) <filename> [-v]: Execute script <filename> [verbose mode]" << std::endl;
       std::cout << std::endl;
     }
  }

  return 0;
}

