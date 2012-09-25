#ifndef QtSlicer_h
#define QtSlicer_h

#include <QtSlicerGUI.h>
#include <itkImage.h>
#include <QtSlicerHelpGUI.h>

class QtSlicer : public Gui
{ 
public:
    
  QtSlicer( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
  ~QtSlicer();
  
  typedef itk::Image<double,3> ImageType;
  typedef itk::Image<unsigned short,3> OverlayType;


  void DisplayPosition(int x,int y ,int z,float value);
  void Help();
  void SetInputImage(ImageType * newImData);
  void SetInputOverlay(OverlayType * newImData);
  void QtSlicer::DisplaySliceNumber(int number);
  void QtSlicer::DisplayIMin(int value);
  void QtSlicer::DisplayIMax(int value);
  void LoadImg();
  void LoadSeg();
};

#endif
