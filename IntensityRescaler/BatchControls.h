#ifndef BatchControls_h
#define BatchControls_h
  
#include "Batch.h"
#include "ImageIntensityNormalizer.h"
#include <qlineedit.h>
#include <qlistview.h>

class BatchControls  : Batch
{
public:

  BatchControls( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
  ~BatchControls();
  void Show();
private:
  void Rescale();
  void SelectSource();
  void SelectTargetImage();
  void SelectTargetSegmentation();
  void SelectOutputDir();
  void AddLabel();
  void RemoveLabel();

  ImageIntensityNormalizer *m_rescaler;
};

#endif

