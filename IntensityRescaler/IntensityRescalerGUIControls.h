#ifndef ImageViewerGUIControls_h
#define ImageViewerGUIControls_h

#include "QtImageWriter.h"

#include "itkImage.h"
#include "QtWindow2D.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "ImageIntensityNormalizer.h"
#include "IntensityRescalerGUI.h"
#include "ImageManager.h"
#include "QtImageLoader.h"
#include "QtImageWriter.h"
#include "IntensityCurveVTK.h"
#include "IntensityCurveBox.h"
#include "BatchControls.h"
#include "HistoGUI.h"

class IntensityRescalerGUIControls : public IntensityRescalerGUI
{ 
    Q_OBJECT
public:
    
  IntensityRescalerGUIControls( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
  ~IntensityRescalerGUIControls();

  typedef unsigned short PixelType;
  typedef unsigned short OverlayPixelType;
  typedef itk::Image<PixelType,3> ImageType;
  typedef itk::Image<OverlayPixelType,3> OverlayType;
  typedef ImageType::Pointer ImagePointer;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;

  void LoadImg();
  void LoadSeg();
  void ChangeSliceX(int value);
  void ChangeSliceY(int value);
  void ChangeSliceZ(int value);
  void UnloadImage();
  void UnloadOverlay();

  void SelectSourceImage(int value);
  void SelectTargetImage(int value);
  void ChangeImageOverlay(int value);
  void ChangeLabelOverlay(int value);
  void SelectTargetOverlay(int value);
  void SelectSourceOverlay(int value);

  void Compute();
  void Batch();

  void ShowHistoSource();
  void ShowHistoTarget();
  int* ComputeHistogram(ImagePointer image,int* _max,int* _maxy);
  int* ComputeDistribution(ImagePointer image,ImagePointer seg,int* _maxx,int* _maxy,int label);
  void SaveImage();

public slots:
   void IntensityMouseMove();

private:
   ImageManager<PixelType,OverlayPixelType>* m_imagemanager;
   IntensityCurveVTK::Pointer m_intensitycurve;

   ImageType::Pointer m_image1;
   ImageType::Pointer m_image2;
   ImageType::Pointer m_overlay1;
   ImageType::Pointer m_overlay2;

   int m_max1;
   int m_max2;
   int m_min1;
   int m_min2;

   std::vector<ImagePointer> m_SourceImage;
   std::vector<ImagePointer> m_SourceOverlay;
   ImageIntensityNormalizer* m_imageIntensityNormalizer;
   std::vector<int> m_labellist;
   int* m_correspondingarray;
   BatchControls*   m_batch;
};

#endif
