#ifndef _ImageManager_H
#define _ImageManager_H

#include "itkImage.h"
#include <qlabel.h>
#include "QtWindow2D.h"

class QtWindow2D;

template<class TPixelImage, class TPixelOverlay>
class ImageManager
{
  public:

  typedef TPixelImage PixelImageType;
  typedef itk::Image<PixelImageType,3> ImageType;
  typedef TPixelOverlay PixelOverlayType;
  typedef itk::Image<PixelOverlayType,3> OverlayType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename OverlayType::Pointer OverlayPointer;
      typedef typename ImageType::SizeType SizeType;
      typedef typename ImageType::SpacingType SpacingType;

  ImageManager();
  ~ImageManager(); 

  void SetSourceImage(ImageType* image);
  void SetTargetImage(ImageType* image);
  void SetSourceOverlay(OverlayType* image);
  void SetTargetOverlay(OverlayType* image);

  void UnSetSourceImage();
  void UnSetTargetImage();
  void UnSetSourceOverlay();
  void UnSetTargetOverlay();

  void SetSourceGrayLabel(QLabel* graylabel);
  void SetSourceSegLabel(QLabel* seglabel);
  void SetTargetGrayLabel(QLabel* graylabel);
  void SetTargetSegLabel(QLabel* seglabel);

  void UpdateLabel(int sourcegrayvalue,int sourcesegvalue,int targetgrayvalue,int targetsegvalue);
  void SetBlendingMode(int mode);

  ImagePointer GetSourceImage();
  ImagePointer GetTargetImage();
  OverlayPointer GetSourceOverlay();
  OverlayPointer GetTargetOverlay();

  void GetCrosshair(int *x,int* y, int* z);
  void SetCrosshair(unsigned int x,unsigned int y,unsigned int z);
  void SetZoomFactor(float zoom);
  float GetZoomFactor();

  void SetViewPosition(float posx,float posy);

  void UpdateCrosshair();
  SizeType   GetImageSize();
  SpacingType  GetImageSpacing();

  bool IsSourceImage();
  bool IsTargetImage();
  bool IsSourceOverlay();
  bool IsTargetOverlay();

  void AddWindow2D(QtWindow2D* window2D);
  void DelWindow2D(QtWindow2D* window2D);
  void ChangeSliceX(int value);
  void ChangeSliceY(int value);
  void ChangeSliceZ(int value);
  void ChangeAlpha(int value);
  void ChangeLabelOverlay(int value);
  void SetIntensityMax(int value);
  void SetIntensityMin(int value);

  void Update();

  protected:
  ImagePointer m_SourceImage;
  OverlayPointer m_SourceOverlay;
  ImagePointer m_TargetImage;
  OverlayPointer m_TargetOverlay;
  bool m_IsSourceOverlay;
  bool m_IsTargetOverlay;
  bool m_IsSourceImage;
  bool m_IsTargetImage;

  private:
  std::vector<QtWindow2D*> m_window2Dlist;
  int m_crosshairx;
  int m_crosshairy;
  int m_crosshairz;
  QLabel* m_sourcegraylabel;
  QLabel* m_sourceseglabel;
  QLabel* m_targetgraylabel;
  QLabel* m_targetseglabel;
  bool m_issourcegraylabel;
  bool m_issourceseglabel;
  bool m_istargetgraylabel;
  bool m_istargetseglabel;
  int m_alpha;
  float m_zoom;



};

#include "ImageManager.txx"

#endif
