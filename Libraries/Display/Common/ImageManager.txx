#ifndef _ImageManager_TXX
#define _ImageManager_TXX

//#include "ImageManager.h"

template< class TPixelImage, class TPixelOverlay >
ImageManager<TPixelImage, TPixelOverlay>
::ImageManager()
{
  m_IsSourceOverlay = false;
  m_IsTargetOverlay = false;
  m_IsSourceImage = false;
  m_IsTargetImage = false;
  m_issourcegraylabel = false;
  m_issourceseglabel = false;
  m_istargetgraylabel = false;
  m_istargetseglabel = false;
}

template< class TPixelImage, class TPixelOverlay >
ImageManager<TPixelImage, TPixelOverlay>
::~ImageManager()
{
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetSourceImage(ImageType* image)
{
  m_SourceImage = image;
  m_IsSourceImage = true;
    for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->SetInputImage(image);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetTargetImage(ImageType* image)
{
  m_TargetImage = image;
  m_IsTargetImage = true;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->SetInputImage2(image);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetSourceOverlay(OverlayType* image)
{
  m_SourceOverlay = image;
  m_IsSourceOverlay = true;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->SetInputOverlay(image);
}


template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetTargetOverlay(OverlayType* image)
{
  m_TargetOverlay = image;
  m_IsTargetOverlay = true;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->SetInputOverlay2(image);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::UnSetSourceImage()
{
  m_IsSourceImage = false;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->UnSetInputImage();
}


template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::UnSetSourceOverlay()
{
  m_IsSourceOverlay = false;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->UnSetInputOverlay();
}


template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::UnSetTargetImage()
{
  m_IsTargetImage = false;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->UnSetInputImage2();
}


template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::UnSetTargetOverlay()
{
  m_IsTargetOverlay = false;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
        (*it)->UnSetInputOverlay2();
}


template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::AddWindow2D(QtWindow2D* window2D)
{
  m_window2Dlist.push_back(window2D);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::DelWindow2D(QtWindow2D* window2D)
{
  std::vector<QtWindow2D*>::iterator it = m_window2Dlist.begin();
  for (int i=0;i<(int)m_window2Dlist.size();i++)
  {
    if (m_window2Dlist[i] == window2D)
      m_window2Dlist.erase(it);
    it++;
  }
}



template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::Update()
{


}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::ChangeSliceX(int value)
{
/*for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
    {
    if ((*it)->GetId() == 0)
      (*it)->ChangeSlice(value);
  }*/
  SetCrosshair(m_crosshairx,m_crosshairy, value);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::ChangeSliceY(int value)
{
 /*for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
    {
    if ((*it)->GetId() == 1)
      (*it)->ChangeSlice(value);
  }*/
  SetCrosshair(m_crosshairx,value, m_crosshairz);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::ChangeSliceZ(int value)
{
 /*for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
    {
    if ((*it)->GetId() == 2)
      (*it)->ChangeSlice(value);
  }*/
   SetCrosshair(value,m_crosshairy, m_crosshairz);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::ChangeAlpha(int value)
{
 m_alpha = value;
 for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
    {
      (*it)->ChangeAlpha(value);

  }
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::ChangeLabelOverlay(int value)
{
 for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
    {
      (*it)->ChangeLabelOverlay(value);
  }
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetCrosshair(unsigned int x,unsigned int y,unsigned int z)
{
  SizeType m_imagesize=GetImageSize();
  if ((m_imagesize[0] == 0) && (m_imagesize[1] == 0) && (m_imagesize[2] == 0))
    return;

  if (x>=m_imagesize[0])
    x = m_imagesize[0]-1;

  if (y>=m_imagesize[1])
    y = m_imagesize[1]-1;

  if (z>=m_imagesize[2])
    z = m_imagesize[2]-1;

  if (x<0) x = 0;
  if (y<0) y = 0;
  if (z<0) z = 0;
  
  if (x>=m_imagesize[0]) x = m_imagesize[0]-1;
  if (y>=m_imagesize[1]) y = m_imagesize[1]-1;
  if (z>=m_imagesize[2]) z = m_imagesize[2]-1;

  m_crosshairx = x;
  m_crosshairy = y;
  m_crosshairz = z;


  //Update Label
  itk::Index<3> m_index;
  m_index[0] = m_crosshairx;
  m_index[1] = m_crosshairy;
  m_index[2] = m_crosshairz;

  int m_sourcegrayvalue=-1;
  int m_sourcesegvalue=-1; 
  int m_targetgrayvalue=-1;
  int m_targetsegvalue=-1; 

  if (m_IsSourceImage)
    m_sourcegrayvalue = m_SourceImage->GetPixel(m_index);

  if (m_IsSourceOverlay)
    m_sourcesegvalue = m_SourceOverlay->GetPixel(m_index);

  if (m_IsTargetImage)
    m_targetgrayvalue = m_TargetImage->GetPixel(m_index);

  if (m_IsTargetOverlay)
    m_targetsegvalue = m_TargetOverlay->GetPixel(m_index);

 for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
    {
    if ((*it)->GetId() == 0)
    {
      (*it)->SetCrosshair(x,y);
      (*it)->ChangeSlice(z);
    }
    if ((*it)->GetId() == 1)
    {
      (*it)->SetCrosshair(x,z);
      (*it)->ChangeSlice(y);
    }
    if ((*it)->GetId() == 2)
    {
      (*it)->SetCrosshair(y,z);
      (*it)->ChangeSlice(x);
    }
  }


  UpdateLabel(m_sourcegrayvalue,m_sourcesegvalue,m_targetgrayvalue,m_targetsegvalue);
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::GetCrosshair(int *x,int* y, int* z)
{
  *x = m_crosshairx;
  *y = m_crosshairy;
  *z = m_crosshairz;

}

template< class TPixelImage, class TPixelOverlay>
void 
ImageManager<TPixelImage, TPixelOverlay>
::SetZoomFactor(float zoom)
{
  m_zoom = zoom;
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
  {
   (*it)->SetZoomFactor(m_zoom);
  }
}
  
template< class TPixelImage, class TPixelOverlay>
float 
ImageManager<TPixelImage, TPixelOverlay>
::GetZoomFactor()
{
  return m_zoom;
}



template< class TPixelImage, class TPixelOverlay>
void 
ImageManager<TPixelImage, TPixelOverlay>
::SetViewPosition(float posx,float posy)
{
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
  {
   (*it)->SetViewPosition(posx,posy);
  }
}

template< class TPixelImage, class TPixelOverlay>
void 
ImageManager<TPixelImage, TPixelOverlay>
::SetIntensityMin(int value)
{
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
  {
   (*it)->SetIntensityMin(value);
  }
}

template< class TPixelImage, class TPixelOverlay>
void 
ImageManager<TPixelImage, TPixelOverlay>
::SetIntensityMax(int value)
{
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
  {
   (*it)->SetIntensityMax(value);
  }
}


template< class TPixelImage, class TPixelOverlay>
void 
ImageManager<TPixelImage, TPixelOverlay>
::SetBlendingMode(int mode)
{
  for(std::vector<QtWindow2D*>::iterator it=m_window2Dlist.begin();it!=m_window2Dlist.end();++it)
  {
   (*it)->SetBlendingMode(mode);
  }
}



template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::UpdateCrosshair()
{
  SetCrosshair(m_crosshairx,m_crosshairy,m_crosshairz);  
}


template< class TPixelImage, class TPixelOverlay>
typename ImageManager<TPixelImage, TPixelOverlay>::SizeType
ImageManager<TPixelImage, TPixelOverlay>
::GetImageSize()
{
  if (m_IsSourceImage)
    return m_SourceImage->GetLargestPossibleRegion().GetSize();

  if (m_IsSourceOverlay)
    return m_SourceOverlay->GetLargestPossibleRegion().GetSize();

  if (m_IsTargetImage)
    return m_TargetImage->GetLargestPossibleRegion().GetSize();

  if (m_IsTargetOverlay)
    return m_TargetOverlay->GetLargestPossibleRegion().GetSize();

  itk::Size<3> m_size;
  m_size[0] = 0;
  m_size[1] = 0;  
  m_size[2] = 0;
  return m_size;
}


template< class TPixelImage, class TPixelOverlay>
typename ImageManager<TPixelImage, TPixelOverlay>::SpacingType
ImageManager<TPixelImage, TPixelOverlay>
::GetImageSpacing()
{
  return m_SourceImage->GetSpacing();
}

template< class TPixelImage, class TPixelOverlay>
typename ImageManager<TPixelImage, TPixelOverlay>::ImagePointer 
ImageManager<TPixelImage, TPixelOverlay>
::GetSourceImage()
{
  return m_SourceImage;
}

template< class TPixelImage, class TPixelOverlay>
typename ImageManager<TPixelImage, TPixelOverlay>::ImagePointer 
ImageManager<TPixelImage, TPixelOverlay>
::GetTargetImage()
{
 return m_TargetImage;
}

template< class TPixelImage, class TPixelOverlay>
typename ImageManager<TPixelImage, TPixelOverlay>::OverlayPointer 
ImageManager<TPixelImage, TPixelOverlay>
::GetSourceOverlay()
{
  return m_SourceOverlay;
}

template< class TPixelImage, class TPixelOverlay>
typename ImageManager<TPixelImage, TPixelOverlay>::OverlayPointer 
ImageManager<TPixelImage, TPixelOverlay>
::GetTargetOverlay()
{
  return m_IsTargetOverlay;
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetSourceGrayLabel(QLabel* graylabel)
{
  m_sourcegraylabel = graylabel;
  m_issourcegraylabel = true;
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetSourceSegLabel(QLabel* seglabel)
{
  m_sourceseglabel = seglabel;
  m_issourceseglabel = true;
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetTargetGrayLabel(QLabel* graylabel)
{
  m_targetgraylabel = graylabel;
  m_istargetgraylabel = true;
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::SetTargetSegLabel(QLabel* seglabel)
{
  m_targetseglabel = seglabel;
  m_istargetseglabel = true;
}

template< class TPixelImage, class TPixelOverlay>
void
ImageManager<TPixelImage, TPixelOverlay>
::UpdateLabel(int sourcegrayvalue,int sourcesegvalue,int targetgrayvalue,int targetsegvalue)
{
  if (m_issourcegraylabel)
  {
    if (sourcegrayvalue != -1)
      m_sourcegraylabel->setText(QString("%1").arg(sourcegrayvalue));
    else
      m_sourcegraylabel->setText("Na");
  }

  if (m_issourceseglabel)
  {
    if (sourcesegvalue != -1)
      m_sourceseglabel->setText(QString("%1").arg(sourcesegvalue));
    else
      m_sourceseglabel->setText("Na");
  }

  if (m_istargetgraylabel)
  {
    if (targetgrayvalue != -1)
      m_targetgraylabel->setText(QString("%1").arg(targetgrayvalue));
    else
      m_targetgraylabel->setText("Na");
  }

  if (m_istargetseglabel)
  {
    if (targetsegvalue != -1)
      m_targetseglabel->setText(QString("%1").arg(targetsegvalue));
    else
      m_targetseglabel->setText("Na");
  }
}

#endif
