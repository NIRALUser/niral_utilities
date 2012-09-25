#include "QtWindow2D.h"
#include "ImageManager.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "Palette.h"
#include "GL/glut.h"

/**
 * This definition is needed to use RGBA pixels for compilation
 */
typedef itk::RGBAPixel<unsigned char> ColorPixel;

namespace itk {
template<>
class NumericTraits<ColorPixel>
{
public:
  typedef ColorPixel ValueType;
  typedef ColorPixel PrintType;
  typedef ColorPixel AbsType;
  typedef ColorPixel AccumulateType;
  static const ColorPixel Zero;
  static const ColorPixel One;

  static ColorPixel NonpositiveMin() { return Zero; }
  static bool IsPositive(ColorPixel val) { return true; }
  static bool IsNonpositive(ColorPixel val) { return false; }
  static bool IsNegative(ColorPixel val) { return false; }
  static bool IsNonnegative(ColorPixel val) {return true; }
private:

  static const unsigned char ZeroArray[4];
  static const unsigned char OneArray[4];
};

} // End of namespace

const unsigned char itk::NumericTraits<ColorPixel>::ZeroArray[4] = {0,0,0,0};
const ColorPixel itk::NumericTraits<ColorPixel>::Zero = 
  ColorPixel(itk::NumericTraits<ColorPixel>::ZeroArray);

const unsigned char itk::NumericTraits<ColorPixel>::OneArray[4] = {1,1,1,1};
const ColorPixel itk::NumericTraits<ColorPixel>::One = 
  ColorPixel(itk::NumericTraits<ColorPixel>::OneArray);


QtWindow2D::QtWindow2D( QWidget *parent, const char *name)
: QGLWidget(parent, name)
{
  cValidOverlayData     = false;
  cViewOverlayData      = false;
  cViewOverlayCallBack  = NULL;
  cOverlayOpacity       = 0.5;
  cWinOverlayData       = NULL;
  cViewAxisLabel = true;
  cViewOverlayData = true;
  cViewDetails = true;
  cViewCrosshairs = true;
  cViewValue = true;
  cClickMode = CM_SELECT;
  cClickSelect[0] = 0;
  cClickSelect[1] = 0;
  cClickSelect[2] = 0;

  cWinOffset = 0;
  view_pos_x = 0;
  view_pos_y = 0;
  view_zoom = 1;
  m_alpha = 100;
  m_overlaymax = 0;
  cColorTable = ColorTableType::New();
  cColorTable->useDiscrete();
  cW = 0;
  cH = 0;
  for(unsigned int i=0;i<3;i++)
  {
    cFlipX[i]=false;
    cFlipY[i]=false;
    cFlipZ[i]=false;
  }
//  cWinImData = NULL;
//  cWinZBuffer = NULL;

  m_isslicer = false;
  m_islabel = false;
  cValidImData = false;
  cValidImData2 = false;
  cValidOverlayData2 = false;

  m_isImage1 = false;
  m_isImage2 = false;
  m_isOverlay1 = false;
  m_isOverlay2 = false;
  m_id=0;
  m_rotation = 0;
  cShowCrosshair = true;
  m_Color = Palette::Generate(Palette::PALETTE_RGB);

  m_drawimagename=false;
  m_drawposition = false;
  m_drawintensity = false;
  m_drawimageinfo = false;
  m_overlayvisible = true;
  m_imagename = "";
  m_imageinfotype = "";
  // Initialize the Grey slice texture
  m_LabelTexture.SetGlComponents(4);
  m_LabelTexture.SetGlFormat(GL_RGBA);
  m_blendingmode = 0;

  m_manualintensity = false;
  m_crosshairvisible = true;
  m_displayoverlayzero = false;
  m_x = 0;
  m_y = 1;

  cWinZoom = 1;
  cWinOrientation = m_id;
  cWinOrder[0] = 2;
  cWinOrder[1] = 1;
  cWinOrder[2] = 0;
  cDimSize[0] = 0;
  cDimSize[1] = 0;
  cDimSize[2] = 0;
  cWinCenter[0] = 0;
  cWinCenter[1] = 0;
  cWinCenter[2] = 0;
  cSpacing[0]= 1.0;
  cSpacing[1]= 1.0;

  m_imageinfosize[0] = 0;
  m_imageinfosize[1] = 0;
  m_imageinfosize[2] = 0;
  
  m_imageinfopixdim[0] = 0.0;
  m_imageinfopixdim[1] = 0.0;
  m_imageinfopixdim[2] = 0.0;
}
  
  
// QtWindow2D::
// QtWindow2D( QGLFormat glf, QWidget *parent, const char *name)
// : QGLWidget(glf,parent, name)
// {    
//   cValidOverlayData     = false;
//   cViewOverlayData      = false;
//   cViewOverlayCallBack  = NULL;
//   cOverlayOpacity       = 0.5;
//   cWinOverlayData       = NULL;
//   cColorTable = ColorTableType::New();
//   cColorTable->useDiscrete();
  
//   // Initialize the Grey slice texture
// }
  

QtWindow2D::~QtWindow2D()
{
  delete[] m_Color;
}
  
void 
QtWindow2D::
SetInputImage(ImageType * newImData)
{
  cImData = newImData;

  typedef MinimumMaximumImageCalculator<ImageType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();

  calculator->SetImage( cImData );
  calculator->Compute();
   
  if (!m_manualintensity)
  {
    cIWMin      = calculator->GetMinimum();
    cIWMax      = calculator->GetMaximum();
  }
    
  if( !newImData )    
  {
    return;
  }

  RegionType region = newImData->GetLargestPossibleRegion();
  if( region.GetNumberOfPixels() == 0 ) 
  {
    return;
  }

  SizeType   size   = region.GetSize();
  if( cValidOverlayData )
  {
    RegionType overlay_region = cOverlayData->GetLargestPossibleRegion();
    SizeType   overlay_size   = overlay_region.GetSize();
      
    for( int i=0; i<3; i++ )
    {
      if( size[i] != overlay_size[i] )
      {
        return;
      }
    }
  } 


  cDimSize[0]=size[0];
  cDimSize[1]=size[1];
  cDimSize[2]=size[2];

  cIWModeMin  = IW_MIN;
  cIWModeMax  = IW_MAX;
    
  cImageMode = IMG_VAL;
    
  cWinCenter[0] = cDimSize[0]/2;
  cWinCenter[1] = cDimSize[1]/2;
  cWinCenter[2] = 0;
  
  cWinZoom = 1;
  cWinOrientation = m_id;
  cWinOrder[0] = 2;
  cWinOrder[1] = 1;
  cWinOrder[2] = 0;

  cWinSizeX = cDimSize[0];
  cWinOffset = cDimSize[0]%2;
  cWinMinX  = 0;
  cWinMinY  = 0;
  if (m_id == 0)
  {
  cWinDataSizeX = cDimSize[0];
  cWinDataSizeY = cDimSize[1];
  cWinMaxX  = cDimSize[0] - 1;  
  cWinMaxY  = cDimSize[1] - 1;
  numslice = cDimSize[2]/2;
  cSpacing[0]= cImData->GetSpacing()[0];
  cSpacing[1]= cImData->GetSpacing()[1];
  }

  if (m_id == 1)
  {
  cWinDataSizeX = cDimSize[0];
  cWinDataSizeY = cDimSize[2];
  cWinMaxX  = cDimSize[0] - 1;  
  cWinMaxY  = cDimSize[2] - 1;
  numslice = cDimSize[1]/2;
  cSpacing[0]= cImData->GetSpacing()[0];
  cSpacing[1]= cImData->GetSpacing()[2];
  }

  if (m_id == 2)
  {
  cWinMaxX  = cDimSize[1] - 1;  
  cWinMaxY  = cDimSize[2] - 1;
  cWinDataSizeX = cDimSize[1];
  cWinDataSizeY = cDimSize[2];
  numslice = cDimSize[0]/2;
  cSpacing[0]= cImData->GetSpacing()[1];
  cSpacing[1]= cImData->GetSpacing()[2];
  }

  cViewImData  = true;
  cValidImData = true;

  m_isImage1 = true;
  
  this->update();

  if (m_isslicer)
  {
    if (m_id == 0)
    {
    cWinMaxZ = newImData->GetLargestPossibleRegion().GetSize()[2]-1;
    m_slicer->setMaxValue(newImData->GetLargestPossibleRegion().GetSize()[2]-1);
    m_slicer->setValue(newImData->GetLargestPossibleRegion().GetSize()[2]/2);
    }

    if (m_id == 1)
    {
    cWinMaxZ = newImData->GetLargestPossibleRegion().GetSize()[1]-1;
    m_slicer->setMaxValue(newImData->GetLargestPossibleRegion().GetSize()[1]-1);
    m_slicer->setValue(newImData->GetLargestPossibleRegion().GetSize()[1]/2);
    }

    if (m_id == 2)
    {
    cWinMaxZ = newImData->GetLargestPossibleRegion().GetSize()[0]-1;
    m_slicer->setMaxValue(newImData->GetLargestPossibleRegion().GetSize()[0]-1);
    m_slicer->setValue(newImData->GetLargestPossibleRegion().GetSize()[0]/2);
    }

  }
   update();
   this->updateGL();
}

void 
QtWindow2D::
ChangeId(int id)
{
  m_id = id;
  cWinOrientation = m_id;

  if (m_id == 0)
  {
  m_x = 0;
  m_y = 1;
  }

  if (m_id == 1)
  {
  m_x = 0;
  m_y = 2;
  }

  if (m_id == 2)
  {
  m_x = 1;
  m_y = 2;
  }


  cWinMaxX  = cDimSize[m_x] - 1;  
  cWinMaxY  = cDimSize[m_y] - 1;
  cWinDataSizeX = cDimSize[m_x];
  cWinDataSizeY = cDimSize[m_y];
  cSpacing[0]= cImData->GetSpacing()[m_x];
  cSpacing[1]= cImData->GetSpacing()[m_y];


 /* if(cWinImData != NULL)
  {
    delete [] cWinImData;
    delete [] cWinRotData;
  }
    
  cWinImData = new unsigned char[ cWinDataSizeX * cWinDataSizeY ];
  cWinRotData = new unsigned char[ cWinDataSizeX * cWinDataSizeY ];


   if(cWinOverlayData != NULL) 
   {
     delete [] cWinOverlayData;
   }
    
   cWinOverlayData = new unsigned char[ cWinDataSizeX * cWinDataSizeY * 4 ];


  if(cWinZBuffer != NULL) 
  {
    delete [] cWinZBuffer;
  }
    
  cWinZBuffer = new unsigned short[ cWinDataSizeX * cWinDataSizeY ];*/

  
  int x,y,z;
  m_manager->GetCrosshair(&x,&y,&z);
  int cross[3];
  cross[0] = x;
  cross[1] = y;
  cross[2] = z;
  //Set Cursor
  SetCrosshair(cross[m_x],cross[m_y]);


  if (m_isslicer)
  {
    if (m_id == 0)
    {
    cWinMaxZ = cImData->GetLargestPossibleRegion().GetSize()[2]-1;
    m_slicer->setMaxValue(cImData->GetLargestPossibleRegion().GetSize()[2]-1);
    ChangeSlice(z);
    }

    if (m_id == 1)
    {
    cWinMaxZ = cImData->GetLargestPossibleRegion().GetSize()[1]-1;
    m_slicer->setMaxValue(cImData->GetLargestPossibleRegion().GetSize()[1]-1);
    ChangeSlice(y);
    //m_slicer->setValue(y);
    }

    if (m_id == 2)
    {
    cWinMaxZ = cImData->GetLargestPossibleRegion().GetSize()[0]-1;
    m_slicer->setMaxValue(cImData->GetLargestPossibleRegion().GetSize()[0]-1);
    ChangeSlice(x);
    //m_slicer->setValue(x);
    }
  }

  //this->update();
  //this->updateGL();  
}

void 
QtWindow2D::
SetInputImage2(ImageType * newImData)
{
  if( !newImData )    
  {
    return;
  }

  RegionType region = newImData->GetLargestPossibleRegion();
  if( region.GetNumberOfPixels() == 0 ) 
  {
    return;
  }

  SizeType   size   = region.GetSize();
  if( cValidOverlayData )
  {
    RegionType overlay_region = cOverlayData->GetLargestPossibleRegion();
    SizeType   overlay_size   = overlay_region.GetSize();
      
    for( int i=0; i<3; i++ )
    {
      if( size[i] != overlay_size[i] )
      {
        return;
      }
    }
  } 

  cImData2 = newImData;
    
  typedef MinimumMaximumImageCalculator<ImageType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( cImData2 );
  calculator->Compute();
  cIWMin2      = calculator->GetMinimum();
  cIWMax2      = calculator->GetMaximum();
  
  //Select the alpha in the middle slice
 /* if (m_id == 0)
  {
    m_alpha = newImData->GetLargestPossibleRegion().GetSize()[1]/2;
  }

  if (m_id == 1)
  {
    m_alpha = newImData->GetLargestPossibleRegion().GetSize()[2]/2;
  }

  if (m_id == 2)
  {
    m_alpha = newImData->GetLargestPossibleRegion().GetSize()[2]/2;
  }*/

 cValidImData2 = true;
 m_isImage2 = true;
  this->update();
  this->updateGL();
}



const QtWindow2D::ImagePointer &
QtWindow2D
::GetInputImage(void) const
{
  return cImData;
}


void 
QtWindow2D
::SetInputOverlay( OverlayType * newOverlayData )
{
  RegionType newoverlay_region = newOverlayData->GetLargestPossibleRegion();
  SizeType   newoverlay_size  = newoverlay_region.GetSize();
 
  SizeType   cImData_size;
  cImData_size[0] = 0;
  cImData_size[1] = 0;
  cImData_size[2] = 0;

  if (m_isImage1)
  {
  RegionType cImData_region =  cImData->GetLargestPossibleRegion();
  cImData_size   = cImData_region.GetSize();
  }

  if ( (!m_isImage1) || (cImData_size == newoverlay_size))
  {
    cOverlayData = newOverlayData;
    
    cViewOverlayData  = true;
    cValidOverlayData = true;
    //cOverlayOpacity   = (float)0.5;

    m_isOverlay1 = true;

  /*  if(cWinOverlayData != NULL) 
    {
      delete [] cWinOverlayData;
    }
    */

    typedef MinimumMaximumImageCalculator<OverlayType> CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
    calculator->SetImage(newOverlayData);
    calculator->Compute();     
    int m_min = calculator->GetMinimum();
    int m_max = calculator->GetMaximum();
    SetOverlayMinMax(m_min,m_max);

    //const unsigned long bufferSize = (cWinDataSizeX) * cWinDataSizeY * 4;
    //cWinOverlayData = new unsigned char[ bufferSize ];
    //this->update();
  }
   update();
   this->updateGL();
}


void 
QtWindow2D
::SetInputOverlay2( OverlayType * newOverlayData )
{
  RegionType newoverlay_region = newOverlayData->GetLargestPossibleRegion();
  SizeType   newoverlay_size  = newoverlay_region.GetSize();
 
  SizeType   cImData_size;
  cImData_size[0] = 0;
  cImData_size[1] = 0;
  cImData_size[2] = 0;

  if (m_isImage2)
  {
    RegionType cImData_region =  cImData2->GetLargestPossibleRegion();
    cImData_size   = cImData_region.GetSize();
  }

  if ( (!m_isImage2) || (cImData_size == newoverlay_size))
  {
    cOverlayData2 = newOverlayData;
    m_isOverlay2 = true;

    typedef MinimumMaximumImageCalculator<OverlayType> CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
    calculator->SetImage(newOverlayData);
    calculator->Compute();     
    int m_min = calculator->GetMinimum();
    int m_max = calculator->GetMaximum();
    SetOverlayMinMax2(m_min,m_max);

  }

  update();
  this->updateGL();
}


void 
QtWindow2D
::UnSetInputImage()
{
  m_isImage1 = false;
  update();
}


void 
QtWindow2D
::UnSetInputImage2()
{
  m_isImage2 = false;
  update();
}

void 
QtWindow2D
::UnSetInputOverlay()
{
  m_isOverlay1 = false;
  cValidOverlayData = false;
  update();
}

void 
QtWindow2D
::UnSetInputOverlay2()
{
  m_isOverlay2 = false;
  cValidOverlayData2 = false;
  update();
}





const QtWindow2D::OverlayType::Pointer &
QtWindow2D::GetInputOverlay( void ) 
const
{
  return cOverlayData;
}


void 
QtWindow2D::
ViewOverlayData( bool newViewOverlayData)
{ 
  cViewOverlayData = newViewOverlayData;
  
  if( cViewOverlayCallBack != NULL )
  {
    cViewOverlayCallBack();
  }
  
  this->paintGL();
}



bool 
QtWindow2D::ViewOverlayData(void)
{
  return cViewOverlayData;
}



void 
QtWindow2D::ViewOverlayCallBack(
void (* newViewOverlayCallBack)(void) 
)
{
  cViewOverlayCallBack = newViewOverlayCallBack;
}


void 
QtWindow2D::OverlayOpacity(float newOverlayOpacity)
{
  cOverlayOpacity = newOverlayOpacity; 
  if(cViewOverlayCallBack != NULL) 
  {
    cViewOverlayCallBack();
  }
}

float 
QtWindow2D::OverlayOpacity(void)
{
  return cOverlayOpacity;
}


QtWindow2D::ColorTableType 
* QtWindow2D::GetColorTable(void)
{
  return cColorTable.GetPointer();
}


void 
QtWindow2D::update()
{
  if( !cValidImData ) 
  {
    return;
  }
  
  IndexType ind;
  
  int l;
  
  float tf  = 0.0;
  float tf2 = 0.0;

  if (m_id == 0) ind[2] = numslice;
  if (m_id == 1) ind[1] = numslice;
  if (m_id == 2) ind[0] = numslice;


  //Create 2D image
  Image2DType::Pointer m_2Dimage = Image2DType::New();
  float values[2];
  values[0]= cSpacing[0]; 
  values[1]= cSpacing[1];

  float origin_x= ((cWinDataSizeX/2)*values[0]*(-1));
  float origin_y=((cWinDataSizeY/2)*values[1]*(-1));
  float origin[2] = {origin_x, origin_y};
 
  Image2DType::RegionType region2D;
  Image2DType::SizeType size2D;
  size2D[0]= cWinDataSizeX;
  size2D[1]= cWinDataSizeY;

  region2D.SetSize( size2D );
  m_2Dimage->SetRegions( region2D );
  m_2Dimage->Allocate();
  m_2Dimage->SetOrigin(origin);
  m_2Dimage->SetSpacing(values);

  Iterator2DType it2DS(m_2Dimage,m_2Dimage->GetLargestPossibleRegion());
  it2DS.GoToBegin();

  if (cValidImData2 == false)
    m_alpha = 100;

  for(int k=cWinMinY; k <= cWinMaxY; k++)
  {
    if (m_id == 0) ind[1] = k;
    if (m_id == 1) ind[2] = k;
    if (m_id == 2) ind[2] = k;

    for(int j=cWinMinX; j <= cWinMaxX; j++) 
    { 
      if (m_id == 0) ind[0] = j;
      if (m_id == 1) ind[0] = j;
      if (m_id == 2) ind[1] = j;

      if (cImData->GetPixel(ind) > cIWMax)
         tf = 255;
      else
        if (cImData->GetPixel(ind) < cIWMin)
         tf = 0;
      else
       tf = (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
   
      l = j + k*cWinDataSizeX;

      if (m_blendingmode == 0)
      {
        if (k<=(cWinMaxY*m_alpha)/100)
        {
          if (m_isImage1)
              it2DS.Set((unsigned char)tf);
          else
               it2DS.Set((unsigned char)0);
      
        }
         else
        {
          if (m_isImage2)
          {  
              tf2 = (float)((cImData2->GetPixel(ind)-cIWMin2)/(cIWMax2-cIWMin2)*255);
              it2DS.Set((unsigned char)tf2);
          }
          else
              it2DS.Set((unsigned char)0);  
        }
      }
      else
      {
          if ((m_isImage1) && (!m_isImage2))
              it2DS.Set((unsigned char)tf);
          else
          {
            if ((m_isImage1) && (m_isImage2))
            {  
              tf2 = (float)((cImData2->GetPixel(ind)-cIWMin2)/(cIWMax2-cIWMin2)*255);
              float m_alphab = ((float)m_alpha)/(float)100;
              it2DS.Set((unsigned char)((tf*m_alphab)+(tf2*(1-m_alphab))));
            } 
        }
      }

      ++it2DS;
    }
  }

  /*if (m_GreyTexture)
    delete m_GreyTexture;*/

  //Image2DType::SizeType m_size = m_2Dimage->GetLargestPossibleRegion().GetSize();

  m_GreyTexture.SetImage(m_2Dimage);


  //int m_texturealpha = 100;
  
  if (m_overlaymax <= 3)
  {
    //Update color label table
    m_Color[1*255].r = 1*65533; m_Color[1*255].g = 0; m_Color[1*255].b = 0;
    m_Color[2*255].r = 0;       m_Color[2*255].g = 1*65533; m_Color[2*255].b = 0;
    m_Color[3*255].r = 0;       m_Color[3*255].g = 0; m_Color[3*255].b = 1*65533;
  }

  if (m_isOverlay1)
  {
    //Create 2D Label image
    Image2DLabelType::Pointer m_2Dlabelimage = Image2DLabelType::New();
    float values[2];
    values[0]= cSpacing[0]; 
    values[1]= cSpacing[1];

    float origin_x= ((cWinDataSizeX/2)*values[0]*(-1));
    float origin_y=((cWinDataSizeY/2)*values[1]*(-1));
    float origin[2] = {origin_x, origin_y};
 
    Image2DType::RegionType region2D;
    Image2DType::SizeType size2D;
    size2D[0]= cWinDataSizeX;
    size2D[1]= cWinDataSizeY;

    region2D.SetSize( size2D );
    m_2Dlabelimage->SetRegions( region2D );
    m_2Dlabelimage->Allocate();
    m_2Dlabelimage->SetOrigin(origin);
    m_2Dlabelimage->SetSpacing(values);

    Iterator2DLabelType it2LabelDS(m_2Dlabelimage,m_2Dlabelimage->GetLargestPossibleRegion());
    it2LabelDS.GoToBegin();

    for(int k=cWinMinY; k <= cWinMaxY; k++)
    {
      if (m_id == 0) ind[1] = k;
      if (m_id == 1) ind[2] = k;
      if (m_id == 2) ind[2] = k;

      for(int j=cWinMinX; j <= cWinMaxX; j++) 
      { 
        if (m_id == 0) ind[0] = j;
        if (m_id == 1) ind[0] = j;
        if (m_id == 2) ind[1] = j;

        if (cOverlayData->GetPixel(ind) < m_overlaymin)
          tf = 255;
        else
          if (cOverlayData->GetPixel(ind) > m_overlaymax)
            tf = 0;
          else
          {
            if (m_overlaymax > 3)
             tf = 255-(float)(((float)cOverlayData->GetPixel(ind)-(float)m_overlaymin)/((float)m_overlaymax-(float)m_overlaymin))*255;
            else
             tf = cOverlayData->GetPixel(ind);
          }
       
        if (m_isOverlay2)
        {
          if (m_overlaymax2 > 3)
            tf2 = 255-(float)(((float)cOverlayData2->GetPixel(ind)-(float)m_overlaymin2)/((float)m_overlaymax2-(float)m_overlaymin2))*255;
          else
            tf2 = cOverlayData2->GetPixel(ind);
        }
      
        l = j + k*cWinDataSizeX;

        if (k<=(cWinMaxY*m_alpha)/100)
        {
          LabelPixelType m_value;
          m_value.SetRed((unsigned char)(m_Color[(int)tf*255].r/256));
          m_value.SetGreen((unsigned char)(m_Color[(int)tf*255].g/256));
          m_value.SetBlue((unsigned char)(m_Color[(int)tf*255].b/256));
          m_value.SetAlpha((unsigned char)(255)); 


          if (m_displayoverlayzero)
            it2LabelDS.Set((LabelPixelType)m_value);
          else
          {
            if (m_overlaymax > 3)
            {
              if (255-tf != 0)
                it2LabelDS.Set((LabelPixelType)m_value);
              else
                it2LabelDS.Set((unsigned char)0);
            }
            else
            {
              if (tf != 0)
                it2LabelDS.Set((LabelPixelType)m_value);
              else
                it2LabelDS.Set((unsigned char)0);
            }
          }

        }
         else
        {
          if (m_isOverlay2)
          {  
            LabelPixelType m_value;
            m_value.SetRed((unsigned char)(m_Color[(int)tf2*255].r/256));
            m_value.SetGreen((unsigned char)(m_Color[(int)tf2*255].g/256));
            m_value.SetBlue((unsigned char)(m_Color[(int)tf2*255].b/256));
            m_value.SetAlpha((unsigned char)(255)); 

            if (m_overlaymax2 > 3)
            {
              if (255-tf2 != 0)
                  it2LabelDS.Set((LabelPixelType)m_value);
              else
                  it2LabelDS.Set((unsigned char)0);
            }
            else
            {
              if (tf2 != 0)
                it2LabelDS.Set((LabelPixelType)m_value);
              else
                it2LabelDS.Set((unsigned char)0);
            }

          }
          else
              it2LabelDS.Set((unsigned char)0);  
        }
        ++it2LabelDS;
      }
    }

    m_LabelTexture.SetImage(m_2Dlabelimage);
  }
  this->updateGL();
}

void QtWindow2D::size(int w, int h)
{
  this->update();
  this->paintGL();
}



/** Set up the OpenGL view port, matrix mode, etc. */
void QtWindow2D::resizeGL( int w, int h )
{
  glViewport( 0, 0, (GLint)w, (GLint)h );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 1.0, 200.0);
}

/** Initialize the OpenGL Window */
void QtWindow2D::initializeGL() 
{
 /* glClearColor((float)0.0, (float)0.0, (float)0.0, (float)0.0);          
  glShadeModel(GL_FLAT);*/
    
  //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  //if you don't include this
    //image size differences distort
    //glPixelStorei(GL_PACK_ALIGNMENT, 1);
 // makeRasterFont();

}

/** Draw */
void QtWindow2D::paintGL(void)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0,width(),0.0,height());
  glViewport(0,0,width(),height());

  // Establish the model view matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

 // Clear the display, using a blue shade when under focus
  glClearColor(0,0,0,1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    

  //Set up lighting attributes
  glPushAttrib(GL_LIGHTING_BIT | GL_DEPTH_BUFFER_BIT | 
               GL_PIXEL_MODE_BIT | GL_TEXTURE_BIT );  
  
  glDisable(GL_LIGHTING);

  if( !cImData ) 
    return;

  // Prepare for overlay drawing.  The model view is set up to correspond
  // to pixel coordinates of the slice
  glPushMatrix();
  glTranslated(0.5 * width(),0.5 * height(),0.0);
  float m_zoom = (float)view_zoom;
  glScalef(m_zoom,m_zoom,1.0);

  // Set the view position (position of the center of the image?)
  float m_ViewPositionX = cWinDataSizeX * cSpacing[0] * 0.5f;
  float m_ViewPositionY = cWinDataSizeY * cSpacing[1] * 0.5f;

  glTranslated(-m_ViewPositionX+view_pos_x,-m_ViewPositionY+view_pos_y,0.0);
  glScalef(cSpacing[0],cSpacing[1],1.0);

  m_GreyTexture.Draw();


  if (m_isOverlay1 && m_overlayvisible)
  {
   m_LabelTexture.DrawTransparent((unsigned char)(cOverlayOpacity*255));
  } 
     
  
 
  if (m_crosshairvisible)
  {
    glPushAttrib(GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
    double clrCrosshair[4];
    clrCrosshair[0] = 1;
    clrCrosshair[1] = 0.1;
    clrCrosshair[2] = 0;
    clrCrosshair[2] = 0.85;

    int m_crosshairx = image_x;
    int m_crosshairy = image_y;

    glColor4dv(clrCrosshair);
    // Refit matrix so that the lines are centered on the current pixel
    glPushMatrix();
  
    glTranslated( m_crosshairx+0.5f, m_crosshairy+0.5f, 0.0 );
 


    glBegin(GL_LINES); 
    glVertex2f(0, 0); glVertex2f(cWinDataSizeX-m_crosshairx, 0);
    glVertex2f(0, 0); glVertex2f(-m_crosshairx, 0);
    glVertex2f(0, 0); glVertex2f(0, cWinDataSizeY-m_crosshairy);
    glVertex2f(0, 0); glVertex2f(0, -m_crosshairy);
    glEnd();

    glPopMatrix();
    glPopAttrib();
  }

  // Clean up the GL state
  glPopMatrix();
  glPopAttrib(); 

  DrawString();

  // Display!
  glFlush();
}

void QtWindow2D::DrawString()
{
  int offset = height()-13;
  glColor4f((float)252/(float)255, (float)163/(float)255, (float)21/(float)255,1);
  if (m_drawimagename && (m_imagename != ""))
  {
    OutputString(0,offset,m_imagename.latin1());
    offset -= 14;
  }

  ImageType::IndexType ind;
  
  ind[0] = (unsigned long)cClickSelect[0];
  ind[1] = (unsigned long)cClickSelect[1];
  ind[2] = (unsigned long)cClickSelect[2];
 
  if (m_drawposition)
  {
    OutputString(0,offset,QString("Pos: [%1,%2,%3]").arg(ind[0]).arg(ind[1]).arg(ind[2]));
    offset -= 14;
  }

 if (m_drawintensity)
  {
    OutputString(0,offset,QString("Value: %1 [%2|%3]").arg(cImData->GetPixel(ind)).arg(m_imagemin).arg(m_imagemax));
    offset -= 14;
    if (m_isOverlay1)
    {
      OutputString(0,offset,QString("Overlay: %1 [%2|%3]").arg(cOverlayData->GetPixel(ind)).arg(m_overlaymin).arg(m_overlaymax));
      offset -= 14;
    }

    if (m_isImage2)
    {
      OutputString(0,offset,QString("Value2: %1 [%2|%3]").arg(cImData2->GetPixel(ind)).arg(cIWMin2).arg(cIWMax2));
      offset -= 14;
    }
  }

 if (m_drawimageinfo && (m_imageinfotype != ""))
  {
    OutputString(0,offset,QString("Size: %1x%2x%3").arg(m_imageinfosize[0]).arg(m_imageinfosize[1]).arg(m_imageinfosize[2]));
    offset -= 14;
    OutputString(0,offset,QString("Pixdim: %1x%2x%3").arg(m_imageinfopixdim[0]).arg(m_imageinfopixdim[1]).arg(m_imageinfopixdim[2]));
    offset -= 14;
    OutputString(0,offset,QString("Type: %1").arg(m_imageinfotype));
    offset -= 14;
  }
}

void QtWindow2D::SetImageName(QString imagename)
{
    m_imagename = imagename;
    m_drawimagename = true;

}


void QtWindow2D::SetManualIntensityOn()
{
  m_manualintensity = true;
}

void QtWindow2D::SetImageMinMax(int min,int max)
{
  m_imagemin = min;
  m_imagemax = max;
}

void QtWindow2D::GetImageMinMax(int& min,int& max)
{
  min = m_imagemin;
  max = m_imagemax;
}

void QtWindow2D::SetDisplayOverlayZero(bool flag)
{
  m_displayoverlayzero = flag;
}

void QtWindow2D::SetOverlayMinMax(int min,int max)
{
  m_overlaymin = min;
  m_overlaymax = max;
}

void QtWindow2D::SetOverlayMinMax2(int min,int max)
{
  m_overlaymin2 = min;
  m_overlaymax2 = max;
}

void QtWindow2D::DrawImageName(bool flag)
{
    m_drawimagename = flag;
    updateGL();
}

void QtWindow2D::DrawPosition(bool flag)
{
    m_drawposition = flag;
    updateGL();
}

void QtWindow2D::DrawIntensity(bool flag)
{
    m_drawintensity = flag;
    updateGL();
}

void QtWindow2D::SetImageInfo(int* size,float* pixdim,QString type)
{
  for (int i=0;i<3;i++)
  {
    m_imageinfosize[i] = size[i];
    m_imageinfopixdim[i] = pixdim[i];
  }
  m_imageinfotype = type;
}

void QtWindow2D::DrawImageInfo(bool flag)
{
    m_drawimageinfo = flag;
    updateGL();
}


void QtWindow2D::OutputString ( float x, float y,const char *string)
{
  int length;
  length = strlen ( string );
  glRasterPos2f ( x, y );
  for (int i = 0; i < length; i++ )
  glutBitmapCharacter ( GLUT_BITMAP_8_BY_13 , string[i] );
} 

void QtWindow2D::SetZoomFactor(float zoom)
{
   view_zoom = zoom;
   updateGL();
}

void QtWindow2D::SetViewPosition(float posx,float posy)
{
   view_pos_x = posx; 
   view_pos_y = posy;

   updateGL();
}

void QtWindow2D::mouseMoveEvent( QMouseEvent *event ) 
{
  QPoint pos = static_cast<QMouseEvent *>( event )->pos();
  if (m_currentbutton == 1)
  {
        shift_x = (float)(pos.x() - start_x);
        shift_y = (float)(pos.y() - start_y);
        view_pos_x += (float)(shift_x/(view_zoom)); //*cSpacing[0]));//*cSpacing[0]);
        view_pos_y += (float)(-shift_y/(view_zoom));//*cSpacing[1])); //*cSpacing[1]);

        if( event->state() & ShiftButton ) 
          m_manager->SetViewPosition(view_pos_x,view_pos_y);

        start_x = pos.x();
        start_y = pos.y();

        updateGL();
  }

  if (m_currentbutton == 3)
  {
        shift_y = (float)(pos.y() - start_y);
        view_zoom += (float)shift_y/20;
        if (view_zoom < 0) view_zoom = 0;
        start_y = pos.y();
        
        if( event->state() & ShiftButton ) 
          m_manager->SetZoomFactor(view_zoom);

        updateGL();
  }

   this->updateGL();
}



void QtWindow2D::mousePressEvent( QMouseEvent *event ) 
{
   if (m_isImage1 || m_isImage2 || m_isOverlay1 || m_isOverlay2)
  {
     
    if(( event->button() & LeftButton ) && !( event->state() & ControlButton ))
    {
      QPoint pos = static_cast<QMouseEvent *>( event )->pos();
      start_x = (float)pos.x();
      start_y = (float)pos.y();
      m_currentbutton = 1;

      if( event->state() & ShiftButton ) 
      {
      }
     }
     else if (( event->button() & MidButton ) || ((event->state() & ControlButton) && ( event->button() & LeftButton )))
     {
       QPoint pos = static_cast<QMouseEvent *>( event )->pos();

       //float shift_x = ((float)view_pos_x*view_zoom*cSpacing[0]) + (float)width()/2.0  - ((float)cWinDataSizeX* (float)view_zoom*cSpacing[0])/2.0;
       //float shift_y = ((float)view_pos_y*view_zoom*cSpacing[1]) + (float)height()/2.0 - ((float)cWinDataSizeY* (float)view_zoom*cSpacing[1])/2.0;
    
       float shift_x = ((float)view_pos_x*view_zoom) + (float)width()/2.0  - ((float)cWinDataSizeX* (float)view_zoom*cSpacing[0])/2.0;
       float shift_y =  ((float)view_pos_y*view_zoom) + (float)height()/2.0 - ((float)cWinDataSizeY* (float)view_zoom*cSpacing[1])/2.0;
   
       float m_image_x = (int)(((float)pos.x() - shift_x)/((float)view_zoom*cSpacing[0]));
       float m_image_y = (int)(((height()-(float)pos.y())- shift_y)/((float)view_zoom*cSpacing[1]));
       if (m_image_x <0) 
         image_x = 0;
       else
         image_x = (unsigned int) m_image_x;


       if (m_image_y <0) 
         image_y = 0;
       else
         image_y = (unsigned int) m_image_y;

      
       if (image_x >= cDimSize[m_x])
         image_x = cDimSize[m_x]-1;

      if (image_y >= cDimSize[m_y])
         image_y = cDimSize[m_y]-1;

      m_currentbutton = 2;
      
      if (m_id == 0)
        m_manager->SetCrosshair(image_x,image_y,numslice);

      if (m_id == 1)
        m_manager->SetCrosshair(image_x,numslice,image_y);

      if (m_id == 2)
        m_manager->SetCrosshair(numslice,image_x,image_y);


      m_currentbutton = 2;
      // middle mouse button
      //this->mouseEventActive = true;
      //QObject::connect( this->stepTimer, SIGNAL(timeout()),
      //                  this->middleButtonFunction );
     }
     else if( event->button() & RightButton ) 
     {
      QPoint pos = static_cast<QMouseEvent *>( event )->pos();
      start_x = (float)pos.x();
      start_y = (float)pos.y();
      m_currentbutton = 3; 

     }
   }
 
   emit Clicked(m_currentbutton);
   this->updateGL();
}


void QtWindow2D::SetCrosshair(int x,int y)
{
  image_x = x;
  image_y = y;


  if (m_id == 0)
  {
    cClickSelect[0] = x;
    cClickSelect[1] = y;
  }

  if (m_id == 1)
  {
    cClickSelect[0] = x;
    cClickSelect[2] = y;
  }

  if (m_id == 2)
  {
    cClickSelect[1] = x;
    cClickSelect[2] = y;
  }

  this->updateGL();
}


void QtWindow2D::GetCrosshair(int& x,int& y)
{
  x = image_x;
  y = image_y;
}




void QtWindow2D::ChangeSlice(unsigned int value)
{
  if(value>=cDimSize[cWinOrder[m_id]])
  return;

  if (m_isImage1 || m_isImage2 || m_isOverlay1 || m_isOverlay2)
  {
  sliceNum(value);
  update();
  SetSliderValue(value);
  if (m_islabel)
  {
    char text[20];
    sprintf(text,"%i/%i",value,cWinMaxZ);
    m_label->setText(text);
  }
  }
}



void QtWindow2D::sliceNum(unsigned int newSliceNum)
{
  if(newSliceNum>=cDimSize[cWinOrder[m_id]])
    newSliceNum = cDimSize[cWinOrder[m_id]]-1;
  cWinCenter[cWinOrder[m_id]] = newSliceNum;
  numslice = newSliceNum;
  cClickSelect[2-m_id] = numslice;
  /*if(cSliceNumCallBack != NULL)
    cSliceNumCallBack();
  if(cSliceNumArgCallBack != NULL)
    cSliceNumArgCallBack(cSliceNumArg);*/
}

unsigned int QtWindow2D::sliceNum()
{
  return cWinCenter[cWinOrder[m_id]];
}

void QtWindow2D::clickSelect(float newX, float newY, float newZ)
  {    
  cClickSelect[0] = newX;
  if(cClickSelect[0]<0)
    cClickSelect[0] = 0;
  if(cClickSelect[0] >= cDimSize[0])
    cClickSelect[0] = cDimSize[0]-1;
  
  cClickSelect[1] = newY;
  if(cClickSelect[1]<0)
    cClickSelect[1] = 0;
  if(cClickSelect[1] >= cDimSize[1])
    cClickSelect[1] = cDimSize[1]-1;
  
  cClickSelect[2] = newZ;
  if(cClickSelect[2]<0)
    cClickSelect[2] = 0;
  if(cClickSelect[2] >= cDimSize[2])
    cClickSelect[2] = cDimSize[2]-1;
  
  ImageType::IndexType ind;
  
  ind[0] = (unsigned long)cClickSelect[0];
  ind[1] = (unsigned long)cClickSelect[1];
  ind[2] = (unsigned long)cClickSelect[2];
  cClickSelectV = cImData->GetPixel(ind);
 
  emit Position(ind[0],ind[1],ind[2],cClickSelectV);
}


void QtWindow2D::SetIntensityMax(int value)
{
  cIWMax = value;
  update();
}

void QtWindow2D::SetIntensityMin(int value)
{
  //OverlayOpacity((float)value/(float)256);
  cIWMin = value;
  update();
}
 
void QtWindow2D::ZoomIn()
{
  cWinZoom += 1;
  update();
  //this->updateGL();
}
 
void QtWindow2D::ZoomOut()
{
  cWinZoom -= 0.1;
  update();
  //this->updateGL();
}


void QtWindow2D::SetManager(ManagerType* manager)
{
  m_manager = manager;
  manager->AddWindow2D(this);
}


void QtWindow2D::SetSlider(QSlider* slicer)
{
  m_slicer = slicer;
  m_isslicer = true;
}

void QtWindow2D::SetLabel(QLabel* label)
{
  m_label = label;
  m_islabel = true;
}


void QtWindow2D::SetSliderValue(int value)
{
  if (m_isslicer)
    m_slicer->setValue(value);
}


void QtWindow2D::SetId(int id)
{
  m_id = id;
  cWinOrientation = m_id;
}

int QtWindow2D::GetId()
{
   return m_id;
}

void QtWindow2D::makeRasterFont(void)
{
 /*  GLuint i, j;
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   fontOffset = glGenLists (128);
   for (i = 0,j = 'A'; i < 26; i++,j++) {
      glNewList(fontOffset + j, GL_COMPILE);
      glBitmap(8,13, 0.0, 2.0, 10.0, 0.0, letters[i]);
      glEndList();
   }
   glNewList(fontOffset + ' ', GL_COMPILE);
   glBitmap(4, 6, 0.0, 2.0, 10.0, 0.0, space); /// 8 13
   glEndList();*/
}

void QtWindow2D::printString(char *s)
{
   glPushAttrib (GL_LIST_BIT);
   glListBase(fontOffset);
   glCallLists(strlen(s), GL_UNSIGNED_BYTE, (GLubyte *) s);
   glPopAttrib ();
}

void QtWindow2D::ChangeAlpha(int value)
{
  m_alpha = value;
  this->update();
}



void QtWindow2D::ChangeLabelOverlay(int value)
{
  cOverlayOpacity = (float)value/(float)100;
  update();

  this->update();
}

void QtWindow2D::SetRotation(float angle)
{
  m_rotation = angle;
  this->update();
}

void QtWindow2D::ShowCrosshair(bool flag)
{
  cShowCrosshair = flag;
  this->update();
}


int QtWindow2D::Bicubic_Interpol(unsigned char *m_image, int newx, int newy, int xsize, int ysize, float a, float b, int f)
{

  /* Perform the bicubic interpolation method */

   int m, n;
   float sum, cubsum;
   int index;
   int imageval;

   sum = 0;

   /* Find pixel value */ 
   index = ((newx + newy*xsize));   

   for (m = -1; m <=2; m++) {
       for (n = -1; n <= 2; n++) {
   if ((newx+m >= 1) && (newx+m < xsize-1)) {
              if ((newy+n >=1) && (newy+n < ysize-1)) {  
    /* Find surrounding pixels values */ 
                 imageval = m_image[index + (m + xsize*n)];
              } else {
     imageval = 0;
              }
         } else {
              imageval=0;
         }
         /* Weight... */
         cubsum = CubicFunc(f, (float)(-a+m))*CubicFunc(f, (float)(-b+n));
         /* Add it to do a complete weighted average */
         sum = sum + imageval*cubsum;
      
       } 
   }

   return((int) sum);
}


/* Interpolation functions */

float QtWindow2D::CubicBSpline(float x) {

  /* Return the value of the cubic B spline
     between -2 .. 2          */

   float y, y2, y3;
   float result;

   result=0;
 
   if (x < 0) {
      y = -x; 
   } else {
      y = x;
   }
   
   y2=y*y;
   y3=y2*y;

   if (y < 1.0) {
      result=(2.0/3.0 + 0.5 * y3 - y2);
   } else {
     if (y <= 2.0) {
         result=(1.0/6.0)*(2-y)*(2-y)*(2-y);
     }
   }
   return(result);
}

float QtWindow2D::CubicInterpolation(float x) {

  /* Returns the value of an interpolated
     sine between -2 .. 2  */

   float a, y, y2, y3, result;
  
   a=0.1;  /* arbitrary value */
   result=0;
 
   if (x < 0) {
      y = -x;
   } else {
      y=x;
   }

   y2=y*y;
   y3=y2*y;

   if (y <= 1.0) {
      result=((a+2)*y3 - (a+3)*y2 + 1);
   } else {
     if (y <= 2.0) {
        result=(a*y3 - 5*a*y2 + 8*a*y - 4*a);
     }
   }
   
   return(result);
}

float QtWindow2D::SquareFunction(float x) {

  /* Returns 1 between -0.5 .. 0.5 */

  float result;

  if ((x >= -0.5) && (x <= 0.5)) {
     result = 1.0;
  } else {
     result = 0;
  }
  return(result);
}

float QtWindow2D::TriangleFunction(float x) {

  /* Returns the triangle function 
     between -1.0 .. 1.0            */ 

   float result;

   result = 0;

   if ((x >= -1) && (x <= 0)){ 
      result = x + 1;
   } else {
      if ((x > 0) && (x <= 1)) result = 1 - x;
   }

   return(result);
}

float QtWindow2D::CubicFunc(int f, float x) {

  /* Returns the different cubic functions */
  /* 0 = B-Spline, 1 = Interpolated sine   */

   float result = 0.0;

   switch (f) {
        case 0: result = CubicBSpline(x);
                break;
        case 1: result = CubicInterpolation(x);
                break;
   }
   
   return(result);
}

float QtWindow2D::BilinearFunc(int f, float x) {
    
  /* Returns the different linear functions */
  /* 0 = Square ; 1 = Triangle              */

   float result = 0.0;

   switch (f) {
        case 0: result = SquareFunction(x);
                break;
        case 1: result = TriangleFunction(x);
                break;
   }
   
   return(result);
}



void QtWindow2D::ResetView()
{
  view_pos_x = 0;
  view_pos_y = 0;
  view_zoom = 1;
  this->updateGL();
}


bool QtWindow2D::OverlayVisible()
{
  return m_overlayvisible;
}

void QtWindow2D::HideOverlay()
{
  if (m_isOverlay1 || m_isOverlay2)
  {
    m_overlayvisible = false;
    this->update();
    this->updateGL();
  }
}

void QtWindow2D::ShowOverlay()
{
  if (m_isOverlay1 || m_isOverlay2)
  {
    m_overlayvisible = true;
    this->update();
    this->updateGL();
  }
}

QtWindow2D::SizeType QtWindow2D::GetImageSize()
{
   return cImData->GetLargestPossibleRegion().GetSize();
}


QPixmap QtWindow2D::RenderPixmap()
{
  makeCurrent();
  paintGL();
  QPixmap pixmap = renderPixmap();
   //QPixmap pixmap2(width(),height());
  // 
  // paintGL();
  // glReadPixels(0,0,width(),height(),GL_RGB,GL_UNSIGNED_BYTE,&pixmap2);
  // bitBlt(&pixmap,0,0,&pixmap2,0,0,width(),height(),CopyROP);
   return pixmap;
}


void QtWindow2D::SetBlendingMode(int mode)
{
  m_blendingmode = mode;
  this->update();
  this->updateGL();
}


bool QtWindow2D::IsCrosshairVisible()
{
  return m_crosshairvisible;
}

void QtWindow2D::SetCrosshairVisible(bool flag)
{
  m_crosshairvisible = flag;
    this->update();
  this->updateGL();
}
