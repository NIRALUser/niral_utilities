#include "IntensityRescalerGUIControls.h"
#include <qslider.h>
#include <qcombobox.h>
#include <qcheckbox.h>
#include <qpushbutton.h>
#include "Vector2D.h"
#include <qmessagebox.h>
#include "FunctionPlot2DBox.h"
  
#define _DEBUG false

IntensityRescalerGUIControls::IntensityRescalerGUIControls( QWidget* parent,  const char* name, WFlags fl )
    : IntensityRescalerGUI( parent, name,  fl )
{
  m_imagemanager = new ImageManager<PixelType,OverlayPixelType>();
  m_imagemanager->SetTargetGrayLabel(g_sourcegrayvalue);
  m_imagemanager->SetTargetSegLabel(g_sourcesegvalue);
  m_imagemanager->SetSourceGrayLabel(g_targetgrayvalue);
  m_imagemanager->SetSourceSegLabel(g_targetsegvalue);
  g_Window2DX->SetManager(m_imagemanager);
  g_Window2DX->SetSlider(g_SliderX);
  g_Window2DX->SetLabel(g_LabelX);
  g_Window2DX->SetId(0);

  g_Window2DY->SetManager(m_imagemanager);
  g_Window2DY->SetSlider(g_SliderY);
  g_Window2DY->SetLabel(g_LabelY);
  g_Window2DY->SetId(1);
  
  g_Window2DZ->SetManager(m_imagemanager);
  g_Window2DZ->SetSlider(g_SliderZ);
  g_Window2DZ->SetLabel(g_LabelZ);
  g_Window2DZ->SetId(2);
  
  m_intensitycurve = IntensityCurveVTK::New();
  m_intensitycurve->Initialize(5);
  g_intensitycurve->SetCurve(m_intensitycurve);
  connect(g_intensitycurve, SIGNAL(MouseMove()), this, SLOT(IntensityMouseMove()));

  m_imageIntensityNormalizer = new ImageIntensityNormalizer();
  m_image1 = 0;
  m_image2 = 0;
  m_overlay1 = 0;
  m_overlay2 = 0;
}


IntensityRescalerGUIControls::~IntensityRescalerGUIControls()
{
  delete m_imagemanager;
  delete m_imageIntensityNormalizer;
}


void IntensityRescalerGUIControls::LoadImg()
{
  QtImageLoader<PixelType> * m_loader = new QtImageLoader<PixelType>();  
  ImageType::Pointer m_inputimage = m_loader->GetOutput();
  g_LoadImg->setDown(false);
  if (m_inputimage.IsNull()) return;
  g_targetimage->insertItem(m_loader->GetFileName());
  g_sourceimage->insertItem(m_loader->GetFileName());
  g_image->insertItem(m_loader->GetFileName());
  //Add Image pointer to the list
  m_SourceImage.push_back(m_inputimage);
  //Load the first image
  if (m_SourceImage.size() == 1)
  {
    g_targetimage->setCurrentItem(1);
    SelectTargetImage(1);
  }
  if (m_SourceImage.size() == 2)
  {
    g_sourceimage->setCurrentItem(2);
    SelectSourceImage(2);
  }

  delete m_loader;
}

void IntensityRescalerGUIControls::LoadSeg()
{
  QtImageLoader<OverlayPixelType> * m_loader = new QtImageLoader<OverlayPixelType>();
  OverlayType::Pointer m_inputimage = m_loader->GetOutput();
  g_LoadSeg->setDown(false);
  if (m_inputimage.IsNull()) return;
  g_targetoverlay->insertItem(m_loader->GetFileName());
  g_sourceoverlay->insertItem(m_loader->GetFileName());
  g_overlay->insertItem(m_loader->GetFileName());
  //Add Image pointer to the list
  m_SourceOverlay.push_back(m_inputimage);
  if (m_SourceOverlay.size() == 1)
  {
    g_targetoverlay->setCurrentItem(1);
    SelectTargetOverlay(1);
  }
  if (m_SourceOverlay.size() == 2)
  {
    g_sourceoverlay->setCurrentItem(2);
    SelectSourceOverlay(2);
  }

  delete m_loader;
}

void IntensityRescalerGUIControls::UnloadImage()
{
  m_SourceImage.erase(m_SourceImage.begin() + g_image->currentItem()); 

  if (g_targetimage->currentItem() == g_image->currentItem()+1);
  {
    SelectTargetImage(0);
    g_targetimage->setCurrentItem(0);
  }

  if (g_sourceimage->currentItem() == g_image->currentItem()+1);
  {
    SelectSourceImage(0);
    g_sourceimage->setCurrentItem(0);
  }

  g_targetimage->removeItem(g_image->currentItem()+1);
  g_sourceimage->removeItem(g_image->currentItem()+1);
  g_image->removeItem(g_image->currentItem());
}     

void IntensityRescalerGUIControls::UnloadOverlay()
{
  m_SourceOverlay.erase(m_SourceOverlay.begin() + g_overlay->currentItem()); 
  
  if (g_targetoverlay->currentItem() == g_overlay->currentItem()+1);
  {
    SelectTargetOverlay(0);
    g_targetoverlay->setCurrentItem(0);
  }

  if (g_sourceoverlay->currentItem() == g_overlay->currentItem()+1);
  {
    SelectSourceOverlay(0);
    g_sourceoverlay->setCurrentItem(0);
  }
  
  g_targetoverlay->removeItem(g_overlay->currentItem()+1);
  g_sourceoverlay->removeItem(g_overlay->currentItem()+1);
  g_overlay->removeItem(g_overlay->currentItem());
}     


void IntensityRescalerGUIControls::ChangeSliceX(int value)
{
  m_imagemanager->ChangeSliceX(value);
}

void IntensityRescalerGUIControls::ChangeSliceY(int value)
{
  m_imagemanager->ChangeSliceY(value);
}

void IntensityRescalerGUIControls::ChangeSliceZ(int value)
{
  m_imagemanager->ChangeSliceZ(value);
}

void IntensityRescalerGUIControls::IntensityMouseMove()
{
  if (m_image2.IsNull())
    return;

  double t = 0.0;
  double tStep = 1.0f / m_max2;
  
  for (int i=0;i<=m_max2;i++)
  {
    m_correspondingarray[i] = (int)(m_intensitycurve->Evaluate(t)*m_max2);
      t+=tStep;
  }
  ImagePointer m_image;
  m_image = m_imageIntensityNormalizer->Normalize(m_image2,m_correspondingarray);
  m_imagemanager->SetTargetImage(m_image);
}


void IntensityRescalerGUIControls::SelectTargetImage(int value)
{
  if (value == 0)
  {
    m_imagemanager->UnSetSourceImage();
    m_image1 = 0;
  }
  else
  {
    ImagePointer m_inputimage = m_SourceImage[value-1];
    m_imagemanager->SetSourceImage(m_inputimage);
    m_image1 = m_inputimage;
    if (m_image2) 
      m_imagemanager->ChangeAlpha(50);
    else
     m_imagemanager->ChangeAlpha(100);
  }
}

void IntensityRescalerGUIControls::SelectSourceImage(int value)
{
  if (value == 0)
  {
    m_imagemanager->UnSetTargetImage();
    m_image2 = 0;
  }
  else
  {
    ImagePointer m_inputimage = m_SourceImage[value-1];
    m_imagemanager->SetTargetImage(m_inputimage);
    m_image2 = m_inputimage;
    typedef MinimumMaximumImageCalculator<ImageType> CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
    calculator->SetImage( m_image2 );
    calculator->Compute();     
    m_min2 = calculator->GetMinimum();
    m_max2 = calculator->GetMaximum();
    if (m_correspondingarray) delete m_correspondingarray;
      m_correspondingarray = new int[m_max2+1];
    if (m_image1) 
      m_imagemanager->ChangeAlpha(50);
    else
     m_imagemanager->ChangeAlpha(0);
  }
}

void IntensityRescalerGUIControls::SelectTargetOverlay(int value)
{
  if (value == 0)
  {
    m_imagemanager->UnSetSourceOverlay();
    g_labellist->clear();
    m_overlay1 = 0;
  }
    else
  {
    ImagePointer m_inputimage = m_SourceOverlay[value-1];
    m_imagemanager->SetSourceOverlay(m_inputimage);

    m_overlay1 = m_inputimage;
    g_labellist->clear();
    std::vector<int> m_labellist = m_imageIntensityNormalizer->ListLabel(m_inputimage);
    for ( unsigned int i = 0; i < m_labellist.size(); i++ )
    (void)new QCheckListItem( g_labellist, QString( "%1" ).arg( m_labellist[i] ), QCheckListItem::CheckBox );


    QCheckListItem* item = (QCheckListItem*)g_labellist->firstChild();
    while (item)
    {
      item->setOn(true);
      item = (QCheckListItem*)item->nextSibling();
    }
  }
}

void IntensityRescalerGUIControls::SelectSourceOverlay(int value)
{
  if (value == 0)
  {
    m_imagemanager->UnSetTargetOverlay();
    m_overlay2 = 0;
  }
  else
  {
    ImagePointer m_inputimage = m_SourceOverlay[value-1];
    m_imagemanager->SetTargetOverlay(m_inputimage);
    m_overlay2 = m_inputimage;
  }
}

void IntensityRescalerGUIControls::ChangeLabelOverlay(int value)
{
  m_imagemanager->ChangeLabelOverlay(value);
}

void IntensityRescalerGUIControls::ChangeImageOverlay(int value)
{
  m_imagemanager->ChangeAlpha(value);
}

void IntensityRescalerGUIControls::Compute()
{
  if (m_image1.IsNull()) 
  {
    QMessageBox::information( this, "IntensityRescaler (Qt)","No target image");  
    g_compute->setDown(false);
    return;
  }
  if (m_image2.IsNull()) 
  {
    QMessageBox::information( this, "IntensityRescaler (Qt)","No source image");  
    g_compute->setDown(false);
    return;
  }  
  if (m_overlay1.IsNull()) 
  {
    QMessageBox::information( this, "IntensityRescaler (Qt)","No target segmentation");  
    g_compute->setDown(false);
    return;
  }
  if (m_overlay2.IsNull()) 
  {
    QMessageBox::information( this, "IntensityRescaler (Qt)","No source segmentation");  
    g_compute->setDown(false);
    return;
  }

  float sigma = 3;
  double newmax = m_max2;
  ImagePointer m_targetrescaled;
  ImagePointer m_sourcerescaled;
  ImageIntensityNormalizer::VectorList* m_vectorlist;
  ImageIntensityNormalizer::VectorList* m_vectorlist2;

  //Compute LabelList
  std::vector<int> m_labellist = m_imageIntensityNormalizer->ListLabel(m_overlay1);


  //Add all label list
  m_imageIntensityNormalizer->SetLabelList(m_labellist);

  //Add considering label(s)
  m_imageIntensityNormalizer->ClearLabel();

  QCheckListItem* item = (QCheckListItem*)g_labellist->firstChild();
  while (item)
  {
    if (item->isOn())
      m_imageIntensityNormalizer->AddLabel(atoi(item->text(0).latin1()));

    item = (QCheckListItem*)item->nextSibling();
  }
  
  if (g_targetwindowing->isChecked())
  {
    //Compute Mean values for target
    m_vectorlist = m_imageIntensityNormalizer->ComputeMax(m_image1,m_overlay1);

    //Intensity Windowing for target
    m_targetrescaled = m_imageIntensityNormalizer->TargetIntensityWindowing(m_image1,*m_vectorlist,sigma,&newmax);
    m_imagemanager->SetSourceImage(m_targetrescaled);
    
    //Update Combo Box
    g_targetimage->insertItem(g_targetimage->currentText()+ "-rescaled");
    g_targetimage->setCurrentItem(g_targetimage->count()-1);
    g_sourceimage->insertItem(g_targetimage->currentText());
    g_image->insertItem(g_targetimage->currentText());
    m_SourceImage.push_back(m_targetrescaled);
    m_image1 = m_targetrescaled;
  }

  if (g_sourcewindowing->isChecked())
  {
    //Compute Mean values for target
    m_vectorlist = m_imageIntensityNormalizer->ComputeMax(m_targetrescaled,m_overlay1);
    m_vectorlist2 = m_imageIntensityNormalizer->ComputeMax(m_image2,m_overlay2);

    //Compute max value for x*sigma
    //double meanmin=999999999;
    //double meanmax=0;
    //double sourcenewmax;
    //Fix maximum value

    m_sourcerescaled = m_imageIntensityNormalizer->IntensityWindowing(m_image2,*m_vectorlist2,sigma,newmax);
    m_imagemanager->SetTargetImage(m_sourcerescaled);

    if (!g_classematching->isChecked())
    {
      //Update Combo Box
      g_targetimage->insertItem(g_sourceimage->currentText()+ "-rescaled");
      g_sourceimage->insertItem(g_sourceimage->currentText());
      g_sourceimage->setCurrentItem(g_sourceimage->count()-1);
      g_image->insertItem(g_sourceimage->currentText());
      m_SourceImage.push_back(m_sourcerescaled);
    }
    m_image2 = m_sourcerescaled;
  }


  //Update IntensityCurve
  if (g_classematching->isChecked())
  {  
    //Update Values
    ImagePointer m_image;
    ImagePointer m_imagetemp;


    double m_error = 99998;
    double m_previouserror = 99999;
    
    m_vectorlist = m_imageIntensityNormalizer->ComputeMean(m_image1,m_overlay1);
    m_vectorlist2 = m_imageIntensityNormalizer->ComputeMean(m_image2,m_overlay2);

    ImageIntensityNormalizer::VectorList::iterator targetit = m_vectorlist->begin();
    ImageIntensityNormalizer::VectorList::iterator sourceit = m_vectorlist2->begin();

    if (_DEBUG)  std::cout << "Vector Size: " << m_vectorlist2->size() <<std::endl;
    m_intensitycurve->Initialize(m_vectorlist2->size()+2);
    
    int offset=1;
    for (int i=0;i<(int) m_vectorlist2->size();i++)
    {
      m_intensitycurve->UpdateControlPoint(offset,(*sourceit)[1]/newmax,(*targetit)[1]/newmax);  
      if (_DEBUG)  std::cout << "Update: " <<  (*sourceit)[1] << " -> " << (*targetit)[1] <<std::endl;
      sourceit++;
      targetit++;
      offset++;
    }  
  
    while (m_error < m_previouserror)
    {

      m_previouserror = m_error;

      //Update Intensity on image
      double t = 0.0;
      double tStep = 1.0f / newmax;
      for (int i=0;i<=(int) newmax;i++)
      {
        m_correspondingarray[i] = (int)(m_intensitycurve->Evaluate(t)*newmax);
        t+=tStep;
      }
      m_imagetemp = m_image;
      m_image = m_imageIntensityNormalizer->Normalize(m_image2,m_correspondingarray);
      m_imagemanager->SetTargetImage(m_image);


      //Compute Mean values for source
      m_vectorlist2 = m_imageIntensityNormalizer->ComputeMean(m_image,m_overlay2);

      //Compute Error
      ImageIntensityNormalizer::VectorList::iterator targetit = m_vectorlist->begin();
      ImageIntensityNormalizer::VectorList::iterator sourceit = m_vectorlist2->begin();
      double m_errorSum = 0;

      int offset=1;
      while (targetit != m_vectorlist->end())
      {
        m_errorSum += (((*targetit)[1]) - ((*sourceit)[1])) * (((*targetit)[1]) - ((*sourceit)[1]));
        m_intensitycurve->UpdateControlPoint(offset,(*sourceit)[1]/newmax,(*targetit)[1]/newmax);  
        offset++;
        targetit++;
        sourceit++;
      }

      g_intensitycurve->repaint();
      m_error = sqrt(m_errorSum);
      if (_DEBUG)  std::cout << "Error: " << m_error <<std::endl;
    }
      
    //Update Combo Box
    g_targetimage->insertItem(g_sourceimage->currentText()+ "-rescaled");
    g_sourceimage->insertItem(g_sourceimage->currentText()+ "-rescaled");
    g_sourceimage->setCurrentItem(g_sourceimage->count()-1);
    g_image->insertItem(g_sourceimage->currentText());
    m_SourceImage.push_back(m_imagetemp);
  }


  

}

void IntensityRescalerGUIControls::Batch()
{
  m_batch = new BatchControls(0);
  m_batch->Show();
}


void IntensityRescalerGUIControls::ShowHistoSource()
{  
  HistoGUI* m_histogui = new HistoGUI(0);
  //Compute Histogramm
  int* histo;
  int maxx,maxy;

  histo = ComputeHistogram(m_image2,&maxx,&maxy);
  

  //int m_nbPoint = 256;
  float* x = new float[maxx];
  float* y = new float[maxx];

  for (int i=0;i<maxx;i++)
  {
    x[i] = i; //(i*m_nbPoint)/maxx;
    y[i] = histo[i]; //((histo[i]*m_nbPoint)/maxy);
  }

    Vector2D<float> plotmin(0,0);
  Vector2D<float> plotmax(maxx,maxy);
  m_histogui->Plot1->GetPlotter().GetSettings().SetPlotRangeMin(plotmin);
  m_histogui->Plot1->GetPlotter().GetSettings().SetPlotRangeMax(plotmax);
  m_histogui->Plot1->GetPlotter().SetDataPoints(x,y,maxx);



  //Add mean
    std::vector<int> m_labellist = m_imageIntensityNormalizer->ListLabel(m_overlay2);
  m_imageIntensityNormalizer->SetLabelList(m_labellist);
  
  ImageIntensityNormalizer::VectorList* m_vectorlist;
  m_vectorlist = m_imageIntensityNormalizer->ComputeMean(m_image2,m_overlay2);
  if (_DEBUG)  std::cout << "Mean ---" <<m_vectorlist->size() << std::endl;

  


  for (ImageIntensityNormalizer::VectorList::iterator j = m_vectorlist->begin(); j != m_vectorlist->end(); j++)
  {
    if (_DEBUG) std::cout << "Mean: " << (*j)[1] << " | " << (*j)[2] << std::endl;
    m_histogui->Plot1->GetPlotter().AddMean((*j)[1]);
        if ((*j)[0] == 250)
    m_histogui->Plot1->GetPlotter().AddMean((*j)[1]+3*(*j)[2]);
  }

  m_histogui->show();
}

void IntensityRescalerGUIControls::ShowHistoTarget()
{
  HistoGUI* m_histogui = new HistoGUI(0);
  //Compute Histogramm
  int* histo;
  int maxx,maxy;

  histo = ComputeHistogram(m_image1,&maxx,&maxy);
  

  //int m_nbPoint = 256;
  float* x = new float[maxx];
  float* y = new float[maxx];

  for (int i=0;i<maxx;i++)
  {
    x[i] = i; //(i*m_nbPoint)/maxx;
    y[i] = histo[i]; //((histo[i]*m_nbPoint)/maxy);
  }

    Vector2D<float> plotmin(0,0);
  Vector2D<float> plotmax(maxx,maxy);
  m_histogui->Plot1->GetPlotter().GetSettings().SetPlotRangeMin(plotmin);
  m_histogui->Plot1->GetPlotter().GetSettings().SetPlotRangeMax(plotmax);
  m_histogui->Plot1->GetPlotter().SetDataPoints(x,y,maxx);

  //Add distribution
  std::vector<int> m_labellist = m_imageIntensityNormalizer->ListLabel(m_overlay1);
  
  for (int k=0;k<(int)m_labellist.size();k++)
  {
    if (_DEBUG) std::cout << "Compute distibution label for k=" << k << std::endl;
    int* histo;
    int maxx,maxy;

    histo = ComputeDistribution(m_image1,m_overlay1,&maxx,&maxy,m_labellist[k]);
    //int m_nbPoint = 256;
    float* x = new float[maxx];
    float* y = new float[maxx];

    for (int i=0;i<maxx;i++)
    {
      x[i] = i;
      y[i] = histo[i];
    }
    m_histogui->Plot1->GetPlotter().SetDataPoints(x,y,maxx,k+1);
  }


  //Add mean
//  std::vector<int> m_labellist = m_imageIntensityNormalizer->ListLabel(m_overlay1);
  m_imageIntensityNormalizer->SetLabelList(m_labellist);
  
  ImageIntensityNormalizer::VectorList* m_vectorlist;
  m_vectorlist = m_imageIntensityNormalizer->ComputeMean(m_image1,m_overlay1);
  //std::cout << "Mean ---" <<m_vectorlist->size() << std::endl;

  


  for (ImageIntensityNormalizer::VectorList::iterator j = m_vectorlist->begin(); j != m_vectorlist->end(); j++)
  {
    if (_DEBUG) std::cout << "Mean: " << (*j)[1] << " | " << (*j)[2] << std::endl;
    m_histogui->Plot1->GetPlotter().AddMean((*j)[1]);
    if ((*j)[0] == 250)
    m_histogui->Plot1->GetPlotter().AddMean((*j)[1]+3*(*j)[2]);
  }

  m_histogui->show();  
}


int* IntensityRescalerGUIControls::ComputeHistogram(ImagePointer image,int* _maxx,int* _maxy)
{

  int* m_histo = new int[65536];

  //Initialize m_histo
  for (int i=0;i<65535;i++)
      m_histo[i] =0; 

  int maxx = 0;
  int maxy = 0;

  IteratorType m_itimage(image,image->GetLargestPossibleRegion());
//  ImageType::SizeType m_imagesize = (image->GetLargestPossibleRegion().GetSize());
  
  m_itimage.GoToBegin();
  int val;
  while (!m_itimage.IsAtEnd())
  {
    val = (int)m_itimage.Get();
    if (val !=0)
    {
      if (val > maxx) maxx = val;
      m_histo[val] += m_itimage.Get();
  
      if(m_histo[val]>maxy) maxy=m_histo[val];
    }
    ++m_itimage;
  }


  *_maxx = maxx;
  *_maxy = maxy;
  return m_histo;
}


int* IntensityRescalerGUIControls::ComputeDistribution(ImagePointer image,ImagePointer seg,int* _maxx,int* _maxy,int label)
{
  int* m_histo = new int[65536];

  //Initialize m_histo
  for (int i=0;i<65535;i++)
      m_histo[i] =0; 

  int maxx = 0;
  int maxy = 0;

  IteratorType m_itimage(image,image->GetLargestPossibleRegion());
  IteratorType m_itseg(seg,seg->GetLargestPossibleRegion());
  
    
  if (_DEBUG)  std::cout << "Compute distibution label for label" << label << std::endl;
  
  m_itimage.GoToBegin();
  m_itseg.GoToBegin();  
  int val;
  while (!m_itimage.IsAtEnd())
  {
    if ((int)m_itseg.Get() == label)
    {
      val = (int)m_itimage.Get();
      if (val > maxx) maxx = val;
      m_histo[val] += m_itimage.Get();
  
      if(m_histo[val]>maxy) maxy=m_histo[val];
    }

    ++m_itimage;
    ++m_itseg;
  }


  *_maxx = maxx;
  *_maxy = maxy;
  return m_histo;
}


void IntensityRescalerGUIControls::SaveImage()
{
  QtImageWriter<PixelType> * m_writer = new QtImageWriter<PixelType>(g_image->currentText()+".gipl");  
  m_writer->SetInput(m_SourceImage[g_image->currentItem()]);  
  m_writer->Write();
}
