#include "QtSlicer.h"
#include <iostream>
#include "QtGlSliceView.h"
#include <qlineedit.h>
#include <qslider.h>


#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <qfiledialog.h>

/**
 *
 */
QtSlicer::QtSlicer( QWidget* parent,  const char* name, bool modal, WFlags fl )
    : Gui( parent, name, modal, fl )
{
}

/**  
 *  Destroys the object and frees any allocated resources
 */
QtSlicer::~QtSlicer()
{
}

void QtSlicer::DisplayPosition(int x,int y ,int z,float value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",x);
  PositionX->setText(tr(tempchar));
  sprintf(tempchar,"%d",y);
  PositionY->setText(tr(tempchar));
  sprintf(tempchar,"%d",z);
  PositionZ->setText(tr(tempchar));
  sprintf(tempchar,"%3.1f",value);
  PixelValue->setText(tr(tempchar));
  delete tempchar;
}

void QtSlicer::Help()
{
  HelpWindow* helpWindow = new HelpWindow(this,"Help ...");
  helpWindow->show();
}

void QtSlicer::SetInputImage(ImageType * newImData)
{
  this->OpenGlWindow->SetInputImage(newImData);
  this->Slider1->setMaxValue(newImData->GetLargestPossibleRegion().GetSize()[2]-1);

  // Set the slice slider at z/2
  this->Slider1->setValue(newImData->GetLargestPossibleRegion().GetSize()[2]/2);
  this->DisplaySliceNumber(newImData->GetLargestPossibleRegion().GetSize()[2]/2);

  this->IntensityMin->setMinValue( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMin->setMaxValue( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  this->IntensityMin->setValue( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMax->setMinValue( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMax->setMaxValue( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  this->IntensityMax->setValue( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  
  char* tempchar = new char[20];
  sprintf(tempchar,"%.0f",this->OpenGlWindow->GetIntensityMin());
  this->IntensityMinDisplay->setText(tempchar);
  sprintf(tempchar,"%.0f",this->OpenGlWindow->GetIntensityMax());
  this->IntensityMaxDisplay->setText(tempchar);
  delete tempchar;

  this->OpenGlWindow->show();
  this->OpenGlWindow->update();
}


void QtSlicer::SetInputOverlay(OverlayType * newImData)
{
  this->OpenGlWindow->SetInputOverlay(newImData);
  this->OpenGlWindow->ViewOverlayData(true);
  this->OpenGlWindow->OverlayOpacity(0.8);
  this->OpenGlWindow->show();
  this->OpenGlWindow->update();
}



void QtSlicer::DisplaySliceNumber(int number)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",number);
  this->SliceValue->setText(tempchar);
  delete tempchar;
}

void QtSlicer::DisplayIMin(int value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",value);
  this->IntensityMinDisplay->setText(tempchar);
  delete tempchar;
}

void QtSlicer::DisplayIMax(int value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",value);
  this->IntensityMaxDisplay->setText(tempchar);
  delete tempchar;
}

void QtSlicer::LoadImg()
{
  typedef double                            PixelType;
  typedef itk::Image<PixelType, 3>          ImageType;
  typedef itk::ImageFileReader<ImageType>   ReaderType;



  ReaderType::Pointer reader = ReaderType::New();
  
  QString s = QFileDialog::getOpenFileName(".","Images (*.*)", 0, "open file dialog","Chose an image filename" );

  reader->SetFileName( s.latin1() );
  
  std::cout << "loading image " << s.latin1() << " ... ";
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return;
    }
 
   SetInputImage( reader->GetOutput() );
}

void QtSlicer::LoadSeg()
{
  typedef itk::Image<unsigned short, 3>      OverlayType;
    typedef itk::ImageFileReader<OverlayType>   OverlayReaderType;
   OverlayReaderType::Pointer reader2 = OverlayReaderType::New();
  
 QString s = QFileDialog::getOpenFileName(".","Images (*.*)", 0, "open file dialog","Chose an image filename" );

  reader2->SetFileName( s.latin1() );
  
  std::cout << "loading image " << s.latin1() << " ... ";
  try
    {
    reader2->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return;
    }

	SetInputOverlay(reader2->GetOutput() );
}