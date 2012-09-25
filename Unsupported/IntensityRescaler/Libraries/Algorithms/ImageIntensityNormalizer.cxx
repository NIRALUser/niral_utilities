#include "ImageIntensityNormalizer.h"

#define NL_DEBUG false

ImageIntensityNormalizer::ImageIntensityNormalizer()
{

}

ImageIntensityNormalizer::~ImageIntensityNormalizer()
{

}

ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::ReadImage(const char* filename)
{
  ReaderType::Pointer m_reader = ReaderType::New();
  m_reader->SetFileName(filename);
  try
  {
      m_reader->Update();
  }
  catch (...)
  {
     return 0;
  }

  return m_reader->GetOutput();
}

void ImageIntensityNormalizer::LoadImage(const char* filename)
{
  ReaderType::Pointer m_reader = ReaderType::New();
  m_reader->SetFileName(filename);
  m_reader->Update();
  m_SourceImage.push_back(m_reader->GetOutput());
}

void ImageIntensityNormalizer::LoadImage(ImagePointer image)
{
  m_SourceImage.push_back(image);
}

void ImageIntensityNormalizer::SaveImage(ImagePointer image,const char* filename)
{
  WriterType::Pointer m_writer = WriterType::New();
  m_writer->SetInput(image);
  m_writer->SetFileName(filename);
  m_writer->Update();
}



ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::GetImage(int index)
{
  return m_SourceImage[index];
}

ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::GetSegmentation(int index)
{
  return m_SegImage[index];
}

void ImageIntensityNormalizer::SaveImage(const char* filename,const char* type)
{
}

void ImageIntensityNormalizer::LoadSegmentation(const char* filename)
{
  ReaderType::Pointer m_reader = ReaderType::New();
  m_reader->SetFileName(filename);
  m_reader->Update();
  m_SegImage.push_back(m_reader->GetOutput());
}

void ImageIntensityNormalizer::UnloadImage(int index)
{
     m_SourceImage.erase(m_SourceImage.begin() + index); 
}
  
void ImageIntensityNormalizer::UnloadSegmentation(int index)
{
     m_SegImage.erase(m_SegImage.begin() + index); 
}


void ImageIntensityNormalizer::AddLabel(int label)
{
  m_labellist.push_back(label);
}

void ImageIntensityNormalizer::ClearLabel()
{
  m_labellist.clear();
}

void ImageIntensityNormalizer::SetLabelList(std::vector<int>& _labellist)
{
  for (std::vector<int>::iterator j = _labellist.begin(); j != _labellist.end(); j++)
  {
    m_labellist.push_back(*j);
  }
}

void ImageIntensityNormalizer::RemoveLabel(int label)
{
  for (std::vector<int>::iterator j = m_labellist.begin(); j != m_labellist.end(); j++)
  {
      if ((*j) == label) m_labellist.erase(j);
  }
}


ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::DuplicateImage(ImagePointer image)
{
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(image);
    duplicator->Update();
  return duplicator->GetOutput();
}


std::vector<int> ImageIntensityNormalizer::ListLabel(ImagePointer segmentation)
{
  std::vector<int> labellist;
  bool label[256];
  double numlabel[256];
  for (int i=0;i<256;i++)
  {
    label[i]=false;
    numlabel[i]=0;
  }

  IteratorType m_itimage(segmentation,segmentation->GetLargestPossibleRegion());
  ImageType::SizeType m_imagesize = (segmentation->GetLargestPossibleRegion().GetSize());
  
  m_itimage.GoToBegin();
  for (unsigned int z=0;z<m_imagesize[2];z++)
    for (unsigned int y=0;y<m_imagesize[1];y++)
      for (unsigned int x=0;x<m_imagesize[0];x++)
      {
        if (m_itimage.Get() != 0)
        {
          label[(unsigned char)m_itimage.Get()] = true;
          numlabel[(unsigned char)m_itimage.Get()]++;
        }
        ++m_itimage;
      }

  for (int i=0;i<256;i++)
  {
    if (label[i]) 
    {
      labellist.push_back(i);
    }
  }
  return labellist;
}



ImageIntensityNormalizer::VectorList* ImageIntensityNormalizer::ComputeMean(ImagePointer image,ImagePointer segmentation)
{
  

  VectorList* m_vectorlist= new VectorList();
  
  double label[256];
  double numlabel[256];

  for (int i=0;i<256;i++)
  {
    label[i]=0;
    numlabel[i]=0;
  }

  IteratorType m_itimage(image,image->GetLargestPossibleRegion());
  ImageType::SizeType m_imagesize = (image->GetLargestPossibleRegion().GetSize());
  
  IteratorType m_itsegmentation(segmentation,segmentation->GetLargestPossibleRegion());

  m_itsegmentation.GoToBegin();
  m_itimage.GoToBegin();
  for (unsigned int z=0;z<m_imagesize[2];z++)
  {
    for (unsigned int y=0;y<m_imagesize[1];y++)
      for (unsigned int x=0;x<m_imagesize[0];x++)
      {
        if (m_itsegmentation.Get() != 0)
        {
        bool ok = false;
        for (std::vector<int>::iterator j = m_labellist.begin(); j != m_labellist.end(); j++)
        {
          if ((*j) == m_itsegmentation.Get()) ok = true;
        }
        
        if (ok)
        {
            label[(unsigned char)m_itsegmentation.Get()] += m_itimage.Get();
            numlabel[(unsigned char)m_itsegmentation.Get()]++;  
        }  
        }
        ++m_itsegmentation;
        ++m_itimage;
      }
  }



  double mean;
  double dist;

  //Compute variance based on distribution function
  for (std::vector<int>::iterator j = m_labellist.begin(); j != m_labellist.end(); j++)
  {
    int* m_histo;
    int maxx,maxy;
    m_histo= ComputeDistribution(image,segmentation,&maxx,&maxy,(*j));
  
    if (NL_DEBUG)  std::cout << " Maxx: " << maxx << " and Maxy: " << maxy << std::endl;

    m_itsegmentation.GoToBegin();
    m_itimage.GoToBegin();
    double variance = 0;
    double numvar = 0;
    for (unsigned int z=0;z<m_imagesize[2];z++)
    {
      for (unsigned int y=0;y<m_imagesize[1];y++)
        for (unsigned int x=0;x<m_imagesize[0];x++)
        {
          if (m_itsegmentation.Get() == (*j))
          {
            mean = label[(*j)]/numlabel[(*j)];
            if (m_itimage.Get()>mean)
            {
              dist = (m_itimage.Get()-mean);
              variance += (dist*dist);
              numvar++;
            }
          }
          ++m_itsegmentation;
          ++m_itimage;
        }
    }

    VectorType m_vector;
    m_vector[0] = (*j);
    m_vector[1] = label[(*j)]/numlabel[(*j)]; 
    m_vector[2] = sqrt(variance/numvar);
    m_vectorlist->push_back(m_vector);
  }

  return m_vectorlist;
}



ImageIntensityNormalizer::VectorList* ImageIntensityNormalizer::ComputeMax(ImagePointer image,ImagePointer segmentation)
{
  

  VectorList* m_vectorlist= new VectorList();
  
  double label[256];
  double numlabel[256];

  for (int i=0;i<256;i++)
  {
    label[i]=0;
    numlabel[i]=0;
  }

  IteratorType m_itimage(image,image->GetLargestPossibleRegion());
  ImageType::SizeType m_imagesize = (image->GetLargestPossibleRegion().GetSize());
  
  IteratorType m_itsegmentation(segmentation,segmentation->GetLargestPossibleRegion());


  m_itsegmentation.GoToBegin();
  m_itimage.GoToBegin();
  for (unsigned int z=0;z<m_imagesize[2];z++)
  {
    for (unsigned int y=0;y<m_imagesize[1];y++)
      for (unsigned int x=0;x<m_imagesize[0];x++)
      {
        if (m_itsegmentation.Get() != 0)
        {
        bool ok = false;
        for (std::vector<int>::iterator j = m_labellist.begin(); j != m_labellist.end(); j++)
        {
          if ((*j) == m_itsegmentation.Get()) ok = true;
        }
        
        if (ok)
        {
            label[(unsigned char)m_itsegmentation.Get()] += m_itimage.Get();
            numlabel[(unsigned char)m_itsegmentation.Get()]++;  
        }  
        }
        ++m_itsegmentation;
        ++m_itimage;
      }
  }



  double mean;
  double dist;

  if (NL_DEBUG) std::cout << "LabelList Size: " << m_labellist.size() <<std::endl;

  //Compute variance based on distribution function
  for (std::vector<int>::iterator j = m_labellist.begin(); j != m_labellist.end(); j++)
  {
    int* m_histo;
    int maxx,maxy;
    m_histo= ComputeDistribution(image,segmentation,&maxx,&maxy,(*j));
  
    if (NL_DEBUG) std::cout << "Maxx: " << maxx << "and maxy: " << maxy << std::endl;

    m_itsegmentation.GoToBegin();
    m_itimage.GoToBegin();
    double variance = 0;
    double numvar = 0;
    for (unsigned int z=0;z<m_imagesize[2];z++)
    {
      for (unsigned int y=0;y<m_imagesize[1];y++)
        for (unsigned int x=0;x<m_imagesize[0];x++)
        {
          if (m_itsegmentation.Get() == (*j))
          {
            mean = label[(*j)]/numlabel[(*j)];
            if (m_itimage.Get()>maxx)
            {
              dist = (m_itimage.Get()-maxx);
              variance += (dist*dist);
              numvar++;
            }
          }
          ++m_itsegmentation;
          ++m_itimage;
        }
    }

    VectorType m_vector;
    m_vector[0] = (*j);
    m_vector[1] = maxx; //label[(*j)]/numlabel[(*j)]; 
    m_vector[2] = sqrt(variance/numvar);
    m_vectorlist->push_back(m_vector);
  }

  return m_vectorlist;
}


ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::IntensityWindowing(ImagePointer image,VectorList& m_vectorlist,double sigma,double _newmax)
{
  typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
  CalculatorType::Pointer calculator  =  CalculatorType::New();
  calculator->SetImage(image );
  calculator->Compute();
  int minbefore = calculator->GetMinimum();
  int maxbefore = calculator->GetMaximum(); 

  if (NL_DEBUG) std::cout << "Min target before: " << minbefore <<std::endl;
  if (NL_DEBUG)  std::cout << "Max target before: " << maxbefore <<std::endl;

  //find min and max mean
  //int min= calculator->GetMinimum();
  //int max = calculator->GetMaximum();
  
  int newmin;
  int newmax;

  double meanmin=99999999;
  double meanmax=0;

  //Fix maximum value
  for (VectorList::iterator j = m_vectorlist.begin(); j != m_vectorlist.end(); j++)
  {
    if ((*j)[1] < meanmin) 
    {
      newmin = (int)((*j)[1]-(*j)[2]);
      meanmin = (*j)[1];
    }

    if ((*j)[1]>meanmax) 
    {
      newmax = (int)((*j)[1]+(sigma*(*j)[2]));
      meanmax = (*j)[1];
    }
  }

  if (maxbefore<newmax) newmax = maxbefore;
  if (NL_DEBUG) std::cout << "_Newmax: " << _newmax << " | Newmax: " << newmax << std::endl;

  WindowType::Pointer filter = WindowType::New();
  filter->SetInput(image);
  filter->SetOutputMinimum(0);
  filter->SetOutputMaximum((int)_newmax);
  filter->SetWindowMinimum(0);
  filter->SetWindowMaximum(newmax);

  filter->Update();

  calculator  =  CalculatorType::New();
  calculator->SetImage(filter->GetOutput());
  calculator->Compute();
  if (NL_DEBUG) std::cout << "Min target after: " << calculator->GetMinimum() <<std::endl;
  if (NL_DEBUG) std::cout << "Max target after: " << calculator->GetMaximum() <<std::endl;


  return filter->GetOutput();
}


ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::IntensityWindowing(ImagePointer image,double inmin,double inmax,double outmin,double outmax)
{
  typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
  CalculatorType::Pointer calculator  =  CalculatorType::New();
  calculator->SetImage(image );
  calculator->Compute();
  int minbefore = calculator->GetMinimum();
  int maxbefore = calculator->GetMaximum(); 

  if (NL_DEBUG) std::cout << "Min target before: " << minbefore <<std::endl;
  if (NL_DEBUG) std::cout << "Max target before: " << maxbefore <<std::endl;

  WindowType::Pointer filter = WindowType::New();
  filter->SetInput(image);
  filter->SetOutputMinimum((int)outmin);
  filter->SetOutputMaximum((int)outmax);
  filter->SetWindowMinimum((int)inmin);
  filter->SetWindowMaximum((int)inmax);

  filter->Update();

  calculator  =  CalculatorType::New();
  calculator->SetImage(filter->GetOutput());
  calculator->Compute();
  if (NL_DEBUG) std::cout << "Min target after: " << calculator->GetMinimum() <<std::endl;
  if (NL_DEBUG) std::cout << "Max target after: " << calculator->GetMaximum() <<std::endl;


  return filter->GetOutput();
}

ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::IntensityWindowing(ImagePointer image,double max,double _newmax)
{
  typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
  CalculatorType::Pointer calculator  =  CalculatorType::New();
  calculator->SetImage(image );
  calculator->Compute();
  int minbefore = calculator->GetMinimum();
  int maxbefore = calculator->GetMaximum(); 

  if (NL_DEBUG) std::cout << "Min target before: " << minbefore <<std::endl;
  if (NL_DEBUG) std::cout << "Max target before: " << maxbefore <<std::endl;


  WindowType::Pointer filter = WindowType::New();
  filter->SetInput(image);
  filter->SetOutputMinimum(0);
  filter->SetOutputMaximum((int)_newmax);
  filter->SetWindowMinimum(0);
  filter->SetWindowMaximum((int)max);

  filter->Update();

  calculator  =  CalculatorType::New();
  calculator->SetImage(filter->GetOutput());
  calculator->Compute();
  if (NL_DEBUG) std::cout << "Min target after: " << calculator->GetMinimum() <<std::endl;
  if (NL_DEBUG) std::cout << "Max target after: " << calculator->GetMaximum() <<std::endl;


  return filter->GetOutput();
}


ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::TargetIntensityWindowing(ImagePointer image,VectorList& m_vectorlist,double sigma,double* _newmax)
{
  typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
  CalculatorType::Pointer calculator  =  CalculatorType::New();
  calculator->SetImage(image);
  calculator->Compute();
  if (NL_DEBUG) std::cout << "Min source before: " << calculator->GetMinimum() <<std::endl;
  if (NL_DEBUG) std::cout << "Max source before: " << calculator->GetMaximum() <<std::endl;

  
  double newmax;
  double meanmax=0;

  //Fix maximum value
  for (VectorList::iterator j = m_vectorlist.begin(); j != m_vectorlist.end(); j++)
  {
    if ((*j)[1]>meanmax) 
    {
      if (NL_DEBUG) std::cout << "Mean target " << (*j)[0] << " : " << sigma << " | " << (*j)[1] << " | " << (*j)[2] << std::endl;
      newmax = (int) ((*j)[1] + (sigma*((*j)[2])));
      meanmax = (*j)[1];
    }
  }

  if (newmax > calculator->GetMaximum()) newmax = calculator->GetMaximum();

  WindowType::Pointer filter = WindowType::New();
  filter->SetInput(image);
  filter->SetWindowMinimum(0);
  filter->SetWindowMaximum((int)newmax);
  filter->SetOutputMinimum(0);
  filter->SetOutputMaximum((int)newmax);

  filter->Update();

  calculator  =  CalculatorType::New();
  calculator->SetImage(filter->GetOutput());
  calculator->Compute();
  if (NL_DEBUG) std::cout << "Min source after: " << calculator->GetMinimum() <<std::endl;
  if (NL_DEBUG) std::cout << "Max source after: " << calculator->GetMaximum() <<std::endl;

  *_newmax = (int)newmax;
  if (NL_DEBUG) std::cout << "_Newmax: " << *_newmax << std::endl;

  return filter->GetOutput();
}


ImageIntensityNormalizer::ImagePointer  ImageIntensityNormalizer::AdjustClasses(ImagePointer target,ImagePointer targetseg, ImagePointer source,  ImagePointer sourceseg)
{
  ImageIntensityNormalizer::VectorList* m_vectorlisttarget;
  ImageIntensityNormalizer::VectorList* m_vectorlistsource;
  
  //Compute Mean values for target
  m_vectorlisttarget = ComputeMean(target,targetseg);

  //Compute Mean values for target
  m_vectorlistsource = ComputeMean(source,sourceseg);
  
  for (ImageIntensityNormalizer::VectorList::iterator j = m_vectorlistsource->begin(); j != m_vectorlistsource->end(); j++)
  {
    if (NL_DEBUG) std::cout <<"Labels: " <<(*j)[0] << " | " << (*j)[1]  << " | " << (*j)[2]  << std::endl;
  }

  //Get number of control points
  int m_nbcontrols = m_vectorlistsource->size();
  
  //Copmute Max of these images
  typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
  CalculatorType::Pointer calculator  =  CalculatorType::New();
  calculator->SetImage(target);
  calculator->Compute();
  double m_imagemax = calculator->GetMaximum();


  //Create VTK curve 
  m_curve = IntensityCurveVTK::New();
  m_curve->Initialize(m_nbcontrols+2);

  VectorList::iterator j = m_vectorlistsource->begin();
  VectorList::iterator k = m_vectorlisttarget->begin();
  for (int i=0;i<m_nbcontrols;i++)
  {
    m_curve->UpdateControlPoint(i+1,(*j)[1]/m_imagemax,(*k)[1]/m_imagemax);  
    j++;
    k++;
  }

  //Update image with corresponding points
  int* m_correspondingarray;
  m_correspondingarray = (int*) malloc((int)(m_imagemax+1)*sizeof(int));

  float t = 0.0;
    float tStep = 1.0f / m_imagemax;
  for (int i=0;i<=m_imagemax;i++)
  {
    m_correspondingarray[i] = (int)(m_curve->Evaluate(t)*m_imagemax);
      t+=tStep;
  }

  ImagePointer m_image;
    m_image = Normalize(source,m_correspondingarray);

  //Compute Mean values for target
  m_vectorlistsource = ComputeMean(m_image,sourceseg);

  for (ImageIntensityNormalizer::VectorList::iterator j = m_vectorlistsource->begin(); j != m_vectorlistsource->end(); j++)
  {
    if (NL_DEBUG) std::cout <<"Labels: " <<(*j)[0] << " | " << (*j)[1]  << " | " << (*j)[2]  << std::endl;
  }

  return m_image;
}



ImageIntensityNormalizer::ImagePointer ImageIntensityNormalizer::Normalize(ImagePointer image,int* m_correspondingarray)
{
  ImagePointer m_image;
  m_image = DuplicateImage(image);

  IteratorType itS(m_image,m_image->GetLargestPossibleRegion());
  IteratorType itSorig(image,image->GetLargestPossibleRegion());


  itS.GoToBegin();
  itSorig.GoToBegin();

  while (!itS.IsAtEnd())
  {
    itS.Set(m_correspondingarray[(int)itSorig.Get()]);
    ++itS;
    ++itSorig;
  }

  return m_image;
}


int* ImageIntensityNormalizer::ComputeHistogram(ImagePointer image,int* _maxx,int* _maxy)
{

  int* m_histo = new int[65536];

  //Initialize m_histo
  for (int i=0;i<65535;i++)
      m_histo[i] =0; 

  int maxx = 0;
  int maxy = 0;

  IteratorType m_itimage(image,image->GetLargestPossibleRegion());
  //ImageType::SizeType m_imagesize = (image->GetLargestPossibleRegion().GetSize());
  
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


int* ImageIntensityNormalizer::ComputeDistribution(ImagePointer image,ImagePointer seg,int* _maxx,int* _maxy,int label)
{
  int* m_histo = new int[65536];

  //Initialize m_histo
  for (int i=0;i<65535;i++)
      m_histo[i] = 0; 

  int maxx = 0;
  int maxy = 0;

  IteratorType m_itimage(image,image->GetLargestPossibleRegion());
  IteratorType m_itseg(seg,seg->GetLargestPossibleRegion());
  
    
  if (NL_DEBUG) std::cout << "Compute distibution label for label" << label << std::endl;
  
  m_itimage.GoToBegin();
  m_itseg.GoToBegin();  
  int val;
  while (!m_itimage.IsAtEnd())
  {
    if ((int)m_itseg.Get() == label)
    {
      val = (int)m_itimage.Get();
      // only positive values
      if (val < 0) val = 0;
      if (val > 65535) val = 65535;
      m_histo[val] += m_itimage.Get();
  
      if(m_histo[val]>maxy) 
      {
        maxx = val;
        maxy = m_histo[val];
      }
    }
  
    ++m_itimage;
    ++m_itseg;
  }

  *_maxx = maxx;
  *_maxy = maxy;
  return m_histo;
}
