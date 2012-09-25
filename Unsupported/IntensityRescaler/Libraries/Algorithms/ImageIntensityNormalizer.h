#ifndef __ImageIntensityNormalizer_h_
#define __ImageIntensityNormalizer_h_

  

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkIntensityWindowingImageFilter.h"
#include "IntensityCurveVTK.h"

//#include "ItkFastIO.h"

class ImageIntensityNormalizer
{
public:
  ImageIntensityNormalizer();
  ~ImageIntensityNormalizer();
  
  typedef unsigned short            PixelType;
  typedef itk::Image< PixelType,  3 >      ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;
  typedef ImageType::Pointer          ImagePointer;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;
  typedef itk::ImageDuplicator<ImageType>    DuplicatorType;
  typedef itk::Vector<double,3>        VectorType;
  typedef std::vector<VectorType>        VectorList;
  typedef itk::IntensityWindowingImageFilter<ImageType,ImageType>   WindowType;


  void LoadImage(const char* filename);
  ImagePointer ReadImage(const char* filename);
  void LoadImage(ImagePointer image);
  void SaveImage(ImagePointer image,const char* filename);
  void SaveImage(const char* filename,const char* type);
  void LoadSegmentation(const char* filename);
  void UnloadImage(int index);
  void UnloadSegmentation(int index);
  std::vector<int> ListLabel(ImagePointer segmentation);
  void AddLabel(int label);
  void ClearLabel();
  void SetLabelList(std::vector<int>& _labellist);
  void RemoveLabel(int label);
  ImagePointer GetImage(int index);
  ImagePointer GetSegmentation(int index);
  ImagePointer DuplicateImage(ImagePointer image);
     VectorList* ComputeMean(ImagePointer image,ImagePointer segmentation);
     VectorList* ComputeMax(ImagePointer image,ImagePointer segmentation);
  ImagePointer IntensityWindowing(ImagePointer image,VectorList& m_vectorlist,double sigma,double newmax);
  ImagePointer IntensityWindowing(ImagePointer image,double max,double _newmax);
  ImagePointer IntensityWindowing(ImagePointer image,double inmin,double inmax, double outmin,double outmax);
  
  ImagePointer TargetIntensityWindowing(ImagePointer image,VectorList& m_vectorlist,double sigma,double* _newmax);
  ImagePointer AdjustClasses(ImagePointer target,ImagePointer targetseg,ImagePointer source, ImagePointer sourceseg);
  ImagePointer Normalize(ImagePointer image,int* m_correspondingarray);
  int* ComputeHistogram(ImagePointer image,int* _maxx,int* _maxy);
  int* ComputeDistribution(ImagePointer image,ImagePointer seg,int* _maxx,int* _maxy,int label);

  std::vector<ImagePointer> m_SourceImage;
  std::vector<ImagePointer> m_SegImage;
  std::vector<int> m_labellist;
  IntensityCurveVTK::Pointer m_curve;
};


#endif // __ImageIntensityNormalizer_h_
