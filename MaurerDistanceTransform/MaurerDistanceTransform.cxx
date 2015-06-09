#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "MaurerDistanceTransformCLP.h"

int main(int argc, char * argv[])
{
    typedef itk::Image<unsigned char , 3> UnsignedCharImageType;
    typedef itk::Image<float , 3> FloatImageType;
    PARSE_ARGS;
    UnsignedCharImageType::Pointer image = UnsignedCharImageType::New();
    typedef itk::ImageFileReader<UnsignedCharImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(input);
    reader->Update();
    image = reader->GetOutput();

    typedef itk::SignedMaurerDistanceMapImageFilter< UnsignedCharImageType, FloatImageType > SignedMaurerDistanceMapImageFilterType;
    SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter =
            SignedMaurerDistanceMapImageFilterType::New();
    distanceMapImageFilter->SetInput( image ) ;
    distanceMapImageFilter->Update();
    typedef itk::ImageFileWriter< FloatImageType > WriterType;
    WriterType::Pointer writer = WriterType::New() ;
    writer->SetInput( distanceMapImageFilter->GetOutput() ) ;
    writer->SetFileName( output ) ;
    writer->Update() ;
    return EXIT_SUCCESS ;
}
