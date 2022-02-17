#include <iostream>
#include <fstream>
#include <iomanip>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkImageRegionIterator.h"
#include "DWI_NiftiNrrdConversionCLP.h"

typedef itk::Image<unsigned short, 4> ImageType;
typedef itk::VectorImage < unsigned short , 3 > VectorImageType;
typedef itk::ImageFileReader <ImageType> ReaderType;
typedef itk::ImageFileReader <VectorImageType> ReaderVectorType;
typedef itk::ImageFileWriter <ImageType> WriterType;
typedef itk::ImageFileWriter <VectorImageType> WriterVectorType;


void ImageToVector(ImageType::Pointer& ImageIn, VectorImageType::Pointer& ImageOut)
{
	ImageType::RegionType RegionIn=ImageIn->GetLargestPossibleRegion();
	VectorImageType::RegionType RegionOut;
	
	ImageType::SizeType RegionSizeIn=RegionIn.GetSize();
	VectorImageType::SizeType RegionSizeOut;
	RegionSizeOut[0]=RegionSizeIn[0];
	RegionSizeOut[1]=RegionSizeIn[1];
	RegionSizeOut[2]=RegionSizeIn[2];
	ImageType::IndexType RegionIndexIn=RegionIn.GetIndex();
	VectorImageType::IndexType RegionIndexOut;
	RegionIndexOut[0]=RegionIndexIn[0];
	RegionIndexOut[1]=RegionIndexIn[1];
	RegionIndexOut[2]=RegionIndexIn[2];
	ImageType::PointType OriginIn=ImageIn->GetOrigin();
	VectorImageType::PointType OriginOut;
	OriginOut[0]=OriginIn[0];
	OriginOut[1]=OriginIn[1];
	OriginOut[2]=OriginIn[2];
	ImageType::SpacingType SpacingIn=ImageIn->GetSpacing();
	VectorImageType::SpacingType SpacingOut;
	SpacingOut[0]=SpacingIn[0];
	SpacingOut[1]=SpacingIn[1];
	SpacingOut[2]=SpacingIn[2];
	ImageType::DirectionType DirectionIn=ImageIn->GetDirection();
	VectorImageType::DirectionType DirectionOut;
	DirectionOut[0][0]=DirectionIn[0][0];
	DirectionOut[0][1]=DirectionIn[0][1];
	DirectionOut[0][2]=DirectionIn[0][2];
	DirectionOut[1][0]=DirectionIn[1][0];
	DirectionOut[1][1]=DirectionIn[1][1];
	DirectionOut[1][2]=DirectionIn[1][2];
	DirectionOut[2][0]=DirectionIn[2][0];
	DirectionOut[2][1]=DirectionIn[2][1];
	DirectionOut[2][2]=DirectionIn[2][2];
	
	RegionOut.SetSize(RegionSizeOut);
	RegionOut.SetIndex(RegionIndexOut);
	ImageOut->SetRegions(RegionOut);
	ImageOut->SetOrigin(OriginOut);
	ImageOut->SetSpacing(SpacingOut);
	ImageOut->SetDirection(DirectionOut);
	ImageOut->SetVectorLength(RegionSizeIn[3]);
	ImageOut->Allocate();
	
	itk::ImageRegionIterator<VectorImageType> Out(ImageOut, ImageOut->GetLargestPossibleRegion());
	itk::ImageRegionIterator<ImageType> In(ImageIn, ImageIn->GetLargestPossibleRegion());

	
	itk::ImageRegionIterator<ImageType>::IndexType Id,Start;
	itk::VariableLengthVector<unsigned short> Value;
	Value.SetSize(RegionSizeIn[3]);
	In.GoToBegin();
	Start=In.GetIndex();
	for(Out.GoToBegin(); !Out.IsAtEnd(); ++Out)
	{
		for(unsigned int i=0; i<RegionSizeIn[3]; i++)
		{
			Id=In.GetIndex();
			Value.SetElement(i,In.Get());
			Id[3]++;
			In.SetIndex(Id);
		}
		In.SetIndex(Start);
		++In;
		Start=In.GetIndex();
		Out.Set(Value);
	}
}

void VectorToImage(VectorImageType::Pointer& ImageIn, ImageType::Pointer& ImageOut)
{
	VectorImageType::RegionType RegionIn=ImageIn->GetLargestPossibleRegion();
	ImageType::RegionType RegionOut;
	
	VectorImageType::SizeType RegionSizeIn=RegionIn.GetSize();
	ImageType::SizeType RegionSizeOut;
	RegionSizeOut[0]=RegionSizeIn[0];
	RegionSizeOut[1]=RegionSizeIn[1];
	RegionSizeOut[2]=RegionSizeIn[2];
	RegionSizeOut[3]=ImageIn->GetVectorLength();
	VectorImageType::IndexType RegionIndexIn=RegionIn.GetIndex();
	ImageType::IndexType RegionIndexOut;
	RegionIndexOut[0]=RegionIndexIn[0];
	RegionIndexOut[1]=RegionIndexIn[1];
	RegionIndexOut[2]=RegionIndexIn[2];
	RegionIndexOut[3]=0;
	VectorImageType::PointType OriginIn=ImageIn->GetOrigin();
	ImageType::PointType OriginOut;
	OriginOut[0]=OriginIn[0];
	OriginOut[1]=OriginIn[1];
	OriginOut[2]=OriginIn[2];
	OriginOut[3]=0;
	VectorImageType::SpacingType SpacingIn=ImageIn->GetSpacing();
	ImageType::SpacingType SpacingOut;
	SpacingOut[0]=SpacingIn[0];
	SpacingOut[1]=SpacingIn[1];
	SpacingOut[2]=SpacingIn[2];
	SpacingOut[3]=RegionSizeIn[0];
	VectorImageType::DirectionType DirectionIn=ImageIn->GetDirection();
	ImageType::DirectionType DirectionOut;
	DirectionOut[0][0]=DirectionIn[0][0];
	DirectionOut[0][1]=DirectionIn[0][1];
	DirectionOut[0][2]=DirectionIn[0][2];
	DirectionOut[1][0]=DirectionIn[1][0];
	DirectionOut[1][1]=DirectionIn[1][1];
	DirectionOut[1][2]=DirectionIn[1][2];
	DirectionOut[2][0]=DirectionIn[2][0];
	DirectionOut[2][1]=DirectionIn[2][1];
	DirectionOut[2][2]=DirectionIn[2][2];
	DirectionOut[3][0]=0;
	DirectionOut[3][1]=0;
	DirectionOut[3][2]=0;
	DirectionOut[3][3]=1;
	
	RegionOut.SetSize(RegionSizeOut);
	RegionOut.SetIndex(RegionIndexOut);
	ImageOut->SetRegions(RegionOut);
	ImageOut->SetOrigin(OriginOut);
	ImageOut->SetSpacing(SpacingOut);
	ImageOut->SetDirection(DirectionOut);
	ImageOut->Allocate();
	
	itk::ImageRegionIterator<ImageType> Out(ImageOut, ImageOut->GetLargestPossibleRegion());
	itk::ImageRegionIterator<VectorImageType> In(ImageIn, ImageIn->GetLargestPossibleRegion());

	itk::ImageRegionIterator<ImageType>::IndexType Id, Start;
	Out.GoToBegin();
	Start=Out.GetIndex();
	for(In.GoToBegin(); !In.IsAtEnd(); ++In)
	{
		itk::VariableLengthVector<unsigned short> Value=In.Get();
		for(unsigned int i=0; i<ImageIn->GetVectorLength(); i++)
		{
			Out.Set(Value[i]);
			Id=Out.GetIndex();
			Id[3]++;
			Out.SetIndex(Id);
		}
		Out.SetIndex(Start);
		++Out;
		Start=Out.GetIndex();
	}
}

double GetMax(std::string Values)
{
	std::istringstream iss(Values);
	double Max;
	iss>>Max;
	while(!iss.eof())
	{
		double Value;
		iss>>Value;
		if(Value>Max)
			Max=Value;
	}
	return Max;
}

void GetFileNameInfo(std::string FullFileName, std::string& Path, std::string& FileName, std::string& Extension)
{
	Path=FullFileName.substr(0,FullFileName.find_last_of("/")+1);
	
	Extension=FullFileName.substr(FullFileName.find_last_of(".")+1,FullFileName.size()-FullFileName.find_last_of(".")+1);
	
	FileName=FullFileName.substr(FullFileName.find_last_of("/")+1,FullFileName.size()-FullFileName.find_last_of("/")-Extension.size()-2);
}
	

int main(int argc, char* argv[])
{
	PARSE_ARGS;
	
	std::string PathIn, FileNameIn, ExtensionIn;
	std::string PathOut, FileNameOut, ExtensionOut;
	GetFileNameInfo(input,PathIn,FileNameIn,ExtensionIn);
	GetFileNameInfo(output,PathOut,FileNameOut,ExtensionOut);
	
	if(ExtensionIn=="nii" && ExtensionOut=="nrrd")
	{
		std::ifstream bvalfile(b_vals.c_str(), std::ios::in);
		std::string bval;
		double BValueMax;
		if(bvalfile)
		{
			getline(bvalfile, bval);
			BValueMax=GetMax(bval);
			bvalfile.close();
		}
		else
		{
			std::cerr<<"Cannot open B-Values file."<<std::endl;
			return 0;
		}
		
		std::ifstream bvectfile(b_vects.c_str(), std::ios::in);
		std::string bvect[3];
		if(bvectfile)
		{
			getline(bvectfile, bvect[0]);
			getline(bvectfile, bvect[1]);
			getline(bvectfile, bvect[2]);
			bvectfile.close();
		}
		else
		{
			std::cerr<<"Cannot open B-Vectors file."<<std::endl;
			return 0;
		}
		
		ReaderType::Pointer reader=ReaderType::New();
		reader->SetFileName(input);
		reader->Update();
		
		ImageType::Pointer ImageIn=reader->GetOutput();
		VectorImageType::Pointer ImageOut=VectorImageType::New();
		ImageToVector(ImageIn,ImageOut);
		
		itk::MetaDataDictionary dict;
		
		std::istringstream bvectiss[3], bvaliss(bval);
		std::ostringstream oss;
		bvectiss[0].str(bvect[0].c_str());
		bvectiss[1].str(bvect[1].c_str());
		bvectiss[2].str(bvect[2].c_str());
		
		float Gradient[3], BValue;
		bvaliss>>BValue;
		bvectiss[0]>>Gradient[0];
		bvectiss[1]>>Gradient[1];
		bvectiss[2]>>Gradient[2];
		
		int Count=0;
		while(!bvectiss[0].eof() && !bvectiss[1].eof() && !bvectiss[2].eof() && !bvaliss.eof())
		{
			double ScaleFactor=sqrt(BValue/BValueMax);
			
			oss<<std::setfill('0')<<std::setw(4)<<Count;
			std::string key="DWMRI_gradient_"+oss.str();
			oss.str("");
			
			oss<<Gradient[0]*ScaleFactor<<" "<<Gradient[1]*ScaleFactor<<" "<<Gradient[2]*ScaleFactor;
			itk::EncapsulateMetaData<std::string>(dict,key,oss.str());
			
			bvaliss>>BValue;
			bvectiss[0]>>Gradient[0];
			bvectiss[1]>>Gradient[1];
			bvectiss[2]>>Gradient[2];
			Count++;
			oss.str("");
		}
		itk::EncapsulateMetaData<std::string>(dict,"modality","DWMRI");
		itk::EncapsulateMetaData<std::string>(dict,"ITK_InputFilterName","NrrdImageIO");
		oss<<BValueMax;
		itk::EncapsulateMetaData<std::string>(dict,"DWMRI_b-value",oss.str());
		ImageOut->SetMetaDataDictionary(dict);
		
		WriterVectorType::Pointer writer=WriterVectorType::New();
		writer->SetFileName(output);
		writer->SetInput( ImageOut );
		
		try{
			writer->Update();
		}
		catch(itk::ExceptionObject & err)
		{
			std::cerr<<"ExceptionObject caught!" <<std::endl;
			std::cerr<<err<<std::endl;
		}
	}
	else if(ExtensionIn=="nrrd" && ExtensionOut=="nii")
	{
		ReaderVectorType::Pointer reader=ReaderVectorType::New();
		reader->SetFileName(input);
		reader->Update();
		
		VectorImageType::Pointer ImageIn=reader->GetOutput();
		ImageType::Pointer ImageOut=ImageType::New();
		VectorToImage(ImageIn,ImageOut);
		
		itk::MetaDataDictionary DictIn=ImageIn->GetMetaDataDictionary();
		
		std::string BValFileName, BVectFileName, BVect[3];
		double XYZ[3], XYZNorm[3], Norm;
		int BValueMax, B;
		
		BValFileName=PathOut+FileNameOut+"_bvals";
		BVectFileName=PathOut+FileNameOut+"_bvects";
		std::ofstream BValFile(BValFileName.c_str(), std::ios::out);
		std::ofstream BVectFile(BVectFileName.c_str(), std::ios::out);
		std::vector<std::string> Keys=DictIn.GetKeys();
		for(unsigned int i=0; i<Keys.size(); i++)
		{
			if(Keys[i].find("DWMRI_b-value") != std::string::npos)
			{
				std::string Value;
				itk::ExposeMetaData<std::string>(DictIn,Keys[i],Value);
				std::istringstream iss(Value.c_str());
				iss>>BValueMax;
			}
		}
		
		for(unsigned int i=0; i<Keys.size(); i++)
		{
			if(Keys[i].find("DWMRI_gradient") != std::string::npos)
			{
				std::string Vect;
				itk::ExposeMetaData<std::string>(DictIn,Keys[i],Vect);
				
				std::istringstream iss(Vect.c_str());
				std::ostringstream oss;
				
				iss>>XYZ[0];
				iss>>XYZ[1];
				iss>>XYZ[2];
				
				Norm=sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2]);
				if(Norm!=0)
				{
					XYZNorm[0]=XYZ[0]/Norm;
					XYZNorm[1]=XYZ[1]/Norm;
					XYZNorm[2]=XYZ[2]/Norm;
					B=((XYZ[0]*XYZ[0])/(XYZNorm[0]*XYZNorm[0]))*BValueMax;
				}
				else
				{
					XYZNorm[0]=0;
					XYZNorm[1]=0;
					XYZNorm[2]=0;
					B=0;
				}
				
				BValFile<<B<<" ";
				
				oss<<XYZNorm[0]<<" ";
				BVect[0]+=oss.str();
				oss.str("");
				oss<<XYZNorm[1]<<" ";
				BVect[1]+=oss.str();
				oss.str("");
				oss<<XYZNorm[2]<<" ";
				BVect[2]+=oss.str();
			}
		}
		
		BVectFile<<BVect[0]<<std::endl;
		BVectFile<<BVect[1]<<std::endl;
		BVectFile<<BVect[2]<<std::endl;
		
		itk::MetaDataDictionary DictOut;
		std::ostringstream Value;
		
		itk::EncapsulateMetaData<std::string>(DictOut,"ITK_InputFilterName","NiftiImageIO");
		
		itk::EncapsulateMetaData<std::string>(DictOut,"dim[0]","4");
		ImageType::SizeType RegionSize=ImageOut->GetLargestPossibleRegion().GetSize();
		Value<<RegionSize[0];
		itk::EncapsulateMetaData<std::string>(DictOut,"dim[1]",Value.str());
		Value.str("");
		Value<<RegionSize[1];
		itk::EncapsulateMetaData<std::string>(DictOut,"dim[2]",Value.str());
		Value.str("");
		Value<<RegionSize[2];
		itk::EncapsulateMetaData<std::string>(DictOut,"dim[3]",Value.str());
		Value.str("");
		Value<<RegionSize[3];
		itk::EncapsulateMetaData<std::string>(DictOut,"dim[4]",Value.str());
		
		itk::EncapsulateMetaData<std::string>(DictOut,"pixdim[0]","1");
		ImageType::SpacingType Spacing=ImageOut->GetSpacing();
		Value.str("");
		Value<<Spacing[0];
		itk::EncapsulateMetaData<std::string>(DictOut,"pixdim[1]",Value.str());
		Value.str("");
		Value<<Spacing[1];
		itk::EncapsulateMetaData<std::string>(DictOut,"pixdim[2]",Value.str());
		Value.str("");
		Value<<Spacing[2];
		itk::EncapsulateMetaData<std::string>(DictOut,"pixdim[3]",Value.str());
		Value.str("");
		Value<<Spacing[3];
		itk::EncapsulateMetaData<std::string>(DictOut,"pixdim[4]",Value.str());
		
		ImageOut->SetMetaDataDictionary(DictOut);
		
		WriterType::Pointer writer=WriterType::New();
		writer->SetFileName(output);
		writer->SetInput( ImageOut );
		
		try{
			writer->Update();
		}
		catch(itk::ExceptionObject & err)
		{
			std::cerr<<"ExceptionObject caught!" <<std::endl;
			std::cerr<<err<<std::endl;
		}
	}
	else
		std::cout<<"Wrong input and output file extension combination"<<std::endl;
	
	return 0;
}
