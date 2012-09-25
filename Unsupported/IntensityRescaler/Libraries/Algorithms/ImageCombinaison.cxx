#include "ImageCombinaison.h"
  
ImageCombinaison::ImageCombinaison()
{
}


ImageCombinaison::~ImageCombinaison()
{
}


ImageCombinaison::ImagePointer ImageCombinaison::Multiply(ImagePointer image1, ImagePointer image2)
{
	//Duplicate first Image
	ImagePointer image = DuplicateImage(image1);

	IteratorType itS(image,image->GetLargestPossibleRegion());
	IteratorType itS1(image1,image1->GetLargestPossibleRegion());
	IteratorType itS2(image2,image2->GetLargestPossibleRegion());


	itS1.GoToBegin();
	itS2.GoToBegin();
	itS.GoToBegin();

	while (!itS1.IsAtEnd())
	{
		itS.Set(itS1.Get()*itS2.Get());
		++itS;
		++itS1;
		++itS2;
	}

	return image;
}


ImageCombinaison::ImagePointer ImageCombinaison::Minus(float val,ImagePointer image1)
{
	//Duplicate first Image
	ImagePointer image = DuplicateImage(image1);

	IteratorType itS(image,image->GetLargestPossibleRegion());
	IteratorType itS1(image1,image1->GetLargestPossibleRegion());

	itS1.GoToBegin();
	itS.GoToBegin();

	while (!itS1.IsAtEnd())
	{
		itS.Set(val-itS1.Get());
		++itS;
		++itS1;
	}

	return image;
}

ImageCombinaison::ImagePointer ImageCombinaison::Minus(ImagePointer image1,float val)
{
	//Duplicate first Image
	ImagePointer image = DuplicateImage(image1);

	IteratorType itS(image,image->GetLargestPossibleRegion());
	IteratorType itS1(image1,image1->GetLargestPossibleRegion());

	itS1.GoToBegin();
	itS.GoToBegin();

	while (!itS1.IsAtEnd())
	{
		itS.Set(itS1.Get()-val);
		++itS;
		++itS1;
	}

	return image;
}

ImageCombinaison::ImagePointer ImageCombinaison::Divide(ImagePointer image1,float div)
{
	//Duplicate first Image
	ImagePointer image = DuplicateImage(image1);
	IteratorType itS(image,image->GetLargestPossibleRegion());
	IteratorType itS1(image1,image1->GetLargestPossibleRegion());
	itS1.GoToBegin();
	itS.GoToBegin();

	while (!itS1.IsAtEnd())
	{
		itS.Set(itS1.Get()/div);
		++itS;
		++itS1;
	}

	return image;
}





ImageCombinaison::ImagePointer ImageCombinaison::DuplicateImage(ImagePointer image)
{
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(image);
    duplicator->Update();
	return duplicator->GetOutput();
}
