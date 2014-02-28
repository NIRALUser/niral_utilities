/*=========================================================================
 *
 * Program :   $Insight Segmentation & Registration Toolkit $
 * Module  :   $DMDBiomarkerTool: DMDBase.h $
 * Purpose :   $The base class of Duchenne Muscle Dystrophy biomarker tools $
 * Language:   $C++ $
 * Date    :   $Date: 2010-06-09 19:07:21 $
 * Version :   $Revision: 0.10 $
 * Authors :   $Jiahui Wang, Martin Styner $
 * Update  :   1. Created the t2 fitting class (06-07-10)
 * Copyright (c) Neuro Image Research and Analysis Lab.@UNC All rights reserved.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even 
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.
 *
=========================================================================*/

#ifndef DMDT2FITTING
#define DMDT2FITTING    // prevent for redefine
#include "DMDData.h"

class DMDT2Fitting : public virtual DMDData
{
    private:       
    public:   
        DMDT2Fitting() {                 
                  
	}	
	typedef itk::RecursiveGaussianImageFilter< DMDData::OrientedImageType, DMDData::OrientedImageType >  FilterType;	
	void smooth( DMDData::OrientedSeriesReaderType::Pointer reader, DMDData::OrientedWriterType::Pointer & writer );
 
}; 
//////////////////////////////////////////////////////////////////////////
void DMDT2Fitting::smooth( DMDData::OrientedSeriesReaderType::Pointer reader, DMDData::OrientedWriterType::Pointer & writer )
{
    // specify the filter type -- gaussian filter
    typedef itk::RecursiveGaussianImageFilter< DMDData::OrientedImageType, DMDData::OrientedImageType >  FilterType;
    FilterType::Pointer                 filterX = FilterType::New();
    FilterType::Pointer                 filterY = FilterType::New();
    FilterType::Pointer                 filterZ = FilterType::New();
    
    filterX->SetDirection( 0 );   // 0 --> X direction
    filterY->SetDirection( 1 );   // 1 --> Y direction
    filterZ->SetDirection( 2 );   // 1 --> Z direction

    filterX->SetOrder( FilterType::ZeroOrder );
    filterY->SetOrder( FilterType::ZeroOrder );
    filterZ->SetOrder( FilterType::ZeroOrder );

    filterX->SetNormalizeAcrossScale( false );
    filterY->SetNormalizeAcrossScale( false );
    filterZ->SetNormalizeAcrossScale( false );

    filterX->SetInput(reader->GetOutput());
    filterY->SetInput(filterX->GetOutput());
    filterZ->SetInput(filterY->GetOutput());

    // Set the Filter width
    const double                        sigma = 3;
    filterX->SetSigma( sigma );
    filterY->SetSigma( sigma );
    filterZ->SetSigma( sigma );

    // run the pipeline
    filterZ->Update();

    //   return filterZ;
    writer->SetInput( filterZ->GetOutput() );

}

#endif
