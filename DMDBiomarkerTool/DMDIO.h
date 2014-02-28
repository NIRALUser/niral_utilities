/*=========================================================================
 *
 * Program :   $Insight Segmentation & Registration Toolkit $
 * Module  :   $DMDBiomarkerTool: DMDBase.h $
 * Purpose :   $The base class of Duchenne Muscle Dystrophy biomarker tools $
 * Language:   $C++ $
 * Date    :   $Date: 2010-06-04 12:36:34 $
 * Version :   $Revision: 0.10 $
 * Authors :   $Jiahui Wang, Martin Styner $
 * Update  :   1. Created the base class (06-04-10)
 * Copyright (c) Nero Image Research and Analysis Lab.@UNC All rights reserved.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even 
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.
 *
=========================================================================*/

#include "DMDBase.h"

class DMDIO : public DMDBase
{
    public:
        // declare the pixel type and image dimension
        typedef itk::OrientedImage< PixelType, 3 >       ImageType;
        DMDData( int d ) : DMDBase(d) {
	  //  typedef float                                           PixelType;
	  //  const unsigned int                                      Dimension = 3;
	  //  typedef itk::OrientedImage< PixelType, 3 >              ImageType;
	}
    
};
