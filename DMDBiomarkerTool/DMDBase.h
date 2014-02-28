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

#include "itkOrientedImage.h"

class DMDBase {
    private:
  
    public:
        const unsigned int                                          Dimension;
	typedef float                                               PixelType;
        DMDBase( int d ) : Dimension(3)    //  use a member initialization list to initialize a data member 
	{                       
	  //    typedef float                                           PixelType;
	    //            typedef itk::OrientedImage< PixelType, Dimension>       ImageType;   
	}
};

      
