/*=========================================================================

  Program:   FiberViewer
  Module:    $RCSfile: FunctionPlot2DBox.h,v $
  Language:  C++
  Date:      $Date: 2007/09/20 16:28:22 $
  Version:   $Revision: 1.3 $
  Author:    Matthieu Jomier

  Copyright (c) 2004 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __FunctionPlot2DBox_h_  
#define __FunctionPlot2DBox_h_  
  
#include "FunctionPlot2D.h"  
#include <qgl.h>   
  
class FunctionPlot2DBox : public QGLWidget  
{  
public:  
  /** The standard FLTK window constructor */  
  FunctionPlot2DBox( QWidget *parent=0, const char *name=0);  
    
  /** Access the plotting object associated with this window */  
  FunctionPlot2D &GetPlotter()  
  {  
    return m_Plotter;  
  }  
  
  
  void SetPlotter(FunctionPlot2D & plotter)  
  {  
    m_Plotter = plotter;  
  }  
  
  /** The draw method */  
  void paintGL();  
  void initializeGL();  
  
private:  
  FunctionPlot2D m_Plotter;  
};  
  
#endif  
