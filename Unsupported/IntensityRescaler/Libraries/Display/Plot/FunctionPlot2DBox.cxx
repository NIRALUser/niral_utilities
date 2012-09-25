/*=========================================================================

  Program:   FiberViewer
  Module:    $RCSfile: FunctionPlot2DBox.cxx,v $
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
#include "FunctionPlot2DBox.h"  
  
FunctionPlot2DBox  
::FunctionPlot2DBox( QWidget *parent, const char *name )  
    : QGLWidget( parent, name )  
{  
}  
  
  
void   
FunctionPlot2DBox  
::paintGL()  
{  
/*  qglClearColor( Qt :: green);  
  makeCurrent();*/  
  
  // The standard 'valid' business  
  //if(!valid())  
    /*{  
    // Set up the basic projection with a small margin  
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity();  
    gluOrtho2D(-0.05,1.05,-0.05,1.05);  
    glViewport(0,0,width(),height());  
  
    // Establish the model view matrix  
    glMatrixMode(GL_MODELVIEW);  
    glLoadIdentity();      
    }*/  
    initializeGL();  
  // Clear the viewport  
  glClearColor(0.75,0.75,0.75,1.0);  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);   
  
  // Push the related attributes  
  glPushAttrib(GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_LINE_BIT);  
  
  // Disable lighting  
  glDisable(GL_LIGHTING);  
  
  // Call the plotting routine  
  m_Plotter.Draw();  
   
  // Pop the attributes  
  glPopAttrib();  
  
  // Done  
  glFlush();  
}  
  
void   
FunctionPlot2DBox  
::initializeGL()    
{  
    // Set up the basic projection with a small margin  
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity();  
    gluOrtho2D(-0.05,1.05,-0.05,1.05);  
    glViewport(0,0,width(),height());  
  
    // Establish the model view matrix  
    glMatrixMode(GL_MODELVIEW);  
    glLoadIdentity();      
}  
