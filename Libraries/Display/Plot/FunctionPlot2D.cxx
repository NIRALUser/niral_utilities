/*=========================================================================

  Program:   FiberViewer
  Module:    $RCSfile: FunctionPlot2D.cxx,v $
  Language:  C++
  Date:      $Date: 2007/09/20 16:28:22 $
  Version:   $Revision: 1.5 $
  Author:    Matthieu Jomier

  Copyright (c) 2004 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "FunctionPlot2D.h"  
  
#include <qgl.h>   
#include <GL/glut.h>  
    
FunctionPlot2D  
::FunctionPlot2D()  
{  
  m_Settings = FunctionPlot2DSettings::GetDefaultSettings();  
}  
  
void  
FunctionPlot2D  
::SetDataPoints(float *xPoints,float *yPoints,unsigned int nPoints,int plot)  
{  
  // Intialize the array  
  m_Points[plot].clear();  
  m_Points[plot].reserve(nPoints);  
    
  
  float xmax=-999999999.0f;  
  float xmin=9999999999.0f;  
  float ymin=999999999.0f;  
  float ymax=-99999999.0f;  
  nb_plot = plot+1;  
  
  
  // Fill the array  
  for(unsigned int i=0;i<nPoints;i++)  
  {  
    m_Points[plot].push_back(Vector2D<float>(xPoints[i],yPoints[i]));      
  if (xPoints[i] > xmax) xmax = xPoints[i];  
  if (xPoints[i] < xmin) xmin = xPoints[i];  
  
  if (yPoints[i] > ymax) ymax = yPoints[i];  
  if (yPoints[i] < ymin) ymin = yPoints[i];  
  }  
  
  //Set Min/Max values  
  //m_Settings.SetPlotRangeMin(Vector2D<float>(xmin,ymin));  
  //m_Settings.SetPlotRangeMax(Vector2D<float>(xmax,ymax));  
}  
  
void  
FunctionPlot2D  
::Draw()  
{  
  // Push the color/blending attributes  
  glPushAttrib(GL_COLOR_BUFFER_BIT);  
  
  // Draw the background rectangle  
  glColor3f(m_Settings.GetBackgroundColor()[0],m_Settings.GetBackgroundColor()[1],m_Settings.GetBackgroundColor()[2]);  
  glBegin(GL_QUADS);  
  glVertex2f(0,0);  
  glVertex2f(0,1);  
  glVertex2f(1,1);  
  glVertex2f(1,0);  
  glEnd();  
  
  // Draw the frame around the plot if requested  
  if(m_Settings.GetShowFrame())  
    {  
  glColor3f(m_Settings.GetFrameColor()[0],m_Settings.GetFrameColor()[1],m_Settings.GetFrameColor()[2]);  
    glBegin(GL_LINE_LOOP);  
    glVertex2f(0,0);  
    glVertex2f(0,1);  
    glVertex2f(1,1);  
    glVertex2f(1,0);  
    glEnd();  
    }  
  
  // These are the extents of the plot region  
  Vector2D<float> xMin = m_Settings.GetPlotRangeMin();  
  Vector2D<float> xMax = m_Settings.GetPlotRangeMax();  
  
  // Make sure these points create a rectangle  
//  assert(xMin(0) < xMax(0) && xMin(1) < xMax(1));  
    
  // Change the view matrix to the passed in coordinates  
  glMatrixMode(GL_MODELVIEW);  
  glPushMatrix();  
  glScalef(1.0f / (xMax[0] - xMin[0]), 1.0f / (xMax[1] - xMin[1]), 1.0f);  
  glTranslatef(-xMin[0],-xMin[1],0.0f);  
  
  // Do we draw the axes?  
  if(m_Settings.GetShowAxes())  
    {  
  
    // First draw the ticks, or it may look ugly  
    if(m_Settings.GetShowMajorTicks())  
      {  
        
      }  
  
    if(m_Settings.GetShowMinorTicks())  
      {  
  
      }  
  
    // Now draw the axes  
    Vector2D<float> xAxes = m_Settings.GetAxesPosition();  
  glColor3f(m_Settings.GetAxesColor()[0],m_Settings.GetAxesColor()[1],m_Settings.GetAxesColor()[2]);  
  
  
    glBegin(GL_LINES);  
    glVertex2f(xAxes[0],xMin[1]);  
    glVertex2f(xAxes[0],xMax[1]);  
    glVertex2f(xMin[0], xAxes[1]);  
    glVertex2f(xMax[0], xAxes[1]);  
    glEnd();  
  
    }  
  
  
    if(m_Settings.GetIntermediateHAxes())  
    {  
         // Now draw the axes  
    Vector2D<float> xAxes = m_Settings.GetAxesPosition();  
    glColor3f(m_Settings.GetAxesColor()[0],m_Settings.GetAxesColor()[1],m_Settings.GetAxesColor()[2]);  
    glEnable(GL_LINE_STIPPLE);  
    glLineStipple(2,0xAAAA);  
    float pas = (xMax[1]-xMin[1])/4;  
    for (int nb_line=0;nb_line<3;nb_line++)  
    {  
      glBegin(GL_LINES);  
      glVertex2f(xMin[0], xMin[1]+(nb_line+1)*pas);  
      glVertex2f(xMax[0], xMin[1]+(nb_line+1)*pas);  
      glEnd();  
    }  
    glDisable(GL_LINE_STIPPLE);  
  
  }  
  
  
    if(m_Settings.GetIntermediateVAxes())  
    {  
         // Now draw the axes  
    Vector2D<float> xAxes = m_Settings.GetAxesPosition();  
    glColor3f(m_Settings.GetAxesColor()[0],m_Settings.GetAxesColor()[1],m_Settings.GetAxesColor()[2]);  
    glEnable(GL_LINE_STIPPLE);  
    glLineStipple(2,0xAAAA);  
    float pas = (xMax[0]-xMin[0])/4;  
    for (int nb_line=0;nb_line<3;nb_line++)  
    {  
      glBegin(GL_LINES);  
      glVertex2f(xMin[0]+(nb_line+1)*pas,xMin[1]);  
      glVertex2f(xMin[0]+(nb_line+1)*pas,xMax[1]);  
      glEnd();  
    }  
    glDisable(GL_LINE_STIPPLE);  
  }  
  
  
    //Display Mean  
        // Now draw the axes  
    Vector2D<float> xAxes = m_Settings.GetAxesPosition();  
    glColor3f(0.0,0.0,1.0);  
  
    for (std::vector<float>::iterator j = m_vplot.begin(); j != m_vplot.end(); j++)  
    {  
      glBegin(GL_LINES);  
      glVertex2f(*j,xMin[1]);  
      glVertex2f(*j,xMax[1]);  
      glEnd();  
    }  
  
  
  
  glutWireTeapot(0.5);  
    glEnable(GL_BLEND);  
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
   glColor4f(1.0, 0, 0, (float)0.75);  
  // glRasterPos2i(0,height()/2);  
  
   glDisable(GL_BLEND);  
    
  
  // Allow line smoothing  
  glEnable(GL_BLEND);  
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
  glLineWidth(1.0);  
  
  // Finally, plot the function  
  glColor3f(m_Settings.GetPlotColor()[0],m_Settings.GetPlotColor()[1],m_Settings.GetPlotColor()[2]);  
  
  for (int nbpoint=0;nbpoint<nb_plot;nbpoint++)  
  {  
  glBegin(GL_LINE_STRIP);  
  // Plot all the points  
  
  for(PointVector::iterator it=m_Points[nbpoint].begin();it!=m_Points[nbpoint].end();++it)  
    {  
        glVertex2f((*it)[0],(*it)[1]);  
    }  
  
  glEnd();  
  }  
  
  // Restore the model matrix  
  glPopMatrix();    
  
  // Restore the attributes  
  glPopAttrib();  
}  
  
  
void  
FunctionPlot2D  
:: OutputString ( float x, float y, char *string )  
{  
  int length;  
  length = strlen ( string );  
  glRasterPos2f ( x, y );  
  for (int i = 0; i < length; i++ )  
  glutBitmapCharacter ( GLUT_BITMAP_9_BY_15, string[i] );  
}   
  
void  
FunctionPlot2D  
::AddMean(float xpoint)  
{  
  m_vplot.push_back(xpoint);  
}  
  
FunctionPlot2DSettings  
FunctionPlot2DSettings  
::GetDefaultSettings()  
{  
  FunctionPlot2DSettings settings;  
  
  settings.m_PlotRangeMin = Vector2D<float>(-1.0f,-1.0f);  
  settings.m_PlotRangeMax = Vector2D<float>(1.0f,1.0f);  
  settings.m_AxesPosition = Vector2D<float>(0.0f,0.0f);  
  settings.m_MajorTickSpacing = Vector2D<float>(0.5f,0.5f);  
  settings.m_MinorTickSpacing = Vector2D<float>(0.1f,0.1f);  
  
  settings.m_ShowAxes = true;  
  settings.m_ShowFrame = true;  
  settings.m_ShowMajorTicks = true;  
  settings.m_ShowMinorTicks = true;  
  settings.m_IntermediateHAxes = true;  
  settings.m_IntermediateVAxes = true;  
     
  settings.m_AxesColor = Vector3D<float>(0.0f,0.0f,0.0f);  
  settings.m_FrameColor = Vector3D<float>(0.0f,0.0f,0.0f);  
  settings.m_PlotColor = Vector3D<float>(1.0f,0.0f,0.0f);  
  settings.m_BackgroundColor = Vector3D<float>(1.0f,1.0f,1.0f);  
  
  return settings;  
}  
  
  
Vector2D<float> FunctionPlot2DSettings::GetPlotRangeMin()  
{  
  return m_PlotRangeMin;  
}  
  
void FunctionPlot2DSettings::SetPlotRangeMin(Vector2D<float> val)  
{  
  m_PlotRangeMin = val;  
}  
    
Vector2D<float> FunctionPlot2DSettings::GetPlotRangeMax()  
{  
  return m_PlotRangeMax;  
}  
  
void FunctionPlot2DSettings::SetPlotRangeMax(Vector2D<float> val)  
{  
  m_PlotRangeMax = val;  
}  
    
Vector2D<float> FunctionPlot2DSettings::GetAxesPosition()  
{  
  return m_AxesPosition;  
}  
  
void FunctionPlot2DSettings::SetAxesPosition(Vector2D<float> val)  
{  
  m_AxesPosition = val;  
}  
    
Vector2D<float> FunctionPlot2DSettings::GetMajorTickSpacing()  
{  
  return m_MajorTickSpacing;  
}  
  
void FunctionPlot2DSettings::SetMajorTickSpacing(Vector2D<float> val)  
{  
  m_MajorTickSpacing = val;  
}  
    
Vector2D<float> FunctionPlot2DSettings::GetMinorTickSpacing()  
{  
  return m_MinorTickSpacing;  
}  
  
void FunctionPlot2DSettings::SetMinorTickSpacing(Vector2D<float> val)  
{  
  m_MinorTickSpacing = val;  
}  
    
bool FunctionPlot2DSettings::GetShowAxes()  
{  
  return m_ShowAxes;  
}  
  
void FunctionPlot2DSettings::SetShowAxes(bool val)  
{  
  m_ShowAxes = val;  
}  
    
bool FunctionPlot2DSettings::GetShowFrame()  
{  
  return m_ShowFrame;  
}  
  
void FunctionPlot2DSettings::SetShowFrame(bool val)  
{  
  m_ShowFrame = val;  
}  
    
bool FunctionPlot2DSettings::GetShowMajorTicks()  
{  
  return m_ShowMajorTicks;  
}  
  
void FunctionPlot2DSettings::SetShowMajorTicks(bool val)  
{  
  m_ShowMajorTicks = val;  
}  
    
bool FunctionPlot2DSettings::GetShowMinorTicks()  
{  
  return m_ShowMinorTicks;  
}  
  
void FunctionPlot2DSettings::SetShowMinorTicks(bool val)  
{  
  m_ShowMinorTicks = val;  
}  
    
Vector3D<float> FunctionPlot2DSettings::GetAxesColor()  
{  
  return m_AxesColor;  
}  
  
void FunctionPlot2DSettings::SetAxesColor(Vector3D<float> val)  
{  
  m_AxesColor = val;  
}  
    
Vector3D<float> FunctionPlot2DSettings::GetFrameColor()  
{  
  return m_FrameColor;  
}  
  
void FunctionPlot2DSettings::SetFrameColor(Vector3D<float> val)  
{  
  m_FrameColor = val;  
}  
    
Vector3D<float> FunctionPlot2DSettings::GetPlotColor()  
{  
  return m_PlotColor;  
}  
  
void FunctionPlot2DSettings::SetPlotColor(Vector3D<float> val)  
{  
  m_PlotColor = val;  
}  
    
Vector3D<float> FunctionPlot2DSettings::GetBackgroundColor()  
{  
  return m_BackgroundColor;  
}  
  
void FunctionPlot2DSettings::SetBackgroundColor(Vector3D<float> val)  
{  
  m_BackgroundColor = val;  
}  
    
bool FunctionPlot2DSettings::GetIntermediateHAxes()  
{  
  return m_IntermediateHAxes;  
}  
  
bool FunctionPlot2DSettings::GetIntermediateVAxes()  
{  
  return m_IntermediateVAxes;  
}  
