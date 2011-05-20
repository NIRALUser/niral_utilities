/*=========================================================================

  Program:   FiberViewer
  Module:    $RCSfile: BarGraph2D.cxx,v $
  Language:  C++
  Date:      $Date: 2007/09/20 16:28:22 $
  Version:   $Revision: 1.9 $
  Author:    Matthieu Jomier

  Copyright (c) 2004 NeuroImaging Lab @ UNC. All rights reserved.
  See NeuroLibCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "BarGraph2D.h"  
  
#include <qgl.h>   
#include <GL/glut.h>  
    
BarGraph2D  
::BarGraph2D()  
{  
  m_Settings = BarGraph2DSettings::GetDefaultSettings();  
  nb_plot = 0;  
  m_XValue = 0;
}  
 
 
void  
BarGraph2D  
::SetDataPoints(float *xPoints,float *yPoints,unsigned int nPoints,int plot)  
{  
  int i;
  
  // Initialize the array  
  for (i=0;i<nb_plot;i++)
    m_Points[i].clear();  
  m_Points[plot].reserve(nPoints);  
    
  float xmax=-999999999.0f;  
  float xmin=9999999999.0f;  
  float ymin=999999999.0f;  
  float ymax=-99.0f;  
  
  if (plot != 0)  
    {  
    Vector2D<float> m_min = m_Settings.GetPlotRangeMin();  
    Vector2D<float> m_max = m_Settings.GetPlotRangeMax();  
    xmax = m_max[0];  
    ymax = m_max[1];  
    xmin = m_min[0];  
    ymin = m_min[1];  
    }  
  
  nb_plot = plot+1;  
  
  // Fill the array  
  for(unsigned int point=0;point<nPoints;point++)  
    {  
    m_Points[plot].push_back(Vector2D<float>(xPoints[point],yPoints[point]));      
    if (xPoints[point] > xmax) xmax = xPoints[point];  
    if (xPoints[point] < xmin) xmin = xPoints[point];  
  
    if (yPoints[point] > ymax) ymax = yPoints[point];  
    if (yPoints[point] < ymin) ymin = yPoints[point];  
    }  

  ymin = 0;  
  
  //Set Min/Max values  
  m_Settings.SetPlotRangeMin(Vector2D<float>(xmin,ymin));  
  m_Settings.SetPlotRangeMax(Vector2D<float>(xmax,ymax));  
} 


void  BarGraph2D::SetDataPoints1(float *xPoints,float *yPoints,float *y1Points, float *y2Points, unsigned int nPoints,int plot)  
{  
  int i;
  
  // Initialize the array  
  for (i=0; i<nb_plot; i++)
    m_Points[i].clear();  
  for (i=0; i<=plot; i++)
    m_Points[i].reserve(nPoints);
  
  float xmax=-999999999.0f;  
  float xmin=9999999999.0f;  
  float ymin=999999999.0f;  
  float ymax=-99.0f;  
  
  if (plot != 2)  
    {  
    Vector2D<float> m_min = m_Settings.GetPlotRangeMin();  
    Vector2D<float> m_max = m_Settings.GetPlotRangeMax();  
    xmax = m_max[0];  
    ymax = m_max[1];  
    xmin = m_min[0];  
    ymin = m_min[1];  
    }  
    
  nb_plot = plot+1;  
   
  // Fill the array  
  for(i=0;i<nPoints;i++)  
    {  
    m_Points[0].push_back(Vector2D<float>(xPoints[i],yPoints[i])); 
    m_Points[1].push_back(Vector2D<float>(xPoints[i],y1Points[i]));
    m_Points[2].push_back(Vector2D<float>(xPoints[i],y2Points[i]));
    
    if (xPoints[i] > xmax) 
      xmax = xPoints[i];  
    if (xPoints[i] < xmin) 
      xmin = xPoints[i];  
  
    if (y1Points[i] > ymax) 
      ymax = y1Points[i];  
    
    if (y2Points[i] < ymin) 
      ymin = y2Points[i];
    }  

  if (ymin >= 0)
    ymin = 0;  

  
  //Set Min/Max values  
  m_Settings.SetPlotRangeMin(Vector2D<float>(xmin,ymin));  
  m_Settings.SetPlotRangeMax(Vector2D<float>(xmax,ymax));  
} 

void  BarGraph2D::SetXValue(float value)
{
  m_XValue = value;
}
 
void  BarGraph2D::Draw(int width,int height)  
{  
  float offsetx = 40;  
  float offsety = 30;  
  float boundaryx = 10;  
  float boundaryy = 5;  
  float gboundaryy = 5;  
  m_Settings.SetTickXNumber(5);  
  m_Settings.SetTickYNumber(4);
  
  // Push the color/blending attributes  
  glPushAttrib(GL_COLOR_BUFFER_BIT);  
  
  // Draw the background rectangle  
  glColor3f(m_Settings.GetBackgroundColor()[0],m_Settings.GetBackgroundColor()[1],m_Settings.GetBackgroundColor()[2]);  
  glBegin(GL_QUADS);  
  glVertex2f(offsetx,offsety);  
  glVertex2f(offsetx,height-boundaryy);  
  glVertex2f(width-boundaryx,height-boundaryy);  
  glVertex2f(width-boundaryx,offsety);  
  glEnd();  
  
  // Draw the frame around the plot if requested  
  if(m_Settings.GetShowFrame())  
    {  
    glColor3f(m_Settings.GetFrameColor()[0],m_Settings.GetFrameColor()[1],m_Settings.GetFrameColor()[2]);  
    glBegin(GL_LINE_LOOP);  
    glVertex2f(offsetx,offsety);  
    glVertex2f(offsetx,height-boundaryy);  
    glVertex2f(width-boundaryx,height-boundaryy);  
    glVertex2f(width-boundaryx,offsety);  
    glEnd();  
    }  
    
  
  // These are the extents of the plot region  
  Vector2D<float> xMin = m_Settings.GetPlotRangeMin();  
  Vector2D<float> xMax = m_Settings.GetPlotRangeMax();  
  
  
  
  if(m_Settings.GetIntermediateHAxes())  
    {  
    // Now draw the axes  
    Vector2D<float> xAxes = m_Settings.GetAxesPosition();  
    glColor3f(m_Settings.GetAxesColor()[0],m_Settings.GetAxesColor()[1],m_Settings.GetAxesColor()[2]);  
    glEnable(GL_LINE_STIPPLE);  
    glLineStipple(2,0xAAAA);  
    float  pas = (height-offsety-boundaryy)/m_Settings.GetTickYNumber();  
    for (int nb_line=0;nb_line<m_Settings.GetTickYNumber()-1;nb_line++)  
      {  
      glBegin(GL_LINES);  
      glVertex2f(offsetx, offsety+(nb_line+1)*pas);  
      glVertex2f(width-boundaryx, offsety+(nb_line+1)*pas);  
      glEnd();  
      }  
    glDisable(GL_LINE_STIPPLE);  
  
    if (nb_plot != 0)  
      {  
      pas = (height-offsety-boundaryy-gboundaryy)/m_Settings.GetTickYNumber();  
      for (int nb_line=0;nb_line<=m_Settings.GetTickYNumber();nb_line++)  
        {  
        char text[20];  
        
        char param[6] = "%.";  
        sprintf(param,"%s%i",param,m_Settings.GetFloatingPointY());  
        strcat(param,"f");  
        sprintf(text,param,(xMax[1]-xMin[1])*nb_line/m_Settings.GetTickYNumber() + xMin[1]);  

        int PosX = offsetx-10-(strlen(text)*6);
        if (PosX < 0)
          PosX=0;
        OutputString(PosX, offsety+gboundaryy+(nb_line)*pas-5,text); 
        //OutputString(offsetx-10-(strlen(text)*6), offsety+gboundaryy+(nb_line)*pas-5,text);  
        }  
      }  
    }  
  
  
  if(m_Settings.GetIntermediateVAxes())  
    {  
    // Now draw the axes  
    Vector2D<float> xAxes = m_Settings.GetAxesPosition();  
    glColor3f(m_Settings.GetAxesColor()[0],m_Settings.GetAxesColor()[1],m_Settings.GetAxesColor()[2]);  
    glEnable(GL_LINE_STIPPLE);  
    glLineStipple(2,0xAAAA);  
    float pas = (width-offsetx-boundaryx)/m_Settings.GetTickXNumber();  
    for (int nb_line=0;nb_line<m_Settings.GetTickXNumber()-1;nb_line++)  
      {  
      glBegin(GL_LINES);  
      glVertex2f(offsetx+(nb_line+1)*pas,offsety);  
      glVertex2f(offsetx+(nb_line+1)*pas,height-boundaryy);  
      glEnd();  
      }  
    glDisable(GL_LINE_STIPPLE);  
    
    if (nb_plot == 0)  
      return;  
    
    if (m_Settings.GetTickXNumber() >  m_Points[0].size())  
      m_Settings.SetTickXNumber(m_Points[0].size());  
  
  
    if (m_Settings.GetPlotType() == 0)  
      {  
      float pasx = (width-offsetx-boundaryx)/(m_Points[0].size());  
      int offset = (int)((m_Points[0].size())/m_Settings.GetTickXNumber());  
      for (unsigned int nb_line=0;nb_line<m_Points[0].size();nb_line+=offset)  
        {  
        char text[20];  
        char param[6] = "%.";  
        sprintf(param,"%s%i",param,m_Settings.GetFloatingPointX());  
        strcat(param,"f");  
        sprintf(text,param,m_Points[0][nb_line][0]);  
        //if (offset > max)  
        {  
        OutputString(offsetx+nb_line*pasx + (pasx/2) - (strlen(text)*5.5)/2,offsety-13,text);  
            
        }  
        }  
      }  
  
    if (m_Settings.GetPlotType() == 1)  
      {  
      for (int nb_line=0;nb_line<=m_Settings.GetTickXNumber();nb_line++)  
        {  
        char text[20];  
        sprintf(text,"%.0f",((xMax[0]-xMin[0])*nb_line)/m_Settings.GetTickXNumber() + xMin[0]);  
        if (nb_line==m_Settings.GetTickXNumber())  
          OutputString(offsetx+(nb_line)*pas-(((strlen(text))*5.5)/2)-6,offsety-13,text);      
        else  
          OutputString(offsetx+(nb_line)*pas-(((strlen(text))*5.5)/2),offsety-13,text);  
        }  
      }  
      
  
  
    }  
  
  float maxlen = 0;  
  for (int m_nbplot=0;m_nbplot<nb_plot;m_nbplot++)  
    if (strlen(m_Settings.GetScaleText(m_nbplot).c_str()) > maxlen)  
      maxlen = strlen(m_Settings.GetScaleText(m_nbplot).c_str());  
  
  //Draw X Text  
  OutputString((width-(strlen(m_Settings.GetXText().c_str())*4))/2,4,(char*)m_Settings.GetXText().c_str());  
  
  //Draw scale text  
  //Plot 0  
/*    glBegin(GL_LINES);  
    glVertex2f(width-maxlen*6.5-boundaryx,height-boundaryy-gboundaryy-6+3); //5.5  
    glVertex2f(width-maxlen*6.5-boundaryx+15,height-boundaryy-gboundaryy-6+3);  
    glEnd();  
    OutputString(width-maxlen*6.5-boundaryx+17,height-boundaryy-gboundaryy-6,(char*)m_Settings.GetScaleText(0).c_str());    
     
   //Plot 1  
   if (nb_plot > 0)  
   {  
       glEnable(GL_LINE_STIPPLE);  
      glLineStipple(2,0xAAAA);  
      glBegin(GL_LINES);  
      glVertex2f(width-maxlen*6.5-boundaryx,height-boundaryy-gboundaryy-6-10);  
      glVertex2f(width-maxlen*6.5-boundaryx+15,height-boundaryy-gboundaryy-6-10);  
      glEnd();  
      glDisable(GL_LINE_STIPPLE);  
      OutputString(width-maxlen*6.5-boundaryx+17,height-boundaryy-gboundaryy-6-13,(char*)m_Settings.GetScaleText(1).c_str());    
   }  
*/  
  
  //Draw Y Text  
//   OutputVerticalString(5,(height-(strlen(m_Settings.GetYText().c_str())*4))/2,(char*)m_Settings.GetYText().c_str());  
  
  
  glEnable(GL_BLEND);  
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
  glColor4f(1.0, 0, 0, (float)0.75);  
  glDisable(GL_BLEND);  
    
  
  // Allow line smoothing  
  glEnable(GL_BLEND);  
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
  glLineWidth(1.0);  
  
  // Finally, plot the function  
  for (int nbpoint=0;nbpoint<nb_plot;nbpoint++)  
    {    
    // Plot all the points  
    float i=offsetx+1;  
    float pasx = (width-offsetx-boundaryx)/(m_Points[nbpoint].size());     
    float pasy = (height-offsety-boundaryy-(gboundaryy*2))/(xMax[1]-xMin[1]);  
        
    char text[20];  
    sprintf(text,"%.0f",xMax[0]);  
      
    float abspas = (float)((width-offsetx-boundaryx-(strlen(text)*13))/(strlen(text)*13));  
    if (abspas > (float)m_Points[nbpoint].size()) 
      abspas = 1;  
    else  
      abspas = ((float)m_Points[nbpoint].size())/(float)abspas;  
  
    float ordpas = (float)((height-offsety-boundaryy-(gboundaryy*2))/20);  
    if (ordpas > (float)30) 
      ordpas = 1;  
    else  
      ordpas = ((float)30)/(float)ordpas;  
  
    if (ordpas == 0) ordpas=30;  
    if (abspas == 0) abspas=(strlen(text)*13);  
  
  
    int absoff = 0;  
    float xold;  
    float yold;  
    bool first = true;  
    float m_distance = 0;  
  
    //DrawColormap(offsetx-17,offsety+gboundaryy,15,height-boundaryy-(gboundaryy*2)-offsety);  
  
  
    for(PointVector::iterator it=m_Points[nbpoint].begin();it!=m_Points[nbpoint].end();++it)  
      {  
      //std:: cout << "bar nbpoint : "<<nbpoint<<" it: "<<(*it)<<std ::endl;
      //std:: cout << " it[0]: "<<(*it)[0]<<std ::endl;
    
    
//    float m_xvalue = m_FaAnalysis.Get_XValue();

      //Draw Bar
      if (m_Settings.GetPlotType() == 0)  
        {  
        glBegin(GL_QUADS);  
        glColor3f(m_Settings.GetPlotColor()[0],m_Settings.GetPlotColor()[1],m_Settings.GetPlotColor()[2]);  
        glVertex2f(i,(float)(offsety+gboundaryy));  
        glVertex2f(i+pasx-2,(float)(offsety+gboundaryy));  
        glVertex2f(i+pasx-2,(float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy);  
        glVertex2f(i,(float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy);  
        glEnd();  
  
        glColor3f(m_Settings.GetFrameColor()[0],m_Settings.GetFrameColor()[1],m_Settings.GetFrameColor()[2]);  
        glBegin(GL_LINE_LOOP);  
        glVertex2f(i,(float)(offsety+gboundaryy));  
        glVertex2f(i+pasx-2,(float)(offsety+gboundaryy));  
        glVertex2f(i+pasx-2,(float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy);  
        glVertex2f(i,(float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy);  
        glEnd();  
 
        /*   if((*it)[0] == m_XValue)
      {
        glPointSize(6.0);
        glBegin(GL_POINTS);
          glColor3f(1.0,0.0,0.0);
          glVertex2f(i-1+pasx/2,(float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy);
        glEnd();
        glPointSize(1.0);
      }
 */
        }  
        
   
      //Draw Line  
      if (m_Settings.GetPlotType() == 1)  
        {  
        if (first)  
          {  
          xold = i+(pasx/2);  
          yold = (float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy;  
          first = false;  
          }  
        else  
          {  
          m_distance += sqrt((xold-(i+(pasx/2)))*(xold-(i+(pasx/2)))+(yold-((float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy))*(yold-((float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy)));  
          }  
                
        //********************************************
        if (nbpoint==0)
          glColor3f(m_Settings.GetPlotColor()[0],m_Settings.GetPlotColor()[1],m_Settings.GetPlotColor()[2]);
        if ((nbpoint==1)||(nbpoint==2))
          {
          glColor3f(1.0,0.0,0.0);
          glEnable(GL_LINE_STIPPLE);  
          glLineStipple(2,0xAAAA); 
          }
        //********************************************
          
        //glColor3f(m_Settings.GetPlotColor()[0],m_Settings.GetPlotColor()[1],m_Settings.GetPlotColor()[2]);  
   
        glBegin(GL_LINES);  
        glVertex2f(xold,yold);  
        glVertex2f(i+(pasx/2),(float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy);  
        glEnd();  
    
        if((*it)[0] == m_XValue)
          {
          glPointSize(6.0);
          glBegin(GL_POINTS);
          glColor3f(1.0,0.0,0.0);
          glVertex2f(i+(pasx/2),(float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy);
          glEnd();
          glPointSize(1.0);
          }
                
                
        if (nbpoint != 0)  
          glDisable(GL_LINE_STIPPLE);  
  
        xold =  i+(pasx/2);  
        yold = (float)(offsety+gboundaryy)+((*it)[1]-xMin[1])*pasy;  
        }  
  
      //Write abscisse  
      
      if (m_distance >= 30)  
        {  
        m_distance = 0;  
        // if (nbpoint ==0)  DrawRectangle(xold,yold,2);  
        //  if (nbpoint ==1)  DrawTriangle(xold,yold,3);  
        //DrawCroix(xold,yold,3);  
        //DrawCircle(xold,yold,1.5);  
        }  
  
  
      /* if ((absoff%abspas == 0) && (absoff+abspas<m_Points[nbpoint].size()))  
        {  
          sprintf(text,"%.0f",(*it)[0]);  
          OutputString(i+((pasx-(strlen(text)-1)*9+8)/2),offsety/2-7,text);  
        }  
        */  
  
      //glVertex2f((*it)[0],(*it)[1]);  
      i += pasx;  
      absoff++;  
      }  
    
    // sprintf(text,"%.0f",m_Points[nbpoint][m_Points[nbpoint].size()-1][0]);  
    // std::cout << "TEXT: " << text << std::endl;  
//  OutputString((width-boundaryx-((strlen(text)-1)*9+8)/2),offsety/2-7,text);  
  
    /*if (m_Points[nbpoint].size() != 0)  
    for (int l=0;l<=30;l++)  
    {  
      //Write ordonnee  
      if (l%ordpas == 0)  
      {  
        sprintf(text,"%.0f",(xMax[1]-xMin[1])*l/30 + xMin[1]);  
        OutputString(5,offsety+gboundaryy+l*ordpasy-3,text);  
      }  
    }*/  
    }  
  // Restore the model matrix  
  //glPopMatrix();    
  
  // Restore the attributes  
  // glPopAttrib();  
}  
  
  
  
void  
BarGraph2D  
:: OutputString ( float x, float y, char *string )  
{  
  glColor3f(m_Settings.GetTextColor()[0],m_Settings.GetTextColor()[1],m_Settings.GetTextColor()[2]);  
  int length;  
  length = strlen ( string );  
  glRasterPos2f ( x, y );  
  for (int i = 0; i < length; i++ )  
    glutBitmapCharacter ( GLUT_BITMAP_HELVETICA_10 , string[i] );// GLUT_BITMAP_TIMES_ROMAN_10  //GLUT_BITMAP_9_BY_15  
   
  //delete string;  
}   
  
void  
BarGraph2D  
:: DrawColormap ( float x, float y, float width,float height,int type)  
{  
  switch (type)  
    {  
    case 0:   
      for (int i=0;i<255;i++)  
        {  
        glBegin(GL_QUADS);  
        glColor3f(((float)i)/255,((float)i)/255,((float)i)/255);  
        glVertex2f(x,y+(i*height)/255);  
        glVertex2f(x+width,y+(i*height)/255);  
        glVertex2f(x+width,y+((i+1)*height)/255);  
        glVertex2f(x,y+((i+1)*height)/255);  
        glEnd();         
        }  
      break;  
  
    }  
}  
  
void  
BarGraph2D  
:: DrawCircle ( float x, float y, float radius)  
{  
  float vectorY1=y;  
  float vectorX1=x;  
  float vectorX;  
  float vectorY;  
  float angle;  
  glBegin(GL_TRIANGLES);    
  for(int i=0;i<=360;i++)  
    {  
    angle=(float)(((double)i)/57.29577957795135);    
    vectorX=x+(radius*(float)sin((double)angle));  
    vectorY=y+(radius*(float)cos((double)angle));      
    glVertex2d(x,y);  
    glVertex2d(vectorX1,vectorY1);  
    glVertex2d(vectorX,vectorY);  
    vectorY1=vectorY;  
    vectorX1=vectorX;    
    }  
  glEnd();  
}   
  
void  
BarGraph2D  
:: DrawRectangle ( float x, float y, float size)  
{  
  glBegin(GL_QUADS);  
  glVertex2f(x-size,y-size);  
  glVertex2f(x-size,y+size);  
  glVertex2f(x+size,y+size);  
  glVertex2f(x+size,y-size);  
  glEnd();  
}   
  
void  
BarGraph2D  
:: DrawCroix ( float x, float y, float size)  
{  
  glBegin(GL_LINES);  
  glVertex2f(x-size,y-size);  
  glVertex2f(x+size,y+size);  
  glEnd();  
  
  glBegin(GL_LINES);  
  glVertex2f(x+size,y-size);  
  glVertex2f(x-size,y+size);  
  glEnd();  
}   
  
void  
BarGraph2D  
:: DrawTriangle ( float x, float y, float size)  
{  
  glBegin(GL_TRIANGLES);    
  for(int i=0;i<=360;i++)  
    {  
    glVertex2d(x-size,y-size);  
    glVertex2d(x+size,y-size);  
    glVertex2d(x,y+size);  
    }  
  glEnd();  
}   
  
  
void  
BarGraph2D  
::AddMean(float xpoint)  
{  
  m_vplot.push_back(xpoint);  
}  
  
BarGraph2DSettings  
BarGraph2DSettings  
::GetDefaultSettings()  
{  
  BarGraph2DSettings settings;  
  
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
  settings.m_ClearColor = Vector3D<float>(0.75f,0.75f,0.75f);  
  settings.m_TextColor = Vector3D<float>(0.0f,0.0f,0.0f);  
  settings.m_floatingpointx = 0;  
  settings.m_floatingpointy = 0;  
  
  return settings;  
}  
  
  
Vector2D<float> BarGraph2DSettings::GetPlotRangeMin()  
{  
  return m_PlotRangeMin;  
}  
  
void BarGraph2DSettings::SetPlotRangeMin(Vector2D<float> val)  
{  
  m_PlotRangeMin = val;  
}  
    
Vector2D<float> BarGraph2DSettings::GetPlotRangeMax()  
{  
  return m_PlotRangeMax;  
}  
  
void BarGraph2DSettings::SetPlotRangeMax(Vector2D<float> val)  
{  
  m_PlotRangeMax = val;  
}  
    
Vector2D<float> BarGraph2DSettings::GetAxesPosition()  
{  
  return m_AxesPosition;  
}  
  
void BarGraph2DSettings::SetAxesPosition(Vector2D<float> val)  
{  
  m_AxesPosition = val;  
}  
    
Vector2D<float> BarGraph2DSettings::GetMajorTickSpacing()  
{  
  return m_MajorTickSpacing;  
}  
  
void BarGraph2DSettings::SetMajorTickSpacing(Vector2D<float> val)  
{  
  m_MajorTickSpacing = val;  
}  
    
Vector2D<float> BarGraph2DSettings::GetMinorTickSpacing()  
{  
  return m_MinorTickSpacing;  
}  
  
void BarGraph2DSettings::SetMinorTickSpacing(Vector2D<float> val)  
{  
  m_MinorTickSpacing = val;  
}  
    
bool BarGraph2DSettings::GetShowAxes()  
{  
  return m_ShowAxes;  
}  
  
void BarGraph2DSettings::SetShowAxes(bool val)  
{  
  m_ShowAxes = val;  
}  
    
bool BarGraph2DSettings::GetShowFrame()  
{  
  return m_ShowFrame;  
}  
  
void BarGraph2DSettings::SetShowFrame(bool val)  
{  
  m_ShowFrame = val;  
}  
    
bool BarGraph2DSettings::GetShowMajorTicks()  
{  
  return m_ShowMajorTicks;  
}  
  
void BarGraph2DSettings::SetShowMajorTicks(bool val)  
{  
  m_ShowMajorTicks = val;  
}  
    
bool BarGraph2DSettings::GetShowMinorTicks()  
{  
  return m_ShowMinorTicks;  
}  
  
void BarGraph2DSettings::SetShowMinorTicks(bool val)  
{  
  m_ShowMinorTicks = val;  
}  
    
Vector3D<float> BarGraph2DSettings::GetAxesColor()  
{  
  return m_AxesColor;  
}  
  
void BarGraph2DSettings::SetAxesColor(Vector3D<float> val)  
{  
  m_AxesColor = val;  
}  
    
Vector3D<float> BarGraph2DSettings::GetFrameColor()  
{  
  return m_FrameColor;  
}  
  
void BarGraph2DSettings::SetFrameColor(Vector3D<float> val)  
{  
  m_FrameColor = val;  
}  
    
Vector3D<float> BarGraph2DSettings::GetPlotColor()  
{  
  return m_PlotColor;  
}  
  
void BarGraph2DSettings::SetPlotColor(Vector3D<float> val)  
{  
  m_PlotColor = val;  
}  
    
Vector3D<float> BarGraph2DSettings::GetBackgroundColor()  
{  
  return m_BackgroundColor;  
}  
  
void BarGraph2DSettings::SetBackgroundColor(Vector3D<float> val)  
{  
  m_BackgroundColor = val;  
}  
  
Vector3D<float> BarGraph2DSettings::GetTextColor()  
{  
  return m_TextColor;  
}  
  
void BarGraph2DSettings::SetTextColor(Vector3D<float> val)  
{  
  m_TextColor = val;  
}  
  
  
Vector3D<float> BarGraph2DSettings::GetClearColor()  
{  
  return m_ClearColor;  
}  
  
void BarGraph2DSettings::SetClearColor(Vector3D<float> val)  
{  
  m_ClearColor = val;  
}    
  
  
bool BarGraph2DSettings::GetIntermediateHAxes()  
{  
  return m_IntermediateHAxes;  
}  
  
bool BarGraph2DSettings::GetIntermediateVAxes()  
{  
  return m_IntermediateVAxes;  
}  
  
  
int BarGraph2DSettings::GetPlotType()  
{  
  return m_PlotType;  
}  
  
void BarGraph2DSettings::SetPlotType(int val)  
{  
  m_PlotType = val;  
}  
  
  
  
float BarGraph2DSettings::GetTickXNumber()  
{  
  return m_TickXnumber;  
}  
  
  
void BarGraph2DSettings::SetTickXNumber(int value)  
{  
  m_TickXnumber = value;  
}  
  
float BarGraph2DSettings::GetTickYNumber()  
{  
  return m_TickYnumber;  
}  
  
void BarGraph2DSettings::SetTickYNumber(int value)  
{  
  m_TickYnumber = value;  
}  
  
  
void BarGraph2DSettings::SetXText(std::string value)  
{  
  m_Xtext = value;  
}  
  
void BarGraph2DSettings::SetYText(std::string value)  
{  
  m_Ytext = value;  
}  
  
  
std::string BarGraph2DSettings::GetXText()  
{  
  return m_Xtext;  
}  
  
std::string BarGraph2DSettings::GetYText()  
{  
  return m_Ytext;  
}  
  
void BarGraph2DSettings::SetScaleText(std::string value,int offset)  
{  
  m_ScaleText[offset] = value;  
}  
  
std::string BarGraph2DSettings::GetScaleText(int offset)  
{  
  return m_ScaleText[offset];  
}  
  
  
void BarGraph2DSettings::SetFloatingPointX(int value)  
{  
  m_floatingpointx = value;  
}  
  
void BarGraph2DSettings::SetFloatingPointY(int value)  
{  
  m_floatingpointy = value;  
}  
  
int BarGraph2DSettings::GetFloatingPointX()  
{  
  return m_floatingpointx;  
}  
  
int BarGraph2DSettings::GetFloatingPointY()  
{  
  return m_floatingpointy;  
}
