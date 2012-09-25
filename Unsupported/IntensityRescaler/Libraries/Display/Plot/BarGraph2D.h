/*=========================================================================

  Program:   FiberViewer
  Module:    $RCSfile: BarGraph2D.h,v $
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
#ifndef __BarGraph2D_h_  
#define __BarGraph2D_h_  
  
//#include "IRISTypes.h"  
#include <vector>  
#include "Vector2D.h"  
#include "Vector3D.h"
#include <qgl.h>  
  
/**  
 * Settings for plotting 2D functions.  Based on Mathematica options  
 * to the Plot command (a small subset, of course).  
 */  
class BarGraph2DSettings  
{  
  public:  
  
  Vector2D<float> GetPlotRangeMin();  
  void SetPlotRangeMin(Vector2D<float>);  
    
  Vector2D<float> GetPlotRangeMax();  
  void SetPlotRangeMax(Vector2D<float>);  
  
  Vector2D<float> GetAxesPosition();  
  void SetAxesPosition(Vector2D<float>);  
    
  Vector2D<float> GetMajorTickSpacing();  
  void SetMajorTickSpacing(Vector2D<float>);  
  
  Vector2D<float> GetMinorTickSpacing();  
  void SetMinorTickSpacing(Vector2D<float>);  
  
  bool GetIntermediateHAxes();  
  bool GetIntermediateVAxes();  
  
  bool GetShowAxes();  
  void SetShowAxes(bool);  
  
  bool GetShowFrame();  
  void SetShowFrame(bool);  
  
  bool GetShowMajorTicks();  
  void SetShowMajorTicks(bool);  
  
  bool GetShowMinorTicks();  
  void SetShowMinorTicks(bool);  
  
  Vector3D<float> GetAxesColor();  
  void SetAxesColor(Vector3D<float>);  
  
  Vector3D<float> GetFrameColor();  
  void SetFrameColor(Vector3D<float>);  
  
  Vector3D<float> GetPlotColor();  
  void SetPlotColor(Vector3D<float>);  
  
  Vector3D<float> GetBackgroundColor();  
  void SetBackgroundColor(Vector3D<float>);  
  
  Vector3D<float> GetTextColor();  
  void SetTextColor(Vector3D<float>);  
  
  Vector3D<float> GetClearColor();  
  void SetClearColor(Vector3D<float>);  
  
  float GetTickXNumber();  
  void SetTickXNumber(int);  
  
  float GetTickYNumber();  
  void SetTickYNumber(int);  
  
  int GetPlotType();  
  void SetPlotType(int);  
  
  void SetXText(std::string);  
  void SetYText(std::string);  
  
  std::string GetXText();  
  std::string GetYText();  
  
  void SetScaleText(std::string value,int offset);  
  std::string GetScaleText(int offset);  
  
  
  void SetFloatingPointX(int value);  
  void SetFloatingPointY(int value);  
  int GetFloatingPointX();  
  int GetFloatingPointY();  
  
  
  
  static BarGraph2DSettings GetDefaultSettings();  
  
private:  
  Vector2D<float> m_PlotRangeMin;  
  Vector2D<float> m_PlotRangeMax;  
  Vector2D<float> m_AxesPosition;  
  Vector2D<float> m_MajorTickSpacing;  
  Vector2D<float> m_MinorTickSpacing;  
  
  bool m_ShowAxes;  
  bool m_ShowFrame;  
  bool m_ShowMajorTicks;  
  bool m_ShowMinorTicks;  
  bool m_IntermediateHAxes;  
  bool m_IntermediateVAxes;  
    
  Vector3D<float> m_AxesColor;  
  Vector3D<float> m_FrameColor;  
  Vector3D<float> m_PlotColor;  
  Vector3D<float> m_BackgroundColor;  
  Vector3D<float> m_TextColor;  
  Vector3D<float> m_ClearColor;  
  
  int m_TickXnumber;  
  int m_TickYnumber;  
  int m_PlotType;  
  std::string m_Xtext;  
  std::string m_Ytext;  
  std::string m_ScaleText[10];  
  int m_floatingpointx;  
  int m_floatingpointy;  
  
};  
  
/**  
 * \class BarGraph2D  
 * \brief A UI component for plotting 2D graphs using GL functions.  
 */  
class BarGraph2D  
{  
public:  
  /** Get a reference to the settings in this object */  
  BarGraph2DSettings &GetSettings()   
  {  
    return m_Settings;  
  }  
  
  /** Set the data to be plotted */  
  void SetDataPoints(float *xPoints,float *yPoints,unsigned int nPoints,int plot=0);
  void SetDataPoints1(float *xPoints,float *yPoints,float *y1Points,float *y2Points,unsigned int nPoints,int plot=2);    //*******
  void AddMean(float xpoint);  
  void SetXValue(float value);
  
  /**   
   * Draw the function in the standard OpenGL context.  The method will place  
   * the plot inside the region ((0,0),(1,1)).    
   */  
  void Draw(int width,int height);  
  
  /** Constructor, takes on default plot settings */  
  BarGraph2D();  
  
  void OutputString ( float x, float y, char *string );  
  void DrawCircle ( float x, float y, float radius);  
  void DrawRectangle ( float x, float y, float size);  
  void DrawCroix ( float x, float y, float size);  
  void DrawColormap ( float x, float y, float width,float height,int type=0);  
  void DrawTriangle ( float x, float y, float size);  
  
private:  
  typedef Vector2D<float> vector2D;  
  typedef std::vector<vector2D> PointVector;  
  GLuint fontOffset;  
  
  float m_XValue;
  BarGraph2DSettings m_Settings;    
  PointVector m_Points[10]; 
  int nb_plot;  
  std::vector<float> m_vplot;  
};  
  
#endif  
