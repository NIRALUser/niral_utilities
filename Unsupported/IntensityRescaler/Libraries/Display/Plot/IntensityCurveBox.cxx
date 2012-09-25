/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: IntensityCurveBox.cxx,v $
  Language:  C++
  Date:      $Date: 2005/07/27 21:08:47 $
  Version:   $Revision: 1.3 $
  Copyright (c) 2003 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "IntensityCurveBox.h"
//#include "IntensityCurveUILogic.h"

#include <cassert>
#include <cstdio>
#include <iostream>

#include <qgl.h> 

#define CURVE_RESOLUTION 256

IntensityCurveBox
::IntensityCurveBox( QWidget *parent, const char *name, threedwindow* threeparent )
    : QGLWidget( parent, name )
{
      // Start with the blank curve
       m_Curve = NULL;
      //  m_Parent = parent;
      m_iscurve = false;
       // Set up the default handler
      //  PushInteractionMode(&m_DefaultHandler);
    installEventFilter(this);
}

void 
IntensityCurveBox
::paintGL() 
{       
 
      if (m_iscurve == false)
            return;

  //if (!valid()) 
    {
    // Set up the basic projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-0.05,1.05,-0.05,1.05);
    glViewport(0,0,width(),height());

    // Establish the model view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    double modelMatrix[16], projMatrix[16];
    int viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
#ifdef __APPLE__
    glGetIntegerv(GL_VIEWPORT, (GLint *) viewport);
#else
    glGetIntegerv(GL_VIEWPORT,viewport);
#endif
    }

  // Clear the viewport
  glClearColor(0.81,0.81,0.81,1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

  // Push the related attributes
  glPushAttrib(GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_LINE_BIT);

  // Disable lighting
  glDisable(GL_LIGHTING);

  // Draw the plot area
  // glColor3d(0.85,0.85,0.85);
  glBegin(GL_QUADS);
  glColor3d(0.,0.,0.);
  glVertex2d(0,0);
  // glColor3d(0.,0.,0.85);
  glColor3d(1.,1.,1.);
  glVertex2d(0,1);
  glColor3d(1.,1.,1.);
  glVertex2d(1,1);
  glColor3d(0.,0.,0.);
  glVertex2d(1,0);
  glEnd();

  // Draw the box around the plot area
  glColor3d(0,0,0);
  glBegin(GL_LINE_LOOP);
  glVertex2d(0,0);
  glVertex2d(0,1);
  glVertex2d(1,1);
  glVertex2d(1,0);
  glEnd();

  // Set up the smooth line drawing style
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glLineWidth(2.0);

  // Draw the curve using linear segments
  glColor3d(1.0,0.0,0.0);
  glBegin(GL_LINE_STRIP);

  float t = 0.0;
  float tStep = 1.0f / (CURVE_RESOLUTION);
  for (unsigned int i=0;i<=CURVE_RESOLUTION;i++) 
    {
    glVertex2f(t,m_Curve->Evaluate(t));
    t+=tStep;
    }

  glEnd();

  // Draw the handles
  for (unsigned int c=0;c<m_Curve->GetControlPointCount();c++) 
    {
    // Get the next control point
    float t,x;
    m_Curve->GetControlPoint(c,t,x);

    // Draw a quad around the control point

    double rx = 5.0 / width();
    double ry = 5.0 / height();

    glColor3d(1,1,0.5);
    glBegin(GL_QUADS);
    glVertex2d(t,x-ry);
    glVertex2d(t+rx,x);
    glVertex2d(t,x+ry);
    glVertex2d(t-rx,x);
    glEnd();       

    glColor3d(0,0,0);
    glLineWidth(1.0);
    glColor3d(1.0,0.0,0.0);
    glBegin(GL_LINE_LOOP);
    glVertex2d(t,x-ry);
    glVertex2d(t+rx,x);
    glVertex2d(t,x+ry);
    glVertex2d(t-rx,x);
    glEnd();

    }

  // Pop the attributes
  glPopAttrib();

  // Done
  glFlush();
}

int 
IntensityCurveBox
::GetControlPointInVincinity(float x, float y, int pixelRadius) 
{
  float rx = pixelRadius * 1.0f / width();
  float ry = pixelRadius * 1.0f / height();
  float fx = 1.0f / (rx * rx);
  float fy = 1.0f / (ry * ry);

  float minDistance = 1.0f;
  int nearestPoint = -1;

  for (unsigned int c=0;c<m_Curve->GetControlPointCount();c++) {
    // Get the next control point
    float cx,cy;
    m_Curve->GetControlPoint(c,cx,cy);

    std::cout << cx << " | " << cy << " and " << x << " | " << y << std::endl;

    // Check the distance to the control point
    float d = (cx - x) * (cx - x) * fx + (cy - y) * (cy - y) * fy;
    if (minDistance >= d) {
      minDistance = d;
      nearestPoint = c;
    }
  }

  // Negative: return -1
  return nearestPoint;
}

IntensityCurveInterface * IntensityCurveBox
::GetCurve()
{
  return m_Curve;
}

void IntensityCurveBox
::SetCurve(IntensityCurveInterface * val)
{
      m_Curve = val;
      m_iscurve = true;
}


void IntensityCurveBox::mousePressEvent( QMouseEvent *e )
{
  // Check the control point affected by the event
      QPoint pos = static_cast<QMouseEvent *>( e )->pos();
      float x,y;
      x = pos.x();
      y = pos.y();
      updatePos(&x,&y);
  m_MovingControlPoint = GetControlPointInVincinity(x,y,5);
  //SetCursor(m_MovingControlPoint);
  std::cout << m_MovingControlPoint << std::endl;
  mousepressed = true;
  //return 1;
}

void IntensityCurveBox::mouseReleaseEvent( QMouseEvent *e)
{
  mousepressed = false;
}

void IntensityCurveBox::mouseMoveEvent(QMouseEvent * e)
{
      QPoint pos = static_cast<QMouseEvent *>( e )->pos();
      float x,y;
      x = pos.x();
      y = pos.y();
      updatePos(&x,&y);
   if (m_MovingControlPoint >= 0) 
    {
    // Update the moving control point
        Vector3D<float> myvector(x,y,0);
     if (UpdateControl(myvector) )
      {
             repaint();
      }
    }      
   
   emit MouseMove();
//   m_Parent->UpdateIntensity();
}




void IntensityCurveBox::updatePos(float* x,float*y)
{
  int XCanvas[3];
  XCanvas[0] = (int)*x;
  XCanvas[1] = height() - 1 - (int)*y;

  // Make this window the current context
//  make_current();
makeCurrent(); 
  // Convert the event coordinates into the model view coordinates
  double modelMatrix[16], projMatrix[16];
  int viewport[4];
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
#ifdef __APPLE__
  glGetIntegerv(GL_VIEWPORT, (GLint *) viewport);
#else
  glGetIntegerv(GL_VIEWPORT,viewport);
#endif

  // Projection works with doubles, event is a float
  Vector3D<double> xProjection;
#ifdef __APPLE__
  gluUnProject(XCanvas[0],XCanvas[1],0,
                 modelMatrix,projMatrix, (GLint *) viewport,
                 &xProjection[0],&xProjection[1],&xProjection[2]);
#else
  gluUnProject(XCanvas[0],XCanvas[1],0,
                 modelMatrix,projMatrix,viewport,
                 &xProjection[0],&xProjection[1],&xProjection[2]);
#endif
  
  *x = xProjection[0];
  *y = xProjection[1];

}



/*IntensityCurveBox::DefaultHandler
::DefaultHandler(IntensityCurveBox *parent) 
{
  this->m_Parent = parent;
}*/

/*int 
IntensityCurveBox::DefaultHandler
::OnMousePress(const FLTKEvent &event)
{
  // Check the control point affected by the event
  m_MovingControlPoint = 
    m_Parent->GetControlPointInVincinity(event.XSpace[0],event.XSpace[1],5);

  SetCursor(m_MovingControlPoint);

  return 1;
}

int 
IntensityCurveBox::DefaultHandler
::OnMouseRelease(const FLTKEvent &event,
                 const FLTKEvent &irisNotUsed(dragEvent))
{
  if (m_MovingControlPoint >= 0) {
    // Update the control point
    if (UpdateControl(event.XSpace))
      {
      // Repaint parent
      m_Parent->redraw();

      // Fire the update event (should this be done on drag?)
      m_Parent->GetParent()->OnCurveChange();
      }

    // Set the cursor back to normal
    SetCursor(-1);
  }

  return 1;
}

int 
IntensityCurveBox::DefaultHandler
::OnMouseDrag(const FLTKEvent &event,
              const FLTKEvent &irisNotUsed(dragEvent))
{
  if (m_MovingControlPoint >= 0) 
    {
    // Update the moving control point
    if (UpdateControl(event.XSpace)) 
      {
      // Repaint parent
      m_Parent->redraw();

      // Fire the update event (should this be done on drag?)
      m_Parent->GetParent()->OnCurveChange();
      }
    }

  return 1;
}

int 
IntensityCurveBox::DefaultHandler
::OnMouseEnter(const FLTKEvent &irisNotUsed(event))
{
  return 1;
}

int 
IntensityCurveBox::DefaultHandler
::OnMouseLeave(const FLTKEvent &irisNotUsed(event))
{
  return 1;
}

int 
IntensityCurveBox::DefaultHandler
::OnMouseMotion(const FLTKEvent &event)
{
  int cp = 
    m_Parent->GetControlPointInVincinity(event.XSpace[0],event.XSpace[1],5);

  SetCursor(cp);

  return 1;
}

void 
IntensityCurveBox::DefaultHandler
::SetCursor(int cp) 
{
  if (cp < 0) 
    {
    fl_cursor(FL_CURSOR_DEFAULT);
    }   
  else if (cp == 0 || 
           cp == (int)(m_Parent->GetCurve()->GetControlPointCount()-1)) 
    {
    fl_cursor(FL_CURSOR_WE);
    } 
  else     
    {
    fl_cursor(FL_CURSOR_MOVE);
    }
}
*/

bool IntensityCurveBox::UpdateControl(Vector3D<float> &p) 
{
  // Take a pointer to the spline for convenience
  IntensityCurveInterface &curve = *m_Curve;
  int last = curve.GetControlPointCount()-1;

  // Check whether the motion is out of range
  //if (p[0] < 0.0 || p[0] > 1.0 || p[1] < 0.0 || p[1] > 1.0 )
   // return false;

      if (p[0] < 0.0) p[0] = 0;
      if (p[0] > 1.0) p[0] = 1.0;
      if (p[1] < 0.0) p[1] = 0.0;
      if (p[1] > 1.0) p[1] = 1.0;

  // First and last control points are treated specially because they
  // provide windowing style behavior
  if (m_MovingControlPoint == 0 || m_MovingControlPoint == last) 
    {
    // Get the current domain
    float xMin,xMax,yMin,yMax;        
    curve.GetControlPoint(0,xMin,yMin);
    curve.GetControlPoint(last,xMax,yMax);

    // Check if the new domain is valid
    float epsilon = 0.02;
    if (m_MovingControlPoint == 0 && p[0] < xMax - epsilon && yMin < 0+epsilon) 
      {
       xMin = p[0];
         curve.ScaleControlPointsToWindow(xMin,xMax);
         yMin = 0;
         curve.ScaleControlPointsYToWindow(yMin,yMax); 
      } 
    else if (m_MovingControlPoint == last && p[0] > xMin + epsilon && yMax > 1-epsilon)  
      {
        xMax = p[0];
          curve.ScaleControlPointsToWindow(xMin,xMax);
            yMax = 1;
          curve.ScaleControlPointsYToWindow(yMin,yMax); 
      }

    if (m_MovingControlPoint == 0 && p[1] < yMax - epsilon  && xMin < 0+ epsilon) 
      {
            yMin = p[1];
            xMin = 0;
          curve.ScaleControlPointsToWindow(xMin,xMax);
            curve.ScaleControlPointsYToWindow(yMin,yMax); 
      } 
      else if (m_MovingControlPoint == last && p[1] > yMin + epsilon  && xMax > 1- epsilon) 
    {
            yMax = p[1];
            xMax = 1;
          curve.ScaleControlPointsToWindow(xMin,xMax);
            curve.ScaleControlPointsYToWindow(yMin,yMax);
    } 


    } 
  else if (m_MovingControlPoint > 0) 
    {
    // Check whether the Y coordinate is in range
    if (p[1] < 0.0 || p[1] > 1.0)
      return false;

    // Record the position of the curve before motion
    float x,y;
    curve.GetControlPoint(m_MovingControlPoint,x,y);

    // Update the control point
    curve.UpdateControlPoint(m_MovingControlPoint,(float) p[0],(float) p[1]);

    // Check the curve for monotonicity
    if (!curve.IsMonotonic()) {
      curve.UpdateControlPoint(m_MovingControlPoint,x,y);
      return false;
    }
    }


  return true;
}

