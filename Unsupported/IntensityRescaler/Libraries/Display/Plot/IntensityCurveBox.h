/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: IntensityCurveBox.h,v $
  Language:  C++
  Date:      $Date: 2004/07/28 13:05:57 $
  Version:   $Revision: 1.2 $
  Copyright (c) 2003 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __IntensityCurveBox_h_
#define __IntensityCurveBox_h_

//#include <FL/Fl_Gl_Window.H>

#include <qgl.h> 
#include "Vector3D.h"
//#include <IRISTypes.h>
#include "IntensityCurveInterface.h"
//#include "threedwindow.h"
//#include <FLTKCanvas.h>
//#include <InteractionMode.h>

/**
 * \class IntensityCurveBox
 * \brief An FLTK Box (Gl_Window) used to paint intensity mapping curves.
 */
class threedwindow;

class IntensityCurveBox : public QGLWidget 
{
	  Q_OBJECT
public:
  IntensityCurveBox( QWidget *parent=0, const char *name=0,threedwindow* threeparent=0); 

  /**
   * Handle displaying the curve
   */
  void paintGL();

  // Get/set the intensity curve
  IntensityCurveInterface * GetCurve();
  void SetCurve(IntensityCurveInterface *);
  void updatePos(float* x,float*y);
  bool UpdateControl(Vector3D<float> &p);
  // Get/set the parent object
/*  IntensityCurveUILogic *GetParent();
  void SetParent(IntensityCurveUILogic *);*/

  /**
   * The resolution of the curve displayed on the screen
   * TODO: Control over curve resolution
   */
  static unsigned int CURVE_RESOLUTION;
bool mousepressed;


signals:
  void MouseMove();

protected:
    virtual void   mousePressEvent( QMouseEvent *e );
    virtual void   mouseReleaseEvent(QMouseEvent *e);
	virtual void   mouseMoveEvent ( QMouseEvent * e );

private:

  /**
   * Check if a control point is close to another point (i.e. mouse position)
   */
  int GetControlPointInVincinity(float x, float y, int pixelRadius); 

  /**
   * The intensity mapping curve
   */
  IntensityCurveInterface *m_Curve;
    int m_MovingControlPoint;
  /** 
   * Parent object
   */
//  IntensityCurveUILogic *m_Parent;

  /**
   * Interaction handler for control point manipulation
   */
/*  class DefaultHandler : public InteractionMode {
  public:
    DefaultHandler(IntensityCurveBox *parent);

    int OnMousePress(const FLTKEvent &event);
    int OnMouseRelease(const FLTKEvent &event, const FLTKEvent &pressEvent);
    int OnMouseDrag(const FLTKEvent &event, const FLTKEvent &pressEvent);
    int OnMouseEnter(const FLTKEvent &event);
    int OnMouseLeave(const FLTKEvent &event);
    int OnMouseMotion(const FLTKEvent &event);

  private:
    // Pointer to the parent canvas
    IntensityCurveBox *m_Parent;

    // Set cursor depending on selected point
    void SetCursor(int iControlPoint);

    // Update a control point to reflect mouse motion to p
    bool UpdateControl(const Vector3f &p);

    // Control point being currently edited
    int m_MovingControlPoint;

  } m_DefaultHandler;

  // Allow access to private data
  friend class DefaultHandler;*/
  threedwindow	*m_Parent;
  bool m_iscurve;
};

#endif // __IntensityCurveBox_h_

