#include "Rotation3D.h"
#define M_PI 3.14159265358
#include <math.h>

Rotation3D::Rotation3D()
{
	m_centerx = 0;
	m_centery = 0;
	m_centerz = 0;
}

Rotation3D::~Rotation3D()
{
}

void Rotation3D::UpdateMatrix()
{
  double cosrx = cos(_rx*(M_PI/180.0));
  double cosry = cos(_ry*(M_PI/180.0));
  double cosrz = cos(_rz*(M_PI/180.0));
  double sinrx = sin(_rx*(M_PI/180.0));
  double sinry = sin(_ry*(M_PI/180.0));
  double sinrz = sin(_rz*(M_PI/180.0));

  // Add other transformation parameters to transformation matrix
  _matrix[0][0] = cosry*cosrz;
  _matrix[0][1] = cosry*sinrz;
  _matrix[0][2] = -sinry;
  _matrix[0][3] = _tx;
  _matrix[1][0] = (sinrx*sinry*cosrz-cosrx*sinrz);
  _matrix[1][1] = (sinrx*sinry*sinrz+cosrx*cosrz);
  _matrix[1][2] = sinrx*cosry;
  _matrix[1][3] = _ty;
  _matrix[2][0] = (cosrx*sinry*cosrz+sinrx*sinrz);
  _matrix[2][1] = (cosrx*sinry*sinrz-sinrx*cosrz);
  _matrix[2][2] = cosrx*cosry;
  _matrix[2][3] = _tz;
  _matrix[3][3] = 1.0;
}

void Rotation3D::UpdateParameter()
{
  double tmp;

  _tx = _matrix[0][3];
  _ty = _matrix[1][3];
  _tz = _matrix[2][3];

  tmp = asin(-_matrix[0][2]);

  _rx = atan2(_matrix[1][2]/cos(tmp), _matrix[2][2]/cos(tmp)) * 180.0/M_PI;
  _ry = tmp * 180.0/M_PI;
  _rz = atan2(_matrix[0][1]/cos(tmp), _matrix[0][0]/cos(tmp)) * 180.0/M_PI;
}

double Rotation3D::Get(int i)
{
  switch (i) {
  case 0:
    return _tx;
    break;
  case 1:
    return _ty;
    break;
  case 2:
    return _tz;
    break;
  case 3:
    return _rx;
    break;
  case 4:
    return _ry;
    break;
  case 5:
    return _rz;
    break;
  default:
    return 0;
  }
}

void Rotation3D::Put(int i, double x)
{
  switch (i) {
  case 0:
    _tx = x;
    break;
  case 1:
    _ty = x;
    break;
  case 2:
    _tz = x;
    break;
  case 3:
    _rx = x;
    break;
  case 4:
    _ry = x;
    break;
  case 5:
    _rz = x;
    break;
  default:
	  return;
  }
  
  // Update transformation matrix
  this->UpdateMatrix();
}


void Rotation3D::SetCenter(double x, double y, double z)
{
	m_centerx = x;
	m_centery = y;
	m_centerz = z;
}

void Rotation3D::Rotate(double& x, double& y, double& z)
{
   double a, b, c;

   // Pre-multiply point with transformation matrix
   a = _matrix[0][0]*(x-m_centerx)+_matrix[0][1]*(y-m_centery)+_matrix[0][2]*(z-m_centerz) + m_centerx;
   b = _matrix[1][0]*(x-m_centerx)+_matrix[1][1]*(y-m_centery)+_matrix[1][2]*(z-m_centerz) + m_centery;
   c = _matrix[2][0]*(x-m_centerx)+_matrix[2][1]*(y-m_centery)+_matrix[2][2]*(z-m_centerz) + m_centerz;

   // Copy result back
   x = a;
   y = b;
   z = c;
}
