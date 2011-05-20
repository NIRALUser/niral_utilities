#ifndef Vector2D_txx
#define Vector2D_txx

#include <iostream>
#include <cmath>
//#include <minmax.h>

template <class T>
inline
Vector2D<T>
::Vector2D()
  : x(), y()
{}

template <class T>
inline
Vector2D<T>
::Vector2D(const T& xIn, 
      const T& yIn)
  : x(xIn), y(yIn)
{}

#ifdef WIN32
template <class T, class U>
#else
template <class T>
template <class U>
#endif
inline
Vector2D<T>
::Vector2D(const Vector2D<U>& rhs)
  : x(rhs.x), y(rhs.y)
{}

#ifdef WIN32
template <class T, class U>
#else
template <class T>
template <class U>
#endif
inline
Vector2D<T>&
Vector2D<T>
::operator=(const Vector2D<U>& rhs)
{
  x = rhs.x;
  y = rhs.y;
  return *this;
}

template <class T>
inline
T&
Vector2D<T>
::operator[](unsigned int index) 
{
  return (&x)[index];
}
  
template <class T>
inline
const T&
Vector2D<T>
::operator[](unsigned int index) const
{
  return (&x)[index];
}

template <class T>
inline
double
Vector2D<T>
::length() const
{
  return sqrt(x*x + y*y);
}

template <class T>
inline
double
Vector2D<T>
::lengthSquared() const
{
  return x*x + y*y;
}

template <class T>
inline
double
Vector2D<T>
::productOfElements() const
{
  return x * y;
}

template <class T>
inline
double
Vector2D<T>
::sumOfElements() const
{
  return x + y;
}

template <class T>
inline
T
Vector2D<T>
::maxElement() const
{
  return std::max(x, y);
}

template <class T>
inline
T
Vector2D<T>
::minElement() const
{
  return std::min(x, y);
}

template <class T>
inline
void
Vector2D<T>
::set(const T& x,
      const T& y)
{ 
  this->x = x;
  this->y = y;
}

template <class T>
inline
void
Vector2D<T>
::scale(const T& sx,
   const T& sy)
{
  x *= sx;
  y *= sy;
}

template <class T>
inline
void
Vector2D<T>
::translate(const T& tx,
       const T& ty)
{
  x += tx;
  y += ty;
}

template <class T>
inline
void
Vector2D<T>
::invert()
{
  x = -x;
  y = -y;
}  

template <class T>
inline
void
Vector2D<T>
::normalize()
{
  double len = 1.0 / length();
  x *= len;
  y *= len;
}

template <class T>
inline
double
Vector2D<T>
::distance(const Vector2D<T>& rhs) const
{
  return sqrt(distanceSquared(rhs));
}

template <class T>
inline
double
Vector2D<T>
::distanceSquared(const Vector2D<T>& rhs) const
{
  return 
    (x - rhs.x) * (x - rhs.x) +
    (y - rhs.y) * (y - rhs.y);
}

template <class T>
inline
double
Vector2D<T>
::dot(const Vector2D<T>& rhs) const
{
  return (x * rhs.x + y * rhs.y);
}

template <class T>
inline
double
Vector2D<T>
::cross(const Vector2D<T>& rhs) const
{
  return (x * rhs.y - y * rhs.x);
}

template <class T>
inline
std::ostream& 
Vector2D<T>
::writeASCII(std::ostream& output) const
{
  output << "(" << x << ", " << y << ")";
  return output;  
}

template <class T>
inline
std::istream& 
Vector2D<T>
::readASCII(std::istream& input)
{
  char paren, comma;
  input >> paren
   >> x >> comma
   >> y >> comma;
  return input;
} 

template <class T>
inline
std::ostream& 
Vector2D<T>
::writeBinary(std::ostream& output) const
{
  output.write((char*)(&x), 2 * sizeof(T));
  return output;
}
  
template <class T>
inline
std::istream& 
Vector2D<T>
::readBinary(std::istream& input)
{
  input.read((char*)(&x), 2 * sizeof(T));
  return input;
}

//
// comparison operators
//
template <class T>
inline
bool 
operator<(const Vector2D<T>& lhs, const Vector2D<T>& rhs)
{
  if (lhs.x < rhs.x) return true;
  if (lhs.x > rhs.x) return false;
  return (lhs.y < rhs.y);
}

template <class T>
inline
bool 
operator==(const Vector2D<T>& lhs, const Vector2D<T>& rhs)
{
  return ((lhs.x == rhs.x) && 
     (lhs.y == rhs.y));
}

template <class T>
inline
bool 
operator!=(const Vector2D<T>& lhs, const Vector2D<T>& rhs)
{
  return ((lhs.x != rhs.x) || 
     (lhs.y != rhs.y));
}

//
// invert operator (unary minus)
//
template <class T>
inline
Vector2D<T> 
operator-(const Vector2D<T>& rhs)
{
  return Vector2D<T>(-rhs.x, -rhs.y);
}

//
// elementwise operators with Vector2D
//
template <class T, class U>
inline
Vector2D<T>& 
operator+=(Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  return lhs;
}

template <class T, class U>
inline
Vector2D<T>& 
operator-=(Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  return lhs;
}

template <class T, class U>
inline
Vector2D<T>& 
operator*=(Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  lhs.x *= rhs.x;
  lhs.y *= rhs.y;
  return lhs;
}

template <class T, class U>
inline
Vector2D<T>& 
operator/=(Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  lhs.x /= rhs.x;
  lhs.y /= rhs.y;
  return lhs;
}

template <class T, class U>
inline
Vector2D<T> 
operator+(const Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  return Vector2D<T>(static_cast<T>(lhs.x + rhs.x), 
           static_cast<T>(lhs.y + rhs.y));
}

template <class T, class U>
inline
Vector2D<T> 
operator-(const Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  return Vector2D<T>(static_cast<T>(lhs.x - rhs.x), 
           static_cast<T>(lhs.y - rhs.y));
}

template <class T, class U>
inline
Vector2D<T> 
operator*(const Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  return Vector2D<T>(static_cast<T>(lhs.x * rhs.x), 
           static_cast<T>(lhs.y * rhs.y));
}

template <class T, class U>
inline
Vector2D<T> 
operator/(const Vector2D<T>& lhs, const Vector2D<U>& rhs)
{
  return Vector2D<T>(static_cast<T>(lhs.x / rhs.x), 
           static_cast<T>(lhs.y / rhs.y));
}

//
// element-wise operators with a double
//
template <class T>
inline
Vector2D<T>& 
operator+=(Vector2D<T>& v, const double& u)
{
  v.x += u;
  v.y += u;
  return v;
}

template <class T>
inline
Vector2D<T>& 
operator-=(Vector2D<T>& v, const double& u)
{
  v.x -= u;
  v.y -= u;
  return v;
}

template <class T>
inline
Vector2D<T>& 
operator*=(Vector2D<T>& v, const double& u)
{
  v.x *= u;
  v.y *= u;
  return v;
}

template <class T>
inline
Vector2D<T>& 
operator/=(Vector2D<T>& v, const double& u)
{
  v.x /= u;
  v.y /= u;
  return v;
}

template <class T>
inline
Vector2D<T> 
operator+(const Vector2D<T>& v, const double& u)
{
  return Vector2D<T>(static_cast<T>(v.x + u), 
           static_cast<T>(v.y + u)); 
}

template <class T>
inline
Vector2D<T> 
operator+(const double& u, const Vector2D<T>& v)
{
  return Vector2D<T>(static_cast<T>(v.x + u),
           static_cast<T>(v.y + u)); 
}

template <class T>
inline
Vector2D<T> 
operator-(const Vector2D<T>& v, const double& u)
{
  return Vector2D<T>(static_cast<T>(v.x - u), 
           static_cast<T>(v.y - u)); 
}

template <class T>
inline
Vector2D<T> 
operator-(const double& u, const Vector2D<T>& v)
{
  return Vector2D<T>(static_cast<T>(u - v.x), 
           static_cast<T>(u - v.y)); 
}

template <class T>
inline
Vector2D<T> 
operator*(const Vector2D<T>& v, const double& u)
{
  return Vector2D<T>(static_cast<T>(v.x * u), 
           static_cast<T>(v.y * u)); 
}

template <class T>
inline
Vector2D<T> 
operator*(const double& u, const Vector2D<T>& v)
{
  return Vector2D<T>(static_cast<T>(v.x * u), 
           static_cast<T>(v.y * u)); 
}

template <class T>
inline
Vector2D<T> 
operator/(const Vector2D<T>& v, const double& u)
{
  return Vector2D<T>(static_cast<T>(v.x / u), 
           static_cast<T>(v.y / u)); 
}
template <class T>
inline
Vector2D<T> 
operator/(const double& u, const Vector2D<T>& v)
{
  return Vector2D<T>(static_cast<T>(u / v.x), 
           static_cast<T>(u / v.y)); 
}

//
// input/output
//
template <class T>
std::ostream& 
operator<<(std::ostream& output, const Vector2D<T>& v)
{
  return v.writeASCII(output);
}

template <class T>
std::istream& 
operator>>(std::istream& input, const Vector2D<T>& v)
{
  return v.readASCII(input);
}

#endif
