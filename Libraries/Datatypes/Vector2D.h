#ifndef Vector2D_h
#define Vector2D_h

#include <iosfwd>
#include <iostream>

//
// This struct is meant to be a fast 2D vector class.
// It is intended for use with build in types (int, float, etc.)
// 
// Optimization is important for this struct
// - almost everything is inlined
// - readability and maintainability tradeoffs have been made 
//   for speed gains
//
// bcd 2003
//

template <class T>
class Vector2D 
{
public:
  //
  // IMPORTANT: DO NOT change the order of x, y.  
  // DO NOT add other data members between x, y.
  // Array style access depends on these being adjacent
  // in memory.  This is guarenteed for standards compliant
  // compilers only if they are declared as they are here.
  // DO NOT TOUCH THIS LINE
  //
  T x, y;
  // IMPORTANT DO NOT TOUCH THE ABOVE LINE //

  // 
  // constructors
  //

  Vector2D();
  Vector2D(const T& xIn, const T& yIn);

  template <class U>
  Vector2D(const Vector2D<U>& rhs);

  template <class U>
  Vector2D<T>& operator=(const Vector2D<U>& rhs);

  //
  // array style coordinate access
  //
  // with optomization turned on, this is just as fast
  // as .x, .yaccess
  //
  T& operator[](unsigned int index);
  const T& operator[](unsigned int index) const;

  //
  // the euclidean length of this vector
  //
  double length() const;
  double lengthSquared() const;

  //
  // element accumulation
  //
  double productOfElements() const;
  double sumOfElements() const;

  //
  // min and max coordinate values
  //
  T maxElement() const;
  T minElement() const;

  //
  // manipulation methods
  // 
  void set(const T& x, const T& y);
  void scale(const T& sx, const T& sy);
  void translate(const T& tx, const T& ty);
  void invert();
  void normalize();

  //
  // the distance between the tip of this and the tip of rhs
  //
  double distance(const Vector2D<T>& rhs) const;
  double distanceSquared(const Vector2D<T>& rhs) const;

  //
  // inner (dot) product
  //
  double dot(const Vector2D<T>& rhs) const;

  //
  // vector (cross) product
  //
  double cross(const Vector2D<T>& rhs) const;

  //
  // output
  //
  std::ostream& writeASCII(std::ostream& output = std::cerr) const;
  std::istream& readASCII(std::istream& input = std::cin);
  std::ostream& writeBinary(std::ostream& output = std::cerr) const;
  std::istream& readBinary(std::istream& input = std::cin);
};

//
// comparison operators
//
template <class T> 
bool operator<(const Vector2D<T>& lhs, const Vector2D<T>& rhs);
template <class T>
bool operator==(const Vector2D<T>& lhs, const Vector2D<T>& rhs);
template <class T>
bool operator!=(const Vector2D<T>& lhs, const Vector2D<T>& rhs);

//
// invert operator (unary minus)
//
template <class T>
Vector2D<T> operator-(const Vector2D<T>& rhs);

//
// elementwise operators with Vector2D
//
template <class T, class U>
Vector2D<T>& operator+=(Vector2D<T>& lhs, const Vector2D<U>& rhs);
template <class T, class U>
Vector2D<T>& operator-=(Vector2D<T>& lhs, const Vector2D<U>& rhs);
template <class T, class U>
Vector2D<T>& operator*=(Vector2D<T>& lhs, const Vector2D<U>& rhs);
template <class T, class U>
Vector2D<T>& operator/=(Vector2D<T>& lhs, const Vector2D<U>& rhs);
template <class T, class U>
Vector2D<T> operator+(const Vector2D<T>& lhs, const Vector2D<U>& rhs);
template <class T, class U>
Vector2D<T> operator-(const Vector2D<T>& lhs, const Vector2D<U>& rhs);
template <class T, class U>
Vector2D<T> operator*(const Vector2D<T>& lhs, const Vector2D<U>& rhs);
template <class T, class U>
Vector2D<T> operator/(const Vector2D<T>& lhs, const Vector2D<U>& rhs);

//
// element-wise operators with a double
//
template <class T>
Vector2D<T>& operator+=(Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T>& operator-=(Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T>& operator*=(Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T>& operator/=(Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T> operator+(const Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T> operator+(const double& u, const Vector2D<T>& v);
template <class T>
Vector2D<T> operator-(const Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T> operator-(const double& u, const Vector2D<T>& v);
template <class T>
Vector2D<T> operator*(const Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T> operator*(const double& u, const Vector2D<T>& v);
template <class T>
Vector2D<T> operator/(const Vector2D<T>& v, const double& u);
template <class T>
Vector2D<T> operator/(const double& u, const Vector2D<T>& v);

//
// input/output
//
template <class T>
std::ostream& operator<<(std::ostream& output, const Vector2D<T>& v);
template <class T>
std::istream& operator>>(std::istream& input, const Vector2D<T>& v);

#include "Vector2D.txx"
#endif
