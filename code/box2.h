#pragma once

#include "random.h"
#include "evector.h"
#include <vector>

class Box2
{
protected:
  Vector2 a, b; //!< Lower and upper vertices of the box.
public:
  //! Empty
  Box2() {}
  explicit Box2(const Vector2&, const Vector2&);
  explicit Box2(const Vector2&, const double&);

  Vector2& operator[] (int);
  Vector2 operator[] (int) const;

  Vector2 Size() const;
  double Width() const;
  double Height() const;
  Vector2 Diagonal() const;

  // Access to vertices
  Vector2 Center() const;
  Vector2 Vertex(int) const;
  Vector2 Vertex(int, int, int, int) const;
  
  double Area() const;
  void Extend(const double&);
  void Extend(const Vector2&);
  Box2 Extended(const double&) const;

  void SetParallelepipedic(int, int&, int&);
  void SetParallelepipedic(const double&, int&, int&);

  bool Inside(const Vector2&) const;
  bool Inside(const Vector2&, const double&) const;

  friend std::ostream& operator<<(std::ostream&, const Box2&);

  Vector2 RandomInside(Random&) const;
  Vector2 RandomOn(Random&) const;
  std::vector<Vector2> Poisson(const double&, int, Random&) const;

public:
  static const double epsilon;
  static const Box2 Infinity;
  static const Box2 Null;
  static const Box2 Unit;
};

/*!
\brief Create a box.

Note that is possible to create a box using Vector as parameters
as the compiler will call the constructor Vector2::Vector2(const Vector&).

\param a,b Points.
*/
inline Box2::Box2(const Vector2& a, const Vector2& b)
{
  Box2::a = a;
  Box2::b = b;
}

/*!
\brief Creates a box.

\param c Center.
\param r Radius.
*/
inline Box2::Box2(const Vector2& c, const double& r)
{
  a = c - Vector2(r);
  b = c + Vector2(r);
}

//! Returns either end vertex of the box.
inline Vector2& Box2::operator[] (int i)
{
  if (i == 0) return a;
  else return b;
}

//! Overloaded.
inline Vector2 Box2::operator[] (int i) const
{
  if (i == 0) return a;
  else return b;
}

//! Compute the surface area of a box.
inline double Box2::Area() const
{
	return Width() * Height();
}

/*!
\brief Compute the size (width and height) of a box.
*/
inline Vector2 Box2::Size() const
{
  return b - a;
}

/*!
\brief Compute the width of a box.

\sa Box2::Size()
*/
inline double Box2::Width() const
{
  return b[0] - a[0];
}

/*!
\brief Compute the height of a box.

\sa Box2::Size()
*/
inline double Box2::Height() const
{
  return b[1] - a[1];
}

/*!
\brief Returns the diagonal of the box.
*/
inline Vector2 Box2::Diagonal() const
{
  return (b - a);
}

/*!
\brief Returns the center of the box.
*/
inline Vector2 Box2::Center() const
{
  return 0.5 * (a + b);
}

/*!
\brief Returns the k-th vertex of the box.

The returned vector is computed by analysing the first two bits of k as follows:
\code
Vector2 vertex=Vector2((k&1)?b[0]:a[0],(k&2)?b[1]:a[1]);
\endcode
*/
inline Vector2 Box2::Vertex(int k) const
{
  return Vector2((k & 1) ? b[0] : a[0], (k & 2) ? b[1] : a[1]);
}
