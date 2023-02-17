// Vector  

#ifndef __IntVector__
#define __IntVector__

#include <ostream>

// Class
class Vec3I
{
protected:
  int c[3] = { 0,0,0 }; //!< Integer coordinates.
public:
  //! Empty 
  Vec3I() {}

  Vec3I(int);
  Vec3I(int, int, int);

  // Access members
  int& operator[] (int);
  int operator[] (int) const;

  // Unary operators
  Vec3I operator+ () const;
  Vec3I operator- () const;

  // Assignment operators
  Vec3I& operator+= (const Vec3I&);
  Vec3I& operator-= (const Vec3I&);
  Vec3I& operator*= (const Vec3I&);
  Vec3I& operator*= (int);

  // Binary operators
  friend Vec3I operator+ (const Vec3I&, const Vec3I&);
  friend Vec3I operator- (const Vec3I&, const Vec3I&);

  // Boolean functions
  friend bool operator==(const Vec3I&, const Vec3I&);
  friend bool operator!=(const Vec3I&, const Vec3I&);

  friend std::ostream& operator<<(std::ostream&, const Vec3I&);
};

/*!
\brief Create a vector with the same coordinates.
\param a Integer.
*/
inline Vec3I::Vec3I(int a)
{
  c[0] = c[1] = c[2] = a;
}

/*!
\brief Create a vector with argument coordinates.
\param a,b,c Coordinates.
*/
inline Vec3I::Vec3I(int a, int b, int c)
{
  Vec3I::c[0] = a;
  Vec3I::c[1] = b;
  Vec3I::c[2] = c;
}

//! Gets the i-th coordinate of vector.
inline int& Vec3I::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of vector.
inline int Vec3I::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vec3I Vec3I::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vec3I Vec3I::operator- () const
{
  return Vec3I(-c[0], -c[1], -c[2]);
}

// Assignment unary operators

//! Destructive addition.
inline Vec3I& Vec3I::operator+= (const Vec3I& u)
{
  c[0] += u.c[0]; c[1] += u.c[1]; c[2] += u.c[2];
  return *this;
}

//! Destructive subtraction.
inline Vec3I& Vec3I::operator-= (const Vec3I& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1]; c[2] -= u.c[2];
  return *this;
}

//! Destructive scalar multiply.
inline Vec3I& Vec3I::operator*= (int a)
{
  c[0] *= a; c[1] *= a; c[2] *= a;
  return *this;
}

/*!
\brief Destructively scale a vector by another vector.
*/
inline Vec3I& Vec3I::operator*= (const Vec3I& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1]; c[2] *= u.c[2];
  return *this;
}

//! Adds up two vectors.
inline Vec3I operator+ (const Vec3I& u, const Vec3I& v)
{
  return Vec3I(u.c[0] + v.c[0], u.c[1] + v.c[1], u.c[2] + v.c[2]);
}

//! Difference between two vectors.
inline Vec3I operator- (const Vec3I& u, const Vec3I& v)
{
  return Vec3I(u.c[0] - v.c[0], u.c[1] - v.c[1], u.c[2] - v.c[2]);
}

// Boolean functions

//! Strong equality test.
inline bool operator== (const Vec3I& u, const Vec3I& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]) && (u.c[2] == v.c[2]));
}

//! Strong difference test.
inline bool operator!= (const Vec3I& u, const Vec3I& v)
{
  return (!(u == v));
}

// Class
class Vec2I
{
protected:
  int c[2] = { 0,0 }; //!< Integer coordinates.
public:
  //! Empty 
  Vec2I() {}

  Vec2I(int);
  Vec2I(int, int);

  // Access members
  int& operator[] (int);
  int operator[] (int) const;

  // Unary operators
  Vec2I operator+ () const;
  Vec2I operator- () const;

  // Assignment operators
  Vec2I& operator+= (const Vec2I&);
  Vec2I& operator-= (const Vec2I&);
  Vec2I& operator*= (const Vec2I&);
  Vec2I& operator*= (int);

  // Binary operators
  friend Vec2I operator+ (const Vec2I&, const Vec2I&);
  friend Vec2I operator- (const Vec2I&, const Vec2I&);

  // Boolean functions
  friend bool operator==(const Vec2I&, const Vec2I&);
  friend bool operator!=(const Vec2I&, const Vec2I&);

  friend std::ostream& operator<<(std::ostream&, const Vec2I&);
};

/*!
\brief Create a vector with the same coordinates.
\param a Integer.
*/
inline Vec2I::Vec2I(int a)
{
  c[0] = c[1] = a;
}

/*!
\brief Create a vector with argument coordinates.
\param a,b Coordinates.
*/
inline Vec2I::Vec2I(int a, int b)
{
  Vec2I::c[0] = a;
  Vec2I::c[1] = b;
}

//! Gets the i-th coordinate of vector.
inline int& Vec2I::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of vector.
inline int Vec2I::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vec2I Vec2I::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vec2I Vec2I::operator- () const
{
  return Vec2I(-c[0], -c[1]);
}

// Assignment unary operators

//! Destructive addition.
inline Vec2I& Vec2I::operator+= (const Vec2I& u)
{
  c[0] += u.c[0]; c[1] += u.c[1];
  return *this;
}

//! Destructive subtraction.
inline Vec2I& Vec2I::operator-= (const Vec2I& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1];
  return *this;
}

//! Destructive scalar multiply.
inline Vec2I& Vec2I::operator*= (int a)
{
  c[0] *= a; c[1] *= a;
  return *this;
}

/*!
\brief Destructively scale a vector by another vector.
*/
inline Vec2I& Vec2I::operator*= (const Vec2I& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1];
  return *this;
}

//! Adds up two vectors.
inline Vec2I operator+ (const Vec2I& u, const Vec2I& v)
{
  return Vec2I(u.c[0] + v.c[0], u.c[1] + v.c[1]);
}

//! Difference between two vectors.
inline Vec2I operator- (const Vec2I& u, const Vec2I& v)
{
  return Vec2I(u.c[0] - v.c[0], u.c[1] - v.c[1]);
}

// Boolean functions

//! Strong equality test.
inline bool operator== (const Vec2I& u, const Vec2I& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]));
}

//! Strong difference test.
inline bool operator!= (const Vec2I& u, const Vec2I& v)
{
  return (!(u == v));
}

#endif

