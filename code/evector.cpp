
// Vector  

#include "random.h"

// Self include
#include "evector.h"

/*!
\class Vector evector.h
\brief Vectors in three dimensions.

Most binary operators have been overloaded as expected,
destructive operators, such as addition and subtraction
have been implemented and behave as one could expect.

<P><I>How do I compute the cross product of two vectors?</I>
<BR>Simply use the overloaded Vector::operator/, for instance
\code
Vector c=a/b; // Cross product
\endcode
computes the cross product of a and b.
<P><I>How do I compute the sine of the angle between two vectors?</I>
<BR>Simply use the Sine(const Vector&,const Vector&) function, which internally computes the norm of the cross product divided by the norm of the argument vectors.
\code
double s=Sine(a,b); // Equivalent to Norm(a/b)/(Norm(a)*Norm(b));
\endcode
<P><I>How can I get access to the x, y and z components of a vector?</I>
<BR>Use v[0], v[1] and v[2] to get access to the x, y and z components of a vector v respectively.
<P><I>How do I compute the normal of a triangle?</I>
<BR>Let a,b,c the vertices of the triangle, simply compute the cross product
\code
Vector n=(a-b)/(a-c);  // Cross product
\endcode
or use the member function of the Triangle class:
\code
Vector n=Triangle(a,b,c).Normal(); // Compute the normal
\endcode
<P><I>How can I sort the three elements in a vector?</I>
<BR>Use Vector::Sort() as follows:
\code
Vector s=Vector(2.0,3.0,-1.0).Sort(); // Sort components in ascending order
\endcode

<P><I>How do I perform bi-linear interpolation on vectors?</I>
<BR>Use Vector::Bilinear() with four vectors and bilinear coordinates.
Alternatively, some geometric classes implement bilinear interpolation,
such as Quadrangle::Vertex().

\ingroup MathGroup
*/

const Vector Vector::Null = Vector(0.0, 0.0, 0.0);
const Vector Vector::X = Vector(1.0, 0.0, 0.0);
const Vector Vector::Y = Vector(0.0, 1.0, 0.0);
const Vector Vector::Z = Vector(0.0, 0.0, 1.0);

/*!
\brief Check if three vectors are coplanar.

Simply compute the cross product of a and b, and the dot product with c.
Compare the result with a given tolerance.

\param a,b,c Vectors.
\param epsilon Tolerance parameter.
*/
bool Vector::Coplanar(const Vector& a, const Vector& b, const Vector& c, const double& epsilon)
{
  double s = Math::Sqr((a / b) * c) / (SquaredNorm(a) * SquaredNorm(b) * SquaredNorm(c));
  return (s < epsilon* epsilon);
}

/*!
\brief Test if two vectors are almost equal.

This function computes the difference between the two argument vectors, then the norm
infinity of this difference and check the result againt the epsilon threshold value.
This is a convenience function which is the same as:
\code
Vector a,b; // Two vectors
double e; // Epsilon value
bool e=NormInfinity(Abs(b-a))<e?true:false;
\endcode
*/
bool Vector::Equal(const Vector& a, const Vector& b, const double& epsilon)
{
  Vector ab = Abs(b - a);
  if (ab[0] > epsilon || ab[1] > epsilon || ab[2] > epsilon)
    return false;
  return true;
}

/*!
\brief Normalize a vector, computing the inverse of its norm and scaling
the components.

This function does not check if the vector is null,
which might resulting in errors.
*/
void Normalize(Vector& u)
{
  u *= 1.0 / Norm(u);
}

/*!
\brief Returns the positive sine of two vectors.

Computes the cross product of the vectors and normalizes the result.
\param u, v Vectors.
*/
double Sine(const Vector& u, const Vector& v)
{
  return Norm(u / v) / sqrt((u * u) * (v * v));
}

/*!
\brief Returns the positive cosine of two vectors.

Basically computes the dot product of the normalized vectors.
\param u, v Vectors.
*/
double Cosine(const Vector& u, const Vector& v)
{
  return (u * v) / sqrt((u * u) * (v * v));
}

/*!
\brief Check if two vectors are aligned.

Computes the cosine of the two vectors, and checks for unity.
\param u, v Vectors.
*/
int Aligned(const Vector& u, const Vector& v)
{
  double c = Cosine(u, v);
  c *= c;
  return (c > (1.0 - 0.0001));
}

/*!
\brief Swap two vectors.
\param a, b Vectors.
*/
void Swap(Vector& a, Vector& b)
{
  Vector t = a;
  a = b;
  b = t;
}

/*!
\brief Swap two pointers to (arrays) vectors.
\param a, b Vectors.
*/
void Swap(Vector*& a, Vector*& b)
{
  Vector* t = a;
  a = b;
  b = t;
}

/*!
\brief Check if four points are coplanar.
\param a,b,c,d Four points.
\param epsilon Tolerance parameter.
*/
bool Vector::Coplanar(const Vector& a, const Vector& b, const Vector& c, const Vector& d, const double& epsilon)
{
  return Coplanar(b - a, c - a, d - a, epsilon);
}

/*!
\brief Compute the angle between two vectors.

The two input vectors implitly define a frame
\param b %Argument vector forming the angle.
*/
double Vector::Angle(const Vector& b) const
{
  double angle = acos(Cosine(*this, b));
  //angle = ((*this / b)[2] < 0.0) ? 2.0 * Math::Pi - angle : angle;
  return angle;
}

/*!
\brief Returns a vector orthogonal to the argument vector.

The returned orthogonal vector is not computed randomly.
First, we find the two coordinates of the argument vector with
maximum absolute value. The orthogonal vector is defined by
swapping those two coordinates and changing one sign, whereas
the third coordinate is set to 0.

The returned orthogonal vector lies in the plane orthogonal
to the first vector.
*/
Vector Vector::Orthogonal() const
{
  Vector a = Abs(*this);
  int i = 0;
  int j = 1;
  if (a[0] > a[1])
  {
    if (a[2] > a[1])
    {
      j = 2;
    }
  }
  else
  {
    i = 1;
    j = 2;
    if (a[0] > a[2])
    {
      j = 0;
    }
  }
  a = Vector::Null;
  a[i] = c[j];
  a[j] = -c[i];
  return a;
}

/*!
\brief Computes two random orthonormal vectors to the argument vector.

This function is expensive as it requires the computation
of normalized vectors and the evaluation of a random number.

\param i,j Orthonormal vectors.
*/
void Vector::RandomOrthonormal(Vector& i, Vector& j) const
{
  Vector x, y;
  Orthonormal(x, y);

  double r = Random::R239.Uniform(2.0 * Math::Pi);
  double c = cos(r);
  double s = sin(r);

  i = c * x + s * y;
  j = -s * x + c * y;
}

/*!
\brief Overloaded output-stream operator.
\param u Vector.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Vector& u)
{
  s << "Vector(" << u.c[0] << ',' << u.c[1] << ',' << u.c[2] << ')';
  return s;
}

/*!
\brief Given a vector, creates two vectors xand y that form an orthogonal basis.

This algorithm pickes the minor axis in order to reduce numerical instability
\param x, y Returned vectors such that (x,y,n) form an orthonormal basis (provided n is normalized).
*/
void Vector::Orthonormal(Vector& x, Vector& y) const
{
  x = Normalized(Orthogonal());
  y = Normalized(*this / x);
}

/*!
\brief Computes octant index of a vector with respect to the vector object.

\sa Box::Octant(const Vector&)
\param p %Vector.
*/
int Vector::Octant(const Vector& p) const
{
  return ((p[0] < c[0]) ? 0 : 1) | ((p[1] < c[1]) ? 0 : 2) | ((p[2] < c[2]) ? 0 : 4);
}

/*!
\brief Sort the terms of the vector into ascending order.
*/
Vector Vector::Sort() const
{
  // 0 > 1
  if (c[0] > c[1])
  {
    // 0 > {1,2}
    if (c[0] > c[2])
    {
      // 0 > 1 > 2
      if (c[1] > c[2])
      {
        return Vector(c[0], c[1], c[2]);
      }
      // 0 > 2 > 1
      else
      {
        return Vector(c[0], c[2], c[1]);
      }
    }
    // 0 > 1 && 2 > 0 so 2 > 0 > 1
    else
    {
      return Vector(c[2], c[0], c[1]);
    }
  }
  // 1 > 0
  else
  {
    // 1 > {0,2}
    if (c[1] > c[2])
    {
      // 1 > 0 > 2
      if (c[0] > c[2])
      {
        return Vector(c[1], c[0], c[2]);
      }
      // 1 > 2 > 0
      else
      {
        return Vector(c[1], c[2], c[0]);
      }
    }
    // 1 > 0 && 2 > 1 so 2 > 1 > 0
    else
    {
      return Vector(c[2], c[1], c[0]);
    }
  }
}

/*!
\brief Update the minimum and maximum values given a vector.
\param x Input vector.
\param a, b Lower and upper vectors that will be updated.
*/
void Vector::SetMinMax(const Vector& x, Vector& a, Vector& b)
{
  for (int i = 0; i < 3; i++)
  {
    if (x[i] < a[i])
    {
      a[i] = x[i];
    }
    else if (x[i] > b[i])
    {
      b[i] = x[i];
    }
  }
}

/*!
\brief Trilinear interpolation between eight vectors.

\sa Math::Bilinear()
\param a,b,c,d,e,f,g,h Interpolated values.
\param u,v,w Interpolation coefficients.
*/
Vector Vector::Trilinear(const Vector& a, const Vector& b, const Vector& c, const Vector& d, const Vector& e, const Vector& f, const Vector& g, const Vector& h, const double& u, const double& v, const double& w)
{
  return (1 - w) * Vector::Bilinear(a, b, c, d, u, v) + w * Vector::Bilinear(e, f, g, h, u, v);
}

