// Mathematics
#include <ostream>
#include <limits>
#include "mathematics.h"

/*!
\defgroup MathGroup Core math classes.

\brief Core math classes include several classes such as Vector, Quadric, Cubic and higher order
polynomials and many others that are useful in many graphic applications.
*/

/*!
\class Math mathematics.h
\brief Core class implementing some useful functions and constants.

<P><I>How can I use the constant Pi?</I>
<BR>Simply use the static constant Math::Pi as follows:
\code
double v = 4.0*Math::Pi*r*r*r/3.0; // Volume of a sphere.
\endcode
Note that in this case, you could also have used:
\code
double v = Sphere::Volume(r);
\endcode

<P><I>How many min/max functions have been implemented?</I>
<BR>Up to four arguments are supported; for a larger number of
arguments, a specific routine operating on an array should be written.

<P><I>Is there a function to compute the square of a real sumber?</I>
<BR>Use the following:
\code
double s = Math::Sqr((1.0-sqrt(5.0))/2.0); // Square of golden ratio
\endcode
For a Vector, use SquaredNorm(const Vector&);

<P><I>Are there predefined square roots constants?</I>
<BR>The sqrt function used not te be constexpr, so square roots of reals are not computed at compilation time.
Some constants are provided, such as the following one:
\code
double s = Math::Sqrt3; // Square root of 3
\endcode

<P><I>What is the relative performance of mathematical functions such as square root or cosine?</I>
<BR>The relative cost of commonly used mathematical functions is summarized in the following table.

<table>
<caption>Cost</caption>
<tr><th>Function<th>Cost
<tr><td>+ - * fabs<td>1
<tr><td>/ sqrt modf <td>~4.8
<tr><td>/ sin cos <td>~8.9
<tr><td>/ exp log <td>~10.5
</table>

Math provides some functions such as Math::Sqrt32 or Math::Sqrt4 to avoid computationally intensive standard functions like pow:
\code
double y = Math::Sqrt32(x); // Faster than double y = pow(x,1.5);
\endcode

\ingroup MathGroup

<P><I>How are implemented the step and smooth-step functions that are often used in procedural modeling?</I>
Different smoothing kernels, such as Cubic::Smooth(), are implented in odd-degree polynomials Cubic, Quintic and Septic.
The corresponding step functions, such as Cubic::SmoothStep(), are also implemented.
\sa Linear::Step, Cubic::Smooth, Quintic::Smooth, Cubic::SmoothStep, Quintic::SmoothStep, Septic::SmoothStep.

*/

const double Math::Pi = 3.14159265358979323846;

const double Math::TwoPi = 6.28318530717958647693;

const double Math::HalfPi = Math::Pi / 2.0;

const double Math::e = 2.7182818284590452354;

const double Math::TwoPiOverThree = 2.0943951023931954923084;

const double Math::FourPiOverThree = 4.1887902047863909846168;

const double Math::Infinity = std::numeric_limits<double>::max();

const double Math::Large = 1.0e20;

const double Math::Sqrt5 = sqrt(5.0);

const double Math::Sqrt3 = sqrt(3.0);

const double Math::Sqrt2 = sqrt(2.0);

const double Math::Golden = (sqrt(5.0) + 1.0) / 2.0;

// B binomial coefficients up to 15
int Math::binomials[16][16] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 3, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 5, 10, 10, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 6, 15, 20, 15, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 7, 21, 35, 35, 21, 7, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 8, 28, 56, 70, 56, 28, 8, 1, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 9, 36, 84, 126, 126, 84, 36, 9, 1, 0, 0, 0, 0, 0, 0 },
  { 1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1, 0, 0, 0, 0, 0 },
  { 1, 11, 55, 165, 330, 462, 462, 330, 165, 55, 11, 1, 0, 0, 0, 0 },
  { 1, 12, 66, 220, 495, 792, 924, 792, 495, 220, 66, 12, 1, 0, 0, 0 },
  { 1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1, 0, 0 },
  { 1, 14, 91, 364, 1001, 2002, 3003, 3432, 3003, 2002, 1001, 364, 91, 14, 1, 0 },
  { 1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1 }
};

/*!
\brief Sine wave over unit interval.
\param x Real value.
*/
double Math::Cycloidal(const double& x)
{
  return sin(x - floor(x) * 2.0 * Math::Pi);
}

/*!
\brief %Triangle wave over unit interval.
\param x Real value.
*/
double Math::Triangle(const double& x)
{
  double offset;

  if (x >= 0.0)
  {
    offset = x - floor(x);
  }
  else
  {
    offset = x - (-1.0 - floor(fabs(x)));
  }
  if (offset >= 0.5)
  {
    return (2.0 * (1.0 - offset));
  }
  else
  {
    return (2.0 * offset);
  }
}

/*!
\brief Check if a real number is not NaN.

The code is simply:
\code
return (x == x);
\endcode
This looks like it should always be true, but it's false if x is a NaN.
\param x Real value.
*/
bool Math::IsNumber(double x)
{
  return (x == x);
}

/*!
\brief Check if a real is not infinite.
\param x Real value.
*/
bool Math::IsFinite(double x)
{
  return (x <= Math::Infinity && x >= -Math::Infinity);
}

/*!
\brief Calculate the %Binomial coefficent.
\param n, r %Binomial coefficient parameters.
*/
long Math::Binomial(int n, int r)
{
  if (n < 0 || r<0 || r>n)
  {
    return 0;
  }
  else if (r == n)
  {
    return 1;
  }
  else if (r < 15 && n < 15)
  {
    return binomials[n][r];
  }
  else
  {
    long b = n;
    for (int i = 1; i < r; i++)
    {
      b *= (n - i);
      b /= (i + 1);
    }

    return b;
  }
}

/*!
\brief %Cubic Hermite interpolation.

\param u Interpolation coefficient.
\param a, b Values for t=0 and t=1.
\param ta, tb Derivatives for t=0 and t=1.

\sa Cubic::Hermite()
*/
double Math::Cubic(const double& u, const double& a, const double& b, const double& ta, const double& tb)
{
  return (((tb + ta + 2.0 * (a - b)) * u + -tb - 2.0 * ta + 3.0 * (b - a)) * u + ta) * u + a;
}

/*!
\brief %Cubic point interpolation.

\param u Interpolation coefficient.
\param a, b,c,d Values for t=-1, t=0, t=1 and t=2. The derivatives at t=0 and t=1 will be derived from those

\sa Math::Cubic()
*/
double Math::CubicPoints(const double& u, const double& a, const double& b, const double& c, const double& d)
{
  return Math::Cubic(u, b, c, 0.5 * (c - a), 0.5 * (d - b));
}

/*!
\brief Bi-cubic interpolation between four values, given partial derivatives.

The values are given in trigonometric order.

\param u,v Interpolation coefficients.
\param a00,a10,a11,a01 Interpolated values.
\param u00,u10,u11,u01,v00,v10,v11,v01 Partial derivatives with respect to u and v.
\param x00,x10,x11,x01 Cross derivatives.
*/
double Math::BiCubic(const double& u, const double& v, const double& a00, const double& a10, const double& a11, const double& a01, const double& u00, const double& u10, const double& u11, const double& u01, const double& v00, const double& v10, const double& v11, const double& v01, const double& x00, const double& x10, const double& x11, const double& x01)
{
  double u2 = u * u;
  double v2 = v * v;

  double hu00 = u2 * (2.0 * u - 3.0) + 1.0;
  double hu01 = u2 * (-2.0 * u + 3.0);
  double hu10 = (u2 - 2.0 * u + 1.0) * u;
  double hu11 = u2 * (u - 1.0);

  double hv00 = v2 * (2.0 * v - 3.0) + 1.0;
  double hv01 = v2 * (-2.0 * v + 3.0);
  double hv10 = (v2 - 2.0 * v + 1.0) * v;
  double hv11 = v2 * (v - 1.0);

  return hu00 * (a00 * hv00 + a01 * hv01 + v00 * hv10 + v01 * hv11) +
    hu01 * (a10 * hv00 + a11 * hv01 + v10 * hv10 + v11 * hv11) +
    hu10 * (u00 * hv00 + u01 * hv01 + x00 * hv10 + x01 * hv11) +
    hu11 * (u10 * hv00 + u11 * hv01 + x10 * hv10 + x11 * hv11);
}

/*!
\brief Bi-cubic interpolation between four values, partial derivatives are implicitly defined as null.

This is an optimized implementation.

\param u,v Interpolation coefficients.
\param a00,a10,a11,a01 Interpolated values.
*/
double Math::BiCubic(const double& u, const double& v, const double& a00, const double& a10, const double& a11, const double& a01)
{
  double u2 = u * u;
  double v2 = v * v;

  double hu00 = u2 * (2.0 * u - 3.0) + 1.0;
  double hu01 = u2 * (-2.0 * u + 3.0);

  double hv00 = v2 * (2.0 * v - 3.0) + 1.0;
  double hv01 = v2 * (-2.0 * v + 3.0);

  return hu00 * (a00 * hv00 + a01 * hv01) + hu01 * (a10 * hv00 + a11 * hv01);
}

/*
\brief Maximum of an array of reals.
\param a %Array.
\param n Size.
*/
double Math::MaxArray(double* a, int n)
{
  double t = a[0];
  for (int i = 0; i < n; i++)
  {
    if (a[i] > t)
    {
      t = a[i];
    }
  }
  return t;
}

/*!
\brief Inline version of the atan2() function.
*/

double Math::ArcTan(const double& y, const double& x)
{
  if (x > 0.0)
  {
    return atan(y / x);
  }
  else if (x < 0.0)
  {
    if (y >= 0.0)
    {
      return atan(y / x) + Math::Pi;
    }
    else
    {  // y < 0
      return atan(y / x) - Math::Pi;
    }
  }
  else
  {  // x == 0
    if (y > 0.0)
    {
      return Math::HalfPi;
    }
    else if (y < 0.0)
    {
      return -Math::HalfPi;
    }
    else {  // y == 0
      return 0.0;
    }
  }
}

/*!
\brief Compute the sum of the terms of a geometric series.

Compute 1+x+x<SUP>2</SUP>+...+x<SUP>n-1</SUP>=(1-x<SUP>n</SUP>)/(1-x).
*/
double Math::Geometric(double x, int n)
{
  return (1.0 - pow(x, double(n))) / (1.0 - x);
}

/*!
\brief Compute the Bernstein terms.
*/
void Math::BernsteinSeries(const double& u, int n, double* b)
{
  static double ui[16];
  static double umi[16];

  Math::Powers(u, n, ui);
  Math::Powers(1.0 - u, n, umi);

  for (int i = 0; i < n; i++)
  {
    b[i] = Math::binomials[n - 1][i] * ui[i] * umi[n - 1 - i];
  }
}

/*!
\brief Compute the i-th Bernstein term.
*/
double Math::Bernstein(const double& u, int n, int i)
{
  if (n == 0)
    return 1.0;

  if (n == 1)
  {
    if (i == 0)
    {
      return u;
    }
    else if (i == n)
    {
      return 1.0 - u;
    }
  }

  double ui = 1.0;
  double uni = 1.0;
  for (int k = 1; k <= i; k++)
  {
    ui *= u;
  }
  for (int k = 1; k <= n - i; k++)
  {
    uni *= 1.0 - u;
  }
  return Math::binomials[n - 1][i] * ui * uni;
}