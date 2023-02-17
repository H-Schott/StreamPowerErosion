// Fundamentals

#ifndef __Mathematics__
#define __Mathematics__

#include <math.h>

#include <iostream>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

/*!
\brief Long 64 bits signed integer, will use int64_t
*/
//typedef long long Int64;

/*!
\brief Minimum of two integers.
*/
inline int min(int a, int b)
{
  return (a < b ? a : b);
}

/*!
\brief Minimum of two unsigned integers.
*/
inline unsigned int min(unsigned int a, unsigned int b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two integers.
*/
inline int max(int a, int b)
{
  return (a > b ? a : b);
}

/*!
\brief Minimum of two doubles.
*/
inline double min(const double& a, const double& b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two doubles.
*/
inline double max(const double& a, const double& b)
{
  return (a > b ? a : b);
}

/*!
\brief Minimum of two doubles.
*/
inline float min(const float& a, const float& b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two doubles.
*/
inline float max(const float& a, const float& b)
{
  return (a > b ? a : b);
}

/*
inline double min(const double& a, const double& b, const double& c)
{
  return min(min(a, b), c);
}

inline double max(const double& a, const double& b, const double& c)
{
  return max(max(a, b), c);
}

inline double min(const double& a, const double& b, const double& c, const double& d)
{
  return min(min(a, b), min(c, d));
}

inline double max(const double& a, const double& b, const double& c, const double& d)
{
  return max(max(a, b), max(c, d));
}
*/

/*!
\brief Minimum of four integers.
*/
inline int min(int a, int b, int c, int d)
{
  return min(min(a, b), min(c, d));
}

/*!
\brief Minimum of four unsigned integers.
*/
inline unsigned int min(unsigned int a, unsigned int b, unsigned int c, unsigned int d)
{
  return min(min(a, b), min(c, d));
}

/*!
\brief Maximum of four integers.
*/
inline int max(int a, int b, int c, int d)
{
  return max(max(a, b), max(c, d));
}

//! Clamps an integer value between two bounds.
inline int Clamp(int x, int a, int b)
{
  return (x < a ? a : (x > b ? b : x));
}

/*!
\brief Compute the next integer as s<SUP>p</SUP> >=n.
\param n Integer.

Set all bits on the right-hand side of the most significant set bit to 1 and then increment the value by 1 to &#8220;rollover&#8221; to two&#8217;s nearest power.
*/
inline int NextPower2(int n)
{
  // Decrement n to handle the case when n itself is a power of 2
  n--;

  // Set all bits after the last set bit
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;

  // Increment n and return
  return ++n;
}

/*!
\brief Compute the logarithm of n in base 2.
\param n Integer.
*/
inline int Log2(int n)
{
  int l = 0;
  while (n >>= 1) { l++; }
  return l;
}

//! Clamps a float value between two bounds.
inline float Clamp(float x, float a, float b)
{
  return (x < a ? a : (x > b ? b : x));
}

//! Clamps a double value between two bounds.
inline double Clamp(double x, double a, double b)
{
  return (x < a ? a : (x > b ? b : x));
}

// Math class for constants
class Math
{
public:
  static const double Pi; //!< &pi;.
  static const double TwoPi; //!< 2&pi;.
  static const double HalfPi; //!< Half of &pi;.
  static const double e; //!< Exponential.
  static const double TwoPiOverThree; //!< 2/3 &pi;.
  static const double FourPiOverThree; //!< 4/3 &pi;.
  static const double Sqrt5; //!< Constant &radic;5.
  static const double Sqrt3; //!< Constant &radic;3.
  static const double Sqrt2; //!< Constant &radic;2.
  static const double Golden; //!< Constant (&radic;5+1)/2.
  static const double Large; //!< Large constant 10<SUP>20</SUP>.
  static const double Infinity; //!< Infinity.

  static int binomials[16][16]; //!< %Array of binomial coefficients.
public:
  static double Clamp(const double&, const double& = 0.0, const double& = 1.0);
  static int Clamp(int, int = 0, int = 255);

  // Squares
  static double Sqr(const double&);
  static  double Cube(const double&);
  static double Sqr4(const double&);
  static double SymmetricSqr(const double&);
  static double Pow(const double&, const double&);
  static double Abs(const double&);
  static double Sqrt32(const double&);
  static double Sqrt4(const double&);

  // Series of powers
  static void Powers(const double&, int, double*);

  // Modulo
  static  double Mod(const double&, const double&);

  static double Floor(const double&);
  static double Fract(const double&);
  static double FractFloor(const double&, double&);
  static double Ceil(const double&);

  // Minimum and maximum
  static double Min(const double&, const double&);
  static double Max(const double&, const double&);
  static double Min(const double&, const double&, const double&);
  static double Max(const double&, const double&, const double&);
  static double Min(const double&, const double&, const double&, const double&);
  static double Max(const double&, const double&, const double&, const double&);
  static double MaxArray(double*, int);

  // Angles
  static  double DegreeToRadian(const double&);
  static  double RadianToDegree(const double&);
  static  double Angle(int, int);
  static  double Angle(const double&);
  static double ArcTan(const double&, const double&);

  static double Step(double, double);
  static double Cycloidal(const double&);
  static double Triangle(const double&);

  // Linear interpolation
  static double Lerp(const double&, const double&, const double&);
  static double Bilinear(const double&, const double&, const double&, const double&, const double&, const double&);
  static double Trilinear(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);

  // Cubic interpolation
  static double Cubic(const double&, const double&, const double&, const double&, const double&);
  static double CubicPoints(const double&, const double&, const double&, const double&, const double&);
  static double BiCubic(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double& = 0.0, const double& = 0.0, const double& = 0.0, const double& = 0.0);
  static double BiCubic(const double&, const double&, const double&, const double&, const double&, const double&);

  // Sigmoid
  static double Sigmoid(const double&, const double& = 1.0);
  static double SigmoidQuadric(const double&);

  // Impulse
  static double Impulse(const double&, const double& = 1.0);

  // Gain warping
  static double Warp(const double&);
  static double Warp(const double&, const double&);

  static bool InRange(const double&, const double&, const double&, const double&);
  static void SetMinMax(const double&, double&, double&);

  static int Integer(const double&);
  static void Sort(double&, double&);
  static void Sort(double&, double&, double&);
  static void Swap(double&, double&);

  static bool IsNumber(double);
  static bool IsFinite(double);

  // Signs
  static int IntegerSign(const double&, const double& = 0.0);
  static bool SameSign(const double&, const double&);
  static bool SameSign(const double&, const double&, const double&);
  static double CopySign(const double&, const double&);

  static long Binomial(int, int);

  static  double Unit(int, int);

  static double Gaussian(const double&, const double&);
  static double Fract(const double&, double&);
  static double Fract(const double&, const double&, double&);
  static double Fract(const double&, const double&, const double&, double&);

  static double Geometric(double, int);

  static void Swap(double*&, double*&);
  static void Swap(int*&, int*&);

  static void BernsteinSeries(const double&, int, double*);
  static double Bernstein(const double&, int, int);
  
  static double LinearStep(const double& x, const double& a, const double& b);
  static double CubicSmooth(double x, double r);
};

/*!
\brief Compute the value of a C<SUP>2</SUP> smooth interpolating function (1-x/r)<SUP>3</SUP>.
\param x Squared distance.
\param r Squared radius.
*/
inline double Math::CubicSmooth(double x, double r)
{
    return (1.0f - x / r) * (1.0f - x / r) * (1.0f - x / r);
}

/*!
\brief Check if a values lies within a prescribed range.
\param x Value.
\param a, b Interval.
\param epsilon Epsilon tolerance.
*/
inline bool Math::InRange(const double& x, const double& a, const double& b, const double& epsilon)
{
  return ((x + epsilon) >= a) && ((x - epsilon) <= b);
}

/*!
\brief Bi-linear interpolation between four values.

The values are given in trigonometric order.

\image html bilinear.png
\param a00, a10, a11, a01 Interpolated values.
\param u,v Interpolation coefficients.
*/
inline double Math::Bilinear(const double& a00, const double& a10, const double& a11, const double& a01, const double& u, const double& v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}

/*!
\brief Trilinear interpolation between eight values.

\sa Math::Bilinear()
\param a,b,c,d,e,f,g,h Interpolated values.
\param u,v,w Interpolation coefficients.
*/
inline double Math::Trilinear(const double& a, const double& b, const double& c, const double& d, const double& e, const double& f, const double& g, const double& h, const double& u, const double& v, const double& w)
{
  return (1 - w) * Math::Bilinear(a, b, c, d, u, v) + w * Math::Bilinear(e, f, g, h, u, v);
}

/*!
\brief Squares a double value.
\sa Math::SymmetricSqr
\param x Real value.
*/
inline double Math::Sqr(const double& x)
{
  return x * x;
}

/*!
\brief Compute x<SUP>3/2</SUP> as the square root of x<SUP>3</SUP>.

This function is more efficient than pow() which is extremely slow.
\code
double y=Math::Sqrt32(7.0); // Equivalent but faster than y=pow(7.0,1.5);
\endcode
\param x Real value.
*/
inline double Math::Sqrt32(const double& x)
{
  return sqrt(x * x * x);
}

/*!
\brief Compute the fourth root of x<SUP>3</SUP>.

This function is more efficient than pow() which is extremely slow.
\code
double y=Math::Sqrt4(7.0); // Equivalent to y=sqrt(sqrt(7.0));
\endcode
\param x Real value.
*/
inline double Math::Sqrt4(const double& x)
{
  return sqrt(sqrt(x));
}

/*!
\brief Cubes a double value.
\param x Real value.
*/
inline  double Math::Cube(const double& x)
{
  return x * x * x;
}

/*!
\brief Fourth power of a double value.
\param x Real value.
*/
inline double Math::Sqr4(const double& x)
{
  double y = x * x;
  return y * y;
}

/*!
\brief Symmetric of square function over unit interval.

Simply compute 1-(1-x)<SUP>2</SUP>.
\image html mathquadrics.png
\sa Math::Sqr
\param x Real value.
*/
inline double Math::SymmetricSqr(const double& x)
{
  return 1.0 - Math::Sqr(1.0 - x);
}

/*!
\brief Sort two reals.

\sa Math::Swap(double&,double&)
\param a, b Real arguments, which will be swapped so that a<=b.
*/
inline void Math::Sort(double& a, double& b)
{
  if (a > b)
  {
    double c = a;
    a = b;
    b = c;
  }
}

/*!
\brief Sort three reals.

\sa Math::Sort(double&,double&)
\param a, b, c Real arguments, which will be swapped so that a<=b<=c.
*/
inline void Math::Sort(double& a, double& b, double& c)
{
  if (a > b)
  {
    Swap(a, b);
  }
  // a b C : do nothing
  if (c >= b)
    return;

  // C a b : move and place
  if (c <= a)
  {
    double t = c;
    c = b;
    b = a;
    a = t;
    return;
  }
  // a C b : swap b and c
  Swap(b, c);
}

/*!
\brief Linear interpolation.

Returns (1-t)a+tb.

\param a,b Interpolated values.
\param t Interpolant.
*/
inline double Math::Lerp(const double& a, const double& b, const double& t)
{
  return a + t * (b - a);
}

/*!
\brief Convert degrees to randians.
\param a Angle in degrees.
*/
inline  double Math::DegreeToRadian(const double& a)
{
  return a * Math::Pi / 180.0;
}

/*!
\brief Convert radian to degrees.
\param a Angle in radian.
*/
inline  double Math::RadianToDegree(const double& a)
{
  return a * 180.0 / Math::Pi;
}

/*!
\brief Compute 2 k &pi; / n.
\param k, n Integers.
*/
inline  double Math::Angle(int k, int n)
{
  return (2.0 * Math::Pi * k) / n;
}

/*!
\brief Swap two reals.
\sa Sort(double&,double&)
\param a, b Arguments.
*/
inline void Math::Swap(double& a, double& b)
{
  double t = a;
  a = b;
  b = t;
}

/*!
\brief Swap two integers.
\sa Swap(double&,double&)
\param a, b Arguments.
*/
inline void Swap(int& a, int& b)
{
  int t = a;
  a = b;
  b = t;
}

/*!
\brief Compute the integer part of a real.

This function handles negative values differently by subtracting 1 from the result.
\param x %Real.
*/
inline int Math::Integer(const double& x)
{
  return (x > 0.0 ? int(x) : int(x) - 1);
}

/*!
\brief Modulus for reals with negative values handled properly.
\param x %Real.
\param a Modulo.
*/
inline  double Math::Mod(const double& x, const double& a)
{
  return (x >= 0.0) ? fmod(x, a) : a - fmod(-x, a);
}

/*!
\brief Angle strictly in [0,2\pi[.
\param a %Angle.
*/
inline  double Math::Angle(const double& a)
{
  return Math::Mod(a, 2.0 * Math::Pi);
}

/*!
\brief Clamp a double value between two bounds.
\param x Input value.
\param a, b Lower and upper bounds.
*/
inline double Math::Clamp(const double& x, const double& a, const double& b)
{
  return (x < a ? a : (x > b ? b : x));
}

/*!
\brief Clamp an integer value between two bounds.
\param x Input value.
\param a, b Lower and upper bounds.
*/
inline int Math::Clamp(int x, int a, int b)
{
  return (x < a ? a : (x > b ? b : x));
}

/*!
\brief Update the minimum and maximum values given a double value.
\param x Input value.
\param a, b Lower and upper bounds that will be updated according to x.
*/
inline void Math::SetMinMax(const double& x, double& a, double& b)
{
  if (x < a)
  {
    a = x;
  }
  else if (x > b)
  {
    b = x;
  }
}

/*!
\brief Minimum of two reals.
\param a, b Real values.
*/
inline double Math::Min(const double& a, const double& b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two reals.
\param a, b Real values.
*/
inline double Math::Max(const double& a, const double& b)
{
  return (a > b ? a : b);
}

/*!
\brief Maximum of three reals.
\param a, b, c Real values.
*/
inline double Math::Max(const double& a, const double& b, const double& c)
{
  return Math::Max(Math::Max(a, b), c);
}

/*!
\brief Minimum of three reals.
\param a, b, c Real values.
*/
inline double Math::Min(const double& a, const double& b, const double& c)
{
  return Math::Min(Math::Min(a, b), c);
}

/*!
\brief Maximum of four reals.
\param a, b, c, d Real values.
*/
inline double Math::Max(const double& a, const double& b, const double& c, const double& d)
{
  return Math::Max(Math::Max(a, b), Math::Max(c, d));
}

/*!
\brief Minimum of four reals.
\param a, b, c, d Real values.
*/
inline double Math::Min(const double& a, const double& b, const double& c, const double& d)
{
  return Math::Min(Math::Min(a, b), Math::Min(c, d));
}

/*!
\brief Fractional part of a real.

Implemented as:
\code
double y=x-Math::Floor(x);
\endcode
\param x Real.
*/
inline double Math::Fract(const double& x)
{
  return x - floor(x);
}

/*!
\brief Fractional part of a real.

Implemented as:
\code
double y=x-Math::Floor(x);
\endcode
\param x Real.
\param y Return floor.
*/
inline double Math::FractFloor(const double& x, double& y)
{
  y = floor(x);
  return x - y;
}

/*!
\brief %Floor function.

While it is easier to use c++ function, this has been implemented to be consistent with Fract().
\param x Real.
*/
inline double Math::Floor(const double& x)
{
  return floor(x);
}

/*!
\brief %Ceil function.

\sa Floor
\param x Real.
*/
inline double Math::Ceil(const double& x)
{
  return ceil(x);
}

/*!
\brief Power.

Returns base x raised to the power exponent e if x>0, and -x raised to the power exponent e otherwise.

\param x Base.
\param e Exponent.
*/
inline double Math::Pow(const double& x, const double& e)
{
  if (x == 0.0)
  {
    return 0.0;
  }
  else if (x > 0.0)
  {
    return pow(x, e);
  }
  else
  {
    return -pow(-x, e);
  }
}

/*!
\brief Absolute value.

\param x Real.
*/
inline double Math::Abs(const double& x)
{
  return fabs(x);
}

/*!
\brief Unit real value in [0,1] from two integers.

Usefull for loops, same as:
\code
for (int i=0;i<n;i++)
{
double t=double(i)/double(n-1); // t=Math::Unit(i,n);
}
\endcode
\param i, n Integers.
*/
inline  double Math::Unit(int i, int n)
{
  return double(i) / double(n - 1);
}

/*!
\brief Sigmoid-like function.

This function is more efficient than the real sigmoid function which requires the computation of the exponential.

\param s Sigma.
\param x Real.
*/
inline double Math::Sigmoid(const double& x, const double& s)
{
  return x / sqrt(s * s + x * x);
}

/*!
\brief Impulse function.

An impulse function that doesn't use exponentials, k controls the
falloff of the function.

Derivative is 2*sqrt(k)*(1-k*x*x)/(k*x*x+1), maximum is reached x = sqrt(1/k), value at maximum is 1.

\param k Control parameter.
\param x Real.
*/
inline double Math::Impulse(const double& x, const double& k)
{
  return 2.0 * sqrt(k) * x / (1.0 + k * x * x);
}

/*!
\brief Compactly supported sigmoid-like function implemented using a C<SUP>1</SUP> piecewise symmetric quadric.

The compact support is [-2,2] and the range [-1,1].

The quadric was obtaied by solving the Hermite Cubic contraints:
\code
Cubic c=2.0*Cubic::Compose(Cubic::Hermite(0, 0.5, 1.0, 0.0), Linear(0.5, 0.0)) <<std::endl;
\endcode

\sa Math::Sigmoid

\param x Real.
*/
inline double Math::SigmoidQuadric(const double& x)
{
  if (x > 0.0)
  {
    if (x > 2.0)
    {
      return 1.0;
    }
    else
    {
      return x * (1.0 - 0.25 * x);
    }
  }
  else
  {
    if (x < -2.0)
    {
      return -1.0;
    }
    else
    {
      return x * (1.0 + 0.25 * x);
    }
  }
}

/*!
\copydoc Quadric::Warp(const double&)
*/
inline double Math::Warp(const double& x)
{
  const double y = 2.0 * x;
  if (x < 0.5)
  {
    return 0.5 * sqrt(y);
  }
  else
  {
    return 1.0 - 0.5 * sqrt(2.0 - y);
  }
}

/*!
\brief Unit interval gain function.

Remaps the unit interval into the unit interval by expanding the sides and compressing the center, keeping 1/2 mapped to 1/2.

This function uses calls to pow and is computationally intensive.

\sa Quadric::Warp(const double&)
\param x Real in [0,1].
\param k Exponent.
*/
inline double Math::Warp(const double& x, const double& k)
{
  const double a = 0.5 * pow(2.0 * ((x < 0.5) ? x : 1.0 - x), k);
  return (x < 0.5) ? a : 1.0 - a;
}

/*!
\brief Check if two reals have the same signs.
\param a,b Reals.
*/
inline bool Math::SameSign(const double& a, const double& b)
{
  return signbit(a) == signbit(b);
}

/*!
\brief Check if three reals have the same signs.
\param a,b,c Reals.
*/
inline bool Math::SameSign(const double& a, const double& b, const double& c)
{
  return (signbit(a) == signbit(b)) && (signbit(a) == signbit(c));
}

/*!
\brief Composes a real with the magnitude of x and the sign of y.
\param x Real.
\param y Real whose sign will be copied.
*/
inline double Math::CopySign(const double& x, const double& y)
{
  return copysign(x, y);
}

/*!
\brief Compute the integer sign of a real.
\param x Real.
\param t Threshold.
*/
inline int Math::IntegerSign(const double& x, const double& t)
{
  if (x < -t)
  {
    return -1;
  }
  else if (x > t)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

/*!
\brief Compute integer and fractionnal parts.
\param x Real.
\param i Returned integer part.
\return Fractional part.
*/
inline double Math::Fract(const double& x, double& i)
{
  i = floor(x);
  return (x - i);
}

/*!
\brief Compute integer and fractionnal parts.
\param x Real.
\param a Modulus.
\param i Returned integer part.
\return Fractional part.
*/
inline double Math::Fract(const double& x, const double& a, double& i)
{
  double y = x / a;
  i = floor(y);
  return (y - i) * a;
}

/*!
\brief Compute integer and fractionnal parts given an interval
\param x Real.
\param a,b Interval.
\param i Returned integer part.
\return Fractional part.
*/
inline double Math::Fract(const double& x, const double& a, const double& b, double& i)
{
  return Math::Fract(x - a, b - a, i);
}

inline void Math::Swap(double*& a, double*& b)
{
  double* t = a;
  a = b;
  b = t;
}

inline void Math::Swap(int*& a, int*& b)
{
  int* t = a;
  a = b;
  b = t;
}

/*!
\brief Gaussian.
\param x Real.
\param s Sigma.
*/
inline double Math::Gaussian(const double& x, const double& s)
{
  return exp(-(x * x) / (s * s));
}

/*!
\brief Compute the powers series of a real value 1, x, ... , x<SUP>n-1</SUP>.
\param x Real.
\param n Maximum power.
\param a Array of real, should be at least of size n.
*/
inline void Math::Powers(const double& x, int n, double* a)
{
  a[0] = 1.0;
  a[1] = x;
  double y = x;
  for (int k = 2; k < n; k++)
  {
    y *= x;
    a[k] = y;
  }
}

/*!
\brief Step function.
\param e Step position.
\param x Variable.
*/
inline double Math::Step(double e, double x)
{
  return (x < e) ? 0.0 : 1.0;
}

/*!
\brief Create a linear step.
\param x Value
\param a, b Interval values.
\return Real in unit inverval.
\sa Cubic::SmoothStep, Quintic::SmoothStep
*/
inline double Math::LinearStep(const double& x, const double& a, const double& b)
{
    if (x < a)
    {
        return 0.0;
    }
    else if (x > b)
    {
        return 1.0;
    }
    else
    {
        return (x - a) / (b - a);
    }
}

#endif
