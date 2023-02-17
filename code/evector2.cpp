// Vector  

// Self include
#include "evector.h"

/*!
\class Vector2 evector.h
\brief Vectors in two dimensions.

This class implements most operators and member functions as for Vector class.

\ingroup MathGroup
*/

const Vector2 Vector2::Null = Vector2(0.0, 0.0);

const Vector2 Vector2::X = Vector2(1.0, 0.0);
const Vector2 Vector2::Y = Vector2(0.0, 1.0);

/*!
\brief Normalize a two dimensional vector.

It computes the inverse of its norm and scaling the components.
This function does not check if the vector is null, which might yield errors.

\param u %Vector.
*/
void Normalize(Vector2& u)
{
    u *= 1.0 / Norm(u);
}

/*!
\brief Computes quadrant index of a vector with respect to the vector object.

\sa Box2::Quadrant(const Vector&)
\param p %Vector.
*/
int Vector2::Quadrant(const Vector2& p) const
{
    return ((p[0] < c[0]) ? 0 : 1) | ((p[1] < c[1]) ? 0 : 2);
}

/*!
\brief Returns the sine of two vectors.

Computes the cross product of the vectors and normalizes the result.
\param u, v Vectors.
*/
double Sine(const Vector2& u, const Vector2& v)
{
    return (u[0] * v[1] - u[1] * v[0]) / sqrt((u * u) * (v * v));
}

/*!
\brief Returns the positive cosine of two vectors.

Computes the dot product of the normalized vectors.
\param u, v Vectors.
*/
double Cosine(const Vector2& u, const Vector2& v)
{
    return (u * v) / sqrt((u * u) * (v * v));
}

/*!
\brief Swap two vectors.
\param a, b Vectors.
*/
void Swap(Vector2& a, Vector2& b)
{
    Vector2 t = a;
    a = b;
    b = t;
}

/*!
\brief Update the minimum and maximum values given a vector.
\param x Input vector.
\param a, b Lower and upper vectors that will be updated.
*/
void Vector2::SetMinMax(const Vector2& x, Vector2& a, Vector2& b)
{
    for (int i = 0; i < 2; i++)
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
\brief Test if two vectors are almost equal.

\sa Vector::Equal
*/
bool Vector2::Equal(const Vector2& a, const Vector2& b, const double& epsilon)
{
    Vector2 ab = Abs(b - a);
    if (ab[0] > epsilon || ab[1] > epsilon)
        return false;
    return true;
}

/*!
\brief Compute the angle of a vector.

The vector need not be normalized.
*/
double Vector2::Angle() const
{
    return atan2(c[1], c[0]);
}

/*!
\brief Compute the angle between two vectors.

Note that both vectors should be unit.

\sa Vector2::Angle() const

\param b %Argument vector forming the angle.
*/
double Vector2::Angle(const Vector2& b) const
{
    Vector2 ba(b * (*this), b * Orthogonal());
    return ba.Angle();
}

/*!
\brief Overloaded output-stream operator.
\param u %Vector.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Vector2& u)
{
    s << "Vector2(" << u.c[0] << ',' << u.c[1] << ')';
    return s;
}

/*!
\brief Sort the terms of the vector into ascending order.
*/
Vector2 Vector2::Sort() const
{
    // 0 > 1
    if (c[0] > c[1])
    {
        return Vector2(c[1], c[0]);
    }
    // 1 > 0
    else
    {
        return *this;
    }
}
