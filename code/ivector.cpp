// Vector  

// Self include
#include "ivector.h"

/*!
\class Vec3I ivector.h
\brief Integer vectors in three dimensions.

\ingroup MathGroup
*/

/*!
\brief Overloaded output-stream operator.
\param u Vector.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Vec3I& u)
{
  s << "Vec3I(" << u.c[0] << ',' << u.c[1] << ',' << u.c[2] << ')';
  return s;
}



/*!
\class Vec2I ivector.h
\brief Integer vectors in two dimensions.

\ingroup MathGroup
*/

/*!
\brief Overloaded output-stream operator.
\param u Vector.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Vec2I& u)
{
  s << "Vec2I(" << u.c[0] << ',' << u.c[1] << ')';
  return s;
}

