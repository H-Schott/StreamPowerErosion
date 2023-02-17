
// Rays

#include "ray.h"

/*!
\class Ray ray.h
\brief A ray characterized by its origin an unit direction vector.

This class should be used in ray-object intersection methods and for ray-tracing.

\ingroup KernelGroup

*/

/*!
\brief Compute the reflected ray from a ray hitting a surface.

The origin of the ray is defined as the argument intersection point,
and direction is set according to the normal of the surface at
intersection point.

\param p Intersection point.
\param normal Surface outgoing normal.
*/
Ray Ray::Reflect(const Vector& p, const Vector& normal)
{
	return Ray(p, n - 2.0 * normal * (n * normal));
}

/*!
\brief Scale the ray.

This function changes both the direction vector and origin.
Note that the direction is kept unit.

\param s Scaling vector.
*/
Ray Ray::Scaled(const Vector& s) const
{
	return Ray(c.Scaled(s), Normalized(n.Scaled(s)));
}

/*!
\brief Translate the ray.

\param t Translation vector.
*/
Ray Ray::Translated(const Vector& t) const
{
	return Ray(c + t, n);
}

/*!
\brief Overloaded.

\param s Stream.
\param ray The ray.
*/
std::ostream& operator<<(std::ostream& s, const Ray& ray)
{
	s << "Ray(" << ray.c << ',' << ray.n << ')';
	return s;
}
