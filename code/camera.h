#pragma once

#include "ray.h"
#include "box2.h"

// Core camera class
class Camera
{
protected:
	Vector eye;       //!< Eye.
	Vector at;        //!< Look at point.
	Vector up;        //!< Up vector.
	Vector view;      //!< View unit vector.

	double width;     //!< Screen width.
	double height;    //!< Screen height.

	double cah;       //!< Camera aperture horizontal. 
	double cav;       //!< Camera aperture vertical.
	double fl;        //!< Focal length.

	double nearplane; //!< Near plane.
	double farplane;  //!< Far plane.

public:
	Camera();
	explicit Camera(const Vector&, const Vector & = Vector::Null, const Vector & = Vector::Z, const double& = 1.0, const double& = 1.0, const double& = 1.0, const double& = 100000.0);
	explicit Camera(const Vector&, const Vector&, const Vector&, const double&, const double& = 1.0, const double& = 100000.0);

	void Translate(const Vector&);

	Vector At() const;
	Vector Eye() const;
	Vector Up() const;
	Vector View() const;

	double GetNear() const;
	double GetFar() const;
	double GetCameraApertureH() const;
	double GetCameraApertureV() const;
	double GetFocalLength() const;
	double GetAngleOfViewH() const;
	double GetAngleOfViewV(double, double) const;

	void Step(const double&);

	bool InsideFrustum(const Vector&) const;

	void LeftRight(const double&);
	void Vertical();
	void UpDown(const double&);
	void SideWay(const double&);
	void BackForth(const double&, bool = false);
	void LeftRightRound(const double&);
	void UpDownRound(const double&);

	void LeftRightFPS(const double&);
	void UpDownRoundFPS(const double&);

	void SlideHorizontal(const double&);
	void SetAt(const Vector&);
	void SetEye(const Vector&);
	void SetPlanes(const double&, const double&);

	// Move camera in a plane
	void UpDownVertical(const double&);
	void LeftRightHorizontal(const double&);

	friend std::ostream& operator<<(std::ostream&, const Camera&);

	// Pixel and sub-pixel sampling
	Ray PixelToRay(int, int, int, int) const;
	Ray PixelToRay(int, int, int, int, int, int, int) const;
	bool VectorToPixel(const Vector&, double&, double&, int, int) const;
};

//! Returns the look-at point.
inline Vector Camera::At() const
{
	return at;
}

//! Returns the eye point.
inline Vector Camera::Eye() const
{
	return eye;
}

//! Returns the up point.
inline Vector Camera::Up() const
{
	return up;
}

/*!
\brief Returns the viewing direction (not normalized).
*/
inline Vector Camera::View() const
{
	return at - eye;
}
