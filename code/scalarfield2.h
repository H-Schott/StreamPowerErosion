#pragma once

#include "array2.h"
#include "texture.h"

class ScalarField2 : public Array2
{
protected:
  std::vector<double> field; //!< Field samples.
public:
  //! Empty.
  ScalarField2() {}
  ScalarField2(const Array2&, const double& = 0.0);
  explicit ScalarField2(const Box2&, int, int, const double& = 0.0);
  explicit ScalarField2(const Box2&, int, int, const std::vector<double>&);
  explicit ScalarField2(const Box2&, const char*, double, double);

  //! Empty
  ~ScalarField2();

  void GetRange(double&, double&) const;
  double Integral() const;
  double Average() const;

  virtual Vector2 Gradient(int, int) const;
  virtual double Value(int, int) const;
  virtual double K() const;

  // Access to elements
  double at(int, int) const;
  double at(const Vec2I&) const;
  double& operator()(int, int);
  double& operator()(const Vec2I&);
  double at(int) const;
  double& operator[](int);

  ScalarField2 SetResolution(int, int, bool = false) const;
  virtual double Value(const Vector2&) const;

  void Fill(const double&);

  // Local editing
  void Gaussian(const Vector2&, const double&, const double&);
  void Smooth();
  void Smooth(int);
  double SmoothPoint(int, int) const;

  // Functions
  void Normalize();
  void Clamp(const double&, const double&);
  void SetRange(const double&, const double&);
  void Sqrt();
  ScalarField2 Sqrted() const;

  // IO
  void Save(const char* filename) const;
  Texture2D CreateImage() const;

  friend std::ostream& operator<<(std::ostream&, const ScalarField2&);
  std::vector<double> GetData() const;
  std::vector<float> GetFloatData() const;
};

inline std::vector<double> ScalarField2::GetData() const
{
	return field;
}

/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline double ScalarField2::at(int i, int j) const
{
  return field.at(VertexIndex(i, j));
}

/*!
\brief Return the field value at a given array vertex.
\param q Point.
*/
inline double ScalarField2::at(const Vec2I& q) const
{
  return field.at(VertexIndex(q[0], q[1]));
}

/*!
\brief Return the field value at a given array vertex.
\param q Point.
*/
inline double& ScalarField2::operator()(const Vec2I& q)
{
  return field[VertexIndex(q[0], q[1])];
}

/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline double& ScalarField2::operator()(int i, int j)
{
  return field[VertexIndex(i, j)];
}

inline double& ScalarField2::operator[](int i)
{
	return field[i];
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline double ScalarField2::at(int c) const
{
  return field.at(c);
}