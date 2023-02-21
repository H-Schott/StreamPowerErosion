#include <immintrin.h>
#include "scalarfield2.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/*!
\class ScalarField2 scalarfield.h
\brief A base two-dimensional field of real values.

\ingroup Structure
*/

/*!
\brief Create the field structure.
\param a Array representing the grid domain.
\param v Constant value of field.
*/
ScalarField2::ScalarField2(const Array2& a, const double& v) :Array2(a)
{
  field.resize(nx * ny, v);
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y Size of the array.
\param v Constant value of field.
*/
ScalarField2::ScalarField2(const Box2& box, int x, int y, const double& v) :ScalarField2(Array2(box, x, y), v)
{
}

/*!
\brief Create the field structure.
\param box The box.
\param x,y Size of the array.
\param v Array of scalar values.
*/
ScalarField2::ScalarField2(const Box2& box, int x, int y, const std::vector<double>& v) :Array2(box, x, y)
{
  field = v;
}

/*!
\brief Create the field structure.
\param box The box.
\param name filename to read
\param minV, maxV min/max range of values
*/
ScalarField2::ScalarField2(const Box2& box, const char* filename, double minV, double maxV) 
{
    a = box[0];
    b = box[1];
    int n;
    unsigned char* rawData = stbi_load(filename, &nx, &ny, &n, 1);
    field.resize(nx * ny);
    for (int i = 0; i < nx * ny; i++)
    {
        double t = double(rawData[i]) / 255.0;
        field[i] = minV * (1.0f - t) + t * maxV;
    }
    stbi_image_free(rawData);

    celldiagonal = Vector2((b[0] - a[0]) / (nx - 1), (b[1] - a[1]) / (ny - 1));
    inversecelldiagonal = celldiagonal.Inverse();
}

/*!
\brief Empty.
*/
ScalarField2::~ScalarField2()
{
}

/*!
\brief Get the integral of the scalar field.

Compute the sum of all the elements and make the product by the area of small cell element.

\sa Average()
*/
double ScalarField2::Integral() const
{
    double s = Average();

    // Scale by the size of a cell element
    s *= Box2(a, b).Area();

    return s;
}

/*!
\brief Compute the average value of the elements in the scalar field.

\sa Integral()
*/
double ScalarField2::Average() const
{
  const int size = nx * ny; 
  double s = 0.0;
  // Process all elements

  for (int i = 0; i < size; i++)
  {
    s += field.at(i);
  }

  // Scale by the size of a cell element
  s /= size;
  return s;
}

/*!
\brief Get the range of the field.

\param a,b Returned minimum and maximum.
*/
void ScalarField2::GetRange(double& a, double& b) const
{
  const int size = nx * ny;
  // Escape
  if (size == 0)
  {
    a = b = 0.0;
    return;
  }

  a = field.at(0);
  b = a;
  
  for (int i = 1; i < field.size(); i++)
  {
    double x = field.at(i);
    if (x < a)
    {
    a = x;
    }
    else if (x > b)
    {
    b = x;
    }
  }
}

/*!
\brief Get the field value with world coordinate system.

This function computes a bi-linear interpolation of the values.

\sa Math::Bilinear
\param p Point position (should be strictly inside the box domain).
*/
double ScalarField2::Value(const Vector2& p) const
{
  double u, v;
  int i, j;
  CellInteger(p, i, j, u, v);

  // Test position
  if (!InsideCellIndex(i, j))
  {
    return 0.0;
  }

  return Math::Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), u, v);
}

/*!
\brief Sets the entire field with a constant value.
\param s Scalar.
*/
void ScalarField2::Fill(const double& s)
{
  field.resize(nx * ny, s);
}

/*!
\brief Change the resolution of the scalar field.

Note that because the box should be strictly inside the original domain, this
function is not the same as:
\code
ScalarField field(Box2(Vector2(-2.0),Vector(3.0),5,5); // Original field
ScalarField s=Sample(field.GetBox(),12,12);
\endcode
\param x,y Sampling size.
\param bicubic Bicubic flag, set to true to use bicubic interpolation.
*/
ScalarField2 ScalarField2::SetResolution(int x, int y, bool bicubic) const
{
  // Sampled scalar field
  ScalarField2 sampled(Box2(a, b), x, y);

  // Corners
  sampled(0, 0) = at(0, 0);
  sampled(0, y - 1) = at(0, ny - 1);
  sampled(x - 1, 0) = at(nx - 1, 0);
  sampled(x - 1, y - 1) = at(nx - 1, ny - 1);

  if (bicubic == false)
  {
    // Borders (use linear interpolation)
    for (int i = 1; i < x - 1; i++)
    {
      double tx = (nx - 1) * (i / double(x - 1));
      int x0 = int(floor(tx));
      int x1 = int(ceil(tx));

      sampled(i, 0) = Math::Lerp(at(x0, 0), at(x1, 0), tx - x0);
      sampled(i, y - 1) = Math::Lerp(at(x0, ny - 1), at(x1, ny - 1), tx - x0);
    }

    for (int j = 1; j < y - 1; j++)
    {
      double ty = (ny - 1) * (j / double(y - 1));
      int y0 = int(floor(ty));
      int y1 = int(ceil(ty));

      sampled(0, j) = Math::Lerp(at(0, y0), at(0, y1), ty - y0);
      sampled(x - 1, j) = Math::Lerp(at(nx - 1, y0), at(nx - 1, y1), ty - y0);
    }

    // Interior
    for (int i = 1; i < x - 1; i++)
    {
      for (int j = 1; j < y - 1; j++)
      {
        sampled(i, j) = Value(sampled.ArrayVertex(i, j));
      }
    }
  }
  else
  {
    // Edges
    for (int i = 1; i < x - 1; i++)
    {
      double tx = (nx - 1) * (i / double(x - 1));
      int x0 = int(floor(tx));
      int x1 = int(ceil(tx));

      sampled(i, 0) = Math::Lerp(at(x0, 0), at(x1, 0), tx - x0);
      sampled(i, y - 1) = Math::Lerp(at(x0, ny - 1), at(x1, ny - 1), tx - x0);
    }

    for (int j = 1; j < y - 1; j++)
    {
      double ty = (ny - 1) * (j / double(y - 1));
      int y0 = int(floor(ty));
      int y1 = int(ceil(ty));

      sampled(0, j) = Math::Lerp(at(0, y0), at(0, y1), ty - y0);
      sampled(x - 1, j) = Math::Lerp(at(nx - 1, y0), at(nx - 1, y1), ty - y0);
    }

    // Interior
    for (int i = 1; i < x - 1; i++)
    {
      for (int j = 1; j < y - 1; j++)
      {
        sampled(i, j) = Value(sampled.ArrayVertex(i, j));
      }
    }
  }
  return sampled;
}

/*!
\brief Return the field gradient at a given array vertex.
\param i,j Integer coordinates of the array vertex.
\sa at(int,int)
*/
double ScalarField2::Value(int i, int j) const
{
  return field.at(VertexIndex(i, j));
}

/*!
\brief Compute the gradient at a given array vertex.

\param i,j Integer coordinates of the array vertex.
*/
Vector2 ScalarField2::Gradient(int i, int j) const
{
  Vector2 n;

  // Gradient along x axis
  if (i == 0)
  {
    n[0] = (at(i + 1, j) - at(i, j)) * inversecelldiagonal[0];
  }
  else if (i == nx - 1)
  {
    n[0] = (at(i, j) - at(i - 1, j)) * inversecelldiagonal[0];
  }
  else
  {
    n[0] = (at(i + 1, j) - at(i - 1, j)) * 0.5 * inversecelldiagonal[0];
  }

  // Gradient along y axis
  if (j == 0)
  {
    n[1] = (at(i, j + 1) - at(i, j)) * inversecelldiagonal[1];
  }
  else if (j == ny - 1)
  {
    n[1] = (at(i, j) - at(i, j - 1)) * inversecelldiagonal[1];
  }
  else
  {
    n[1] = (at(i, j + 1) - at(i, j - 1)) * 0.5 * inversecelldiagonal[1];
  }

  return n;
}

/*!
\brief Compute the Lipschitz constant of the elevation function.

Note that this is equivalent to the following code, with the
difference that this function does not store the norm of the gradient in memory.
\code
ScalarField2 n=f.GradientNorm();
double k=0.0;
for (int i=0;i<f.VertexSize();i++)
{
k=Math::Max(k,f.at(i));
}
\endcode
*/
double ScalarField2::K() const
{
  double k = 0.0;

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      k = Math::Max(k, Norm(Gradient(i, j)));
    }
  }
  return k;
}

/*!
\brief Normalize the values of a scalar field to unit interval.

Compute range using GetRange() and apply the affine transformation to map values to [0,1].

\sa Unitize()
*/
void ScalarField2::Normalize()
{
  const int size = nx * ny;

  // Escape
  if (size == 0)
    return;

  double a, b;
  GetRange(a, b);
  if (a == b)
  {
    field.resize(size, 1.0);
    return;
  }
 
  for (int i = 0; i < field.size(); i++)
  {
    field[i] = (field[i] - a) / (b - a);
  }
}

/*!
\brief Clamp the values of a scalar field.
\param a,b Interval.
*/
void ScalarField2::Clamp(const double& a, const double& b)
{
  for (int i = 0; i < field.size(); i++)
  {
    field[i] = Math::Clamp(field.at(i), a, b);
  }
}

/*!
\brief Scales the scalar field to a given range interval.

\param a,b Interval.
*/
void ScalarField2::SetRange(const double& a, const double& b)
{
  double x, y;
  GetRange(x, y);
  if (x == y)
  {
    field.resize(field.size(), a);
  }
  else
  {
    const double c = (b - a) / (y - x);
    for (int i = 0; i < field.size(); i++)
    {
      field[i] = a + c * (field.at(i) - x);
    }
  }
}

void ScalarField2::Sqrt() {
    for (int i = 0; i < field.size(); i++) {
        field[i] = std::sqrt(field.at(i));
    }
}

ScalarField2 ScalarField2::Sqrted() const {
    ScalarField2 sf = *this;
    sf.Sqrt();
    return sf;
}


/*!
\brief Add material with gaussian distribution.
\param center Center of the distribution.
\param radius Radius
\param height Maximum height of the distribution.
*/
void ScalarField2::Gaussian(const Vector2& center, const double& radius, const double& height)
{
    // Compute modification Area
    Vec2I pa, pb;
    VertexIntegerArea(Box2(center, radius), pa, pb);

    // Compute thickness
    for (int y = pa[1]; y <= pb[1]; y++)
    {
        for (int x = pa[0]; x <= pb[0]; x++)
        {
            // Distance between central point and current point
            double u = SquaredNorm(center - ArrayVertex(x, y));
            if (u < radius * radius)
                field[VertexIndex(x, y)] += height * Math::CubicSmooth(u, radius * radius);
        }
    }
}

/*!
\brief Overloaded.
\param s Stream.
\param scalar The scalar field.
*/
std::ostream& operator<<(std::ostream& s, const ScalarField2& scalar)
{
  s << Array2(scalar) << std::endl;

  for (int j = scalar.ny - 1; j > -1; j--)
  {
    for (int i = 0; i < scalar.nx; i++)
    {
      s << scalar.at(i, j) << ' ';
    }
    s << std::endl;
  }
  return s;
}

/*
\brief Compute a smoothed value at a given point in the scalar field.

The function uses a 3<SUP>2</SUP> approximation of the Gaussian kernel.
\param x, y Integer coordinates of the point.
*/
double ScalarField2::SmoothPoint(int x, int y) const
{
  int k = VertexIndex(x, y);

  if (x == 0)
  {
    // Corner
    if (y == 0)
    {
      return (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k + ny)) / 8.0;
    }
    // Corner
    else if (y == ny - 1)
    {
      return (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k - ny)) / 8.0;
    }
    else
    {
      return (2.0 * at(k - ny) + 4.0 * at(k) + 2.0 * at(k + ny) + at(k + 1 - ny) + 2.0 * at(k + 1) + at(k + 1 + ny)) / 12.0;
    }
  }
  else if (x == nx - 1)
  {
    // Corner
    if (y == 0)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + ny)) / 8.0;
    }
    // Corner
    else if (y == ny - 1)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k - ny)) / 8.0;
    }
    else
    {
      return (2.0 * at(k - ny) + 4.0 * at(k) + 2.0 * at(k + ny) + at(k - 1 - ny) + 2.0 * at(k - 1) + at(k - 1 + ny)) / 12.0;
    }
  }
  else
  {
    if (y == 0)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k + nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 12.0;

    }
    else if (y == ny - 1)
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + at(k - 1 - nx) + at(k + 1 - nx)) / 12.0;
    }
    // Center
    else
    {
      return (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - ny) + 2.0 * at(k + ny) + at(k - 1 - ny) + at(k + 1 - ny) + at(k - 1 + ny) + at(k + 1 + ny)) / 16.0;
    }
  }
}

/*!
\brief Applies several smoothing steps to the scalar field.

\sa Smooth()
\param n Number of smoothing steps.
*/
void ScalarField2::Smooth(int n)
{
  for (int i = 0; i < n; i++)
  {
    Smooth();
  }
}

/*!
\brief Smooth the scalar field using a discrete gaussian kernel.

The function uses a 3<SUP>2</SUP> approximation of the Gaussian kernel.
*/
void ScalarField2::Smooth()
{
    std::vector<double> smoothed;
    smoothed.resize(nx * ny);

    int k;

    // Smooth center
    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {
            k = VertexIndex(i, j);
            smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + 2.0 * at(k + nx) + at(k - 1 - nx) + at(k + 1 - nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 16.0;
        }
    }

    // Smooth edges
    for (int i = 1; i < nx - 1; i++)
    {
        k = VertexIndex(i, 0);
        smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k + nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 12.0;

        k = VertexIndex(i, ny - 1);
        smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + at(k - 1 - nx) + at(k + 1 - nx)) / 12.0;
    }

    for (int j = 1; j < ny - 1; j++)
    {
        k = VertexIndex(0, j);
        smoothed[k] = (2.0 * at(k - nx) + 4.0 * at(k) + 2.0 * at(k + nx) + at(k + 1 - nx) + 2.0 * at(k + 1) + at(k + 1 + nx)) / 12.0;

        k = VertexIndex(nx - 1, j);
        smoothed[k] = (2.0 * at(k - nx) + 4.0 * at(k) + 2.0 * at(k + nx) + at(k - 1 - nx) + 2.0 * at(k - 1) + at(k - 1 + nx)) / 12.0;
    }

    // Corners
    k = VertexIndex(0, 0);
    smoothed[k] = (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k + nx) + 1.0 * at(k + nx + 1)) / 9.0;

    k = VertexIndex(nx - 1, 0);
    smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k + nx) + 1.0 * at(k + nx - 1)) / 9.0;

    k = VertexIndex(0, ny - 1);
    smoothed[k] = (4.0 * at(k) + 2.0 * at(k + 1) + 2.0 * at(k - nx) + 1.0 * at(k - nx + 1)) / 9.0;

    k = VertexIndex(nx - 1, ny - 1);
    smoothed[k] = (4.0 * at(k) + 2.0 * at(k - 1) + 2.0 * at(k - nx) + 1.0 * at(k - nx - 1)) / 9.0;

    // Center
    field = smoothed;
}

/*!
\brief Export the scalarfield as an jpeg file.
\param filename file name
*/
void ScalarField2::Save(const char* filename) const
{
    ScalarField2 ff = *this;
    ff.Normalize();
    unsigned char* rawData = new unsigned char[ff.VertexSize() * 3];
    for (int i = 0; i < ff.VertexSize(); i++)
    {
        double t = ff.at(i) * 255.0;
        rawData[(i * 3) + 0] = (unsigned char)(int)t;
        rawData[(i * 3) + 1] = (unsigned char)(int)t;
        rawData[(i * 3) + 2] = (unsigned char)(int)t;
    }
    stbi_write_jpg(filename, nx, ny, 3, rawData, 98);
    delete[] rawData;
}

Texture2D ScalarField2::CreateImage() const {
    Color8 Cool(97, 130, 234, 255);
    Color8 White(221, 221, 221, 255);
    Color8 Warm(220, 94, 75, 255);

    double low, high;
    GetRange(low, high);
    double norm_coeff = 1. / (high - low);
    std::vector<Color8> colors(nx * ny, Color8());

    for (int i = 0; i < nx * ny; i++) {
        double u = norm_coeff * (field[i] - low);

        if (u < 0.5) {
            colors[i] = Color8::Lerp(u / 0.5, Cool, White);
        } else {
            colors[i] = Color8::Lerp((u - 0.5) / 0.5, White, Warm);
        }
    }

    return Texture2D(colors, nx, ny);
}

std::vector<float> ScalarField2::GetFloatData() const {
    std::vector<float> res(nx * ny, 0);
    for (int i = 0; i < nx * ny; i++) res[i] = float(field[i]);
    return res;
}
