#pragma once

#include "box2.h"
#include "ivector.h"

class Array2 : protected Box2
{
protected:
  int nx=0, ny=0; //!< Sizes.
  Vector2 celldiagonal= Vector2::Null; //!< Cell diagonal.
  Vector2 inversecelldiagonal= Vector2::Null; //!< Inverse cell diagonal.
public:
  Array2();
  explicit Array2(const Box2&, int, int);
  explicit Array2(const Box2&, int);

  //! Empty
  ~Array2() {}

  int VertexSize() const;
  int GetSizeX() const;
  int GetSizeY() const;

  int CellSize() const;
  int CellSizeX() const;
  int CellSizeY() const;
  Vector2 CellCenter(int, int) const;

  Box2 Cell(int i, int j) const;
  Box2 Cell(int c) const;
  Box2 GetBox() const;
  Vector2 CellDiagonal() const;

  Vector2 ArrayVertex(int, int) const;
  Vector2 ArrayVertex(const Vec2I&) const;

  Vector2 Size() const;

  // Subdivision
  void Subdivide();

  // Empty
  bool IsEmpty() const;

  // Vertex queries
  void VertexInteger(const Vector2&, int&, int&) const;
  Vec2I VertexInteger(const Vector2&) const;
  Vec2I VertexInteger(const Vector2&, double&, double&) const;
  void VertexIntegerArea(const Box2& box, Vec2I& pa, Vec2I& pb) const;

  // Cell queries
  void CellInteger(const Vector2&, int&, int&) const;
  void CellInteger(const Vector2&, int&, int&, double&, double&) const;
  Vec2I CellInteger(const Vector2&) const;
  Vec2I CellInteger(const Vector2&, double&, double&) const;

  friend std::ostream& operator<<(std::ostream&, const Array2&);

  // Domain queries
  constexpr bool InsideVertexIndex(int, int) const;
  constexpr bool OutsideVertexIndex(int, int) const;
  bool InsideVertexIndex(const Vec2I&) const;
  constexpr bool InsideVertexIndex(int, int, int) const;

  bool Inside(const Vector2&) const;

  constexpr bool BorderVertexIndex(int, int) const;
  bool BorderVertexIndex(const Vec2I&) const;

  // Indexes for storing elements at vertices
  constexpr int VertexIndex(int, int) const;
  int VertexIndex(const Vec2I&) const;

  Vec2I Next(const Vec2I&, int) const;
  void ClampVertexIndex(int&, int&) const;
  
protected:
  constexpr void InverseVertexIndex(int, int&, int&) const;
  Vec2I InverseVertexIndex(int) const;

  // Indexes for storing elements at cells
  constexpr int CellIndex(int, int) const;
  constexpr bool InsideCellIndex(int, int) const;
  bool InsideCellIndex(const Vec2I&) const;
  void InverseCellIndex(int, int&, int&) const;
  
public:
  static int NeighborCode(int, int);
  static Vec2I codeToDir(int);
protected:
  static const Vec2I next[8]; //!< Array of points in the 1-ring neighborhood.
  static const double length[8]; //!< Length to the i-th neighbor.
  static const double inverselength[8]; //!< Inverse length.
};

/*!
\brief Return the size of the array.
*/
inline Vector2 Array2::Size() const
{
  return Box2::Size();
}

/*!
\brief Detect if the array is empty, i.e., any dimension equal to zero.
*/
inline bool Array2::IsEmpty() const
{
  return (nx <= 0) || (ny <= 0);
}

/*!
\brief Get the vertex size of the array for x axis.
*/
inline int Array2::GetSizeX() const
{
  return nx;
}

/*!
\brief Get the vertex size of the array for y axis.
*/
inline int Array2::GetSizeY() const
{
  return ny;
}

/*!
\brief Get the cell size of the array for x axis.
*/
inline int Array2::CellSizeX() const
{
  return nx - 1;
}

/*!
\brief Get the cell size of the array for y axis.
*/
inline int Array2::CellSizeY() const
{
  return ny - 1;
}

/*!
\brief Return the size of the vertex array.
*/
inline int Array2::VertexSize() const
{
  return nx * ny;
}

/*!
\brief Return the size of the cell array.
*/
inline int Array2::CellSize() const
{
  return (nx - 1) * (ny - 1);
}

/*!
\brief Compute the coordinates of a point on the grid.
\param i,j Integer coordinates.
*/
inline Vector2 Array2::ArrayVertex(int i, int j) const
{
  return Vector2(a[0] + i * celldiagonal[0], a[1] + j * celldiagonal[1]);
}

/*!
\brief Compute the coordinates of a point on the grid.
\param p Point.
*/
inline Vector2 Array2::ArrayVertex(const Vec2I& p) const
{
  return Vector2(a[0] + p[0] * celldiagonal[0], a[1] + p[1] * celldiagonal[1]);
}

/*!
\brief Get the box of the array.
*/
inline Box2 Array2::GetBox() const
{
  return Box2(a, b);
}

/*!
\brief Compute the index of a given cell.
\param i,j Integer coordinates of the cell.
*/
inline constexpr int Array2::VertexIndex(int i, int j) const
{
  return i + nx * j;
}
/*!
\brief Compute the index of a given cell.
\param p Point.
*/
inline int Array2::VertexIndex(const Vec2I& p) const
{
  return p[0] + nx * p[1];
}

/*!
\brief Compute the coordinates of a given cell.
\param c Index of the cell.
\param i,j Integer coordinates of the cell.
*/
inline constexpr void Array2::InverseVertexIndex(int c, int& i, int& j) const
{
  i = c % nx;
  j = c / nx;
}

/*!
\brief Compute the coordinates of a given cell.
\param c Index of the cell.
*/
inline Vec2I Array2::InverseVertexIndex(int c) const
{
  return Vec2I(c % nx, c / nx);
}

/*!
\brief Check if the indexes are within range.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::InsideCellIndex(int i, int j) const
{
  return (i >= 0) && (i < nx - 1) && (j >= 0) && (j < ny - 1);
}

/*!
\brief Check if the indexes are within range.
\param p Point.
*/
inline bool Array2::InsideCellIndex(const Vec2I& p) const
{
  return (p[0] >= 0) && (p[0] < nx - 1) && (p[1] >= 0) && (p[1] < ny - 1);
}

/*!
\brief Check if the indexes are within range.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::InsideVertexIndex(int i, int j) const
{
  return (i >= 0) && (i < nx) && (j >= 0) && (j < ny);
}

/*!
\brief Check if the indexes are outside or on the border.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::OutsideVertexIndex(int i, int j) const
{
  if (i <= 0 || j <= 0 || i >= nx - 1 || j >= ny - 1) return false;
  return true;
}

/*!
\brief Check if the indexes are within range.
\param p Point.
*/
inline bool Array2::InsideVertexIndex(const Vec2I& p) const
{
  return (p[0] >= 0) && (p[0] < nx) && (p[1] >= 0) && (p[1] < ny);
}

/*!
\brief Check if the indexes are on the border.
\param i,j Integer coordinates of the vertex.
*/
inline constexpr bool Array2::BorderVertexIndex(int i, int j) const
{
  return (i == 0) || (i == nx - 1) || (j == 0) || (j == ny - 1);
}

/*!
\brief Check if the indexes are on the border.
\param p Point.
*/
inline bool Array2::BorderVertexIndex(const Vec2I& p) const
{
  return (p[0] == 0) || (p[0] == nx - 1) || (p[1] == 0) || (p[1] == ny - 1);
}

/*!
\brief Check if a point is in the rectangular domain.
\param p Point.
*/
inline bool Array2::Inside(const Vector2& p) const
{
  return Box2::Inside(p);
}

/*!
\brief Check if the indexes are within k-range.
\param i,j Integer coordinates of the vertex.
\param k Thickness of the boundary around the domain
*/
inline constexpr bool Array2::InsideVertexIndex(int i, int j, int k) const
{
  return (i >= 0 + k) && (i < nx - k) && (j >= 0 + k) && (j < ny - k);
}

/*!
\brief Compute the coordinates of a given cell.
\param c Index of the cell.
\param i,j Integer coordinates of the cell.
*/
inline void Array2::InverseCellIndex(int c, int& i, int& j) const
{
  i = c % (nx - 1);
  j = c / (nx - 1);
}

/*!
\brief Compute the point on the grid given an input point.
\param p Point.
\param i,j Integer coordinates of the cell.
\param u,v Coordinates of the point in the corresponding cell.
*/
inline void Array2::CellInteger(const Vector2& p, int& i, int& j, double& u, double& v) const
{
  Vector2 q = p - a;

  /*
    Vector2 d = b - a;

    u = q[0] / d[0];
    v = q[1] / d[1];

    // Scale
    u *= (nx - 1);
    v *= (ny - 1);
  */
  u = q[0] * inversecelldiagonal[0];
  v = q[1] * inversecelldiagonal[1];

  // Integer coordinates
  i = int(u);
  j = int(v);

  // Local coordinates within cell
  u -= i;
  v -= j;
}

/*!
\brief Compute the point on the grid given an input point.
\param p Point.
\param u,v Coordinates of the point in the corresponding cell.
*/
inline Vec2I Array2::CellInteger(const Vector2& p, double& u, double& v) const
{
  Vector2 q = p - a;

  /*
    Vector2 d = b - a;

    u = q[0] / d[0];
    v = q[1] / d[1];

    // Scale
    u *= (nx - 1);
    v *= (ny - 1);
  */
  u = q[0] * inversecelldiagonal[0];
  v = q[1] * inversecelldiagonal[1];

  // Integer coordinates
  int i = int(u);
  int j = int(v);

  // Local coordinates within cell
  u -= i;
  v -= j;

  return Vec2I(i, j);
}

/*!
\brief Compute the point on the grid given an input point.
\param p Point.
\param u, v Coordinates of the point in the corresponding cell.
*/
inline Vec2I Array2::VertexInteger(const Vector2& p, double& u, double& v) const
{
  Vector2 q = p - a;

  /*
    Vector2 d = b - a;

    u = q[0] / d[0];
    v = q[1] / d[1];

    // Scale
    u *= (nx - 1);
    v *= (ny - 1);
  */
  u = q[0] * inversecelldiagonal[0];
  v = q[1] * inversecelldiagonal[1];

  // Integer coordinates
  int i = int(u);
  int j = int(v);

  // Local coordinates within cell
  u -= i;
  v -= j;

  return Vec2I(i, j);
}

/*!
\brief Clamp vertex indexes to the size of the array.
\param i,j %Vertex indexes
*/
inline void Array2::ClampVertexIndex(int& i, int& j) const
{
  if (i < 0) { i = 0; }
  if (i > nx - 1) { i = nx - 1; }
  if (j < 0) { j = 0; }
  if (j > nx - 1) { j = nx - 1; }
}

/*!
\brief Compute the index of a given cell.
\param i,j Integer coordinates of the cell.
*/
inline constexpr int Array2::CellIndex(int i, int j) const
{
  return i + (nx - 1) * j;
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
\param i,j Integer coordinates of the cell.
*/
inline void Array2::VertexInteger(const Vector2& p, int& i, int& j) const
{
  Vector2 q = p - a;
  /*
  Vector2 d = b - a;

  double u = q[0] / d[0];
  double v = q[1] / d[1];

  i = int(u * (nx - 1));
  j = int(v * (ny - 1));
  */
  i = int(q[0] * inversecelldiagonal[0]);
  j = int(q[1] * inversecelldiagonal[1]);
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
*/
inline Vec2I Array2::VertexInteger(const Vector2& p) const
{
  int i, j;
  VertexInteger(p, i, j);
  return Vec2I(i, j);
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
*/
inline Vec2I Array2::CellInteger(const Vector2& p) const
{
  int i, j;
  CellInteger(p, i, j);
  return Vec2I(i, j);
}

/*!
\brief Compute the coordinates of a vertex inside a cell.
\param p Point.
\param i,j Integer coordinates of the cell.
*/
inline void Array2::CellInteger(const Vector2& p, int& i, int& j) const
{
  Vector2 q = p - a;
  /*
  Vector2 d = b - a;
  double u = q[0] / d[0];
  double v = q[1] / d[1];
  i = int(u * (nx - 1));
  j = int(v * (ny - 1));
  */
  i = int(q[0] * inversecelldiagonal[0]);
  j = int(q[1] * inversecelldiagonal[1]);
}

/*!
\brief Compute the point next to another one.
\param p Point.
\param n Next neighbor, should be in [0,7].
*/
inline Vec2I Array2::Next(const Vec2I& p, int n) const
{
  return p + next[n];
}