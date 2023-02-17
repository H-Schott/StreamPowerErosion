#include "texture.h"

/*!
\brief A simple class for handling a texture
*/

/*!
\brief Constructor.
\param w, h width and height
*/
Texture2D::Texture2D(int w, int h) : Array2(Box2::Unit, w, h)
{
	colors.resize(nx * ny);
}

/*!
\brief Constructor.
\param cols color array
\param w, h width and height
*/
Texture2D::Texture2D(const std::vector<Color8>& cols, int w, int h) : Array2(Box2::Unit, w, h), colors(cols)
{
}

/*!
\brief Set the texture to a single color.
\param c color
*/
void Texture2D::Fill(const Color8& c)
{
	for (int i = 0; i < colors.size(); i++)
		colors[i] = c;
}

/*!
\brief Get a constant pointer to the data.
*/
const Color8* Texture2D::Data() const
{
	return colors.data();
}
