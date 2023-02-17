#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include "evector.h"
#include "camera.h"

class ScalarField2;
class Window;
class Texture2D;

class TerrainRaytracingWidget
{
protected:
	// Scene
	double x0, y0;					//!< Reference mouse coordinates.
	Camera camera;					//!< %Camera.

	// Heightfield data
	ScalarField2* hf;				//!< Pointer to heightfield.
	float zMin, zMax;				//!< Min/max elevation of the heighfield.
	float K;						//!< Global Lipschitz constant of the heightfield.
	int shadingMode;				//!< Heightfield shading mode.

	// GL data
	GLuint shaderProgram;			//!< GL program for shader.
	GLuint raytraceVAO;				//!< Raytracer VAO.
	GLuint hfTexture;				//!< Heightfield elevation texture.
	GLuint albedoTexture;			//!< Albedo texture for the heightfield.
	std::vector<float> tmpData;		//!< Temporary float vector.
	float cameraAngleOfViewV;		//!< Precomputed vertical camera angle of view.

	// Parent window
	Window* parent;					//!< Parent %Window class.

public:
	TerrainRaytracingWidget();
	~TerrainRaytracingWidget();

	int width() const;
	int height() const;
	void SetWindowPtr(Window* ptr);
	void ScrollCallback(GLFWwindow* w, double x, double y);
	Camera GetCamera() const;
	void SetCamera(const Camera& cam);
	void SetShadingMode(int shadingMode);
	void SetAlbedo(const Texture2D& tex);

	virtual void Update();
	virtual void SetHeightField(ScalarField2* hfPtr);
	virtual void UpdateInternal();
	virtual void initializeGL();
	virtual void paintGL();
	virtual void ReloadShaders();
};

/*!
\brief Set the camera for the widget.
*/
inline void TerrainRaytracingWidget::SetCamera(const Camera& cam)
{
	camera = cam;
}

/*!
\brief Get the current camera.
*/
inline Camera TerrainRaytracingWidget::GetCamera() const
{
	return camera;
}

/*!
\brief Set the parent window pointer.
*/
inline void TerrainRaytracingWidget::SetWindowPtr(Window* ptr)
{
	parent = ptr;
}

/*!
\brief Set the shading mode used in the shader.
*/
inline void TerrainRaytracingWidget::SetShadingMode(int mode)
{
	shadingMode = mode;
}
