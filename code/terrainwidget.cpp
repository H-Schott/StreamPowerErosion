#include "terrainwidget.h"
#include "shader-api.h"
#include "box2.h"
#include "scalarfield2.h"
#include "window.h"
#include "texture.h"

#include <imgui.h>

/*!
\class TerrainRaytracingWidget terrainwidget.h
\brief %Heightfield rendering widget based on sphere-tracing. The API is done so as to mimick Qt6 QOpenGLWidget API as close as possible.
*/

/*!
\brief Default constructor. No call to OpenGL here.
*/
TerrainRaytracingWidget::TerrainRaytracingWidget()
{
	parent = nullptr;
	hf = nullptr;
	hfTexture = albedoTexture = 0;
	raytraceVAO = shaderProgram = 0;
	shadingMode = 0;
	K = 1.0f;
	zMin = zMax = 0.0f;
	cameraAngleOfViewV = 0.0;
	x0 = y0 = 0.0;
}

/*!
\brief Destructor. Release shader program, vao and textures.
*/
TerrainRaytracingWidget::~TerrainRaytracingWidget()
{
	release_program(shaderProgram);
	glDeleteVertexArrays(1, &raytraceVAO);
	glDeleteTextures(1, &hfTexture);
	glDeleteTextures(1, &albedoTexture);
}

/*!
\brief Initialize OpenGL, shaders, texture buffers, and a camera centered at origin.
*/
void TerrainRaytracingWidget::initializeGL()
{
	ReloadShaders();

	camera = Camera(Vector(-1500.0, -1500.0, 250.0));
	cameraAngleOfViewV = float(camera.GetAngleOfViewV(width(), height()));

	glGenVertexArrays(1, &raytraceVAO);
	glGenTextures(1, &hfTexture);
	glGenTextures(1, &albedoTexture);
}

/*!
\brief Renders the scene.
*/
void TerrainRaytracingWidget::paintGL()
{
	// Clear
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Draw
	glDisable(GL_DEPTH_TEST);
	glUseProgram(shaderProgram);

	// Uniforms - Camera
	glUniform3f(glGetUniformLocation(shaderProgram, "CamPos"), float(camera.Eye()[0]), float(camera.Eye()[1]), float(camera.Eye()[2]));
	glUniform3f(glGetUniformLocation(shaderProgram, "CamLookAt"), float(camera.At()[0]), float(camera.At()[1]), float(camera.At()[2]));
	glUniform3f(glGetUniformLocation(shaderProgram, "CamUp"), float(camera.Up()[0]), float(camera.Up()[1]), float(camera.Up()[2]));
	glUniform1f(glGetUniformLocation(shaderProgram, "CamAngleOfViewV"), cameraAngleOfViewV);
	glUniform2f(glGetUniformLocation(shaderProgram, "iResolution"), float(width()), float(height()));

	if (hf)
	{
		// Uniforms - Heightfield
		Box2 box = hf->GetBox();
		glUniform2f(glGetUniformLocation(shaderProgram, "a"), float(box[0][0]), float(box[0][1]));
		glUniform2f(glGetUniformLocation(shaderProgram, "b"), float(box[1][0]), float(box[1][1]));
		glUniform2f(glGetUniformLocation(shaderProgram, "zRange"), zMin, zMax);
		glUniform1f(glGetUniformLocation(shaderProgram, "K"), K);
		glUniform2i(glGetUniformLocation(shaderProgram, "texSize"), hf->GetSizeX(), hf->GetSizeY());
		glUniform1i(glGetUniformLocation(shaderProgram, "shadingMode"), shadingMode);
	}

	// Draw heightfield
	glBindVertexArray(raytraceVAO);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	glBindVertexArray(0);
	glUseProgram(0);
}

/*!
\brief Custom event callback for mouse scrolling.
\param w raw window pointer
\param x, y scroll values
*/
void TerrainRaytracingWidget::ScrollCallback(GLFWwindow* w, double x, double y)
{
	// Don't move the camera is the mouse is over the GUI
	if (parent->MouseOverGUI())
		return;
	const double MoveScale = Norm(camera.View()) * 0.025;
	if (y > 0)
		camera.BackForth(MoveScale);
	else
		camera.BackForth(-MoveScale);
}

/*!
\brief Internal update functions for camera movements.
*/
void TerrainRaytracingWidget::Update()
{
	// Don't move the camera is the mouse is over the GUI
	Vector2 mousePos = parent->GetMousePosition();
	if (parent->MouseOverGUI())
	{
		x0 = mousePos[0];
		y0 = mousePos[1];
		return;
	}

	GLFWwindow* windowPtr = parent->getPointer();
	const double MoveScale = Norm(camera.View()) * 0.015 * 0.05;
	if (parent->GetMousePressed(GLFW_MOUSE_BUTTON_LEFT) && !parent->GetKey(GLFW_KEY_LEFT_CONTROL))
	{
		camera.LeftRightRound((x0 - mousePos[0]) * 0.01);
		camera.UpDownRound((y0 - mousePos[1]) * 0.005);
	}
	if (parent->GetMousePressed(GLFW_MOUSE_BUTTON_MIDDLE) && !parent->GetKey(GLFW_KEY_LEFT_CONTROL))
	{
		camera.LeftRightHorizontal((mousePos[0] - x0) * MoveScale);
		camera.UpDownVertical((mousePos[1] - y0) * MoveScale);
	}
	x0 = mousePos[0];
	y0 = mousePos[1];
}

/*!
\brief Internal update function of the widget.
Updates the internal height texture buffer from the heightfield pointer.
*/
void TerrainRaytracingWidget::SetHeightField(ScalarField2* hfPtr)
{
	hf = hfPtr;
	tmpData.resize(hf->VertexSize());
	UpdateInternal();
}

/*!
\brief Update the internal texture buffer from the heightfield pointer.
Also recomputes the Lipschitz constant and min/max elevation range.
*/
void TerrainRaytracingWidget::UpdateInternal()
{
	// Min/Max elevation
	double zMinDouble, zMaxDouble;
	hf->GetRange(zMinDouble, zMaxDouble);
	zMin = float(zMinDouble);
	zMax = float(zMaxDouble);

	// This avoids problem with rendering a flat heightfield.
	if (zMin == zMax)
		zMax += 10.0f;

	// Global Lipschitz constant
	K = float(hf->K());
	if (K == 0.0f)
		K = 1.0f;

	// Texture buffer: single channel
	for (int i = 0; i < hf->GetSizeX(); i++) {
		for (int j = 0; j < hf->GetSizeY(); j++) {
			float unitZ = float(Math::LinearStep(hf->at(i, j), zMinDouble, zMaxDouble));
			tmpData[hf->VertexIndex(i, j)] = unitZ;
		}
	}
	glUseProgram(shaderProgram);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, hfTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, hf->GetSizeX(), hf->GetSizeY(), 0, GL_RED, GL_FLOAT, tmpData.data());
	glProgramUniform1i(shaderProgram, glGetUniformLocation(shaderProgram, "heightfield"), 0);
	glUseProgram(0);
}

/*!
\brief Send the albedo texture data to the GPU.
*/
void TerrainRaytracingWidget::SetAlbedo(const Texture2D& tex)
{
	glUseProgram(shaderProgram);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, albedoTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex.GetSizeX(), tex.GetSizeY(), 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.Data());
	glProgramUniform1i(shaderProgram, glGetUniformLocation(shaderProgram, "albedo"), 1);
	glUseProgram(0);
}

/*!
\brief Reload the shaders of the widget. Useful for realtime editing and fine tuning of the rendering.
*/
void TerrainRaytracingWidget::ReloadShaders()
{
	shaderProgram = read_program("./data/shaders/heightfield_raytrace_330.glsl");
}

/*!
\brief Get the width of the associated window.
*/
int TerrainRaytracingWidget::width() const
{
	return parent->width();
}

/*!
\brief Get the height of the associated window.
*/
int TerrainRaytracingWidget::height() const
{
	return parent->height();
}
