#pragma once

#include <GL/glew.h>
#include "GLFW/glfw3.h"
#include "evector.h"

class TerrainRaytracingWidget;

class Window
{
protected:
	GLFWwindow* windowPtr;
	int width_internal, height_internal;
	TerrainRaytracingWidget* widget;
	void (*uiUserFunPtr)();

public:
	Window(const char* windowName, int w, int h);
	~Window();

	void SetWidget(TerrainRaytracingWidget* w);
	bool Exit() const;
	void Update();
	void SetUICallback(void (*funPtr)());
	bool GetKey(int key) const;
	bool GetMousePressed(int mouse) const;
	Vector2 GetMousePosition() const;
	bool MouseOverGUI() const;

	int width() const;
	int height() const;
	GLFWwindow* getPointer();
};
