#include "window.h"
#include "terrainwidget.h"

#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

/*!
\brief Class for handling a window with GLFW3 library.
*/

/*!
\brief Constructor.
\param windowName name of the window.
\param w, h width and height of the window.
*/
Window::Window(const char* windowName, int w, int h)
{
	widget = nullptr;
	uiUserFunPtr = nullptr;

	// Window
	width_internal = w;
	height_internal = h;
	if (!glfwInit())
	{
		std::cout << "GLFW failed to initialize" << std::endl;
		glfwTerminate();
		return;
	}
	GLFWmonitor* monitor = glfwGetPrimaryMonitor();
	if (!monitor)
	{
		std::cout << "GLFW failed to get primary monitor" << std::endl;
		glfwTerminate();
		return;
	}
	const GLFWvidmode* mode = glfwGetVideoMode(monitor);
	width_internal = width_internal <= 0 ? mode->width : width_internal;
	height_internal = height_internal <= 0 ? mode->height : height_internal;

	glfwDefaultWindowHints();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
	glfwWindowHint(GLFW_SAMPLES, 4);
	windowPtr = glfwCreateWindow(width_internal, height_internal, windowName, NULL, NULL);
	if (windowPtr == NULL)
	{
		std::cout << "GLFW failed to create window" << std::endl;
		glfwTerminate();
		return;
	}
	glfwMakeContextCurrent(windowPtr);
	glfwSetInputMode(windowPtr, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

	// OpenGL
	glewInit();
	glEnable(GL_DEPTH_TEST);
	GLenum err = glGetError();
	if (err != GL_NO_ERROR)
	{
		std::cout << "OpenGL failed to initialize" << std::endl;
		glfwTerminate();
		return;
	}

	// Dear ImGui
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(windowPtr, true);
	ImGui_ImplOpenGL3_Init("#version 330");

	glfwGetFramebufferSize(windowPtr, &width_internal, &height_internal);
	glViewport(0, 0, width_internal, height_internal);

	// vsync
	glfwSwapInterval(1);

	// TODO(Axel): log this into the console.
	//spdlog::info("All systems were correctly initialized");
	//spdlog::info("OpenGL device information: Vendor: {}", (const char*)glGetString(GL_VENDOR));
	//spdlog::info("OpenGL device information: Renderer: {}", (const char*)glGetString(GL_RENDERER));
	//spdlog::info("OpenGL device information: Version: {}", (const char*)glGetString(GL_VERSION));
	//spdlog::info("OpenGL device information: GLSL: {}", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
	//spdlog::info("Framebuffer size: {}x{}", width_internal, height_internal);
	//spdlog::info("Dear ImGui: {}", ImGui::GetVersion());
}

/*!
\brief Destructor. Releases Imgui, GLFW, and destroys the internal rendering widget.
*/
Window::~Window()
{
	delete widget;
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	glfwTerminate();
}

/*!
\brief Set the rendering widget.
\param w widget
*/
void Window::SetWidget(TerrainRaytracingWidget* w)
{
	widget = w;
	widget->SetWindowPtr(this);
	widget->initializeGL();
	glfwSetWindowUserPointer(windowPtr, widget);
	glfwSetScrollCallback(windowPtr, [](GLFWwindow* win, double x, double y)
	{
		TerrainRaytracingWidget* ptr = (TerrainRaytracingWidget*) glfwGetWindowUserPointer(win);
		ptr->ScrollCallback(win, x, y);
	});
}

/*!
\brief Set the custom UI user callback, which is called in the Update function.
\param funPtr user function pointer
*/
void Window::SetUICallback(void (*funPtr)())
{
	uiUserFunPtr = funPtr;
}

/*!
\brief Returns true if the mouse is currently over a GUI item, false otherwise.
*/
bool Window::MouseOverGUI() const
{
	if (!ImGui::GetCurrentContext())
		return false;
	return ImGui::IsWindowFocused(ImGuiFocusedFlags_AnyWindow)
		|| ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)
		|| ImGui::IsAnyItemHovered();
}

/*!
\brief Base update function for each frame.
*/
void Window::Update()
{
	// Widget rendering
	if (widget)
	{
		widget->Update();
		widget->paintGL();
	}

	// UI Rendering
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
	{
		if (uiUserFunPtr)
			uiUserFunPtr();
	}
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	// Swap buffers
	glfwSwapBuffers(windowPtr);
	glfwPollEvents();
}

/*
\brief Check if a given key has been pressed in the last frame.
\param key the key
*/
bool Window::GetKey(int key) const
{
	return bool(glfwGetKey(windowPtr, key));
}

/*!
\brief Check if a given mouse button has been pressed in the last frame.
\param mouse mouse button
*/
bool Window::GetMousePressed(int mouse) const
{
	return bool(glfwGetMouseButton(windowPtr, mouse));
}

/*!
\brief Get the current mouse position.
*/
Vector2 Window::GetMousePosition() const
{
	Vector2 ret;
	glfwGetCursorPos(windowPtr, &ret[0], &ret[1]);
	return ret;
}

/*!
\brief Function for knowing if the program has finished (user exited or pressed escape).
*/
bool Window::Exit() const
{
	return glfwWindowShouldClose(windowPtr) || glfwGetKey(windowPtr, GLFW_KEY_ESCAPE);
}

/*!
\brief Get the width of the window.
*/
int Window::width() const
{
	return width_internal;
}

/*!
\brief Get the height of the window.
*/
int Window::height() const
{
	return height_internal;
}

/*!
\brief Get the raw window pointer.
*/
GLFWwindow* Window::getPointer()
{
	return windowPtr;
}
