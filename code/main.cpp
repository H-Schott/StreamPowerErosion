#include "window.h"
#include "terrainwidget.h"
#include "scalarfield2.h"
#include "texture.h"
#include <imgui.h>

static Window* window;
static TerrainRaytracingWidget* widget;
static ScalarField2 hf;
static Texture2D albedoTexture;
static int shadingMode;
static double brushRadius = 250.0;
static double brushStrength = 10.0;

/*!
\brief Compute the intersection between a plane and a ray.
The intersection depth is returned if intersection occurs.
\param ray The ray.
\param t Intersection depth.
*/
static bool PlaneIntersect(const Ray& ray, double& t)
{
	double epsilon = 1e-4f;
	double x = Vector::Z * ray.Direction();
	if ((x < epsilon) && (x > -epsilon))
		return false;
	double c = Vector::Null * Vector::Z;
	double y = c - (Vector::Z * ray.Origin());
	t = y / x;
	return true;
}

/*!
\brief Reset the camera around a given 3D box.
*/
static void ResetCamera()
{
	Box2 box = hf.GetBox();
	Vector2 v = 0.5 * box.Diagonal();
	Camera cam = Camera(box.Center().ToVector(0.0) - (2.0 * v).ToVector(-Norm(v)), box.Center().ToVector(0), Vector(0.0, 0.0, 1.0), 1.0, 1.0, 5.0, 10000);
	widget->SetCamera(cam);
}

/*!
\brief User interface for the application.
*/
static void GUI()
{
	// Menu bar
	static bool newScene = false;
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Save heightfield"))
			{
				hf.Save("saved.png");
			}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}

	ImGui::Begin("Menu");
	{
		// Hardcoded examples
		{
			ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Scenes");
			if (ImGui::Button("Example 1"))
			{
				hf = ScalarField2(Box2(Vector2::Null, 30000.0), "../data/heightfields/hfTest1.png", 0.0, 2500.0);
				widget->SetHeightField(&hf);
			}
			ImGui::SameLine();
			if (ImGui::Button("Example 2"))
			{
				hf = ScalarField2(Box2(Vector2::Null, 10000.0), "../data/heightfields/hfTest2.png", 0.0, 2500.0);
				widget->SetHeightField(&hf);
			}
			ImGui::Separator();
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
		}
		
		// Shading
		{
			ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Shading");
			ImGui::RadioButton("Normal", &shadingMode, 0);
			ImGui::RadioButton("Shaded", &shadingMode, 1);
			widget->SetShadingMode(shadingMode);
			ImGui::Separator();
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
		}

		// Actions
		{
			// Reset camera
			if (ImGui::Button("Reset Camera"))
				ResetCamera();

			// Brushes
			ImGui::InputDouble("Brush radius", &brushRadius);
			ImGui::InputDouble("Brush strength", &brushStrength);
			ImGui::Separator();
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
		}
		
		// Simulation statistics
		{
			ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Statistics");
			ImGui::Text("%.3f ms/frame (%.1f FPS)", 1000.0f / float(ImGui::GetIO().Framerate), float(ImGui::GetIO().Framerate));
		}
	}
	ImGui::End();

	ImGui::Begin("Help panel");
	{
		ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Controls");
		ImGui::TextWrapped("Camera: left click + mouse movements, and scrolling for zoom.");
		ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();

		ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Examples");
		ImGui::TextWrapped("Simple examples are provided. Simply click on the associated buttons.");
		ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();

		ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Editing brushes");
		ImGui::TextWrapped("Using ctrl + left/right click, you can add or remove elevation on the terrain. Brush size and strength can be modified using the left panel.");
		ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();

		ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Other functionalities");
		ImGui::TextWrapped("A new scene can be created using the File > New Scene menu.");
		ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();

		ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Issues");
		ImGui::TextWrapped("Contact me at axel(dot)paris69(at)gmail(dot)com, or report an issue on github.");
	}
	ImGui::End();
}

int main()
{
	// Init
	window = new Window("Arches example app", 1920, 1080);
	widget = new TerrainRaytracingWidget();
	window->SetWidget(widget);
	hf = ScalarField2(Box2(Vector2::Null, 1000), 256, 256, 0.0);
	widget->SetHeightField(&hf);
	window->SetUICallback(GUI);

	albedoTexture = Texture2D(hf.GetSizeX(), hf.GetSizeY());
	albedoTexture.Fill(Color8(225, 225, 225, 255));
	widget->SetAlbedo(albedoTexture);

	ResetCamera();

	// Main loop
	while (!window->Exit())
	{
		// Heightfield editing
		bool leftMouse = window->GetMousePressed(GLFW_MOUSE_BUTTON_LEFT);
		bool rightMouse = window->GetMousePressed(GLFW_MOUSE_BUTTON_RIGHT);
		bool mouseOverGUI = window->MouseOverGUI();
		if (!mouseOverGUI && (leftMouse || rightMouse) && window->GetKey(GLFW_KEY_LEFT_CONTROL))
		{
			Camera cam = widget->GetCamera();
			double xpos, ypos;
			glfwGetCursorPos(window->getPointer(), &xpos, &ypos);
			Ray ray = cam.PixelToRay(int(xpos), int(ypos), window->width(), window->height());
			double t;
			if (PlaneIntersect(ray, t))
			{
				hf.Gaussian(ray(t), brushRadius, leftMouse ? brushStrength : -brushStrength);
				widget->UpdateInternal();
			}
		}

		window->Update();
	}
	return 0;
}
