#include "window.h"
#include "terrainwidget.h"
#include "scalarfield2.h"
#include "spe_shader.h"
#include "texture.h"
#include <imgui.h>

static Window* window;
static TerrainRaytracingWidget* widget;
static ScalarField2 hf;
static ScalarField2 uplift;
static ScalarField2 gpu_drainage;
static GPU_SPE gpu_spe;
static Texture2D albedoTexture;
static int shadingMode;
static float brushRadius = 30;
static bool brushRadius_changed = false;
static float brushStrength = 10.0;
static bool brushStrength_changed = false;
static bool ongoing_gpu_spe = false;
static float delta_time = 100;
static bool delta_time_changed = false;

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
			ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Premade Uplifts");
			if (ImGui::Button("Uplift 0")) {
				uplift = ScalarField2(Box2(Vector2::Null, 150 * 1000), 256, 256, 0.2);
				gpu_spe.SetUplift(uplift);
			}
			ImGui::SameLine();
			if (ImGui::Button("Uplift 1")) {
				uplift = ScalarField2(Box2(Vector2::Null, 150*1000), "../data/uplifts/lambda.png", 0.4, 10.0);
				gpu_spe.SetUplift(uplift);
			}
			ImGui::SameLine();
			if (ImGui::Button("Uplift 2")) {
				uplift = ScalarField2(Box2(Vector2::Null, 150*1000), "../data/uplifts/alpes_noise.png", 0.4, 10.0);
				gpu_spe.SetUplift(uplift);
			}
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
			ImGui::Separator();
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
		}
		
		// Shading
		{
			ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Shading");
			ImGui::RadioButton("Normal", &shadingMode, 0);
			ImGui::RadioButton("Shaded", &shadingMode, 1);
			ImGui::RadioButton("Uplift", &shadingMode, 2);
			if (shadingMode == 1) {
				Texture2D texture = Texture2D(hf.GetSizeX(), hf.GetSizeY());
				texture.Fill(Color8(225, 225, 225, 255));
				widget->SetAlbedo(texture);
				widget->SetShadingMode(1);
			}
			else if (shadingMode == 2) {
				Texture2D texture = uplift.CreateImage();
				widget->SetAlbedo(texture);
				widget->SetShadingMode(1);
			} else widget->SetShadingMode(shadingMode);
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
			ImGui::Separator();
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
		}

		// Actions
		{
			// Reset camera
			if (ImGui::Button("Reset Camera"))
				ResetCamera();

			// Brushes
			ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "CTRL + Left Click to draw mountains");
			brushRadius, brushRadius_changed = ImGui::SliderFloat("radius", &brushRadius, 10, 100);
			brushStrength, brushStrength_changed = ImGui::SliderFloat("strength", &brushStrength, 1, 30);
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
			ImGui::Separator();
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
		}

		{
			ImGui::Checkbox("Ongoing simulation", &ongoing_gpu_spe);
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
			ImGui::Separator();
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
			delta_time, delta_time_changed = ImGui::SliderFloat("dt", &delta_time, 1, 100);
			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
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

	/*ImGui::Begin("Help panel");
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
	ImGui::End();*/
}

int main()
{
	// Init
	window = new Window("Stream Power Erosion", 1920, 1080);
	widget = new TerrainRaytracingWidget();
	window->SetWidget(widget);
	hf = ScalarField2(Box2(Vector2::Null, 150*1000), "../data/heightfields/hfTest2.png", 0.0, 1000.0);
	uplift = ScalarField2(Box2(Vector2::Null, 150*1000), 256, 256, 1.0);
	widget->SetHeightField(&hf);
	window->SetUICallback(GUI);

	// gpu_spe init
	gpu_spe.Init(hf);
	gpu_spe.SetUplift(uplift);
	gpu_drainage = hf;

	albedoTexture = Texture2D(hf.GetSizeX(), hf.GetSizeY());
	albedoTexture.Fill(Color8(225, 225, 225, 255));
	widget->SetAlbedo(albedoTexture);

	ResetCamera();

	// Main loop
	while (!window->Exit()) {
		// Heightfield editing
		bool leftMouse = window->GetMousePressed(GLFW_MOUSE_BUTTON_LEFT);
		bool rightMouse = window->GetMousePressed(GLFW_MOUSE_BUTTON_RIGHT);
		bool mouseOverGUI = window->MouseOverGUI();
		if (!mouseOverGUI && (leftMouse) && window->GetKey(GLFW_KEY_LEFT_CONTROL)) {
			Camera cam = widget->GetCamera();
			double xpos, ypos;
			glfwGetCursorPos(window->getPointer(), &xpos, &ypos);
			Ray ray = cam.PixelToRay(int(xpos), int(ypos), window->width(), window->height());
			double t;
			if (PlaneIntersect(ray, t)) {
				uplift.Gaussian(ray(t), brushRadius * 1000, brushStrength);
				gpu_spe.SetUplift(uplift);
			}
		}
		if (ongoing_gpu_spe) {
			// parameters changes
			if (delta_time_changed) gpu_spe.SetDt(delta_time);

			// simulation step
			gpu_spe.Step(200);
			gpu_spe.GetData(hf);
			widget->UpdateInternal();
		}

		window->Update();
	}
	return 0;
}
