#define NOMINMAX
#include "head.h"
#include "Shader.h"
#include "MyWindow.h"
#include "Vertex.h"
#include "Model.h"

int WINAPI WinMain(HINSTANCE hinstance, HINSTANCE hprevinstance, LPSTR lpcmdline, int ncmdshow)
{
	constexpr int width = 640;
	constexpr int height = 480;
	
	MSG msg;
	MyWindow mwin(hinstance, width, height);
	ShowWindow(mwin.hwnd, SW_SHOW);
	UpdateWindow(mwin.hwnd);

	Shader shader;

	ColorBuffer buffer(width, height);
	buffer.clear(1.f);
	shader.setBuffer(buffer);

	DepthBuffer depthBuffer(width, height);
	depthBuffer.clear(2.f);
	shader.setDepthBuffer(depthBuffer);

	Matrix4x4 model = identityMatrix4x4();
	matScale(model, vec3(3.5f));
	matRotate(model, 135.f, vec3(0, 1, 0));
	matTranslate(model, vec3(-1, -0.8, -105));
	matScale(model, vec3(0.035f));
	shader.setModelMatrix(model);

	Camera camera;
	float z_near = 0.1f, z_far = 100.f;

	vec3 ori(0, 0.0, -4.3), look_at(0, 0, 0), vup(0, 1, 0);
	float v_fov = 60., aspect = float(width) / float(height);
	camera = Camera(ori, look_at, vup, v_fov, aspect);
	shader.setVPMatrix(camera, z_near, z_far);

	PointLight defualtLight1(vec3(-1.8f, 0.f, 0.2f), 0.7 * vec3(0.8f, 0.65f, 0.4f));;
	shader.lightList.push_back(defualtLight1);

	shader.test_mode = TEST_TEXTURE;

	Model myModel("E:/self_program/soft_renderer/sr/sr/resource/backpack/backpack.obj");
	myModel.draw(shader);


	bool quit = false;
	while (!quit)
	{
		if (PeekMessage(&msg, NULL, NULL, NULL, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
				quit = true;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		//buffer.clear(0.f);
		//depthBuffer.clear(2.f);
		//shader.setVPMatrix(camera, z_near, z_far);
		//shader.use(data, indices);//draw to buffer
		mwin.draw(buffer);
	}

	mwin.release();

	return 0;
}