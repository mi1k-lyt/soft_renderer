#pragma once
#ifndef __CAMERA_H__
#define __CAMERA_H__
#include "head.h"
#include "math.h"

class Camera {
public:
	vec3 origin;
	vec3 look_at;
	vec3 up;
	vec3 U;
	vec3 V;
	vec3 N;

	float fov;
	float aspect;
	float half_height;
	float half_width;

	Camera() {}
	Camera(vec3& ori, vec3& look_at_, vec3& vup, float v_fov, float aspect_) {
		origin = ori;
		look_at = look_at_;
		up = vup;
		fov = v_fov;
		aspect = aspect_;
		update();
	}

	void update() {
		N = (look_at - origin).unit();
		U = cross(N, up).unit();
		V = cross(U, N);

		float theta = fov * M_PI / 180.;
		half_height = tan(theta / 2.);
		half_width = aspect * half_height;
	}
};

Matrix4x4 calViewMatrix(Camera& c) {//aspect = width/height;	
	Matrix4x4 view(0.);
	for (int i = 0; i < 3; ++i) {
		view.e[0][i] = c.U[i];
		view.e[1][i] = c.V[i];
		view.e[2][i] = c.N[i];
		view.e[i][3] = -c.origin[i];
	}
	view.e[3][3] = 1;
	return std::move(view);
}

Matrix4x4 calProjectionMatrix(Camera& c, float z_near, float z_far) {

	Matrix4x4 proj(0.f);

	proj.e[0][0] = 1. / c.half_width;
	proj.e[1][1] = 1. / c.half_height;
	proj.e[2][2] = 1. / (z_far - z_near);
	proj.e[2][3] = -z_near / (z_far - z_near);
	proj.e[3][2] = 1;

	return std::move(proj);

}


#endif // !__CAMERA_H__

