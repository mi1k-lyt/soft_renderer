#pragma once
#ifndef __LIGHT_H__
#define __LIGHT_H__
#include "head.h"
#include "math.h"

class PointLight {
public:
	vec4 pos;
	vec3 color;

	PointLight() {}
	PointLight(const vec3& p, const vec3& c) {
		pos = vec4(p, 1.f);
		color = c;
	}

};


inline vec3 reflect(const vec3& in, const vec3& n)
{
	vec3 normal = n.unit();
	return in - 2. * dot(in, normal) * normal;
}

bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted)
{
	vec3 unit_v = v.unit();
	float dt = dot(unit_v, n);
	float dis = 1. - ni_over_nt * ni_over_nt * (1. - dt * dt);
	if (dis > 0)
	{
		refracted = ni_over_nt * (unit_v - n * dt) - n * sqrt(dis);
		return true;
	}
	else
	{
		if (dt < 0)
		{
			refracted = v - 2. * dt * n;
			return true;
		}
		else
		{
			return false;
			//refracted = v + 2.*dt*n;
		}
		return true;
	}
}

float schlick(float cosine, float ref_idx) {
	float r0 = (1. - ref_idx) / (1. + ref_idx);
	r0 *= r0;
	return r0 + (1 - r0) * pow((1. - cosine), 5);
}

#endif // !__LIGHT_H__

