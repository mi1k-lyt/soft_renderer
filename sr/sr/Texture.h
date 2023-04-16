#pragma once
#ifndef __TEXTURE_H__
#define __TEXTURE_H__
#include "head.h"
#include "math.h"
class Texture {
public:
	int width, height, nrChannels;
	std::string fileName;
	std::vector<vec4> _data;

	Texture() {}
	Texture(const char* path) {
		getTexture(path);
		fileName = path;
	}
	
	vec4 sample(float u, float v, int mode) {
		unsigned int x, y;
		while (u < 0.f) { u += 1.f; }
		while (u > 1.f) { u -= 1.f; }
		while (v < 0.f) { v += 1.f; }
		while (v > 1.f) { v -= 1.f; }
		u *= (width - 1);
		v *= (height - 1);
		float next_u = u - floor(u), next_v = v - floor(v);
		float current_u = 1.f - next_u, current_v = 1.f - next_v;

		x = static_cast<int>(u);
		y = static_cast<int>(v);
		if (x == width - 1 || y == height - 1) {
			return (_data)[y * width + x];
		}
		else {
			return current_u * current_v * (_data)[y * width + x]
				+ current_u * next_v * (_data)[(y + 1) * width + x]
				+ next_u * current_v * (_data)[y * width + x + 1]
				+ next_u * next_v * (_data)[(y + 1) * width + x + 1];
		}
	}

private:
	void getTexture(const char* name) {
		stbi_set_flip_vertically_on_load(true);

		unsigned char* image = NULL;
		image = stbi_load(name, &width, &height, &nrChannels, 0);

		if (image == NULL) {
			std::cout << name << "\nError: Texture load failed!" << std::endl;
			return;
		}

		
		_data = std::vector<vec4>(width * height);


		unsigned char* char_input;
		int width_temp = width * 4;
		for (int x = 0; x < width; ++x) {
			for (int y = 0; y < height; ++y) {
				int offset = y * width + x;
				char_input = image + nrChannels * offset;
				for (int c = 0; c < nrChannels; ++c) {
					(_data)[offset][c] = static_cast<float>(*(char_input + c)) / 255.f;
				}
				if (nrChannels < 4) {
					(_data)[offset][3] = 1.f;
				}
			}
		}

		stbi_image_free(image);
	}
};
#endif // !__TEXTURE_H__

