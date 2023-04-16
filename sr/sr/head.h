#pragma once

#include <iostream>
#include <vector>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.142592657

enum {
	SR_TRIANGLE,
	SR_LINE,
	SR_POINT
};

enum {
	TEST_COLOR,
	TEST_LIGHT,
	TEST_TEXTURE
};

struct coord {
	int x, y;
};


struct Color {
	float R, G, B;
};
