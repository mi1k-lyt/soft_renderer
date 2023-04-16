#pragma once

#ifndef __BUFFER_H__
#define __BUFFER_H__
#include "head.h"

class ColorBuffer {
public:
	Color* ptr;
	int width, height;
	Color operator[](int i) const { return ptr[i]; }
	Color& operator[](int i) { return ptr[i]; }

	ColorBuffer(int _width, int _height) {
		ptr = NULL;
		try {
			ptr = new Color[_width * _height];
		}
		catch (std::bad_alloc) {
			std::cout << "Bad alloc: memory alloc failed!" << std::endl;
			ptr = NULL;
		}
		if (ptr) {
			memset(ptr, 0, sizeof(Color) * _width * _height);
			width = _width;
			height = _height;
		}
	}

	void resize(int _width, int _height) {
		Color* temp = ptr;
		try {
			ptr = new Color[_width * _height];
		}
		catch (std::bad_alloc) {
			std::cout << "Bad alloc: resize failed!" << std::endl;
			ptr = temp;
			return;
		}
		if (!ptr) {
			ptr = temp;
		}
		else {
			memset(ptr, 0, sizeof(Color) * _width * _height);
			width = _width;
			height = _height;
			delete[] temp;
		}
	}

	~ColorBuffer() {
		memset(ptr, 0, sizeof(Color) * width * height);
		delete[] ptr;
		ptr = NULL;
	}

	void clear(float v) {
		int* t = (int*)&v;
		memset(ptr, *t, sizeof(Color) * width * height);
	}
};

class DepthBuffer {
public:
	float* ptr;//Depth in float
	int width, height;
	float operator[](int i) const { return ptr[i]; }
	float& operator[](int i) { return ptr[i]; }

	DepthBuffer() { ptr = NULL; width = height = -1; }
	DepthBuffer(int _width, int _height) {
		ptr = NULL;
		try {
			ptr = new float[_width * _height];
		}
		catch (std::bad_alloc) {
			std::cout << "Bad alloc: memory alloc failed!" << std::endl;
			ptr = NULL;

		}
		if (ptr) {
			std::fill(ptr, ptr + _width * _height, -2.f);
			//memset(ptr, 0, sizeof(float)*_width*_height);
			width = _width;
			height = _height;
		}
	}

	void resize(int _width, int _height) {
		float* temp = ptr;
		try {
			ptr = new float[_width * _height];
		}
		catch (std::bad_alloc) {
			std::cout << "Bad alloc: resize failed!" << std::endl;
			ptr = temp;
			return;
		}
		if (!ptr) {
			ptr = temp;
		}
		else {
			std::fill(ptr, ptr + _width * _height, -2.f);
			width = _width;
			height = _height;
			delete[] temp;
		}
	}

	~DepthBuffer() {
		memset(ptr, 0, sizeof(float) * width * height);
		delete[] ptr;
		ptr = NULL;
	}

	void clear(float v) {
		std::fill(ptr, ptr + width * height, v);
	}
};


#endif // !__BUFFER_H__

