#pragma once

#ifndef __MYWINDOW_H__
#define __MYWINDOW_H__
#include <windows.h>
#include <WinUser.h>
#include "Buffer.h"
#define WNDCLASSNAME L"WINCLASS"
LRESULT CALLBACK WinProc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{

	switch (msg)
	{
	case WM_CLOSE: {
		PostQuitMessage(0);
		break;
	}

	}
	return DefWindowProc(hwnd, msg, wparam, lparam);
}

class MyWindow {
public:
	HDC hdc;
	HWND hwnd;
	WNDCLASSEX wndclass;

	int width;
	int height;

	MyWindow(HINSTANCE& hinstance, int w, int h) {
		width = w;
		height = h;

		wndclass.cbSize = sizeof(WNDCLASSEX);
		wndclass.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
		wndclass.lpfnWndProc = WinProc;
		wndclass.cbClsExtra = 0;
		wndclass.cbWndExtra = 0;
		wndclass.hInstance = hinstance;
		wndclass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
		wndclass.hCursor = LoadCursor(NULL, IDC_ARROW);
		wndclass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
		wndclass.lpszMenuName = NULL;
		wndclass.lpszClassName = WNDCLASSNAME;
		wndclass.hIconSm = NULL;

		RegisterClassEx(&wndclass);

		// Create the window 

		hwnd = CreateWindowEx(NULL,
			WNDCLASSNAME,
			L"sr",
			WS_OVERLAPPEDWINDOW,
			0, 0,
			w, h,
			NULL,
			NULL,
			hinstance,
			NULL);

		hdc = GetDC(hwnd);
	}

	void release() {
		ReleaseDC(hwnd, hdc);
	}

	void draw(const ColorBuffer& buffer) {
		COLORREF* arr = new COLORREF[buffer.width * buffer.height];
		Color* c;
		for (int y = 0; y < buffer.height; ++y) {
			int offset_y = (buffer.height - 1 - y) * buffer.width;
			int offset_y1 = y * buffer.width;
			for (int x = 0; x < buffer.width; ++x) {
				c = buffer.ptr + offset_y + x;
				if (x == 799) {
					int a = 0;
				}
				*(arr + offset_y1 + x) = RGB(int(c->B), int(c->G), int(c->R));
			}

		}
		HBITMAP map = CreateBitmap(buffer.width, buffer.height, 1, 4 * 8, (void*)arr);
		HDC src = CreateCompatibleDC(hdc);
		SelectObject(src, map);
		BitBlt(hdc, 0, 0, buffer.width, buffer.height, src, 0, 0, SRCCOPY);
		DeleteObject(map);
		DeleteDC(src);
		delete[] arr;
	}
};

#endif // !MYWINDOW

