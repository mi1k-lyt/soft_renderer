#pragma once

#ifndef __VERTEX_H__
#define __VERTEX_H__
#include "head.h"
#include <string>
#include <iostream>
#include "math.h"

class VertexData {
public:

	std::vector<vec4> positions;
	std::vector<vec4> normals;
	std::vector<vec4> colors;
	std::vector<vec4> fragPosInWorlds;
	std::vector<vec4> tangents;
	std::vector<vec4> bitTangents;
	std::vector<vec4> tangentLightPos;
	std::vector<vec4> tangentViewPos;
	std::vector<vec2> texCoords;


	std::vector<vec4>* attribute[100];
	std::string attributeName[100];
	int attributeCount = 0;

	VertexData() {
		attributeCount = 0;
	}
	VertexData(const float* vertexs, int dataAmount, int dataLength, int offset, int pointLength) {
		attributeCount = 0;
		setData(positions, vertexs, dataAmount, dataLength, offset, pointLength, 4);
	}

	void setNormal(const float* vertexs, int dataAmount, int dataLength, int offset, int pointLength) {
		setData(normals, vertexs, dataAmount, dataLength, offset, pointLength, 4);
	}

	void setTexCoords(const float* vertexs, int dataAmount, int dataLength, int offset, int pointLength) {
		setData(texCoords, vertexs, dataAmount, dataLength, offset, pointLength, 2);
	}
	void setColor(const float* vertexs, int dataAmount, int dataLength, int offset, int pointLength) {
		setData(colors, vertexs, dataAmount, dataLength, offset, pointLength, 3);
	}

	void setTBN(const std::vector<int>& indices) {
		const int size = positions.size();
		tangents = std::vector<vec4>(size);
		bitTangents = std::vector<vec4>(size);

		const int idx_size = indices.size();
		for (int i = 2; i < idx_size; i += 3) {
			vec4 v0 = positions[indices[i - 2]],
				v1 = positions[indices[i - 1]],
				v2 = positions[indices[i]];
			vec3 edge1 = (v1 - v0).toVec3(), edge2 = (v2 - v0).toVec3();
			vec2 deltaUV1 = texCoords[i - 1] - texCoords[i - 2], deltaUV2 = texCoords[i] - texCoords[i - 2];
			vec4 tang, bitTang;
			float delta = deltaUV1.x() * deltaUV2.y() - deltaUV2.x() * deltaUV1.y();
			if (abs(delta) < 1e-6) {
				std::cout << "Calculate TBN failed." << std::endl;
				return;
			}
			else {
				delta = 1.f / delta;
				tang = vec4(delta * (deltaUV2.y() * edge1 - deltaUV1.y() * edge2), 0.f).unit();
				bitTang = vec4(delta * (deltaUV1.x() * edge2 - deltaUV2.x() * edge1), 0.f).unit();
				for (int j = 0; j < 3; ++j) {
					tangents[indices[i - j]] = tang;
					bitTangents[indices[i - j]] = bitTang;
				}
			}
		}
		return;
	}

private:

	template <typename T>
	void setData(std::vector<T>& dataContainer, const float* vertexs, int dataAmount, int dataLength, int offset, int pointLength, int selfDataLength) {
		dataContainer = std::vector<T>(dataAmount);
		for (int i = 0; i < dataAmount; ++i) {
			for (int j = 0; j < pointLength; ++j) {
				dataContainer[i][j] = vertexs[i * dataLength + offset + j];
			}
			for (int j = pointLength; j < dataContainer[i].size(); ++j) {
				dataContainer[i][j] = 1.0f;
			}
		}
	}
};

#endif // !__VERTEX_H__

