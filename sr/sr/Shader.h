#pragma once

#ifndef __SHADER_H__
#define __SHADER_H__
#include "head.h"
#include "Camera.h"
#include "Vertex.h"
#include "MyWindow.h"
#include "Light.h"
#include "Buffer.h"
#include "Texture.h"

// 裁剪CVV六个面
#define X_MINUS1 0
#define X_1 1
#define Y_MINUS1 2
#define Y_1 3
#define Z_0 4
#define Z_1 5
#define W_0 6

void drawPoint(ColorBuffer& buffer, int x, int y, const Color& c) {
	*(buffer.ptr + y * buffer.width + x) = c;
}

void drawLine(ColorBuffer& buffer, int x1, int y1, int x2, int y2, const Color& c) {
	int dx = abs(x2 - x1), dy = abs(y2 - y1);
	if (dx == 0 && dy == 0) {
		drawPoint(buffer, x1, y1, c);
	}
	else if (dx == 0) {
		int dir = (y2 - y1) / dy;
		for (int i = y1; i != y2; i += dir) {
			drawPoint(buffer, x1, i, c);
		}
	}
	else if (dy == 0) {
		int dir = ((x2 - x1) > 0) ? 1 : -1;
		for (int i = x1; i != x2; i += dir) {
			drawPoint(buffer, i, y1, c);
		}
		drawPoint(buffer, x1, y1, c);
		drawPoint(buffer, x2, y1, c);
	}
	else if (dx == dy) {
		int dir_x = (x2 - x1) / dx;
		int dir_y = (y2 - y1) / dy;
		for (int i = x1, j = y1; (i != x2 && j != y2); i += dir_x, j += dir_y) {
			drawPoint(buffer, i, j, c);
		}
	}
	else {
		float k = (float)(y2 - y1) / (float)(x2 - x1);
		int d = (k > 0) ? 1 : -1;
		if (abs(k) > 1) {
			int twoDx = dx * 2, twoDxMinusDy = (dx - dy) * 2;
			int p = twoDx - dy;
			int x, y;
			if (y1 > y2) {
				x = x2; y = y2; y2 = y1;
			}
			else {
				x = x1; y = y1;
			}
			drawPoint(buffer, x, y, c);
			while (y < y2) {
				y++;
				if (p < 0) {
					p += twoDx;
				}
				else {
					x += d;
					p += twoDxMinusDy;
				}
				drawPoint(buffer, x, y, c);
			}
		}
		else {
			int twoDy = dy * 2, twoDyMinusDx = (dy - dx) * 2;
			int p = twoDy - dx;
			int x, y;
			if (x1 > x2) {
				x = x2; y = y2; x2 = x1;
			}
			else {
				x = x1; y = y1;
			}
			drawPoint(buffer, x, y, c);
			//if (k < 0) { twoDy *= -1; }
			while (x < x2) {
				x++;
				if (p < 0) {
					p += twoDy;
				}
				else {
					y += d;
					p += twoDyMinusDx;
				}
				drawPoint(buffer, x, y, c);
			}
		}

	}
}

void drawElement(ColorBuffer& buffer, VertexData& v_array, std::vector<int>& indices, int mode) {
	if (mode == SR_POINT) {
		for (int i = 0; i < indices.size(); ++i) {
			//drawPoint()
		}
	}
	for (int i = 2; i < indices.size(); ++i) {

	}
}


class Frag {
public:

	vec4 position;
	vec4 normal;
	vec4 color;
	vec4 fragPosInWorld;
	vec4 tangentLightPos;
	vec4 tangentViewPos;
	vec2 texCoords;
	coord screenPos;

	Frag() {}

	Frag(const VertexData& va, const std::vector<coord>& sp, const int index, int test_mode) {
		position = va.positions[index];
		normal = va.normals[index];
		if (test_mode == TEST_COLOR) {
			color = va.colors[index];
		}
		if (test_mode == TEST_LIGHT) {
			fragPosInWorld = va.fragPosInWorlds[index];
			tangentLightPos = va.tangentLightPos[index];
			tangentViewPos = va.tangentViewPos[index];
		}
		texCoords = va.texCoords[index];
		screenPos = sp[index];
	}

};

struct bbox {
	int max_x, max_y, min_x, min_y;
};

bbox getBbox(Frag& v0, Frag& v1, Frag& v2) {
	bbox temp;
	temp.max_x = maxTri(v0.screenPos.x, v1.screenPos.x, v2.screenPos.x);
	temp.max_y = maxTri(v0.screenPos.y, v1.screenPos.y, v2.screenPos.y);
	temp.min_x = minTri(v0.screenPos.x, v1.screenPos.x, v2.screenPos.x);
	temp.min_y = minTri(v0.screenPos.y, v1.screenPos.y, v2.screenPos.y);
	return temp;
}


class Shader {
public:
	Matrix4x4 model, view, projection;
	Camera camera;
	float z_near, z_far;

	VertexData vertexData;
	std::vector<int> indices;
	std::vector<coord> screenPos;

	ColorBuffer* buffer;
	DepthBuffer* depthBuffer;
	std::vector<Texture*> textureList;
	std::vector<Texture*> diffuseMaps;
	std::vector<Texture*> specularMaps;
	std::vector<Texture*> normalMaps;
	std::vector<Texture*> heightMaps;

	std::vector<Frag> fragments;
	std::vector<PointLight> lightList;

	int prim_mode = SR_TRIANGLE;
	int test_mode = TEST_LIGHT;

	//Shader() {}
	Shader() {
		model = identityMatrix4x4();
	}


	void setModelMatrix(Matrix4x4& m) {
		model = m;
	}

	void setVPMatrix(Camera& c, float z_n, float z_f) {
		camera = c;
		view = calViewMatrix(c);
		z_near = z_n, z_far = z_f;
		projection = calProjectionMatrix(c, z_near, z_far);
	}

	void setBuffer(ColorBuffer& b) {
		buffer = &b;
	}

	void setDepthBuffer(DepthBuffer& b) {
		depthBuffer = &b;
	}

	void use(VertexData& v, std::vector<int>& i) {
		vertexData = v;
		indices = i;
		vertexProgram();
		clip_w(vertexData, indices);
		perspectiveMap(vertexData);
		clip(vertexData, indices);
		toScreen(vertexData);
		rasterizer(vertexData, indices);
		fragmentProgram();
		drawToBuffer();
	}


	void vertexProgram() {
		// 旋转矩阵不变，缩放矩阵求逆
		Matrix4x4 normalCorrect = model.inverse().tranpose();
		applyTrans(model, vertexData.positions);
		applyTrans(normalCorrect, vertexData.tangents);
		applyTrans(normalCorrect, vertexData.bitTangents);
		applyTrans(normalCorrect, vertexData.normals);
		vertexData.fragPosInWorlds = vertexData.positions;

		vertexData.tangentLightPos = std::vector<vec4>(vertexData.positions.size());
		vertexData.tangentViewPos = std::vector<vec4>(vertexData.positions.size());

		if (test_mode == TEST_LIGHT) {
			const unsigned int data_size = vertexData.positions.size();
			for (int i = 0; i < data_size; ++i) {
				vec4 normal_t = vertexData.normals[i];
				normal_t[3] = 0;
				Matrix4x4 TBN(vertexData.tangents[i], vertexData.bitTangents[i], normal_t.unit(), vec4(0, 0, 0, 1));
				TBN = TBN.tranpose();//求逆
				//所有FS中用到的全部转换到切线空间
				vertexData.tangentLightPos[i] = TBN.mul_vec(lightList[0].pos);
				vertexData.tangentViewPos[i] = TBN.mul_vec(vec4(camera.origin, 0.f));
				vertexData.fragPosInWorlds[i] = TBN.mul_vec(vertexData.fragPosInWorlds[i]);
			}
		}
		applyTrans(view, vertexData.positions);
		applyTrans(projection, vertexData.positions);
	}

	void backFaceClip() {
		if (prim_mode == SR_TRIANGLE) {
			const unsigned int size = indices.size();
			for (int i = 2; i < size; i += 3) {
				vec3 normal = vertexData.normals[indices[i]].toVec3();
				if (dot(normal, camera.N) > 0.f) {
					indices[i] = indices[i - 1] = indices[i - 2] = -1;
				}
			}
		}
		return;
	}

	void clip_w(VertexData& data, std::vector<int>& indices) {

		if (prim_mode == SR_TRIANGLE) {
			/*MVP变换之后z与w呈线性关系且z一般小于w（因为对z变换时减去了一个正的常数项）
			因此先对z=0的平面进行裁剪，再对w=0平面进行裁剪（理论上不是必须的）*/

			//clip 的过程中indices的size会改变（只增加）
			//考虑到新生成的三角形也是符合条件的，就不再进行裁剪判断，size就按照原来的来
			int size = indices.size();
			for (int i = 2; i < size; i += 3) {
				int v0 = indices[i - 2], v1 = indices[i - 1], v2 = indices[i];
				if (v2 == -1) { continue; }
				clip_by_plane(data, indices, i, Z_0);
			}
			size = indices.size();
			for (int i = 2; i < size; i += 3) {
				int v0 = indices[i - 2], v1 = indices[i - 1], v2 = indices[i];
				if (v2 == -1) { continue; }
				clip_by_plane(data, indices, i, W_0);
			}
		}
		return;
	}

	void perspectiveMap(VertexData& data) {
		const unsigned int data_size = data.positions.size();
		for (int i = 0; i < data_size; i++) {
			float w_inv = 1.f / data.positions[i][3];

			for (int j = 0; j < 2; ++j) {
				data.positions[i][j] *= w_inv;
			}
			//nonlinear depth
			data.positions[i][2] = (w_inv - 1. / z_near) / (1. / z_far - 1. / z_near);

			//perspective map
			data.positions[i][3] = w_inv;

			data.normals[i] *= w_inv;
			if (test_mode == TEST_COLOR) { data.colors[i] *= w_inv; }
			data.texCoords[i] *= w_inv;
			if (test_mode == TEST_LIGHT) {
				data.fragPosInWorlds[i] *= w_inv;
				data.tangentLightPos[i] *= w_inv;
				data.tangentViewPos[i] *= w_inv;
			}
		}
		return;
	}

	void clip(VertexData& data, std::vector<int>& indices) {
		if (prim_mode == SR_TRIANGLE) {
			//clip 的过程中indices的size会改变（只增加）
			//考虑到新生成的三角形也是符合条件的，就不再进行裁剪判断，size就按照原来的来
			for (int mode = X_MINUS1; mode < Z_0; ++mode) {
				const int size = indices.size();
				for (int i = 2; i < size; i += 3) {
					int v0 = indices[i - 2], v1 = indices[i - 1], v2 = indices[i];
					if (v2 == -1) { continue; }
					clip_by_plane(data, indices, i, mode);
				}
			}
		}
		return;
	}

	void toScreen(VertexData& data) {
		//to screen
		screenPos = std::vector<coord>(data.positions.size());
		int width = buffer->width - 1;
		int height = buffer->height - 1;
		const int size = data.positions.size();
		for (int i = 0; i < size; i++) {
			float x = (data.positions[i][0] + 1.f) / 2.f;
			float y = (data.positions[i][1] + 1.f) / 2.f;
			screenPos[i].x = int(x * width);
			screenPos[i].y = int(y * height);
		}
		return;

	}

	void rasterizer(VertexData& data, std::vector<int>& indices) {
		fragments.clear();
		if (prim_mode == SR_TRIANGLE) {
			const int size = indices.size();
			for (int i = 2; i < size; i += 3) {
				if (indices[i] == -1) { continue; }
				Frag v2(data, screenPos, indices[i], test_mode);
				Frag v1(data, screenPos, indices[i - 1], test_mode);
				Frag v0(data, screenPos, indices[i - 2], test_mode);

				vec2 l1 = vec2(v1.screenPos.x - v0.screenPos.x, v1.screenPos.y - v0.screenPos.y);
				vec2 l2 = vec2(v2.screenPos.x - v0.screenPos.x, v2.screenPos.y - v0.screenPos.y);

				float l1l1 = dot(l1, l1);
				float l1l2 = dot(l1, l2);
				float l2l2 = dot(l2, l2);
				float delta = l1l1 * l2l2 - l1l2 * l1l2;
				if (abs(delta) < 1e-6) { continue; }

				bbox bb = getBbox(v0, v1, v2);
				float u, v;
				for (int x = bb.min_x; x < bb.max_x; ++x) {
					for (int y = bb.min_y; y < bb.max_y; ++y) {
						if (x < 0 || y < 0 || x>(buffer->width - 1) || y>(buffer->height - 1)) { continue; }
						vec2 l0 = vec2(x - v0.screenPos.x, y - v0.screenPos.y);

						float l0l1 = dot(l0, l1);
						float l0l2 = dot(l0, l2);
						u = (l2l2 * l0l1 - l1l2 * l0l2) / delta;
						v = (l1l1 * l0l2 - l0l1 * l1l2) / delta;
						bool in = (u >= 0.f && v >= 0.f && (u + v) <= 1.f);
						if (!in) { continue; }

						Frag temp = interpolate(v0, v1, v2, u, v);

						float depth = temp.position.z();

						int offset = depthBuffer->width * y + x;
						if (depth > *(depthBuffer->ptr + offset)) {
							continue;
						}


						temp.screenPos.x = x;
						temp.screenPos.y = y;

						float w_inv = 1.f / temp.position.w();
						/*for (int j = 1; j < data.attributeCount; ++j) {
							(*(data.attribute[j]))[i] *= w_inv;
						}*/
						temp.normal *= w_inv;
						temp.texCoords *= w_inv;
						if (test_mode == TEST_COLOR) {
							temp.color *= w_inv;
						}
						if (test_mode == TEST_LIGHT) {
							temp.fragPosInWorld *= w_inv;
							temp.tangentLightPos *= w_inv;
							temp.tangentViewPos *= w_inv;
						}
						fragments.push_back(temp);
					}
				}
			}
		}
		return;
	}

	void fragmentProgram() {
		//fragments[i].color 里就是顶点插值得来的颜色
		if (test_mode == TEST_LIGHT) {
			const int size = fragments.size();
			for (int i = 0; i < size; ++i) {

				vec3 lightPos = fragments[i].tangentLightPos.toVec3();
				vec3 viewPos = fragments[i].tangentViewPos.toVec3();
				vec3 fragPos = fragments[i].fragPosInWorld.toVec3();
				vec3 viewDir = (viewPos - fragPos).unit();

				float numLayers = 10.f;
				float layerDepth = 1.f / numLayers;
				float curLayerDepth = 0.f;

				vec2 p(viewDir.x(), viewDir.y());
				p *= 0.1f;
				vec2 deltaTexCoords = p / numLayers;
				vec2 curTexCoords(fragments[i].texCoords);
				float curDepth = textureList[2]->sample(curTexCoords[0], curTexCoords[1], 0).r();

				while (curLayerDepth < curDepth) {
					curTexCoords -= deltaTexCoords;
					curDepth = textureList[2]->sample(curTexCoords[0], curTexCoords[1], 0).r();
					curLayerDepth += layerDepth;
				}

				float afterDepth = curDepth - layerDepth;

				vec2 prevTexCoords = curTexCoords + deltaTexCoords;
				float beforeDepth = textureList[2]->sample(prevTexCoords[0], prevTexCoords[1], 0).r();
				float weight = afterDepth / (afterDepth - beforeDepth);
				curTexCoords = prevTexCoords * weight + curTexCoords * (1.f - weight);


				vec3 normal = textureList[1]->sample(curTexCoords[0], curTexCoords[1], 0).toVec3();
				normal = (normal * 2.f - 1.f).unit();
				vec3 color = textureList[0]->sample(curTexCoords[0], curTexCoords[1], 0).toVec3();

				vec3 lightDir = (lightPos - fragPos).unit();
				float ambient = 0.05f;
				float diffuse = (std::max)(dot(normal, lightDir), 0.f);
				vec3 halfDir = (viewDir + lightDir).unit();
				float specular = pow((std::max)(dot(halfDir, normal), 0.f), 64);
				color *= (ambient + diffuse + specular) * vec3(0.8f);

				
				fragments[i].color = vec4(color, 1.f);
			}
		}
		else if (test_mode == TEST_TEXTURE) {
			const int size = fragments.size();
			for (int i = 0; i < size; ++i) {
				float u = fragments[i].texCoords.x();
				float v = fragments[i].texCoords.y();
				fragments[i].color = diffuseMaps[0]->sample(u, v, 0) * 1.2f;
				
			}
		}
		else if (test_mode == TEST_COLOR) {
			//do nothing
		}
	}

	void drawToBuffer() {
		const int size = fragments.size();
		for (int i = 0; i < size; ++i) {
			int x = fragments[i].screenPos.x;
			int y = fragments[i].screenPos.y;
			float depth = fragments[i].position.z();

			int offset = depthBuffer->width * y + x;
			if (depth < *(depthBuffer->ptr + offset)) {
				*(depthBuffer->ptr + offset) = depth;
				for (int j = 0; j < 3; ++j) {
					if (fragments[i].color[j] > 1.f) { fragments[i].color[j] = 1.f; }
				}
				*(buffer->ptr + offset) = { fragments[i].color.r() * 255, fragments[i].color.g() * 255, fragments[i].color.b() * 255 };
			}
		}
	}

private:
	Frag interpolate(Frag& v0, Frag& v1, Frag& v2, float u, float v) {
		Frag temp;
		float w = 1.f - u - v;
		temp.position = u * v1.position + v * v2.position + w * v0.position;
		temp.normal = u * v1.normal + v * v2.normal + w * v0.normal;
		if (test_mode == TEST_COLOR) {
			temp.color = u * v1.color + v * v2.color + w * v0.color;
		}
		temp.texCoords = u * v1.texCoords + v * v2.texCoords + w * v0.texCoords;
		if (test_mode == TEST_LIGHT) {
			temp.fragPosInWorld = u * v1.fragPosInWorld + v * v2.fragPosInWorld + w * v0.fragPosInWorld;
			temp.tangentLightPos = u * v1.tangentLightPos + v * v2.tangentLightPos + w * v0.tangentLightPos;
			temp.tangentViewPos = u * v1.tangentViewPos + v * v2.tangentViewPos + w * v0.tangentViewPos;
		}

		return std::move(temp);
	}

	int addVertex(VertexData& data, int cliped) {
		data.positions.push_back(data.positions[cliped]);
		if (test_mode == TEST_COLOR) {
			data.colors.push_back(data.colors[cliped]);
		}
		data.normals.push_back(data.normals[cliped]);
		data.texCoords.push_back(data.texCoords[cliped]);
		if (test_mode == TEST_LIGHT) {
			data.fragPosInWorlds.push_back(data.fragPosInWorlds[cliped]);
			data.tangentLightPos.push_back(data.tangentLightPos[cliped]);
			data.tangentViewPos.push_back(data.tangentViewPos[cliped]);
		}
		return data.positions.size() - 1;
	}

	void clip_along_coord(VertexData& va, int kept, int replaced, int mode) {
		//默认 v0 在平面内 v1 在平面外, 插值之后的结果直接覆盖v1
		int coord = mode / 2;
		float delta = va.positions[replaced][coord] - va.positions[kept][coord];
		float plane;
		if (mode == X_1 || mode == Y_1 || mode == Z_1) {
			plane = 1.f;
		}
		else if (mode == X_MINUS1 || mode == Y_MINUS1) {
			plane = -1.f;
		}
		else if (mode == Z_0) {
			plane = 0.f;
		}
		float t = (plane - va.positions[kept][coord]) / delta;
		va.positions[replaced] = va.positions[kept] + t * (va.positions[replaced] - va.positions[kept]);
		if (test_mode == TEST_COLOR) {
			va.colors[replaced] = va.colors[kept] + t * (va.colors[replaced] - va.colors[kept]);
		}
		va.normals[replaced] = va.normals[kept] + t * (va.normals[replaced] - va.normals[kept]);
		if (test_mode == TEST_LIGHT) {
			va.fragPosInWorlds[replaced] = va.fragPosInWorlds[kept] + t * (va.fragPosInWorlds[replaced] - va.fragPosInWorlds[kept]);
		}
		va.texCoords[replaced] = va.texCoords[kept] + t * (va.texCoords[replaced] - va.texCoords[kept]);

	}

	void clip_one_vertex(VertexData& data, std::vector<int>& indices, int index_cliped, int index_v1, int index_v2, int mode) {
		//create temp vertex
		int v_temp1 = addVertex(data, indices[index_cliped]);
		int v_temp2 = addVertex(data, indices[index_cliped]);

		//clip one vertex and rearrange triangles
		clip_along_coord(data, indices[index_v1], v_temp1, mode);
		clip_along_coord(data, indices[index_v2], v_temp2, mode);

		indices[index_cliped] = v_temp2;
		indices.push_back(v_temp2);
		indices.push_back(v_temp1);
		indices.push_back(indices[index_v1]);
	}

	void clip_by_plane(VertexData& data, std::vector<int>& indices, int index, int mode) {
		int v0 = indices[index - 2], v1 = indices[index - 1], v2 = indices[index];

		bool in_v0;
		bool in_v1;
		bool in_v2;
		int coord = mode / 2;
		float threshold;
		if (mode == X_1 || mode == Y_1 || mode == Z_1) {
			threshold = 1.f;
			in_v0 = data.positions[v0][coord] < threshold;
			in_v1 = data.positions[v1][coord] < threshold;
			in_v2 = data.positions[v2][coord] < threshold;
		}
		else if (mode == X_MINUS1 || mode == Y_MINUS1) {
			threshold = -1.f;
			in_v0 = data.positions[v0][coord] > threshold;
			in_v1 = data.positions[v1][coord] > threshold;
			in_v2 = data.positions[v2][coord] > threshold;
		}
		else if (mode == Z_0 || mode == W_0) {
			threshold = 0.f;
			in_v0 = data.positions[v0][coord] > threshold;
			in_v1 = data.positions[v1][coord] > threshold;
			in_v2 = data.positions[v2][coord] > threshold;
		}

		int result = in_v0 + in_v1 + in_v2;
		if (result == 0) {
			//clip all vertexs
			indices[index - 2] = indices[index - 1] = indices[index] = -1;
			return;
		}
		else if (result >= 3) {
			return;
		}
		else if (result == 1) {
			//only one vertex are qualtified
			if (in_v0) {
				int new_v1 = addVertex(data, v1);
				int new_v2 = addVertex(data, v2);
				clip_along_coord(data, v0, new_v1, mode);
				clip_along_coord(data, v0, new_v2, mode);
				indices[index - 1] = new_v1;
				indices[index] = new_v2;
			}
			else if (in_v1) {
				int new_v0 = addVertex(data, v0);
				int new_v2 = addVertex(data, v2);
				clip_along_coord(data, v1, new_v0, mode);
				clip_along_coord(data, v1, new_v2, mode);
				indices[index - 2] = new_v0;
				indices[index] = new_v2;

			}
			else {
				int new_v0 = addVertex(data, v0);
				int new_v1 = addVertex(data, v1);
				clip_along_coord(data, v2, new_v0, mode);
				clip_along_coord(data, v2, new_v1, mode);
				indices[index - 2] = new_v0;
				indices[index - 1] = new_v1;
			}
		}
		else {
			//only one vertex are not qualtified
			if (!in_v0) {
				clip_one_vertex(data, indices, index - 2, index - 1, index, mode);
			}
			else if (!in_v1) {
				clip_one_vertex(data, indices, index - 1, index - 2, index, mode);
			}
			else {
				clip_one_vertex(data, indices, index, index - 2, index - 1, mode);
			}
		}
	}


};



#endif // !__SHADER_H__

