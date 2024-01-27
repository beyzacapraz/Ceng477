#ifndef SCENE_H
#define SCENE_H
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Matrix4.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double> > depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);
	void deepcopy(Matrix4 &m1, Matrix4 m2);
	void deepcopyVec4(Vec4 &v1, Vec4 v2);
	void deepcopyVec3(Vec3 &v1, Vec3 v2);
	void translation(Translation t, Matrix4 & m_t);
	void rotation(Rotation r, Matrix4 & m_r);
	void scaling(Scaling s, Matrix4 & m_s);
	Matrix4 modelingTransformation(Mesh *mesh);
	Matrix4 cameraTransformation(Camera* cam);
	void Transformation_Rasterization(Camera* cam, Mesh* mesh);
	Matrix4 projectionTransformation(Camera*cam);
	bool clipping(Vec4 & vertex1, Vec4 &  vertex2, Color &c0, Color &c1);
	bool visible(double den, double num, double & t_e, double & t_l);
	Matrix4 viewportTransformation(Camera * cam);
	void lineRasterization(Vec4 vertex1, Vec4 vertex2, Color c1, Color c2, Camera *cam);
	Vec4 vec3tovec4(Vec3 v);
	Vec3 vec4tovec3(Vec4 v);
	void perspectiveDivide(Vec4& v);
	void triangleRasterization(Vec4 vertex1, Vec4 vertex2, Vec4 vertex3, int index1, int index2, int index3, Camera *cam);
};

#endif