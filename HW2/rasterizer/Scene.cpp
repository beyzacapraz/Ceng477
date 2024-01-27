#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::deepcopy(Matrix4 &m1, Matrix4 m2){
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            m1.values[i][j] = m2.values[i][j];
        }
    }
}

void Scene::deepcopyVec4(Vec4 &v1, Vec4 v2){
    v1.colorId = v2.colorId;
	v1.x = v2.x;
	v1.y = v2.y;
	v1.z = v2.z;
	v1.t = v2.t;
}

void Scene::deepcopyVec3(Vec3 &v1, Vec3 v2){
    v1.colorId = v2.colorId;
	v1.x = v2.x;
	v1.y = v2.y;
	v1.z = v2.z;
}


void Scene::translation(Translation t, Matrix4 & m){
	Matrix4 m_t;
	deepcopy(m_t, getIdentityMatrix());
	m_t.values[0][3] = t.tx;
	m_t.values[1][3] = t.ty;
	m_t.values[2][3] = t.tz;

	deepcopy(m, multiplyMatrixWithMatrix(m_t, m));
}

void Scene::scaling(Scaling s, Matrix4 & m){
	Matrix4 m_s;
	deepcopy(m_s, getIdentityMatrix());
	m_s.values[0][0] = s.sx;
	m_s.values[1][1] = s.sy;
	m_s.values[2][2] = s.sz;

	deepcopy(m, multiplyMatrixWithMatrix(m_s, m));
}

void Scene::rotation(Rotation r, Matrix4 & m){
	Vec3 u(r.ux, r.uy, r.uz);
	deepcopyVec3(u, normalizeVec3(u));
	double x = std::abs(u.x);
	double y = std::abs(u.y);
	double z = std::abs(u.z);
	double minn = std::min(x,std::min(y,z));
	Vec3 v, w;

	if (minn == x){
		v.x = 0;
		v.y = -u.z;
		v.z = u.y;
	}
	else if (minn == y){
		v.y = 0;
		v.x = -u.z;
		v.z = u.x;
	}
	else {
		v.z = 0;
		v.x = -u.y;
		v.y = u.x;
	}
	deepcopyVec3(v, normalizeVec3(v));
	deepcopyVec3(w, normalizeVec3(crossProductVec3(u,v)));
	Matrix4 M;
	deepcopy(M, getIdentityMatrix());
	M.values[0][0] = u.x;
	M.values[0][1] = u.y;
	M.values[0][2] = u.z;
	M.values[1][0] = v.x;
	M.values[1][1] = v.y;
	M.values[1][2] = v.z;
	M.values[2][0] = w.x;
	M.values[2][1] = w.y;
	M.values[2][2] = w.z;
	
	Matrix4 R;
	deepcopy(R, getIdentityMatrix());
	double radian = (r.angle * M_PI) / 180.0;
	R.values[1][1] = std::cos(radian);
	R.values[1][2] = (-1) * std::sin(radian);
	R.values[2][1] = std::sin(radian);
	R.values[2][2] = std::cos(radian);

	Matrix4 M_T;
	deepcopy(M_T, getIdentityMatrix());
	M_T.values[0][0] = u.x;
	M_T.values[1][0] = u.y;
	M_T.values[2][0] = u.z;
	M_T.values[0][1] = v.x;
	M_T.values[1][1] = v.y;
	M_T.values[2][1] = v.z;
	M_T.values[0][2] = w.x;
	M_T.values[1][2] = w.y;
	M_T.values[2][2] = w.z;

	deepcopy(m, multiplyMatrixWithMatrix(multiplyMatrixWithMatrix(M_T, multiplyMatrixWithMatrix(R, M)), m));
}

Matrix4 Scene::modelingTransformation(Mesh *mesh){
	Matrix4 M_model;
	deepcopy(M_model, getIdentityMatrix());
	for (int i = 0; i < mesh->numberOfTransformations; i++){
		int index = mesh->transformationIds[i] - 1;
		if (mesh->transformationTypes[i] == 's'){ //scaling
			Scaling *s = scalings[index];
			scaling(*s, M_model);
		}
		else if (mesh->transformationTypes[i] == 't'){ //translation
			Translation *t = translations[index];
			translation(*t, M_model);
		}
		else { //rotation
			Rotation *r = rotations[index];
			rotation(*r, M_model);
		}
	}
	return M_model;
	
}
Matrix4 Scene::cameraTransformation(Camera* cam){
	Matrix4 M_cam;
	deepcopy(M_cam, getIdentityMatrix());
	Vec3 u(cam -> u);
	Vec3 v(cam -> v);
	Vec3 w(cam -> w);
	Vec3 e(cam -> position);
	M_cam.values[0][0] = u.x;
	M_cam.values[0][1] = u.y;
	M_cam.values[0][2] = u.z;
	M_cam.values[0][3] = (-1)*(u.x * e.x + u.y * e.y + u.z * e.z);
	M_cam.values[1][0] = v.x;
	M_cam.values[1][1] = v.y;
	M_cam.values[1][2] = v.z;
	M_cam.values[1][3] = (-1)*(v.x * e.x + v.y * e.y + v.z * e.z);
	M_cam.values[2][0] = w.x;
	M_cam.values[2][1] = w.y;
	M_cam.values[2][2] = w.z;
	M_cam.values[2][3] = (-1)*(w.x * e.x + w.y * e.y + w.z * e.z);
	return M_cam;


}
Matrix4 Scene::projectionTransformation(Camera*cam){
	Matrix4 M_orth;
	deepcopy(M_orth, getIdentityMatrix());
	double r = cam -> right;
	double l = cam -> left;
	double t = cam -> top;
	double b = cam -> bottom;
	double n = cam -> near;
	double f = cam -> far;
	M_orth.values[0][0] = 2/(r-l);
	M_orth.values[0][3] = -(r+l)/(r-l);
	M_orth.values[1][1] = 2/(t-b);
	M_orth.values[1][3] = -(t+b)/(t-b);
	M_orth.values[2][2] = -2/(f-n);
	M_orth.values[2][3] = -(f+n)/(f-n);
	if(cam -> projectionType){
		Matrix4 M_p2o, M_per;
		M_p2o.values[0][0] = n;
		M_p2o.values[1][1] = n;
		M_p2o.values[2][2] = f + n;
		M_p2o.values[2][3] = f * n;
		M_p2o.values[3][2] = -1;
		deepcopy(M_per, multiplyMatrixWithMatrix(M_orth, M_p2o));
		return M_per;

	}
	return M_orth;
}
bool Scene::visible(double den, double num, double & t_e, double & t_l){
	if(den > 0){
		double t = num / den;
		if(t > t_l) return false;
		if(t > t_e) t_e = t;

	}
	else if(den < 0){
		double t = num / den;
		if(t < t_e) return false;
		if(t < t_l) t_l = t;
	}
	else if(num > 0) return false;
	return true;
}
bool Scene::clipping(Vec4 & vertex0, Vec4 & vertex1, Color &c0, Color &c1){ //??????
	double t_e = 0; double t_l = 1;
	bool vis = false;
	double d_x, d_y, d_z;
	double x_0 = vertex0.x;
	double y_0 = vertex0.y;
	double x_1 = vertex1.x;
	double y_1 = vertex1.y;
	Color new_c0 = c0;
	Color new_c1 = c1;
	double m0, m1;
	d_x = vertex1.x - vertex0.x;
	d_y = vertex1.y - vertex0.y;
	d_z = vertex1.z - vertex0.z;
	if (visible(d_x, -1 - vertex0.x, t_e, t_l)){//?????
		if(visible(-(d_x), vertex0.x - 1, t_e, t_l)){
			if(visible(d_y, -1 - vertex0.y, t_e, t_l)){
				if(visible(-(d_y), vertex0.y - 1, t_e, t_l)){
					if(visible(d_z, -1 - vertex0.z, t_e, t_l)){
						if(visible(-(d_z), vertex0.z - 1, t_e, t_l)){
							vis = true;
							if(t_l < 1){
								vertex1.x = vertex0.x + d_x * t_l;
								vertex1.y = vertex0.y + d_y * t_l;
								vertex1.z = vertex0.z + d_z * t_l;
							}
							if(t_e > 0){
								vertex0.x = vertex0.x + d_x * t_e;
								vertex0.y = vertex0.y + d_y * t_e;
								vertex0.z = vertex0.z + d_z * t_e;
							}
							if (vertex0.x != x_0){
								m0 = (vertex0.x - x_0) / (x_1 - x_0);
								new_c0.r = (1-m0)*c0.r + m0*c1.r;
								new_c0.g = (1-m0)*c0.g + m0*c1.g;
								new_c0.b = (1-m0)*c0.b + m0*c1.b;
								m1 = (vertex1.x - x_0) / (x_1 - x_0);
								new_c1.r = (1-m1)*c0.r + m1*c1.r;
								new_c1.g = (1-m1)*c0.g + m1*c1.g;
								new_c1.b = (1-m1)*c0.b + m1*c1.b;
							}
							else if (vertex0.y != y_0) {
								m0 = (vertex0.y - y_0) / (y_1 - y_0);
								new_c0.r = (1-m0)*c0.r + m0*c1.r;
								new_c0.g = (1-m0)*c0.g + m0*c1.g;
								new_c0.b = (1-m0)*c0.b + m0*c1.b;
								m1 = (vertex1.y - y_0) / (y_1 - y_0);
								new_c1.r = (1-m1)*c0.r + m1*c1.r;
								new_c1.g = (1-m1)*c0.g + m1*c1.g;
								new_c1.b = (1-m1)*c0.b + m1*c1.b;
							}
						}
					}
				}
			}
		}
	}
	c0 = new_c0;
	c1 = new_c1;
	return vis;
}
Matrix4 Scene::viewportTransformation(Camera * cam){ //x_min,y_min???
	int n_x = cam -> horRes;
	int n_y = cam -> verRes;
	Matrix4 M_vp;
	M_vp.values[0][0] = (double) n_x / 2;
	M_vp.values[1][1] = (double) n_y / 2;
	M_vp.values[0][3] = (double) (n_x - 1)/ 2;
	M_vp.values[1][3] = (double) (n_y - 1) / 2;
	M_vp.values[2][2] = 0.5;
	M_vp.values[2][3] = 0.5;
	return M_vp;



}
Vec4 Scene::vec3tovec4(Vec3 v){
	Vec4 result;
	result.x = v.x;
	result.y = v.y;
	result.z = v.z;
	result.t = 1;
	result.colorId = v.colorId;
	return result;
}
Vec3 Scene::vec4tovec3(Vec4 v){
	Vec3 result;
	result.x = v.x;
	result.y = v.y;
	result.z = v.z;
	result.colorId = v.colorId;
	return result;
}
void Scene::perspectiveDivide(Vec4& v){
	v.x /= v.t;
	v.y /= v.t;
	v.z /= v.t;
	v.t = 1;

}
void Scene::Transformation_Rasterization(Camera* cam, Mesh* mesh){
	Matrix4 M_model, M_proj, M_cam, M_vp;
	M_model = modelingTransformation(mesh);
	M_cam = cameraTransformation(cam);
	M_proj = projectionTransformation(cam);
	M_vp = viewportTransformation(cam);

	for (int i = 0; i < mesh->numberOfTriangles; i++){
		Triangle triangle = mesh->triangles[i];
		int index1 = triangle.vertexIds[0] - 1;
		int index2 = triangle.vertexIds[1] - 1;
		int index3 = triangle.vertexIds[2] - 1;
		Vec3 v1 = *(vertices[index1]);
		Vec3 v2 = *(vertices[index2]);
		Vec3 v3 = *(vertices[index3]);
		Vec4 vertex1, vertex2, vertex3;
		vertex1 = vec3tovec4(v1);
		vertex2 = vec3tovec4(v2);
		vertex3 = vec3tovec4(v3);
		vertex1 = multiplyMatrixWithVec4(M_proj, multiplyMatrixWithVec4(M_cam,multiplyMatrixWithVec4(M_model,vertex1)));
		vertex2 = multiplyMatrixWithVec4(M_proj, multiplyMatrixWithVec4(M_cam,multiplyMatrixWithVec4(M_model,vertex2)));
		vertex3 = multiplyMatrixWithVec4(M_proj, multiplyMatrixWithVec4(M_cam,multiplyMatrixWithVec4(M_model,vertex3)));
		v1 = vec4tovec3(vertex1);
		v2 = vec4tovec3(vertex2);
		v3 = vec4tovec3(vertex3);
		if(cullingEnabled){
			Vec3 normal = (crossProductVec3((subtractVec3(v2,v1)),(subtractVec3(v3,v1)))); 
			Vec3 center,v;
			center.x = (v1.x + v2.x + v3.x) / 3.0;
			center.y = (v1.y + v2.y + v3.y) / 3.0;
			center.z = (v1.z + v2.z + v3.z) / 3.0;
			
			if(dotProductVec3(normal,center) < 0) continue;
		}
		if(cam -> projectionType){
			perspectiveDivide(vertex1);
			perspectiveDivide(vertex2);
			perspectiveDivide(vertex3);
		}
		if(mesh -> type){ //triangle
			vertex1 = multiplyMatrixWithVec4(M_vp,vertex1);
			vertex2 = multiplyMatrixWithVec4(M_vp,vertex2);
			vertex3 = multiplyMatrixWithVec4(M_vp,vertex3);

			triangleRasterization(vertex1, vertex2, vertex3, index1, index2, index3, cam);
		}
		else{
			bool vis1, vis2, vis3;

			Vec4 vertex11(vertex1);
			Vec4 vertex21(vertex2);
			Vec4 vertex31(vertex3);

			Color c1 = *(colorsOfVertices[index1]);
			Color c2 = *(colorsOfVertices[index2]);
			Color c3 = *(colorsOfVertices[index3]);

			vis1 = clipping(vertex1, vertex2, c1, c2);
			vis2 = clipping(vertex21, vertex3, c2, c3);
			vis3 = clipping(vertex31, vertex11, c3, c1);

			vertex1 = multiplyMatrixWithVec4(M_vp,vertex1);
			vertex2 = multiplyMatrixWithVec4(M_vp,vertex2);
			vertex3 = multiplyMatrixWithVec4(M_vp,vertex3);
			vertex11 = multiplyMatrixWithVec4(M_vp,vertex11);
			vertex21 = multiplyMatrixWithVec4(M_vp,vertex21);
			vertex31 = multiplyMatrixWithVec4(M_vp,vertex31);

			if(vis1) lineRasterization(vertex1, vertex2, c1, c2, cam);
			if(vis2) lineRasterization(vertex21, vertex3, c2, c3, cam);
			if(vis3) lineRasterization(vertex31, vertex11, c3, c1, cam);
		}
		
	}

}


void Scene::lineRasterization(Vec4 vertex1, Vec4 vertex2, Color c1, Color c2, Camera *cam){
	double x_0 = vertex1.x;
	double y_0 = vertex1.y;
	double x_1 = vertex2.x;
	double y_1 = vertex2.y;
	double z_0 = vertex1.z;
	double z_1 = vertex2.z;
	double d, depthh, d_x, d_y;
	d_y = fabs(y_1-y_0);
	d_x = fabs(x_0-x_1);
	Color c, d_c;

	if (d_x == 0 && d_y == 0) {
        depthh = z_0;
		if (y_0 >= 0 && y_0 < cam->verRes && x_0 >= 0 && x_0 < cam->horRes && depth[x_0][y_0] >= depthh){
			assignColorToPixel(x_0, y_0, c1);
			depth[x_0][y_0] = depthh;
		}
    }
	else if (d_x == 0) {
		if (y_0 < y_1){
			c = c1;
			d_c.r = (c2.r - c.r) / d_y;
			d_c.g = (c2.g - c.g) / d_y;
			d_c.b = (c2.b - c.b) / d_y;
		}
		else {
			c = c2;
			d_c.r = (c1.r - c.r) / d_y;
			d_c.g = (c1.g - c.g) / d_y;
			d_c.b = (c1.b - c.b) / d_y;
			swap(z_0,z_1);
		}
		
		depthh = z_0;
        for (int y = min(y_0, y_1); y <= max(y_0, y_1) && y >= 0 && y < cam->verRes && x_0 >= 0 && x_0 < cam->horRes; y++) {
            depthh += ((z_1 - z_0)/ d_y);
			if (depth[x_0][y] >= depthh){
				assignColorToPixel(x_0, y, c);
            	depth[x_0][y] = depthh;
			}
			c.r += d_c.r;
			c.g += d_c.g;
			c.b += d_c.b;
        }
    } 
	else if (d_y == 0) {
		if (x_0 < x_1){
			c = c1;
			d_c.r = (c2.r - c.r) / d_x;
			d_c.g = (c2.g - c.g) / d_x;
			d_c.b = (c2.b - c.b) / d_x;
		}
		else {
			c = c2;
			d_c.r = (c1.r - c.r) / d_x;
			d_c.g = (c1.g - c.g) / d_x;
			d_c.b = (c1.b - c.b) / d_x;
			swap(z_0, z_1);
		}
		depthh = z_0;
        for (int x = min(x_0, x_1); x <= max(x_0, x_1) && y_0 >= 0 && y_0 < cam->verRes && x >= 0 && x < cam->horRes; x++) {
            depthh += (z_1 - z_0)/(d_x);
			if (depth[x][y_0] >= depthh){
				assignColorToPixel(x, y_0, c);
				depth[x][y_0] = depthh;
			}
			c.r += d_c.r;
			c.g += d_c.g;
			c.b += d_c.b;
        }
    }
	double m = d_y / d_x;
	if (m <= 1.0){
		if (x_0 < x_1){
			c = c1;
			d_c.r = (c2.r - c.r) / d_x;
			d_c.g = (c2.g - c.g) / d_x;
			d_c.b = (c2.b - c.b) / d_x;
		}
		else {
			c = c2;
			d_c.r = (c1.r - c.r) / d_x;
			d_c.g = (c1.g - c.g) / d_x;
			d_c.b = (c1.b - c.b) / d_x;
		}
		if (x_0 > x_1) {
            swap(x_0, x_1);
            swap(y_0, y_1);
			swap(z_0, z_1);
        }
		int x = x_0;
        int y = y_0;
		d = 2*d_y - d_x;
		depthh = z_0;
		for (; x <= x_1 && y >= 0 && y < cam->verRes && x >= 0 && x < cam->horRes; x++){
			depthh += ((z_1 - z_0)/d_x);
			if (depth[x][y] >= depthh){
				depth[x][y] = depthh;
				assignColorToPixel(x,y,c);
			}
			if(d >= 0){
				y += (y_1 > y_0) ? 1 : -1;
				d += 2*d_y - 2*d_x;
			}
			else{
				d +=2*d_y;
			}
			c.r += d_c.r;
			c.g += d_c.g;
			c.b += d_c.b;
		}
	}
	
	else if (m > 1.0){
		if (y_0 < y_1){
			c = c1;
			d_c.r = (c2.r - c.r) / d_y;
			d_c.g = (c2.g - c.g) / d_y;
			d_c.b = (c2.b - c.b) / d_y;
		}
		else {
			c = c2;
			d_c.r = (c1.r - c.r) / d_y;
			d_c.g = (c1.g - c.g) / d_y;
			d_c.b = (c1.b - c.b) / d_y;
		}
		if (y_0 > y_1) {
            swap(x_0, x_1);
            swap(y_0, y_1);
			swap(z_0, z_1);
        }
		int x = x_0;
        int y = y_0;
		d = 2*d_x - d_y;
		depthh = z_0;
		for (; y <= y_1 && y >= 0 && y < cam->verRes && x >= 0 && x < cam->horRes; y++){
			depthh += ((z_1 - z_0)/d_y);
			if (depth[x][y] >= depthh){
				depth[x][y] = depthh;
				assignColorToPixel(x,y,c);
			}
			if(d >= 0){
				x += (x_1 > x_0) ? 1 : -1;
				d += 2*d_x - 2*d_y;
			}
			else{
				d += 2*d_x;
			}
			c.r += d_c.r;
			c.g += d_c.g;
			c.b += d_c.b;
		}
	}
}
void Scene::triangleRasterization(Vec4 vertex1, Vec4 vertex2, Vec4 vertex3, int index1, int index2, int index3, Camera *cam){
	double x_0 = vertex1.x;
	double y_0 = vertex1.y;
	double x_1 = vertex2.x;
	double y_1 = vertex2.y;
	double x_2 = vertex3.x;
	double y_2 = vertex3.y;
	double depthh;
	double f_01, f_12, f_20;
	double alpha,beta, gamma;
	double den_12,den_20, den_01;
	double x_min = std::min(std::min(x_0, x_1), x_2);
	double x_max = std::max(std::max(x_0, x_1), x_2);
	double y_min = std::min(std::min(y_0, y_1), y_2);
	double y_max = std::max(std::max(y_0, y_1), y_2);
	den_01 = x_2 * (y_0 - y_1) + y_2 * (x_1 - x_0) + x_0 * y_1 - y_0 * x_1;
	den_12 = x_0 * (y_1 - y_2) + y_0 * (x_2 - x_1) + x_1* y_2 - y_1 * x_2 ;
	den_20 = x_1 * (y_2 - y_0) + y_1 * (x_0 - x_2) + x_2* y_0 - y_2 * x_0;
	Color c;
	Color c_0(*colorsOfVertices[index1]);
	Color c_1(*colorsOfVertices[index2]);
	Color c_2(*colorsOfVertices[index3]);
	for(int y = y_min; y <= y_max && y >= 0 && y < cam->verRes; y++){
		for(int x = x_min; x <= x_max && x >= 0 && x < cam->horRes; x++){
			f_01 = x * (y_0 - y_1) + y * (x_1 - x_0) + x_0 * y_1 - y_0 * x_1;
			f_12 = x * (y_1 - y_2) + y * (x_2 - x_1) + x_1 * y_2 - y_1 * x_2;
			f_20 = x * (y_2 - y_0) + y * (x_0 - x_2) + x_2 * y_0 - y_2 * x_0;
			alpha = f_12 / den_12;
			beta = f_20 / den_20;
			gamma = f_01 / den_01;
			if(alpha >= 0 && beta >= 0 && gamma >= 0){
				depthh = vertex1.z * alpha + vertex2.z * beta + vertex3.z * gamma;
				if(depth[x][y] >= depthh){
					depth[x][y] = depthh;
					c.r = alpha * c_0.r + beta * c_1.r + gamma * c_2.r;
					c.g = alpha * c_0.g + beta * c_1.g + gamma * c_2.g;
					c.b = alpha * c_0.b + beta * c_1.b + gamma * c_2.b;
					assignColorToPixel(x,y,c);
				}
				
			}
		}

	}
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	for (int i = 0; i < depth.size(); i++){
		for (int j = 0; j < depth[i].size(); j++){
			depth[i][j] = 1;
		}
	}
	for(int i = 0; i < meshes.size(); i++){
		Transformation_Rasterization(camera, meshes[i]);
	}

}