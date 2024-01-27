#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <ft2build.h>
#include <time.h>
#include FT_FREETYPE_H
#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[6];
int gWidth = 640, gHeight = 480;
GLuint vao1, vao;
struct Character {
	GLuint TextureID;   // ID handle of the glyph texture
	glm::ivec2 Size;    // Size of glyph
	glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
	GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;

glm::mat4 projectionMatrix; //bunny
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix;
glm::vec3 eyePos(0, 0, 0);

glm::mat4 QuadmodelingMatrix; //quad
glm::vec3 QuadeyePos(0, 0, 0);


glm::mat4 Cube1modelingMatrix; //cube1
glm::vec3 Cube1eyePos(0, 0, 0);

glm::mat4 Cube2modelingMatrix; //cube2
glm::vec3 Cube2eyePos(0, 0, 0);

glm::mat4 Cube3modelingMatrix; //cube3
glm::vec3 Cube3eyePos(0, 0, 0);



bool restart = false;
int activeProgramIndex = 0;
float left_right_offset = 0.f;
float time2 = 50.f;
struct Vertex
{
	Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};

struct Texture
{
	Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
	GLfloat u, v;
};

struct Normal
{
	Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
	GLuint vIndex[3], tIndex[3], nIndex[3];
};

vector<Vertex> gBunnyVertices;
vector<Texture> gBunnyTextures;
vector<Normal> gBunnyNormals;
vector<Face> gBunnyFaces;

vector<Vertex> gQuadVertices;
vector<Texture> gQuadTextures;
vector<Normal> gQuadNormals;
vector<Face> gQuadFaces;

vector<Vertex> gCubeVertices; //????
vector<Texture> gCubeTextures;
vector<Normal> gCubeNormals;
vector<Face> gCubeFaces;

GLuint gVertexAttribBuffer, gIndexBuffer; //Bunny
GLint gInVertexLoc, gInNormalLoc;

GLuint gQuadVertexAttribBuffer, gQuadIndexBuffer;//Quad
GLint gQuadInVertexLoc, gQuadInNormalLoc;

GLuint gCube1VertexAttribBuffer, gCube1IndexBuffer;//Cube1
GLint gCube1InVertexLoc, gCube1InNormalLoc;

GLuint gCube2VertexAttribBuffer, gCube2IndexBuffer;//Cube2
GLint gCube2InVertexLoc, gCube2InNormalLoc;

GLuint gCube3VertexAttribBuffer, gCube3IndexBuffer;//Cube2
GLint gCube3InVertexLoc, gCube3InNormalLoc;

GLuint gTextVBO; // Text

int gQuadVertexDataSizeInBytes, gQuadNormalDataSizeInBytes;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;
int gCubeVertexDataSizeInBytes, gCubeNormalDataSizeInBytes;


bool ParseObj(const string& fileName, vector<Vertex>& gVertices, vector<Texture>& gTextures, vector<Face>& gFaces, vector<Normal>& gNormals)
{
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			stringstream str(curLine);
			GLfloat c1, c2, c3;
			GLuint index[9];
			string tmp;

			if (curLine.length() >= 2)
			{
				if (curLine[0] == 'v')
				{
					if (curLine[1] == 't') // texture
					{
						str >> tmp; // consume "vt"
						str >> c1 >> c2;
						gTextures.push_back(Texture(c1, c2));
					}
					else if (curLine[1] == 'n') // normal
					{
						str >> tmp; // consume "vn"
						str >> c1 >> c2 >> c3;
						gNormals.push_back(Normal(c1, c2, c3));
					}
					else // vertex
					{
						str >> tmp; // consume "v"
						str >> c1 >> c2 >> c3;
						gVertices.push_back(Vertex(c1, c2, c3));
					}
				}
				else if (curLine[0] == 'f') // face
				{
					str >> tmp; // consume "f"
					char c;
					int vIndex[3], nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0];
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1];
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2];

					assert(vIndex[0] == nIndex[0] &&
						vIndex[1] == nIndex[1] &&
						vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

					gFaces.push_back(Face(vIndex, tIndex, nIndex));
				}
				else
				{
					cout << "Ignoring unidentified line in obj file: " << curLine << endl;
				}
			}

			//data += curLine;
			if (!myfile.eof())
			{
				//data += "\n";
			}
		}

		myfile.close();
	}
	else
	{
		return false;
	}

	/*
	for (int i = 0; i < gVertices.size(); ++i)
	{
		Vector3 n;

		for (int j = 0; j < gFaces.size(); ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				if (gFaces[j].vIndex[k] == i)
				{
					// face j contains vertex i
					Vector3 a(gVertices[gFaces[j].vIndex[0]].x,
							  gVertices[gFaces[j].vIndex[0]].y,
							  gVertices[gFaces[j].vIndex[0]].z);

					Vector3 b(gVertices[gFaces[j].vIndex[1]].x,
							  gVertices[gFaces[j].vIndex[1]].y,
							  gVertices[gFaces[j].vIndex[1]].z);

					Vector3 c(gVertices[gFaces[j].vIndex[2]].x,
							  gVertices[gFaces[j].vIndex[2]].y,
							  gVertices[gFaces[j].vIndex[2]].z);

					Vector3 ab = b - a;
					Vector3 ac = c - a;
					Vector3 normalFromThisFace = (ab.cross(ac)).getNormalized();
					n += normalFromThisFace;
				}

			}
		}

		n.normalize();

		gNormals.push_back(Normal(n.x, n.y, n.z));
	}
	*/

	assert(gVertices.size() == gNormals.size());

	return true;
}

bool ReadDataFromFile(
	const string& fileName, ///< [in]  Name of the shader file
	string& data)     ///< [out] The contents of the file
{
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			data += curLine;
			if (!myfile.eof())
			{
				data += "\n";
			}
		}

		myfile.close();
	}
	else
	{
		return false;
	}

	return true;
}

GLuint createVS(const char* shaderName)
{
	string shaderSource;

	string filename(shaderName);
	if (!ReadDataFromFile(filename, shaderSource))
	{
		cout << "Cannot find file name: " + filename << endl;
		exit(-1);
	}

	GLint length = shaderSource.length();
	const GLchar* shader = (const GLchar*)shaderSource.c_str();

	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vs, 1, &shader, &length);
	glCompileShader(vs);

	char output[1024] = { 0 };
	glGetShaderInfoLog(vs, 1024, &length, output);
	printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
	string shaderSource;

	string filename(shaderName);
	if (!ReadDataFromFile(filename, shaderSource))
	{
		cout << "Cannot find file name: " + filename << endl;
		exit(-1);
	}

	GLint length = shaderSource.length();
	const GLchar* shader = (const GLchar*)shaderSource.c_str();

	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fs, 1, &shader, &length);
	glCompileShader(fs);

	char output[1024] = { 0 };
	glGetShaderInfoLog(fs, 1024, &length, output);
	printf("FS compile log: %s\n", output);

	return fs;
}

void initShaders()
{
	// Create the programs

	gProgram[0] = glCreateProgram();
	gProgram[1] = glCreateProgram();
	gProgram[2] = glCreateProgram();
	gProgram[3] = glCreateProgram();
	gProgram[4] = glCreateProgram();
	gProgram[5] = glCreateProgram();
	// Create the shaders for both programs
	glBindAttribLocation(gProgram[5], 2, "vertex");
	GLuint vs1 = createVS("quadvert.glsl");
	GLuint fs1 = createFS("quadfrag.glsl");

	GLuint vs2 = createVS("vert2.glsl");
	GLuint fs2 = createFS("frag2.glsl");

	GLuint vs3 = createVS("vert2.glsl");
	GLuint fs3 = createFS("frag2.glsl");

	GLuint vs4 = createVS("cube1vert.glsl");
	GLuint fs4 = createFS("cube1frag.glsl");

	GLuint vs5 = createVS("cube1vert.glsl");
	GLuint fs5 = createFS("cube1frag.glsl");

	GLuint vs6 = createVS("vert_text.glsl");
	GLuint fs6 = createFS("frag_text.glsl");



	// Attach the shaders to the programs

	glAttachShader(gProgram[0], vs2); //Bunny
	glAttachShader(gProgram[0], fs2);

	glAttachShader(gProgram[1], vs1); // Quad
	glAttachShader(gProgram[1], fs1);

	glAttachShader(gProgram[2], vs3); // Cube1
	glAttachShader(gProgram[2], fs3);

	glAttachShader(gProgram[3], vs4); // Cube2
	glAttachShader(gProgram[3], fs4);

	glAttachShader(gProgram[4], vs5); // Cube3
	glAttachShader(gProgram[4], fs5);

	glAttachShader(gProgram[5], vs6); // Text
	glAttachShader(gProgram[5], fs6);


	// Link the programs

	glLinkProgram(gProgram[0]);
	GLint status;
	glGetProgramiv(gProgram[0], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}

	glLinkProgram(gProgram[1]);
	glGetProgramiv(gProgram[1], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}
	// glUseProgram(gProgram[0]); // ??

	// Get the locations of the uniform variables from both programs
	glLinkProgram(gProgram[2]);
	glGetProgramiv(gProgram[2], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}
	glLinkProgram(gProgram[3]);
	glGetProgramiv(gProgram[3], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}
	glLinkProgram(gProgram[4]);
	glGetProgramiv(gProgram[4], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}
	glLinkProgram(gProgram[5]);
	glGetProgramiv(gProgram[5], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}


}
void initFonts(int windowWidth, int windowHeight)
{
	// Set OpenGL options
	//glEnable(GL_CULL_FACE);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
	glUseProgram(gProgram[5]);
	glUniformMatrix4fv(glGetUniformLocation(gProgram[5], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

	// FreeType
	FT_Library ft;
	// All functions return a value different than 0 whenever an error occurred
	if (FT_Init_FreeType(&ft))
	{
		std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
	}

	// Load font as face
	FT_Face face;
	if (FT_New_Face(ft, "D:/VisualStudio/Projects/opengl_hw/LiberationSerif-Italic.ttf", 0, &face))
	{
		std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
	}

	// Set size to load glyphs as
	FT_Set_Pixel_Sizes(face, 0, 48);

	// Disable byte-alignment restriction
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// Load first 128 characters of ASCII set
	for (GLubyte c = 0; c < 128; c++)
	{
		// Load character glyph 
		if (FT_Load_Char(face, c, FT_LOAD_RENDER))
		{
			std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
			continue;
		}
		// Generate texture
		GLuint texture;
		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexImage2D(
			GL_TEXTURE_2D,
			0,
			GL_RED,
			face->glyph->bitmap.width,
			face->glyph->bitmap.rows,
			0,
			GL_RED,
			GL_UNSIGNED_BYTE,
			face->glyph->bitmap.buffer
		);
		// Set texture options
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		// Now store character for later use
		Character character = {
			texture,
			glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
			glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
			face->glyph->advance.x
		};
		Characters.insert(std::pair<GLchar, Character>(c, character));
	}

	glBindTexture(GL_TEXTURE_2D, 0);
	// Destroy FreeType once we're finished
	FT_Done_Face(face);
	FT_Done_FreeType(ft);

	//
	// Configure VBO for texture quads
	//

	glGenBuffers(1, &gTextVBO);
	glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);



	glGenVertexArrays(1, &vao1);
	assert(vao1 > 0);
	glBindVertexArray(vao1);

	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

}

void initVBO(GLuint& gVertexAttribBuffer, GLuint& gIndexBuffer, vector<Vertex>& gVertices, vector<Texture>& gTextures, vector<Normal>& gNormals, vector<Face>& gFaces, int& gVertexDataSizeInBytes, int& gNormalDataSizeInBytes)
{
	glGenVertexArrays(1, &vao);
	assert(vao > 0);
	glBindVertexArray(vao);
	cout << "vao = " << vao << endl;

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer);
	glGenBuffers(1, &gIndexBuffer);

	assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
	int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat[gVertices.size() * 3];
	GLfloat* normalData = new GLfloat[gNormals.size() * 3];
	GLuint* indexData = new GLuint[gFaces.size() * 3];

	float minX = 1e6, maxX = -1e6;
	float minY = 1e6, maxY = -1e6;
	float minZ = 1e6, maxZ = -1e6;

	for (int i = 0; i < gVertices.size(); ++i)
	{
		vertexData[3 * i] = gVertices[i].x;
		vertexData[3 * i + 1] = gVertices[i].y;
		vertexData[3 * i + 2] = gVertices[i].z;

		minX = std::min(minX, gVertices[i].x);
		maxX = std::max(maxX, gVertices[i].x);
		minY = std::min(minY, gVertices[i].y);
		maxY = std::max(maxY, gVertices[i].y);
		minZ = std::min(minZ, gVertices[i].z);
		maxZ = std::max(maxZ, gVertices[i].z);
	}

	std::cout << "minX = " << minX << std::endl;
	std::cout << "maxX = " << maxX << std::endl;
	std::cout << "minY = " << minY << std::endl;
	std::cout << "maxY = " << maxY << std::endl;
	std::cout << "minZ = " << minZ << std::endl;
	std::cout << "maxZ = " << maxZ << std::endl;

	for (int i = 0; i < gNormals.size(); ++i)
	{
		normalData[3 * i] = gNormals[i].x;
		normalData[3 * i + 1] = gNormals[i].y;
		normalData[3 * i + 2] = gNormals[i].z;
	}

	for (int i = 0; i < gFaces.size(); ++i)
	{
		indexData[3 * i] = gFaces[i].vIndex[0];
		indexData[3 * i + 1] = gFaces[i].vIndex[1];
		indexData[3 * i + 2] = gFaces[i].vIndex[2];
	}


	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

	// done copying to GPU memory; can free now from CPU memory
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

}
void renderText(const std::string& text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
	// Activate corresponding render state	
	glUseProgram(gProgram[5]);
	glBindVertexArray(vao1);
	glUniform3f(glGetUniformLocation(gProgram[5], "textColor"), color.x, color.y, color.z);
	glActiveTexture(GL_TEXTURE0);

	// Iterate through all characters
	std::string::const_iterator c;
	for (c = text.begin(); c != text.end(); c++)
	{
		Character ch = Characters[*c];

		GLfloat xpos = x + ch.Bearing.x * scale;
		GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

		GLfloat w = ch.Size.x * scale;
		GLfloat h = ch.Size.y * scale;

		// Update VBO for each character
		GLfloat vertices[6][4] = {
			{ xpos,     ypos + h,   0.0, 0.0 },
			{ xpos,     ypos,       0.0, 1.0 },
			{ xpos + w, ypos,       1.0, 1.0 },

			{ xpos,     ypos + h,   0.0, 0.0 },
			{ xpos + w, ypos,       1.0, 1.0 },
			{ xpos + w, ypos + h,   1.0, 0.0 }
		};

		// Render glyph texture over quad
		glBindTexture(GL_TEXTURE_2D, ch.TextureID);

		// Update content of VBO memory
		glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

		//glBindBuffer(GL_ARRAY_BUFFER, 0);

		// Render quad
		glDrawArrays(GL_TRIANGLES, 0, 6);
		// Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

		x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
	}

	glBindTexture(GL_TEXTURE_2D, 0);
	glBindVertexArray(vao);
}
void initBunny()
{
	//ParseObj("armadillo.obj");
	ParseObj("bunny.obj", gBunnyVertices, gBunnyTextures, gBunnyFaces, gBunnyNormals);
	glEnable(GL_DEPTH_TEST);
	initShaders();
	initFonts(gWidth, gHeight);
	initVBO(gVertexAttribBuffer, gIndexBuffer, gBunnyVertices, gBunnyTextures, gBunnyNormals, gBunnyFaces, gVertexDataSizeInBytes, gNormalDataSizeInBytes);
}
void initQuad()
{
	//ParseObj("armadillo.obj");
	ParseObj("quad.obj", gQuadVertices, gQuadTextures, gQuadFaces, gQuadNormals);
	glEnable(GL_DEPTH_TEST);
	initShaders();
	initFonts(gWidth, gHeight);
	initVBO(gQuadVertexAttribBuffer, gQuadIndexBuffer, gQuadVertices, gQuadTextures, gQuadNormals, gQuadFaces, gQuadVertexDataSizeInBytes, gQuadNormalDataSizeInBytes);
}
void initCube()
{
	//ParseObj("armadillo.obj");
	ParseObj("cube.obj", gCubeVertices, gCubeTextures, gCubeFaces, gCubeNormals);
	glEnable(GL_DEPTH_TEST);
	initShaders();
	initFonts(gWidth, gHeight);
	initVBO(gCube1VertexAttribBuffer, gCube1IndexBuffer, gCubeVertices, gCubeTextures, gCubeNormals, gCubeFaces, gCubeVertexDataSizeInBytes, gCubeNormalDataSizeInBytes);
	initVBO(gCube2VertexAttribBuffer, gCube2IndexBuffer, gCubeVertices, gCubeTextures, gCubeNormals, gCubeFaces, gCubeVertexDataSizeInBytes, gCubeNormalDataSizeInBytes);
	initVBO(gCube3VertexAttribBuffer, gCube3IndexBuffer, gCubeVertices, gCubeTextures, gCubeNormals, gCubeFaces, gCubeVertexDataSizeInBytes, gCubeNormalDataSizeInBytes);
}


void drawBunnyModel()
{
	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
	glDrawElements(GL_TRIANGLES, gBunnyFaces.size() * 3, GL_UNSIGNED_INT, 0);
}
void drawQuadModel()
{
	glBindBuffer(GL_ARRAY_BUFFER, gQuadVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gQuadIndexBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gQuadVertexDataSizeInBytes));
	glDrawElements(GL_TRIANGLES, gQuadFaces.size() * 3, GL_UNSIGNED_INT, 0);
}
void drawCube1Model()
{
	glBindBuffer(GL_ARRAY_BUFFER, gCube1VertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gCube1IndexBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gCubeVertexDataSizeInBytes));
	glDrawElements(GL_TRIANGLES, gCubeFaces.size() * 3, GL_UNSIGNED_INT, 0);
}
void drawCube2Model()
{
	glBindBuffer(GL_ARRAY_BUFFER, gCube2VertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gCube2IndexBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gCubeVertexDataSizeInBytes));
	glDrawElements(GL_TRIANGLES, gCubeFaces.size() * 3, GL_UNSIGNED_INT, 0);
}
void drawCube3Model()
{
	glBindBuffer(GL_ARRAY_BUFFER, gCube3VertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gCube3IndexBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gCubeVertexDataSizeInBytes));
	glDrawElements(GL_TRIANGLES, gCubeFaces.size() * 3, GL_UNSIGNED_INT, 0);
}
float speed = 0.3f;
bool collision_red = false;
bool collision_yellow = false;
void display()
{
	glClearColor(0, 0, 0, 1);
	glClearDepth(1.0f);
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	static float angle = 0;
	static float posX = 0.f;
	static float posY = -2.2f;
	static float posZ = -3.f;
	const float threshold = 0.1f;
	static bool up = true;
	static float Cube3posX = 5.f;
	static float Cube3posZ = -40.f;
	static float Cube2posX = 0.f;
	static float Cube2posZ = -40.f;
	static float Quadangle = -90;
	static float QuadposX = 0.f;
	float angleRad = (float)(angle / 180.0) * M_PI;
	float QuadposY = -2.8f;
	static float QuadposZ = -2.4f;
	float QuadangleRad = (float)(Quadangle / 180.0) * M_PI;
	static float Cube1posX = -5.f;
	float CubeposY = -4.0f;
	static float Cube1posZ = -40.f;
	static bool hit = false;
	static int score = 0;
	glm::vec3 blackColor = glm::vec3(0.0f, 0.0f, 0.0f);  // Set your desired color for black squares
	glm::vec3 whiteColor = glm::vec3(0.78f, 0.63f, 0.78f);  // Set your desired color for white squares
	float offsetValue = 200.f;  // Set the offset value
	float scaleValue = 0.35f;   // Set the scale value
	static float scaleZ = 0.0;

	if (collision_red) {
		glUseProgram(gProgram[0]);
		// Compute the modeling matrix 
		glm::mat4 matT = glm::translate(glm::mat4(1.0), glm::vec3(posX, posY, posZ));
		glm::mat4 matS = glm::scale(glm::mat4(1.0), glm::vec3(0.6, 0.6, 0.6));
		glm::mat4 matR = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
		glm::mat4 matZ = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(0.0, 0.0, 1.0));
		glm::mat4 matRz = glm::rotate(glm::mat4(1.0), angleRad, glm::vec3(0.0, 0.0, 1.0));
		modelingMatrix = matT * matZ * matRz * matS * matR; // starting from right side, rotate around Y to turn back, then rotate around Z some more at each frame, then translate.

		glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[0], "eyePos"), 1, glm::value_ptr(eyePos));
		// Draw the scene
		drawBunnyModel();

		glUseProgram(gProgram[1]);
		glm::mat4 matTQ = glm::translate(glm::mat4(1.0), glm::vec3(QuadposX, QuadposY, QuadposZ));
		glm::mat4 matSQ = glm::scale(glm::mat4(1.0), glm::vec3(7, 100, 0));
		glm::mat4 matRQ = glm::rotate<float>(glm::mat4(1.0), (-180. / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
		glm::mat4 matRzQ = glm::rotate(glm::mat4(1.0), QuadangleRad, glm::vec3(1.0, 0.0, 0.0));

		QuadmodelingMatrix = matTQ * matRzQ * matSQ * matRQ;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(QuadmodelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[1], "eyePos"), 1, glm::value_ptr(QuadeyePos));


		glUniform3fv(glGetUniformLocation(gProgram[1], "black"), 1, glm::value_ptr(blackColor));
		glUniform3fv(glGetUniformLocation(gProgram[1], "white"), 1, glm::value_ptr(whiteColor));
		glUniform1f(glGetUniformLocation(gProgram[1], "offset"), offsetValue);
		glUniform1f(glGetUniformLocation(gProgram[1], "scaleZ"), scaleZ);
		glUniform1f(glGetUniformLocation(gProgram[1], "scale"), scaleValue);
		drawQuadModel();


		glUseProgram(gProgram[2]);
		// Compute the modeling matrix 
		glm::mat4 matTC = glm::translate(glm::mat4(1.0), glm::vec3(Cube1posX, CubeposY, Cube1posZ));
		glm::mat4 matSC = glm::scale(glm::mat4(1.0), glm::vec3(1.15f, 5.f, 0.5f));
		Cube1modelingMatrix = matTC * matSC;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(Cube1modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[2], "eyePos"), 1, glm::value_ptr(Cube1eyePos));
		drawCube1Model();

		glUseProgram(gProgram[3]);
		// Compute the modeling matrix 
		glm::mat4 matTC1 = glm::translate(glm::mat4(1.0), glm::vec3(Cube2posX, CubeposY, Cube2posZ));
		glm::mat4 matSC1 = glm::scale(glm::mat4(1.0), glm::vec3(1.15f, 5.f, 0.5f));
		Cube2modelingMatrix = matTC1 * matSC1;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(Cube2modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[3], "eyePos"), 1, glm::value_ptr(Cube2eyePos));
		drawCube2Model();
		glUseProgram(gProgram[4]);
		// Compute the modeling matrix 
		glm::mat4 matTC2 = glm::translate(glm::mat4(1.0), glm::vec3(Cube3posX, CubeposY, Cube3posZ));
		glm::mat4 matSC2 = glm::scale(glm::mat4(1.0), glm::vec3(1.15f, 5.f, 0.5f));
		Cube3modelingMatrix = matTC2 * matSC2;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[4], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[4], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[4], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(Cube3modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[4], "eyePos"), 1, glm::value_ptr(Cube3eyePos));
		drawCube3Model();
	}
	else {
		float CubePos[4];
		CubePos[1] = -5;
		CubePos[2] = 0;
		CubePos[3] = 5;

		int index1;
		int index2;
		int index3;

		if (restart) {
			speed = 0.3;
			scaleZ = 0;
			score = 0;
			posX = 0.f;
			posY = -2.2f;
			posZ = -3.f;
			Cube3posX = 5.f;
			Cube3posZ = -40.f;
			Cube2posX = 0.f;
			Cube2posZ = -40.f;
			QuadposX = 0.f;
			QuadposY = -2.8f;
			QuadposZ = -2.4f;
			Cube1posX = -5.f;
			Cube1posZ = -40.f;
			time2 = 50.f;
			collision_red = false;
			angle = 0;
			collision_yellow = false;

			index1 = (rand() % 3) + 1;
			index2 = index1;
			while (index1 == index2) {
				index2 = (rand() % 3) + 1;
			}
			index3 = index2;
			while (index3 == index2 || index3 == index1) {
				index3 = (rand() % 3) + 1;
			}
			Cube1posX = CubePos[index1];
			Cube2posX = CubePos[index2];
			Cube3posX = CubePos[index3];
		}

		if (hit) {
			score += 1000;
			hit = false;

		}
		if (collision_yellow) {
			angle += speed * 25;
			if (angle >= 360) {
				angle = 0;
				collision_yellow = false;
			}
		}
		posX += left_right_offset;
		viewingMatrix = glm::lookAt(glm::vec3(0, 0.5, posZ + 3.5), glm::vec3(0, 0.5, posZ + 3.5) + glm::vec3(0, 0, -1), glm::vec3(0, 1, 0));
		if (up) {
			posY += speed * 0.1;
		}
		else if (!up) {
			posY -= speed * 0.1;
		}
		if (posY >= -2.f) {
			posY = -2.f;
			up = false;
		}
		else if (posY <= -2.4f) {
			posY = -2.4f;
			up = true;
		}

		if (posX > 4.5f) {
			posX = 4.5f;
		}
		else if (posX < -4.5f) {
			posX = -4.5f;
		}


		glUseProgram(gProgram[0]);
		// Compute the modeling matrix 
		glm::mat4 matT = glm::translate(glm::mat4(1.0), glm::vec3(posX, posY, posZ));
		glm::mat4 matS = glm::scale(glm::mat4(1.0), glm::vec3(0.6, 0.6, 0.6));
		glm::mat4 matR = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
		glm::mat4 matRz = glm::rotate(glm::mat4(1.0), angleRad, glm::vec3(0.0, 1.0, 0.0));

		modelingMatrix = matT * matS * matRz * matR; // starting from right side, rotate around Y to turn back, then rotate around Z some more at each frame, then translate.

		glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[0], "eyePos"), 1, glm::value_ptr(eyePos));
		// Draw the scene
		drawBunnyModel();
		//posZ -= abs(offset);

		glUseProgram(gProgram[1]);
		// Compute the modeling matrix 
		glm::mat4 matTQ = glm::translate(glm::mat4(1.0), glm::vec3(QuadposX, QuadposY, QuadposZ));
		glm::mat4 matSQ = glm::scale(glm::mat4(1.0), glm::vec3(7, 100, 0));
		glm::mat4 matRQ = glm::rotate<float>(glm::mat4(1.0), (-180. / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
		glm::mat4 matRzQ = glm::rotate(glm::mat4(1.0), QuadangleRad, glm::vec3(1.0, 0.0, 0.0));
		QuadmodelingMatrix = matTQ * matRzQ * matSQ * matRQ;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(QuadmodelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[1], "eyePos"), 1, glm::value_ptr(QuadeyePos));
		glUniform3fv(glGetUniformLocation(gProgram[1], "black"), 1, glm::value_ptr(blackColor));
		glUniform3fv(glGetUniformLocation(gProgram[1], "white"), 1, glm::value_ptr(whiteColor));
		glUniform1f(glGetUniformLocation(gProgram[1], "offset"), offsetValue);
		glUniform1f(glGetUniformLocation(gProgram[1], "scaleZ"), scaleZ);
		glUniform1f(glGetUniformLocation(gProgram[1], "scale"), scaleValue);
		drawQuadModel();


		glUseProgram(gProgram[2]);
		glUniform3fv(glGetUniformLocation(gProgram[2], "eyePos"), 1, glm::value_ptr(Cube1eyePos));
		// Compute the modeling matrix 
		glm::mat4 matTC = glm::translate(glm::mat4(1.0), glm::vec3(Cube1posX, CubeposY, Cube1posZ));
		glm::mat4 matSC = glm::scale(glm::mat4(1.0), glm::vec3(1.15f, 5.f, 0.5f));
		Cube1modelingMatrix = matTC * matSC;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(Cube1modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[2], "eyePos"), 1, glm::value_ptr(Cube1eyePos));
		drawCube1Model();

		glUseProgram(gProgram[3]);
		// Compute the modeling matrix 
		glm::mat4 matTC1 = glm::translate(glm::mat4(1.0), glm::vec3(Cube2posX, CubeposY, Cube2posZ));
		glm::mat4 matSC1 = glm::scale(glm::mat4(1.0), glm::vec3(1.15f, 5.f, 0.5f));
		Cube2modelingMatrix = matTC1 * matSC1;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(Cube2modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[3], "eyePos"), 1, glm::value_ptr(Cube2eyePos));
		drawCube2Model();

		glUseProgram(gProgram[4]);
		// Compute the modeling matrix 
		glm::mat4 matTC2 = glm::translate(glm::mat4(1.0), glm::vec3(Cube3posX, CubeposY, Cube3posZ));
		glm::mat4 matSC2 = glm::scale(glm::mat4(1.0), glm::vec3(1.15f, 5.f, 0.5f));
		Cube3modelingMatrix = matTC2 * matSC2;

		glUniformMatrix4fv(glGetUniformLocation(gProgram[4], "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[4], "viewingMatrix"), 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(glGetUniformLocation(gProgram[4], "modelingMatrix"), 1, GL_FALSE, glm::value_ptr(Cube3modelingMatrix));
		glUniform3fv(glGetUniformLocation(gProgram[4], "eyePos"), 1, glm::value_ptr(Cube3eyePos));
		drawCube3Model();

		if (abs(posX - Cube2posX) <= 1.3 && abs(posZ - (Cube2posZ + 0.5)) <= 0.5) {
			collision_red = true;
			Cube2posZ = 10;
		}
		else if (abs(posX - Cube3posX) <= 1.3 && abs(posZ - (Cube3posZ + 0.5)) <= 0.5) {
			collision_red = true;
			Cube3posZ = 10;
		}
		else if (abs(posX - Cube1posX) <= 1.3 && abs(posZ - Cube1posZ) <= 0.5) {
			collision_yellow = true;
			hit = true;
			Cube1posX = 10;
		}
		index1 = (rand() % 3) + 1;
		index2 = index1;
		while (index1 == index2) {
			index2 = (rand() % 3) + 1;
		}
		index3 = index2;
		while (index3 == index2 || index3 == index1) {
			index3 = (rand() % 3) + 1;
		}

		if (!collision_red) {
			if (posZ + 6 <= Cube1posZ - threshold) {
				Cube1posZ -= 60;
				Cube1posX = CubePos[index1];


			}
			if (posZ + 6 <= Cube2posZ - threshold) {
				Cube2posZ -= 60;
				Cube2posX = CubePos[index2];


			}
			if (posZ + 6 <= Cube3posZ - threshold) {
				Cube3posZ -= 60;
				Cube3posX = CubePos[index3];

			}

		}
		posZ -= speed;
		QuadposZ -= speed;
		scaleZ += speed / 2.5;
		speed += 0.0003;


	}
	renderText("Score: " + to_string(score), 0, 455, 0.5, glm::vec3(1, 1, 1));
	assert(glGetError() == GL_NO_ERROR);
	if (!collision_red) score++;
}


void reshape(GLFWwindow* window, int w, int h)
{
	w = w < 1 ? 1 : w;
	h = h < 1 ? 1 : h;

	gWidth = w;
	gHeight = h;

	glViewport(0, 0, w, h);

	// Use perspective projection
	float fovyRad = (float)(90.0 / 180.0) * M_PI;

	projectionMatrix = glm::perspective(fovyRad, w / (float)h, 1.0f, 200.0f);


	// Assume default camera position and orientation (camera is at
	// (0, 0, 0) with looking at -z direction and its up vector pointing
	// at +y direction)
	// 
	//viewingMatrix = glm::mat4(1);

	viewingMatrix = glm::lookAt(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0) + glm::vec3(0, 0, -1), glm::vec3(0, 1, 0));
	//QuadviewingMatrix = glm::lookAt(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0) + glm::vec3(0, 0, -1), glm::vec3(0, 1, 0));

}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{

	if (key == GLFW_KEY_Q && action == GLFW_PRESS)
	{
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}
	else if (key == GLFW_KEY_G && action == GLFW_PRESS)
	{
		activeProgramIndex = 0;
	}
	else if (key == GLFW_KEY_P && action == GLFW_PRESS)
	{
		activeProgramIndex = 1;
	}
	else if (key == GLFW_KEY_F && action == GLFW_PRESS)
	{
		glShadeModel(GL_FLAT);
	}
	else if (key == GLFW_KEY_S && action == GLFW_PRESS)
	{
		glShadeModel(GL_SMOOTH);
	}

	else if (key == GLFW_KEY_E && action == GLFW_PRESS)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else if (key == GLFW_KEY_E && action == GLFW_PRESS)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else if (key == GLFW_KEY_A && action == GLFW_PRESS)
	{
		left_right_offset = -speed * 0.5;
		if (left_right_offset <= -0.35) {
			left_right_offset = -0.35;
		}
	}
	else if (key == GLFW_KEY_D && action == GLFW_PRESS)
	{
		left_right_offset = speed * 0.5;
		if (left_right_offset >= 0.35) {
			left_right_offset = 0.35;
		}
	}
	else if (action == GLFW_RELEASE) {
		left_right_offset = 0.0f;
		restart = false;

	}
	else if (key == GLFW_KEY_R && action == GLFW_PRESS)
	{
		restart = true;
		collision_red = false;
	}

}

void mainLoop(GLFWwindow* window)
{
	while (!glfwWindowShouldClose(window))
	{

		display();
		glfwSwapBuffers(window);
		glfwPollEvents();



	}
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
	GLFWwindow* window;
	if (!glfwInit())
	{
		exit(-1);
	}

	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this if on MacOS

	int width = 1000, height = 800;
	window = glfwCreateWindow(width, height, "Simple Example", NULL, NULL);

	if (!window)
	{
		glfwTerminate();
		exit(-1);
	}

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	// Initialize GLEW to setup the OpenGL Function pointers
	if (GLEW_OK != glewInit())
	{
		std::cout << "Failed to initialize GLEW" << std::endl;
		return EXIT_FAILURE;
	}

	char rendererInfo[512] = { 0 };
	strcpy_s(rendererInfo, (const char*)glGetString(GL_RENDERER)); // Use strcpy_s on Windows, strcpy on Linux
	strcat_s(rendererInfo, " - "); // Use strcpy_s on Windows, strcpy on Linux
	strcat_s(rendererInfo, (const char*)glGetString(GL_VERSION)); // Use strcpy_s on Windows, strcpy on Linux
	glfwSetWindowTitle(window, rendererInfo);


	initBunny();
	initQuad();
	initCube();
	glfwSetKeyCallback(window, keyboard);
	glfwSetWindowSizeCallback(window, reshape);

	reshape(window, width, height); // need to call this once ourselves
	mainLoop(window); // this does not return unless the window is closed

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
