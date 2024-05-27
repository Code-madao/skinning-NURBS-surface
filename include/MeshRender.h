#ifndef MESHRENDER_H
#define MESHRENDER_H

#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "Shader.h"
//#include "SurfacePlotter.h"
#include "Camera.h"
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
using namespace std;
class MeshRender {

private:
    Shader shader;

    // SurfacePlotter surfacePlotter;
    float* vertices;
    uint numElements;
    uint* indices;
    uint numIndices;

    // int mode;

public:
    uint VAO, VBO, EBO;
    MeshRender();

    // transformation matrices
    glm::mat4 getViewMatrix(Camera camera);
    glm::mat4 getProjectionMatrix(Camera camera, int widowWidth, int windowHeight);
    glm::mat4 getDefaultModelMatrix(void);

    uint generateBuffer(void);
    uint generateVAO(void);

    void generateVerticesIndices(vector<vector<glm::vec3>> Vertices);

    float* getVertices(void);
    uint getNumElements(void);
    uint* getIndices(void);
    uint getNumIndices(void);

    void Initial(const char* vertexPath, const char* fragmentPath, vector<vector<glm::vec3>> Vertices);
    void Draw(Camera camera, glm::mat4 modelMatrix, glm::vec3 lightPos, int windowWidth, int windowHeight);

};



#endif