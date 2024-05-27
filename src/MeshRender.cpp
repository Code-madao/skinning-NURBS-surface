#include "../include/MeshRender.h"
#include "glm/ext.hpp"
#include "../include/tinynurbs/tinynurbs.h"
MeshRender::MeshRender() : vertices(NULL), numElements(0), indices(NULL), numIndices(0), VAO(0), VBO(0), EBO(0) {
}

// transformation matrices
glm::mat4 MeshRender::getViewMatrix(Camera camera) {
    return camera.getViewMatrix();
}
glm::mat4 MeshRender::getProjectionMatrix(Camera camera, int windowWidth, int windowHeight) {
    return glm::perspective(glm::radians(camera.zoom), (float)windowWidth / (float)windowHeight, 0.1f, 99999.0f);
}
glm::mat4 MeshRender::getDefaultModelMatrix(void) {
    //return glm::mat4(1.0f); // identity
    return glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
}

void MeshRender::generateVerticesIndices(vector<vector<glm::vec3>> Vertices) {
    int numX = Vertices.size();
    int numY = Vertices[0].size();

    // vertices:
    // deallocate old data
    if (this->vertices)
        delete[] this->vertices;

    // allocate memory for new data
    this->numElements = 3 * numX * numY;
    this->vertices = new float[this->numElements];

    for (int x = 0; x < numX; x++)
    {
        for (int y = 0; y < numY; y++) {
            this->vertices[(x * numY + y) * 3 + 0] = Vertices[x][y].x; // x
            this->vertices[(x * numY + y) * 3 + 1] = Vertices[x][y].y; // y
            this->vertices[(x * numY + y) * 3 + 2] = Vertices[x][y].z; // z time-dependent
        }
    }

    // indices:
    // deallocte old data

    if (this->indices)
        delete[] this->indices;

    // determine number of indices
    this->numIndices = (numX * (numY - 1) + numY * (numX - 1)) * 2;

    // allocate memory for new data
    this->indices = new uint[this->numIndices];

    int i = 0;

    for (int x = 0; x < numX; ++x) {
        for (int y = 0; y < numY - 1; ++y) {
            this->indices[i++] = x * numY + y;
            this->indices[i++] = x * numY + y + 1;
        }
    }

    for (int y = 0; y < numY; ++y) {
        for (int x = 0; x < numX - 1; ++x) {
            this->indices[i++] = x * numY + y;
            this->indices[i++] = (x + 1) * numY + y;
        }
    }


}

float* MeshRender::getVertices(void) {
    return this->vertices;
}
uint MeshRender::getNumElements(void) {
    return this->numElements;
}
uint* MeshRender::getIndices(void) {
    return this->indices;
}
uint MeshRender::getNumIndices(void) {
    return this->numIndices;
}

uint MeshRender::generateBuffer(void) {
    uint buf;
    glGenBuffers(1, &buf);
    return buf;
}
uint MeshRender::generateVAO(void) {
    uint vao;
    glGenVertexArrays(1, &vao);
    return vao;
}

// initiate shader and set buffer data
void MeshRender::Initial(const char* vertexPath, const char* fragmentPath, vector<vector<glm::vec3>> Vertices) {

    //Initial shader
    this->shader = Shader(vertexPath, fragmentPath);

    //Generate default data
    generateVerticesIndices(Vertices);

    // generate surface plot VAO and VBO and EBO
    this->VAO = generateVAO();
    this->VBO = generateBuffer();
    this->EBO = generateBuffer();

    glBindVertexArray(this->VAO);

    // set VBO data
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO);
    glBufferData(GL_ARRAY_BUFFER, getNumElements() * sizeof(float), getVertices(), GL_DYNAMIC_DRAW);

    // set EBO data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, getNumIndices() * sizeof(uint), getIndices(), GL_DYNAMIC_DRAW);

    // vertices attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

}

void MeshRender::Draw(Camera camera, glm::mat4 modelMatrix, glm::vec3 lightPos, int windowWidth, int windowHeight) {

    glm::mat4 viewMatrix = getViewMatrix(camera);
    glm::mat4 projectionMatrix = getProjectionMatrix(camera, windowWidth, windowHeight);
    // int zRange = this->surfacePlotter.getZRange();

    this->shader.use();
    this->shader.setMat4Uniform("view", viewMatrix);
    this->shader.setMat4Uniform("projection", projectionMatrix);
    glm::mat4 model = getDefaultModelMatrix() * modelMatrix;
    this->shader.setMat4Uniform("model", model);

    this->shader.use();
    glBindVertexArray(this->VAO);
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO);
    glBufferData(GL_ARRAY_BUFFER, getNumElements() * sizeof(float), getVertices(), GL_STATIC_DRAW);
    glDrawElements(GL_LINES, getNumIndices(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


}