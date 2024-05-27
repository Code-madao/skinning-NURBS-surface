#include "../include/CurveRender.h"
#include "glm/ext.hpp"
#include "../include/tinynurbs/tinynurbs.h"
CurveRender::CurveRender() : vertices(NULL), numElements(0), indices(NULL), numIndices(0), VAO(0), VBO(0), EBO(0) {
}

// transformation matrices
glm::mat4 CurveRender::getViewMatrix(Camera camera) {
    return camera.getViewMatrix();
}
glm::mat4 CurveRender::getProjectionMatrix(Camera camera, int windowWidth, int windowHeight) {
    return glm::perspective(glm::radians(camera.zoom), (float)windowWidth / (float)windowHeight, 0.1f, 99999.0f);
}
glm::mat4 CurveRender::getDefaultModelMatrix(void) {
    //return glm::mat4(1.0f); // identity
    return glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
}

void CurveRender::generateVerticesIndices(vector<glm::vec3> Vertices) {
    this->numElements = Vertices.size() * 3;
    // deallocate old data
    if (this->vertices)
        delete[] this->vertices;

    this->vertices = new float[this->numElements];

    for (int i = 0; i < this->numElements/3; i++)
    {
        this->vertices[i * 3 + 0] = Vertices[i].x;
        this->vertices[i * 3 + 1] = Vertices[i].y;
        this->vertices[i * 3 + 2] = Vertices[i].z;
    }

    // indices:
    // deallocte old data

    if (this->indices)
        delete[] this->indices;

    this->numIndices = (this->numElements/3 - 1) * 2;
    this->indices = new uint[this->numIndices];

    int j = 0;
    for (int i = 0; i < this->numIndices / 2; i++)
    {
        this->indices[j++] = i;
        this->indices[j++] = i + 1;
    }

}

float* CurveRender::getVertices(void) {
    return this->vertices;
}
uint CurveRender::getNumElements(void) {
    return this->numElements;
}
uint* CurveRender::getIndices(void) {
    return this->indices;
}
uint CurveRender::getNumIndices(void) {
    return this->numIndices;
}

uint CurveRender::generateBuffer(void) {
    uint buf;
    glGenBuffers(1, &buf);
    return buf;
}
uint CurveRender::generateVAO(void) {
    uint vao;
    glGenVertexArrays(1, &vao);
    return vao;
}

void CurveRender::Initial(const char* vertexPath, const char* fragmentPath, vector<glm::vec3> Vertices) {

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

void CurveRender::Draw(Camera camera, glm::mat4 modelMatrix, glm::vec3 lightPos, int windowWidth, int windowHeight) {

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