#include "../include/NormalRender.h"
#include "glm/ext.hpp"
#include "../include/tinynurbs/tinynurbs.h"
NormalRender::NormalRender() : vertices(NULL), numElements(0), indices(NULL), numIndices(0), VAO(0), VBO(0), EBO(0) {
}

// transformation matrices
glm::mat4 NormalRender::getViewMatrix(Camera camera) {
    return camera.getViewMatrix();
}
glm::mat4 NormalRender::getProjectionMatrix(Camera camera, int windowWidth, int windowHeight) {
    return glm::perspective(glm::radians(camera.zoom), (float)windowWidth / (float)windowHeight, 0.1f, 99999.0f);
}
glm::mat4 NormalRender::getDefaultModelMatrix(void) {
    //return glm::mat4(1.0f); // identity
    return glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
}

void NormalRender::generateVerticesIndices(vector<vector<glm::vec3>> Vertices, vector<vector<glm::vec3>> OffsetVertex) {

    int numX = Vertices.size();
    int numY = Vertices.size();

    this->numElements = numX * numY * 6;
    // deallocate old data
    if (this->vertices)
        delete[] this->vertices;

    this->vertices = new float[this->numElements];

    for (int x = 0; x < numX; x++)
    {
        for (int y = 0; y < numY; y++)
        {
            this->vertices[(x * numY + y) * 6 + 0] = Vertices[x][y].x;
            this->vertices[(x * numY + y) * 6 + 1] = Vertices[x][y].y;
            this->vertices[(x * numY + y) * 6 + 2] = Vertices[x][y].z;

            this->vertices[(x * numY + y) * 6 + 3] = OffsetVertex[x][y].x;
            this->vertices[(x * numY + y) * 6 + 4] = OffsetVertex[x][y].y;
            this->vertices[(x * numY + y) * 6 + 5] = OffsetVertex[x][y].z;
        }
    }
     
    // indices:
    // deallocte old data

    if (this->indices)
        delete[] this->indices;

    this->numIndices = this->numElements / 6 * 2;
    this->indices = new uint[this->numIndices];

    int j = 0;
    for (int i = 0; i < this->numIndices / 2; i++)
    {
        this->indices[j++] = i*2;
        this->indices[j++] = i*2 + 1;
    }

}

float* NormalRender::getVertices(void) {
    return this->vertices;
}
uint NormalRender::getNumElements(void) {
    return this->numElements;
}
uint* NormalRender::getIndices(void) {
    return this->indices;
}
uint NormalRender::getNumIndices(void) {
    return this->numIndices;
}

uint NormalRender::generateBuffer(void) {
    uint buf;
    glGenBuffers(1, &buf);
    return buf;
}
uint NormalRender::generateVAO(void) {
    uint vao;
    glGenVertexArrays(1, &vao);
    return vao;
}

void NormalRender::Initial(const char* vertexPath, const char* fragmentPath, vector<vector<glm::vec3>> Vertices, vector<vector<glm::vec3>> OffsetVertex) {

    //Initial shader
    this->shader = Shader(vertexPath, fragmentPath);

    //Generate default data
    generateVerticesIndices(Vertices, OffsetVertex);

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

void NormalRender::Draw(Camera camera, glm::mat4 modelMatrix, glm::vec3 lightPos, int windowWidth, int windowHeight) {

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