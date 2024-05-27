#include "../include/SurfaceRender.h"
#include "glm/ext.hpp"
#include "../include/tinynurbs/tinynurbs.h"
SurfaceRender::SurfaceRender(): vertices(NULL), numElements(0), indices(NULL), numIndices(0), VAO(0), VBO(0), EBO(0){


}

// transformation matrices
glm::mat4 SurfaceRender::getViewMatrix(Camera camera) {
    return camera.getViewMatrix();
}
glm::mat4 SurfaceRender::getProjectionMatrix(Camera camera, int windowWidth, int windowHeight) {
    return glm::perspective(glm::radians(camera.zoom), (float)windowWidth / (float)windowHeight, 0.1f, 99999.0f);
}
glm::mat4 SurfaceRender::getDefaultModelMatrix(void) {
    //return glm::mat4(1.0f); // identity
    return glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
}


void SurfaceRender::generateVerticesIndices(vector<vector<glm::vec3>> Vertices, vector<vector<glm::vec3>> Normal)
{
    int numX = Vertices.size();
    int numY = Vertices[0].size();
        
    // vertices:
    // deallocate old data
    if (this->vertices)
        delete[] this->vertices;

    // allocate memory for new data
    this->numElements = 6 * numX * numY;
    this->vertices = new float[this->numElements];

    for (int x = 0; x < numX; x++)
    {
        for (int y = 0; y < numY; y++) {
            this->vertices[(x * numY + y) * 6 + 0] = Vertices[x][y].x; // x
            this->vertices[(x * numY + y) * 6 + 1] = Vertices[x][y].y; // y
            this->vertices[(x * numY + y) * 6 + 2] = Vertices[x][y].z; // z time-dependent

            this->vertices[(x * numY + y) * 6 + 3] = Normal[x][y].x; // x
            this->vertices[(x * numY + y) * 6 + 4] = Normal[x][y].y; // y
            this->vertices[(x * numY + y) * 6 + 5] = Normal[x][y].z; // z time-dependent
        }
    }

    /*
    for (int x = 0; x < numX; x++) {
        for (int y = 0; y < numY; y++) {
            glm::vec3 a, b, c, u1, u2, n;
            glm::vec3 n1, n2, n3, n4;

            if (x < numX - 1 && y < numY - 1) {
                a.x = this->vertices[(x * numY + y) * 6 + 0]; // x
                a.y = this->vertices[(x * numY + y) * 6 + 1]; // y
                a.z = this->vertices[(x * numY + y) * 6 + 2]; // z time-dependent

                b.x = this->vertices[(x * numY + y + 1) * 6 + 0]; // x
                b.y = this->vertices[(x * numY + y + 1) * 6 + 1]; // y
                b.z = this->vertices[(x * numY + y + 1) * 6 + 2]; // z time-dependent

                c.x = this->vertices[((x + 1) * numY + y) * 6 + 0]; // x
                c.y = this->vertices[((x + 1) * numY + y) * 6 + 1]; // y
                c.z = this->vertices[((x + 1) * numY + y) * 6 + 2]; // z time-dependent

                u1.x = b.x - a.x;
                u1.y = b.y - a.y;
                u1.z = b.z - a.z;

                u2.x = c.x - a.x;
                u2.y = c.y - a.y;
                u2.z = c.z - a.z;
            }
            else if (x == numX - 1 && y != numY - 1)
            {

                a.x = this->vertices[(x * numY + y) * 6 + 0]; // x
                a.y = this->vertices[(x * numY + y) * 6 + 1]; // y
                a.z = this->vertices[(x * numY + y) * 6 + 2]; // z time-dependent

                b.x = this->vertices[(x * numY + y + 1) * 6 + 0]; // x
                b.y = this->vertices[(x * numY + y + 1) * 6 + 1]; // y
                b.z = this->vertices[(x * numY + y + 1) * 6 + 2]; // z time-dependent

                c.x = this->vertices[((x - 1) * numY + y) * 6 + 0]; // x
                c.y = this->vertices[((x - 1) * numY + y) * 6 + 1]; // y
                c.z = this->vertices[((x - 1) * numY + y) * 6 + 2]; // z time-dependent

                u1.x = b.x - a.x;
                u1.y = b.y - a.y;
                u1.z = b.z - a.z;

                u2.x = a.x - c.x;
                u2.y = a.y - c.y;
                u2.z = a.z - c.z;
            }
            else if (x != numX - 1 && y == numY - 1) {
                a.x = this->vertices[(x * numY + y) * 6 + 0]; // x
                a.y = this->vertices[(x * numY + y) * 6 + 1]; // y
                a.z = this->vertices[(x * numY + y) * 6 + 2]; // z time-dependent

                b.x = this->vertices[(x * numY + y - 1) * 6 + 0]; // x
                b.y = this->vertices[(x * numY + y - 1) * 6 + 1]; // y
                b.z = this->vertices[(x * numY + y - 1) * 6 + 2]; // z time-dependent

                c.x = this->vertices[((x + 1) * numY + y) * 6 + 0]; // x
                c.y = this->vertices[((x + 1) * numY + y) * 6 + 1]; // y
                c.z = this->vertices[((x + 1) * numY + y) * 6 + 2]; // z time-dependent

                u1.x = a.x - b.x;
                u1.y = a.y - b.y;
                u1.z = a.z - b.z;

                u2.x = c.x - a.x;
                u2.y = c.y - a.y;
                u2.z = c.z - a.z;
            }
            else if (x == numX - 1 && y == numY - 1) {
                a.x = this->vertices[(x * numY + y) * 6 + 0]; // x
                a.y = this->vertices[(x * numY + y) * 6 + 1]; // y
                a.z = this->vertices[(x * numY + y) * 6 + 2]; // z time-dependent

                b.x = this->vertices[(x * numY + y - 1) * 6 + 0]; // x
                b.y = this->vertices[(x * numY + y - 1) * 6 + 1]; // y
                b.z = this->vertices[(x * numY + y - 1) * 6 + 2]; // z time-dependent

                c.x = this->vertices[((x - 1) * numY + y) * 6 + 0]; // x
                c.y = this->vertices[((x - 1) * numY + y) * 6 + 1]; // y
                c.z = this->vertices[((x - 1) * numY + y) * 6 + 2]; // z time-dependent

                u1.x = a.x - b.x;
                u1.y = a.y - b.y;
                u1.z = a.z - b.z;

                u2.x = a.x - c.x;
                u2.y = a.y - c.y;
                u2.z = a.z - c.z;
            }

            n.x = u1.y * u2.z - u2.y * u1.z;
            n.y = u2.x * u1.z - u2.z * u1.x;
            n.z = u1.x * u2.y - u2.x * u1.y;

            // add vertex
            this->vertices[(x * numY + y) * 6 + 3] = n.x; // x
            this->vertices[(x * numY + y) * 6 + 4] = n.y; // y
            this->vertices[(x * numY + y) * 6 + 5] = n.z; // z time-dependent

        }

    }
    */

    // indices:
    // deallocte old data
    if (this->indices)
        delete[] this->indices;

    this->numIndices = (numX - 1) * (numY - 1) * 2 * 3;

    // allocate memory for new data
    this->indices = new uint[this->numIndices];

    int i = 0;
    for (int x = 0; x < numX - 1; ++x) {
        for (int y = 0; y < numY - 1; ++y) {
            // ABC
            this->indices[i++] = x * (numY) + y;
            this->indices[i++] = x * (numY) + (y + 1);
            this->indices[i++] = (x + 1) * (numY ) + (y + 1);
            // CDA
            this->indices[i++] = (x + 1) * (numY) + (y + 1);
            this->indices[i++] = (x + 1) * (numY) + y;
            this->indices[i++] = x * (numY) + y;
        }
    }
}

float* SurfaceRender::getVertices(void) {
    return this->vertices;
}
uint SurfaceRender::getNumElements(void) {
    return this->numElements;
}
uint* SurfaceRender::getIndices(void) {
    return this->indices;
}
uint SurfaceRender::getNumIndices(void) {
    return this->numIndices;
}

uint SurfaceRender::generateBuffer(void) {
    uint buf;
    glGenBuffers(1, &buf);
    return buf;
}
uint SurfaceRender::generateVAO(void) {
    uint vao;
    glGenVertexArrays(1, &vao);
    return vao;
}

// initiate shader and set buffer data
void SurfaceRender::Initial(const char* vertexPath, const char* fragmentPath, vector<vector<glm::vec3>> Vertices, vector<vector<glm::vec3>> Normal)
{
	//Initial shader
	this->shader = Shader(vertexPath, fragmentPath);

	//Generate default data
	generateVerticesIndices(Vertices, Normal);

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
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //normal attributes
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

void SurfaceRender::Draw(Camera camera, glm::mat4 modelMatrix, glm::vec3 lightPos, int windowWidth, int windowHeight)
{

    glm::mat4 viewMatrix = getViewMatrix(camera);
    glm::mat4 projectionMatrix = getProjectionMatrix(camera, windowWidth, windowHeight);
   // int zRange = this->surfacePlotter.getZRange();
    this->shader.use();

    this->shader.setVec3Uniform("objectColor", glm::vec3(1.0f, 1.0f, 0.7f));
    this->shader.setVec3Uniform("lightColor", glm::vec3(0.7f, 0.7f, 0.7f));
    this->shader.setVec3Uniform("lightPos", lightPos);
    this->shader.setVec3Uniform("viewPos", camera.position);
   // this->shader.setFloatUniform("zRange", (zRange == 0) ? 1.0f : zRange);
   // this->shader.setFloatUniform("zMin", this->surfacePlotter.getZMin());
    this->shader.setMat4Uniform("view", viewMatrix);
    this->shader.setMat4Uniform("projection", projectionMatrix);
    this->shader.setMat4Uniform("model", getDefaultModelMatrix() * modelMatrix);

    this->shader.use();
    glBindVertexArray(this->VAO);
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO);
    glBufferData(GL_ARRAY_BUFFER, getNumElements() * sizeof(float), getVertices(), GL_STATIC_DRAW);
    glDrawElements(GL_TRIANGLES, getNumIndices(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}