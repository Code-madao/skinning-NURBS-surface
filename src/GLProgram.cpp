#include "../include/GLProgram.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "../include/mynurbs.h"

#define MYVERSION

const char* vertexShaderPath = "../shaders/vertexShader.vs";
const char* whiteFragmentShaderPath = "../shaders/whiteFragmentShader.fs";
const char* redFragmentShaderPath = "../shaders/redFragmentShader.fs";
const char* greenFragmentShaderPath = "../shaders/greenFragmentShader.fs";
const char* blueFragmentShaderPath = "../shaders/blueFragmentShader.fs";

const char* lightvertexShaderPath = "../shaders/lightvertexShader.vs";
const char* fragmentShaderPath = "../shaders/fragmentShader.fs";

const char* imguiFontPath = "../imgui/misc/fonts/DroidSans.ttf";

void normalizeVec3(glm::vec3& input)
{
    double length = sqrt(pow(input.x, 2) + pow(input.y, 2) + pow(input.z, 2));
    input.x = input.x / length;
    input.y = input.y / length;
    input.z = input.z / length;
}

// 计算沿vec方向延伸一定距离的点  len: 延伸参数
glm::vec3 compExtendedPnt(const glm::vec3& point, const glm::vec3& vec, float len)
{
    glm::vec3 res; //保存延伸后的数据点
    res.x = point.x - vec.x * len;
    res.y = point.y - vec.y * len;
    res.z = point.z - vec.z * len;
    return res;
}

GLProgram::GLProgram() :
    deltaTime(0.0f), prevTime(0.0f) {
    surfaceRender.clear();
    curveRender.clear();
    meshRender.clear();
    normalRender.clear();
}

void GLProgram::init(vector < tinynurbs::RationalCurve<double>> curves, vector<tinynurbs::RationalSurface<double>> surfaces) {

    // lighting
    lightPos = glm::vec3(10.f, 4.0f, 10.0f);

    // initialize window system
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    this->window = glfwCreateWindow(this->windowWidth,this-> windowHeight, "3D Surface Plotter", NULL, NULL);
    if (this->window == NULL)
    {
        std::cout << "FAILED TO CREATE GLFW WINDOW" << std::endl;
        glfwTerminate();
        exit(-1);
    }

    glfwSetWindowUserPointer(window, this);

    glfwMakeContextCurrent(this->window);
    glfwSetFramebufferSizeCallback(this->window, _framebufferSizeCallback);
    glfwSetScrollCallback(window, _scrollCallback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR);
    glfwSetMouseButtonCallback(window, _mouseButtonCallback);
    glfwSetCursorPosCallback(window, _cursorPosCallback);

    // initialize GLAD before making OpenGL calls
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "FAILED TO INITIALIZE GLAD" << std::endl;
        exit(-1);
    }

    // GL calls
    glViewport(0, 0, this->windowWidth, this->windowHeight);
    glEnable(GL_DEPTH_TEST);

    // surfaces init
    this->surfaceRender.resize(surfaces.size());
    this->meshRender.resize(surfaceRender.size());
    this->normalRender.resize(surfaceRender.size());
    this->surfaceControlRender.resize(surfaceRender.size());
    this->UsurfaceDerivateRender.resize(surfaceRender.size());
    this->VsurfaceDerivateRender.resize(surfaceRender.size());
    for (int k = 0; k < this->surfaceRender.size(); k++)
    {
        double u, v, delta;
        u = 0;
        v = 0;
        int numX = 40;
        int numY = 40;
        vector<vector<glm::vec3>> Vertices;
        vector<vector<glm::vec3>> Normal;
        vector<vector<glm::vec3>> OffsetVertex;

        /********surface and surface mesh render init********/
        Normal.resize(numX);
        Vertices.resize(numX);
        OffsetVertex.resize(numX);
        delta = 1.0 / (numX - 1);
        for (int x = 0; x < numX; x++) 
        {
            v = 0;
            Vertices[x].resize(numY);
            OffsetVertex[x].resize(numY);
            Normal[x].resize(numY);

            for (int y = 0; y < numY; y++) 
            {
                // add vertex
    #ifdef MYVERSION
                glm::vec3 tmp=myNurbs::myRationalSurfPoint(surfaces[k],u,v);
                glm::vec3 tmp1=myNurbs::myRationalSurfNormal(surfaces[k],u,v);
    #else 
                glm::vec3 tmp = tinynurbs::surfacePoint(surfaces[k], u, v);
                glm::vec3 tmp1 = tinynurbs::surfaceNormal(surfaces[k], u, v);
    #endif
                normalizeVec3(tmp1);
                glm::vec3 tmp2 = compExtendedPnt(tmp, tmp1, 0.2f); //保存延伸后的数据点

                Vertices[x][y] = tmp;
                OffsetVertex[x][y] = tmp2; // not be used
                Normal[x][y] = tmp1;
                v += delta;
            }
            u += delta;
        }
        this->surfaceRender[k].Initial(lightvertexShaderPath, fragmentShaderPath, Vertices, Normal);
        this->meshRender[k].Initial(vertexShaderPath, whiteFragmentShaderPath, Vertices);
 
        /********surface control points render init********/
        vector<vector<glm::vec3>> MeshVertex;
        MeshVertex.resize(surfaces[k].control_points.rows());
        for (int i = 0; i < surfaces[k].control_points.rows(); i++)
        {
            MeshVertex[i].resize(surfaces[k].control_points.cols());
            for (int j = 0; j < surfaces[k].control_points.cols(); j++)
            {
                MeshVertex[i][j] = surfaces[k].control_points[surfaces[k].control_points.cols() * i + j];
            }

        }
        this->surfaceControlRender[k].Initial(vertexShaderPath, whiteFragmentShaderPath, MeshVertex);

        Vertices.clear();
        OffsetVertex.clear();
        Vertices.resize(1);
        OffsetVertex.resize(1);
        Vertices[0].resize(1);
        OffsetVertex[0].resize(1);

        /********derivatives extension line of one point on surface render init********/
        double up, vp;
        up = 0.5;
        vp = 0.5;
#ifdef MYVERSION
        glm::vec3 tmp = myNurbs::myRationalSurfPoint(surfaces[k],up,vp);
        tinynurbs::array2<glm::vec<3,double>> res = myNurbs::myRationalSurfDerivative(surfaces[k],1,up,vp);
#else 
        glm::vec3 tmp = tinynurbs::surfacePoint(surfaces[k], up, vp);
        tinynurbs::array2<glm::vec<3,float>> res = tinynurbs::surfaceDerivatives(surfaces[k], 1, up, vp);
#endif

        glm::vec3 tmp1 = res(1, 0);  //u方向导数
        normalizeVec3(tmp1);

        glm::vec3 tmp2 = compExtendedPnt(tmp, tmp1, 0.2f); //保存延伸后的数据点

        Vertices[0][0] = tmp;
        OffsetVertex[0][0] = tmp2;

        this->UsurfaceDerivateRender[k].Initial(vertexShaderPath, blueFragmentShaderPath, Vertices, OffsetVertex);

        tmp1 = res(0, 1);  //v方向导数
        normalizeVec3(tmp1);
        tmp2=compExtendedPnt(tmp, tmp1, 0.2f);

        Vertices[0][0] = tmp;
        OffsetVertex[0][0] = tmp2;
        this->VsurfaceDerivateRender[k].Initial(vertexShaderPath, redFragmentShaderPath, Vertices, OffsetVertex);

        /********normal extension line of one point on surface render init********/
    #ifdef MYVERSION
        tmp1 = myNurbs::myRationalSurfNormal(surfaces[k],up,vp);
    #else 
        tmp1 = tinynurbs::surfaceNormal(surfaces[k], up, vp);
    #endif
        normalizeVec3(tmp1);
        tmp2=compExtendedPnt(tmp, tmp1, 0.2f);

        Vertices[0][0] = tmp;
        OffsetVertex[0][0] = tmp2;
        this->normalRender[k].Initial(vertexShaderPath, greenFragmentShaderPath, Vertices, OffsetVertex);
    }

    // curves init
    this->curveRender.resize(curves.size());
    this->curveDerivateRender.resize(curves.size());
    this->curveControlPointLineRender.resize(curves.size());
    for (int k = 0; k < this->curveRender.size(); k++)
    {
        double u, delta;
        vector<glm::vec3> crvvertex;
        crvvertex.resize(1000);
        u = 0;
        delta = 1.0 / (3000 / 3 - 1); 
        for (int i = 0; i < 3000 / 3; i++)
        {
    #ifdef MYVERSION
            glm::vec3 tmp = myNurbs::myRationalCurvePoint(curves[k], u);
    #else 
            glm::vec3 tmp = tinynurbs::curvePoint(curves[k], u);
    #endif
            crvvertex[i] = tmp;
            u += delta;
        }
        this->curveRender[k].Initial(vertexShaderPath, whiteFragmentShaderPath, crvvertex);

        vector<glm::vec3> crvcontrol_points;
        for (int i = 0; i < curves[k].control_points.size(); i++)
        {
            crvcontrol_points.push_back(curves[k].control_points[i]);
        }
        this->curveControlPointLineRender[k].Initial(vertexShaderPath, greenFragmentShaderPath, crvcontrol_points);

        double p = 0.5;
        vector<vector<glm::vec3>> Vertices;
        vector<vector<glm::vec3>> OffsetVertex;
    #ifdef MYVERSION
        glm::vec3 tmp = myNurbs::myRationalCurvePoint(curves[k], p);
        std::vector<glm::vec<3,double>> derivateData = myNurbs::myRationalCurveDerivative(curves[k], 1, p);
    #else 
        glm::vec3 tmp = tinynurbs::curvePoint(curves[k], p);
        std::vector<glm::vec3> derivateData = tinynurbs::curveDerivatives(curves[k], 1, p);
    #endif
        glm::vec3 tmp1 = derivateData[1];
        normalizeVec3(tmp1);
        glm::vec3 tmp2 = compExtendedPnt(tmp, tmp1, 1.f); //保存延伸后的数据点

        Vertices.resize(1);
        OffsetVertex.resize(1);
        Vertices[0].push_back(tmp);
        OffsetVertex[0].push_back(tmp2);
        this->curveDerivateRender[k].Initial(vertexShaderPath, redFragmentShaderPath, Vertices, OffsetVertex);
    }

    // initialize IMGUI
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
    ImGui::StyleColorsLight();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 450");
    io.Fonts->AddFontFromFileTTF(imguiFontPath, 32.0f);
}

void GLProgram::run(vector<tinynurbs::RationalCurve<double>> curves, vector<tinynurbs::RationalSurface<double>> surfaces)
{
    ImGuiIO& io = ImGui::GetIO();
    // imgui state
    bool show_grid = false;
    bool show_tangent = false;
    bool show_normal = false;
    float tangent_u = 0.5f, tangent_v = 0.5f;
    bool show_control_mesh = false;
    bool show_shade_suf = false;
    bool show_curves = false;
    bool show_curves_derivate = false;
    float curves_derivate_u = 0.5f;
    bool show_curve_control_points = false;
    bool show_knot_insertion = false;
    bool show_knot_removal = false;
    bool show_degree_elevation = false;
    bool show_degree_reduction = false;
 
    int taskIndex = 1;

    // main loop
    while (!glfwWindowShouldClose(this->window))
    {
        // per-frame time logic
        float currTime = glfwGetTime();
        this->deltaTime = currTime - this->prevTime;
        this->prevTime = currTime;

        // input
        enableControl = (!io.WantCaptureMouse && !io.WantCaptureKeyboard);

        processInput();

        glClearColor(this->clearColor.r, this->clearColor.g, this->clearColor.b, this->clearColor.alpha);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // --- imgui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Settings");
        ImGui::InputInt("task index(1-6)", &taskIndex);
        ImGui::Checkbox("show grid", &show_grid);
        ImGui::Checkbox("show tangent", &show_tangent);
        if (show_tangent || show_normal)
        {
            ImGui::SliderFloat("u##1", &tangent_u, 0.0f, 1.0f);
            ImGui::SliderFloat("v##1", &tangent_v, 0.0f, 1.0f);
        }
        ImGui::Checkbox("show normal", &show_normal);
        ImGui::Checkbox("show shade surface", &show_shade_suf);
        ImGui::Checkbox("show control mesh", &show_control_mesh);
        ImGui::Checkbox("show curves", &show_curves);
        ImGui::Checkbox("show curves derivate", &show_curves_derivate);
        ImGui::Checkbox("show curve control points line", &show_curve_control_points);
        if (show_curves_derivate)
        {
            ImGui::SliderFloat("u##2", &curves_derivate_u, 0.0f, 1.0f);
        }

        /*
        ImGui::Checkbox("show knot insertion", &show_knot_insertion);
        ImGui::Checkbox("show degree elevation", &show_degree_elevation);
        */
        ImGui::End();
 
        //Draw Surfaces
        for (int k = 0; k < this->surfaceRender.size(); k++)
        {
            if (k != taskIndex-1)continue;

            vector<vector<glm::vec3>> Vertices;
            vector<vector<glm::vec3>> OffsetVertex;
            Vertices.resize(1);
            OffsetVertex.resize(1);
            Vertices[0].resize(1);
            OffsetVertex[0].resize(1);

            double up, vp;
            up = tangent_u;
            vp = tangent_v;
    #ifdef MYVERSION
            glm::vec3 tmp = myNurbs::myRationalSurfPoint(surfaces[k], up, vp);
            tinynurbs::array2<glm::vec<3, double>> res = myNurbs::myRationalSurfDerivative(surfaces[k], 1, up, vp);
    #else
            glm::vec3 tmp = tinynurbs::surfacePoint(surfaces[k], up, vp);
            tinynurbs::array2<glm::vec<3, float>> res = tinynurbs::surfaceDerivatives(surfaces[k], 1, up, vp);
    #endif

            glm::vec3 tmp1 = res(1, 0);//derivative in u direction
            normalizeVec3(tmp1);
            glm::vec3 tmp2 = compExtendedPnt(tmp, tmp1, 1.f); //保存延伸后的数据点
            Vertices[0][0] = tmp;
            OffsetVertex[0][0] = tmp2;
            this->UsurfaceDerivateRender[k].Initial(vertexShaderPath, blueFragmentShaderPath, Vertices, OffsetVertex);

            tmp1 = res(0, 1);//derivative in v direction
            normalizeVec3(tmp1);
            tmp2 = compExtendedPnt(tmp, tmp1, 1.f);
            Vertices[0][0] = tmp;
            OffsetVertex[0][0] = tmp2;
            this->VsurfaceDerivateRender[k].Initial(vertexShaderPath, redFragmentShaderPath, Vertices, OffsetVertex);

            if (show_tangent)
            {
                this->UsurfaceDerivateRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
                this->VsurfaceDerivateRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
            }
    #ifdef MYVERSION
            tmp1 = myNurbs::myRationalSurfNormal(surfaces[k], up, vp);
    #else
            tmp1 = tinynurbs::surfaceNormal(surfaces[k], up, vp);
    #endif
            normalizeVec3(tmp1);
            tmp2 = compExtendedPnt(tmp, tmp1, -1.f);
            Vertices[0][0] = tmp;
            OffsetVertex[0][0] = tmp2;
            this->normalRender[k].Initial(vertexShaderPath, greenFragmentShaderPath, Vertices, OffsetVertex);
            if (show_normal)
                this->normalRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);

            if (show_shade_suf)
                this->surfaceRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        }

        //Draw Curves
        for (int k = 0; k < this->curveRender.size(); k++)
        {
           if (k / 5 != taskIndex - 1)continue;

            if(show_curves)
                this->curveRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
            if(show_curve_control_points)
                this->curveControlPointLineRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
            if (show_curves_derivate) {
                double p = curves_derivate_u;
                vector<vector<glm::vec3>> Vertices;
                vector<vector<glm::vec3>> OffsetVertex;
    #ifdef MYVERSION
                glm::vec3 tmp = myNurbs::myRationalCurvePoint(curves[k], p);
                std::vector<glm::vec<3, double>> derivateData = myNurbs::myRationalCurveDerivative(curves[k], 1, p);
    #else
                glm::vec3 tmp = tinynurbs::curvePoint(curves[k], p);
                std::vector<glm::vec3> derivateData = tinynurbs::curveDerivatives(curves[k], 1, p);
    #endif
                glm::vec3 tmp1 = derivateData[1];
                normalizeVec3(tmp1);
                glm::vec3 tmp2 = compExtendedPnt(tmp, tmp1, 1.f); //保存延伸后的数据点

                Vertices.resize(1);
                OffsetVertex.resize(1);
                Vertices[0].push_back(tmp);
                OffsetVertex[0].push_back(tmp2);

                this->curveDerivateRender[k].Initial(vertexShaderPath, redFragmentShaderPath, Vertices, OffsetVertex);
                this->curveDerivateRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
            }
        }
        //Draw Mesh
        for (int k = 0; k < this->surfaceRender.size(); k++)
        {
            if (k != taskIndex - 1)continue;

            if (show_grid)
                this->meshRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
            if(show_control_mesh)
                this->surfaceControlRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // check and call events and swap buffers
        glfwSwapBuffers((this->window));
        glfwPollEvents();
    }
}


void GLProgram::cleanup(void)
{
    // clean up gl resources
    for (int k = 0; k < this->surfaceRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(this->surfaceRender[k].VAO));
        glDeleteBuffers(1, &(this->surfaceRender[k].VBO));
        glDeleteBuffers(1, &this->surfaceRender[k].EBO);

        glDeleteVertexArrays(1, &(this->VsurfaceDerivateRender[k].VAO));
        glDeleteBuffers(1, &(this->VsurfaceDerivateRender[k].VBO));
        glDeleteBuffers(1, &this->VsurfaceDerivateRender[k].EBO);

        glDeleteVertexArrays(1, &(this->UsurfaceDerivateRender[k].VAO));
        glDeleteBuffers(1, &(this->UsurfaceDerivateRender[k].VBO));
        glDeleteBuffers(1, &this->UsurfaceDerivateRender[k].EBO);
    }

    for (int k = 0; k < this->curveRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(this->curveRender[k].VAO));
        glDeleteBuffers(1, &(this->curveRender[k].VBO));
        glDeleteBuffers(1, &this->curveRender[k].EBO);
    }

    for (int k = 0; k < this->surfaceRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(this->meshRender[k].VAO));
        glDeleteBuffers(1, &(this->meshRender[k].VBO));
        glDeleteBuffers(1, &this->meshRender[k].EBO);
    }

    for (int k = 0; k < this->normalRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(this->normalRender[k].VAO));
        glDeleteBuffers(1, &(this->normalRender[k].VBO));
        glDeleteBuffers(1, &this->normalRender[k].EBO);
    }
    // clean up glfw
    glfwTerminate();
}

void GLProgram::setClearColor(float r, float g, float b, float alpha) {
    this->clearColor = {r, g, b, alpha};
}

void GLProgram::framebufferSizeCallback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
    windowWidth = width;
    windowHeight = height;
}

void GLProgram::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if (this->enableControl) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            if (action == GLFW_PRESS) {
                glfwGetCursorPos(window, &prevMouseX, &prevMouseY);
                mousePressed = true;
            }
            else if (action == GLFW_RELEASE) {
                mousePressed = false;
            }
        }
    }
}

void GLProgram::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    if(this->enableControl)
        camera.processMouseScroll(yoffset);
}

void GLProgram::cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
    if (mousePressed && enableControl)
    {
        // get current cursor coordinates
        double currMouseX, currMouseY;
        glfwGetCursorPos(window, &currMouseX, &currMouseY);

        // get points on arcball
        glm::vec3 va = getArcballVector(prevMouseX, prevMouseY);
        glm::vec3 vb = getArcballVector(currMouseX, currMouseY);

        float speedFactor = 0.1f;
        float angleOfRotation = speedFactor * acos(MIN(1.0f, glm::dot(va, vb)));

        // to get the axis of rotation, need to convert from camera coordinates to world coordinates
        glm::vec3 axisCamera = glm::cross(va, vb);
        glm::mat3 cameraToModel = glm::inverse(glm::mat3(camera.getViewMatrix() * modelMatrix));
        glm::vec3 axisModel = cameraToModel * axisCamera;

        // update model rotation matrix
        float tolerance = 1e-4;
        if (angleOfRotation > tolerance)
            modelMatrix = glm::rotate(modelMatrix, glm::degrees(angleOfRotation), axisModel);

        // update cursor position
        prevMouseX = currMouseX;
        prevMouseY = currMouseY;
    }
}

glm::vec3 GLProgram::getArcballVector(float x, float y)
{
    // get normalized vector from center of the virtual arcball to a point P on the arcball's surface
    // if (x,y) is too far away from the arcball, return the nearest point on the arcball's surface
    glm::vec3 P(x/windowWidth * 2 - 1.0, y/windowHeight * 2 - 1.0, 0.0f);
    P.y = -P.y;

    float radius = 1.0f;
    float OP_squared = P.x * P.x + P.y * P.y;

    if (OP_squared <= radius)
        P.z = sqrt(radius - OP_squared); // apply pythagorean theorem to find z
    else
        P = glm::normalize(P); // nearest point

    return P;
}

void GLProgram::processInput(void)
{
    // close window with 'ESC' key
    if (glfwGetKey(this->window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(this->window, true);

    // camera movement
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS && enableControl)
        camera.processKeyboard(UP, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS && enableControl)
        camera.processKeyboard(DOWN, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS && enableControl)
        camera.processKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS && enableControl)
        camera.processKeyboard(RIGHT, deltaTime);
}

void apiOfCreateLoftSurface(const vector<tinynurbs::RationalCurve3d>& inputSections, tinynurbs::RationalSurface3d& result)
{
   myNurbs::createLoftSurface(inputSections, result);
}
