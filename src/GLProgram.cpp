#include "../include/GLProgram.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "../include/AppPaths.h"
#include "../include/AxisOverlay.h"
#include "../include/CurveDebugRender.h"
#include "../include/GeometryDebugUtils.h"
#include "../include/mynurbs.h"
#include "../include/RulerOverlay.h"
#include "../include/SettingsPanel.h"
#include "../include/SurfaceDebugRender.h"
#include "../include/ViewGizmo.h"
#include "GLProgram.h"

#define MYVERSION
#define AXIS_LENGTH 5

namespace
{
    glm::mat4 buildDefaultModelMatrix()
    {
        return glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    }

    glm::mat4 buildProjectionMatrix(const Camera &camera, int windowWidth, int windowHeight)
    {
        float aspect = static_cast<float>(windowWidth) / static_cast<float>(windowHeight > 0 ? windowHeight : 1);
        float halfHeight = camera.zoom * 0.5f;
        float halfWidth = halfHeight * aspect;
        return glm::ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, -99999.0f, 99999.0f);
    }
}

GLProgram::GLProgram() : deltaTime(0.0f), prevTime(0.0f)
{
    surfaceRender.clear();
    curveRender.clear();
    meshRender.clear();
    normalRender.clear();
}

void GLProgram::computeSceneBounds(const vector<tinynurbs::RationalSurface<double>> &surfaces)
{
    glm::vec3 bboxMin(0.0f);
    glm::vec3 bboxMax(0.0f);
    bool hasScenePoint = false;

    auto expandBBox = [&](const glm::vec3 &p)
    {
        if (!hasScenePoint)
        {
            bboxMin = p;
            bboxMax = p;
            hasScenePoint = true;
            return;
        }

        bboxMin = glm::min(bboxMin, p);
        bboxMax = glm::max(bboxMax, p);
    };

    // for (const auto &curve : curves)
    // {
    //     for (const auto &controlPoint : curve.control_points)
    //         expandBBox(glm::vec3(controlPoint));
    // }

    for (const auto &surface : surfaces)
    {
        for (int i = 0; i < surface.control_points.rows(); ++i)
        {
            for (int j = 0; j < surface.control_points.cols(); ++j)
                expandBBox(glm::vec3(surface.control_points(i, j)));
        }
    }

    sceneBBox = glm::mat2x3(bboxMin, bboxMax);
}

void GLProgram::initializeWindow(void)
{
    lightPos = glm::vec3(10.f, 4.0f, 10.0f);

    // initialize window system
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    // get primary monitor info to set window size and position (in monitor center)
    int windowPosX = 100;
    int windowPosY = 100;
    GLFWmonitor *primaryMonitor = glfwGetPrimaryMonitor();
    const GLFWvidmode *videoMode = primaryMonitor ? glfwGetVideoMode(primaryMonitor) : nullptr;
    if (primaryMonitor && videoMode)
    {
        int monitorPosX = 0;
        int monitorPosY = 0;
        glfwGetMonitorPos(primaryMonitor, &monitorPosX, &monitorPosY);

        windowWidth = videoMode->width / 2;
        windowHeight = videoMode->height / 2;
        windowPosX = monitorPosX + (videoMode->width - windowWidth) / 2;
        windowPosY = monitorPosY + (videoMode->height - windowHeight) / 2;
    }

    window = glfwCreateWindow(windowWidth, windowHeight, "3D Surface Plotter", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "FAILED TO CREATE GLFW WINDOW" << std::endl;
        glfwTerminate();
        exit(-1);
    }

    glfwSetWindowUserPointer(window, this);
    glfwSetWindowPos(window, windowPosX, windowPosY);

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, _framebufferSizeCallback);
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
    glViewport(0, 0, windowWidth, windowHeight);
    glEnable(GL_DEPTH_TEST);
}

void GLProgram::initializeAxes(void)
{
    float axisLength = AXIS_LENGTH;
    objectAxisLength = axisLength;
    vector<glm::vec3> xAxisVertices = {glm::vec3(0.0f), glm::vec3(axisLength, 0.0f, 0.0f)};
    vector<glm::vec3> yAxisVertices = {glm::vec3(0.0f), glm::vec3(0.0f, axisLength, 0.0f)};
    vector<glm::vec3> zAxisVertices = {glm::vec3(0.0f), glm::vec3(0.0f, 0.0f, axisLength)};

    // coord sys
    objectXAxisRender.Initial(AppPaths::kVertexShader, AppPaths::kRedFragmentShader, xAxisVertices);
    objectYAxisRender.Initial(AppPaths::kVertexShader, AppPaths::kGreenFragmentShader, yAxisVertices);
    objectZAxisRender.Initial(AppPaths::kVertexShader, AppPaths::kBlueFragmentShader, zAxisVertices);
}

void GLProgram::initializeSurfaceRenderers(const vector<tinynurbs::RationalSurface<double>> &surfaces)
{
    surfaceRender.resize(surfaces.size());
    meshRender.resize(surfaceRender.size());
    normalRender.resize(surfaceRender.size());
    surfaceControlRender.resize(surfaceRender.size());
    UsurfaceDerivateRender.resize(surfaceRender.size());
    VsurfaceDerivateRender.resize(surfaceRender.size());
    for (int k = 0; k < surfaceRender.size(); k++)
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
                glm::vec3 tmp = myNurbs::myRationalSurfPoint(surfaces[k], u, v);
                glm::vec3 tmp1 = myNurbs::myRationalSurfNormal(surfaces[k], u, v);
#else
                glm::vec3 tmp = tinynurbs::surfacePoint(surfaces[k], u, v);
                glm::vec3 tmp1 = tinynurbs::surfaceNormal(surfaces[k], u, v);
#endif
                GeometryDebugUtils::NormalizeVec3(tmp1);
                glm::vec3 tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 0.2f);

                Vertices[x][y] = tmp;
                OffsetVertex[x][y] = tmp2; // not be used
                Normal[x][y] = tmp1;
                v += delta;
            }
            u += delta;
        }
        surfaceRender[k].Initial(AppPaths::kLightVertexShader, AppPaths::kSurfaceFragmentShader, Vertices, Normal);
        meshRender[k].Initial(AppPaths::kVertexShader, AppPaths::kWhiteFragmentShader, Vertices);

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
        surfaceControlRender[k].Initial(AppPaths::kVertexShader, AppPaths::kWhiteFragmentShader, MeshVertex);

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
        glm::vec3 tmp = myNurbs::myRationalSurfPoint(surfaces[k], up, vp);
        tinynurbs::array2<glm::vec<3, double>> res = myNurbs::myRationalSurfDerivative(surfaces[k], 1, up, vp);
#else
        glm::vec3 tmp = tinynurbs::surfacePoint(surfaces[k], up, vp);
        tinynurbs::array2<glm::vec<3, float>> res = tinynurbs::surfaceDerivatives(surfaces[k], 1, up, vp);
#endif

        glm::vec3 tmp1 = res(1, 0); // u方向导数
        GeometryDebugUtils::NormalizeVec3(tmp1);

        glm::vec3 tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 0.2f); // 保存延伸后的数据点

        Vertices[0][0] = tmp;
        OffsetVertex[0][0] = tmp2;

        UsurfaceDerivateRender[k].Initial(AppPaths::kVertexShader, AppPaths::kBlueFragmentShader, Vertices, OffsetVertex);

        tmp1 = res(0, 1); // v方向导数
        GeometryDebugUtils::NormalizeVec3(tmp1);
        tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 0.2f);

        Vertices[0][0] = tmp;
        OffsetVertex[0][0] = tmp2;
        VsurfaceDerivateRender[k].Initial(AppPaths::kVertexShader, AppPaths::kRedFragmentShader, Vertices, OffsetVertex);

        /********normal extension line of one point on surface render init********/
#ifdef MYVERSION
        tmp1 = myNurbs::myRationalSurfNormal(surfaces[k], up, vp);
#else
        tmp1 = tinynurbs::surfaceNormal(surfaces[k], up, vp);
#endif
        GeometryDebugUtils::NormalizeVec3(tmp1);
        tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 0.2f);

        Vertices[0][0] = tmp;
        OffsetVertex[0][0] = tmp2;
        normalRender[k].Initial(AppPaths::kVertexShader, AppPaths::kGreenFragmentShader, Vertices, OffsetVertex);
    }
}

void GLProgram::initializeCurveRenderers(const vector<tinynurbs::RationalCurve<double>> &curves)
{
    curveRender.resize(curves.size());
    curveDerivateRender.resize(curves.size());
    curveControlPointLineRender.resize(curves.size());
    for (int k = 0; k < curveRender.size(); k++)
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
        curveRender[k].Initial(AppPaths::kVertexShader, AppPaths::kGreenFragmentShader, crvvertex);
        curveRender[k].setLineWidth(2.0f);
        curveRender[k].setDrawOnTop(true);

        vector<glm::vec3> crvcontrol_points;
        for (int i = 0; i < curves[k].control_points.size(); i++)
        {
            crvcontrol_points.push_back(curves[k].control_points[i]);
        }
        curveControlPointLineRender[k].Initial(AppPaths::kVertexShader, AppPaths::kGreenFragmentShader, crvcontrol_points);

        double p = 0.5;
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
        GeometryDebugUtils::NormalizeVec3(tmp1);
        glm::vec3 tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 1.f); // 保存延伸后的数据点

        Vertices.resize(1);
        OffsetVertex.resize(1);
        Vertices[0].push_back(tmp);
        OffsetVertex[0].push_back(tmp2);
        curveDerivateRender[k].Initial(AppPaths::kVertexShader, AppPaths::kRedFragmentShader, Vertices, OffsetVertex);
    }
}

void GLProgram::initializeImGui(void)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
    ImGui::StyleColorsLight();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 450");
    io.Fonts->AddFontFromFileTTF(AppPaths::kImGuiFont, 16.0f);
}

void GLProgram::init(vector<tinynurbs::RationalCurve<double>> curves, vector<tinynurbs::RationalSurface<double>> surfaces)
{
    computeSceneBounds(surfaces);
    initializeWindow();
    initializeAxes();
    initializeSurfaceRenderers(surfaces);
    initializeCurveRenderers(curves);
    initializeImGui();
}

void GLProgram::drawSurfaces(const vector<tinynurbs::RationalSurface<double>> &surfaces)
{
    for (int k = 0; k < surfaceRender.size(); k++)
    {
        // vector<vector<glm::vec3>> Vertices;
        // vector<vector<glm::vec3>> OffsetVertex;
        static vector<vector<glm::vec3>> Vertices;
        static vector<vector<glm::vec3>> OffsetVertex;
        Vertices.resize(1);
        OffsetVertex.resize(1);
        Vertices[0].resize(1);
        OffsetVertex[0].resize(1);

        double up = uiState.tangentU;
        double vp = uiState.tangentV;
#ifdef MYVERSION
        glm::vec3 tmp = myNurbs::myRationalSurfPoint(surfaces[k], up, vp);
        tinynurbs::array2<glm::vec<3, double>> res = myNurbs::myRationalSurfDerivative(surfaces[k], 1, up, vp);
#else
        glm::vec3 tmp = tinynurbs::surfacePoint(surfaces[k], up, vp);
        tinynurbs::array2<glm::vec<3, float>> res = tinynurbs::surfaceDerivatives(surfaces[k], 1, up, vp);
#endif

        glm::vec3 tmp1 = res(1, 0);
        GeometryDebugUtils::NormalizeVec3(tmp1);
        glm::vec3 tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 1.f);
        if (uiState.showTangent)
        {
            Vertices[0][0] = tmp;
            OffsetVertex[0][0] = tmp2;
            UsurfaceDerivateRender[k].updateData(Vertices, OffsetVertex);
        }

        tmp1 = res(0, 1);
        GeometryDebugUtils::NormalizeVec3(tmp1);
        tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 1.f);
        if (uiState.showTangent)
        {
            Vertices[0][0] = tmp;
            OffsetVertex[0][0] = tmp2;
            VsurfaceDerivateRender[k].updateData(Vertices, OffsetVertex);

            UsurfaceDerivateRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
            VsurfaceDerivateRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        }

        if (uiState.showNormal)
        {
#ifdef MYVERSION
            tmp1 = myNurbs::myRationalSurfNormal(surfaces[k], up, vp);
#else
            tmp1 = tinynurbs::surfaceNormal(surfaces[k], up, vp);
#endif
            GeometryDebugUtils::NormalizeVec3(tmp1);
            tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, -1.f);
            Vertices[0][0] = tmp;
            OffsetVertex[0][0] = tmp2;
            normalRender[k].updateData(Vertices, OffsetVertex);
            normalRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        }

        if (uiState.showShadeSurface)
            surfaceRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
    }
}

void GLProgram::drawCurves(const vector<tinynurbs::RationalCurve<double>> &curves)
{
    for (int k = 0; k < curveRender.size(); k++)
    {
        if (uiState.showCurves)
            curveRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        if (uiState.showCurveControlPoints)
            curveControlPointLineRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        if (uiState.showCurvesDerivate)
        {
            double p = uiState.curvesDerivateU;
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
            GeometryDebugUtils::NormalizeVec3(tmp1);
            glm::vec3 tmp2 = GeometryDebugUtils::ComputeExtendedPoint(tmp, tmp1, 1.f);

            Vertices.resize(1);
            OffsetVertex.resize(1);
            Vertices[0].push_back(tmp);
            OffsetVertex[0].push_back(tmp2);

            // curveDerivateRender[k].Initial(AppPaths::kVertexShader, AppPaths::kRedFragmentShader, Vertices, OffsetVertex);
            curveDerivateRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        }
    }
}

void GLProgram::drawMeshes(void)
{
    for (int k = 0; k < surfaceRender.size(); k++)
    {
        if (uiState.showGrid)
            meshRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        if (uiState.showControlMesh)
            surfaceControlRender[k].Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
    }
}

void GLProgram::drawOverlays(void)
{
    if (uiState.showObjectAxis)
    {
        objectXAxisRender.Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        objectYAxisRender.Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
        objectZAxisRender.Draw(camera, modelMatrix, lightPos, windowWidth, windowHeight);
    }

    if (uiState.showRuler || uiState.showObjectAxis || uiState.showViewGizmo)
    {
        glm::mat4 viewMatrix = camera.getViewMatrix();
        glm::mat4 projectionMatrix = buildProjectionMatrix(camera, windowWidth, windowHeight);
        glm::mat4 objectModelMatrix = buildDefaultModelMatrix() * modelMatrix;
        glm::vec3 xAxisEndWorld = glm::vec3(objectModelMatrix * glm::vec4(objectAxisLength, 0.0f, 0.0f, 1.0f));
        glm::vec3 yAxisEndWorld = glm::vec3(objectModelMatrix * glm::vec4(0.0f, objectAxisLength, 0.0f, 1.0f));
        glm::vec3 zAxisEndWorld = glm::vec3(objectModelMatrix * glm::vec4(0.0f, 0.0f, objectAxisLength, 1.0f));

        ImDrawList *overlay = ImGui::GetForegroundDrawList();

        if (uiState.showRuler)
            RulerOverlay::Draw(overlay, windowWidth, windowHeight, camera.zoom);

        if (uiState.showObjectAxis)
        {
            AxisOverlay::DrawAxisLabel(overlay, "X", xAxisEndWorld, viewMatrix, projectionMatrix, windowWidth, windowHeight, IM_COL32(220, 40, 40, 255));
            AxisOverlay::DrawAxisLabel(overlay, "Y", yAxisEndWorld, viewMatrix, projectionMatrix, windowWidth, windowHeight, IM_COL32(40, 170, 40, 255));
            AxisOverlay::DrawAxisLabel(overlay, "Z", zAxisEndWorld, viewMatrix, projectionMatrix, windowWidth, windowHeight, IM_COL32(40, 90, 220, 255));
        }

        if (uiState.showViewGizmo)
            ViewGizmo::Draw(overlay, viewMatrix, objectModelMatrix, windowWidth);
    }
}

void GLProgram::run(vector<tinynurbs::RationalCurve<double>> curves, vector<tinynurbs::RationalSurface<double>> surfaces)
{
    ImGuiIO &io = ImGui::GetIO();
    fittingToObjects();
    while (!glfwWindowShouldClose(window))
    {
        float currTime = glfwGetTime();
        deltaTime = currTime - prevTime;
        prevTime = currTime;

        enableControl = (!io.WantCaptureMouse && !io.WantCaptureKeyboard);

        processKeyboardInput();
        if (glfwWindowShouldClose(window))
        {
            break;
        }

        float lightOffset = glm::max(6.0f, objectAxisLength * 0.35f);
        lightPos = camera.position + camera.front * lightOffset + camera.up * lightOffset;

        glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.alpha);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // --- imgui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        SettingsPanel::Draw(uiState);
        drawSurfaces(surfaces);
        drawCurves(curves);
        drawMeshes();
        drawOverlays();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers((window));
        glfwPollEvents();
    }
}

void GLProgram::cleanup(void)
{
    // clean up gl resources
    for (int k = 0; k < surfaceRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(surfaceRender[k].VAO));
        glDeleteBuffers(1, &(surfaceRender[k].VBO));
        glDeleteBuffers(1, &surfaceRender[k].EBO);

        glDeleteVertexArrays(1, &(VsurfaceDerivateRender[k].VAO));
        glDeleteBuffers(1, &(VsurfaceDerivateRender[k].VBO));
        glDeleteBuffers(1, &VsurfaceDerivateRender[k].EBO);

        glDeleteVertexArrays(1, &(UsurfaceDerivateRender[k].VAO));
        glDeleteBuffers(1, &(UsurfaceDerivateRender[k].VBO));
        glDeleteBuffers(1, &UsurfaceDerivateRender[k].EBO);
    }

    for (int k = 0; k < curveRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(curveRender[k].VAO));
        glDeleteBuffers(1, &(curveRender[k].VBO));
        glDeleteBuffers(1, &curveRender[k].EBO);
    }

    for (int k = 0; k < surfaceRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(meshRender[k].VAO));
        glDeleteBuffers(1, &(meshRender[k].VBO));
        glDeleteBuffers(1, &meshRender[k].EBO);
    }

    for (int k = 0; k < normalRender.size(); k++)
    {
        glDeleteVertexArrays(1, &(normalRender[k].VAO));
        glDeleteBuffers(1, &(normalRender[k].VBO));
        glDeleteBuffers(1, &normalRender[k].EBO);
    }

    glDeleteVertexArrays(1, &(objectXAxisRender.VAO));
    glDeleteBuffers(1, &(objectXAxisRender.VBO));
    glDeleteBuffers(1, &objectXAxisRender.EBO);

    glDeleteVertexArrays(1, &(objectYAxisRender.VAO));
    glDeleteBuffers(1, &(objectYAxisRender.VBO));
    glDeleteBuffers(1, &objectYAxisRender.EBO);

    glDeleteVertexArrays(1, &(objectZAxisRender.VAO));
    glDeleteBuffers(1, &(objectZAxisRender.VBO));
    glDeleteBuffers(1, &objectZAxisRender.EBO);

    // clean up glfw
    glfwTerminate();
}

void GLProgram::setClearColor(float r, float g, float b, float alpha)
{
    clearColor = {r, g, b, alpha};
}

void GLProgram::framebufferSizeCallback(GLFWwindow *window, int width, int height)
{
    glViewport(0, 0, width, height);
    windowWidth = width;
    windowHeight = height;
}

void GLProgram::mouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
    if (enableControl)
    {
        if (button == GLFW_MOUSE_BUTTON_LEFT)
        {
            if (action == GLFW_PRESS)
            {
                glfwGetCursorPos(window, &prevMouseX, &prevMouseY);
                mousePressed = true;
            }
            else if (action == GLFW_RELEASE)
            {
                mousePressed = false;
            }
        }
    }
}

void GLProgram::scrollCallback(GLFWwindow *window, double xoffset, double yoffset)
{
    if (enableControl)
        camera.processMouseScroll(yoffset);
}

void GLProgram::cursorPosCallback(GLFWwindow *window, double xpos, double ypos)
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

        // to get the axis of rotation, need to convert from camera coordinates to object local coordinates
        glm::vec3 axisCamera = glm::cross(va, vb);
        glm::mat4 objectModelMatrix = buildDefaultModelMatrix() * modelMatrix;
        glm::mat3 cameraToModel = glm::inverse(glm::mat3(camera.getViewMatrix() * objectModelMatrix));
        glm::vec3 axisModel = cameraToModel * axisCamera;

        // update model rotation matrix
        float tolerance = 1e-4;
        if (angleOfRotation > tolerance)
        {
            glm::vec3 bboxCenter = 0.5f * (sceneBBox[0] + sceneBBox[1]);
            glm::mat4 rotateAroundCenter = glm::translate(glm::mat4(1.0f), bboxCenter) *
                                           glm::rotate(glm::mat4(1.0f), glm::degrees(angleOfRotation), axisModel) *
                                           glm::translate(glm::mat4(1.0f), -bboxCenter);
            modelMatrix = modelMatrix * rotateAroundCenter;
        }

        // update cursor position
        prevMouseX = currMouseX;
        prevMouseY = currMouseY;
    }
}

glm::vec3 GLProgram::getArcballVector(float x, float y)
{
    // get normalized vector from center of the virtual arcball to a point P on the arcball's surface
    // if (x,y) is too far away from the arcball, return the nearest point on the arcball's surface
    glm::vec3 P(x / windowWidth * 2 - 1.0, y / windowHeight * 2 - 1.0, 0.0f);
    P.y = -P.y;

    float radius = 1.0f;
    float OP_squared = P.x * P.x + P.y * P.y;

    if (OP_squared <= radius)
        P.z = sqrt(radius - OP_squared); // apply pythagorean theorem to find z
    else
        P = glm::normalize(P); // nearest point

    return P;
}

void GLProgram::fittingToObjects()
{
    glm::mat4 defaultModelMatrix = glm::rotate(
        glm::mat4(1.0f),
        glm::radians(45.0f),
        glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 worldModelMatrix = defaultModelMatrix * modelMatrix;

    glm::vec3 localMin = sceneBBox[0];
    glm::vec3 localMax = sceneBBox[1];

    glm::vec3 worldMin(0.0f);
    glm::vec3 worldMax(0.0f);
    bool hasWorldPoint = false;

    for (int x = 0; x < 2; ++x)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int z = 0; z < 2; ++z)
            {
                glm::vec3 corner(
                    x == 0 ? localMin.x : localMax.x,
                    y == 0 ? localMin.y : localMax.y,
                    z == 0 ? localMin.z : localMax.z);
                glm::vec3 worldCorner = glm::vec3(worldModelMatrix * glm::vec4(corner, 1.0f));

                if (!hasWorldPoint)
                {
                    worldMin = worldCorner;
                    worldMax = worldCorner;
                    hasWorldPoint = true;
                }
                else
                {
                    worldMin = glm::min(worldMin, worldCorner);
                    worldMax = glm::max(worldMax, worldCorner);
                }
            }
        }
    }

    glm::vec3 target = 0.5f * (worldMin + worldMax);
    glm::vec3 worldSize = worldMax - worldMin;
    float aspect = static_cast<float>(windowWidth) / static_cast<float>(windowHeight > 0 ? windowHeight : 1);
    float halfHeight = worldSize.z * 0.5f;
    float halfWidth = glm::max(worldSize.x, worldSize.y) * 0.5f;
    float fitHalfHeight = glm::max(halfHeight, halfWidth / glm::max(aspect, 1e-4f)) * 1.1f;
    if (fitHalfHeight < 1.0f)
        fitHalfHeight = 1.0f;
    camera.zoom = fitHalfHeight * 2.0f;

    glm::vec3 viewDir = camera.front;
    if (glm::length(viewDir) < 1e-4f)
        viewDir = glm::vec3(0.0f, 0.0f, -1.0f);
    else
        viewDir = glm::normalize(viewDir);

    float fitDistance = glm::max(glm::length(worldSize), camera.zoom) * 2.0f;
    glm::vec3 fitPosition = target - viewDir * fitDistance;
    camera.relocate(target, fitPosition);
}

void GLProgram::processKeyboardInput(void)
{
    static bool vKeyHandled = false;

    // close window with 'ESC' key
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (enableControl)
    {
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            camera.processKeyboard(UP, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            camera.processKeyboard(DOWN, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            camera.processKeyboard(LEFT, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            camera.processKeyboard(RIGHT, deltaTime);

        bool vKeyPressed = (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS);
        if (vKeyPressed && !vKeyHandled)
        {
            fittingToObjects();
        }

        vKeyHandled = vKeyPressed;
    }
    else
        vKeyHandled = false;
}

void apiOfCreateLoftSurface(const vector<tinynurbs::RationalCurve3d> &inputSections, tinynurbs::RationalSurface3d &result)
{
    myNurbs::createLoftSurface(inputSections, result);
}
