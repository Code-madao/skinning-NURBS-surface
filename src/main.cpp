#include <exception>
#include <iostream>

#include "../include/AppPaths.h"
#include "../include/GLProgram.h"
#include "../include/LoftSceneBuilder.h"
#include "../include/ProjectRoot.h"
#include "../include/TestCaseLoader.h"
#include "../include/TestCaseSelector.h"

#define WINDOW_WIDTH 1600
#define WINDOW_HEIGHT 1200

int GLProgram::windowWidth = WINDOW_WIDTH;
int GLProgram::windowHeight = WINDOW_HEIGHT;
Camera GLProgram::camera;
bool GLProgram::mousePressed = false;
double GLProgram::prevMouseX, GLProgram::prevMouseY;
glm::mat4 GLProgram::modelMatrix = glm::mat4(1.0f);
glm::mat2x3 GLProgram::sceneBBox = glm::mat2x3(glm::vec3(0.0f), glm::vec3(0.0f));

int main(int argc, char **argv)
{
    ProjectRoot::SetWorkingDirectoryToProjectRoot();

    try
    {
        const TestCaseLibrary library = loadTestCaseLibrary(AppPaths::kTestCasesJson);
        const TestCaseData &selectedCase = TestCaseSelector::Select(library, argc, argv);

        std::cout << "Loaded test case: " << selectedCase.name << std::endl;
        if (!selectedCase.notes.empty())
            std::cout << "Notes: " << selectedCase.notes << std::endl;

        LoftScene scene = LoftSceneBuilder::Build(selectedCase);
        LoftSceneBuilder::PrintSurfaceInfo(scene.surfaces.front());

        GLProgram program;
        program.init(scene.curves, scene.surfaces);
        program.setClearColor(0.05f, 0.18f, 0.25f, 1.0f);
        program.run(scene.curves, scene.surfaces);
        program.cleanup();

        return 0;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Failed to start test case: " << ex.what() << std::endl;
        try
        {
            const TestCaseLibrary library = loadTestCaseLibrary(AppPaths::kTestCasesJson);
            TestCaseSelector::PrintAvailableCases(library);
        }
        catch (...)
        {
        }
        return 1;
    }
}
