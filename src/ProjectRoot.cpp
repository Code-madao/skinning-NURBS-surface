#include "../include/ProjectRoot.h"

#include <filesystem>

void ProjectRoot::SetWorkingDirectoryToProjectRoot()
{
    namespace fs = std::filesystem;

    fs::path current = fs::current_path();
    while (!current.empty())
    {
        if (fs::exists(current / "CMakeLists.txt") &&
            fs::exists(current / "shaders") &&
            fs::exists(current / "imgui"))
        {
            fs::current_path(current);
            return;
        }

        fs::path parent = current.parent_path();
        if (parent == current)
            break;

        current = parent;
    }
}
