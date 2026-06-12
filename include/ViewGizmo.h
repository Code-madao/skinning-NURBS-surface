#pragma once

#include <glad/glad.h>
#include <glm/glm.hpp>

struct ImDrawList;

namespace ViewGizmo
{
    void Draw(ImDrawList *drawList, const glm::mat4 &viewMatrix, const glm::mat4 &modelMatrix, int windowWidth);
}
