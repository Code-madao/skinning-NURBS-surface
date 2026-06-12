#pragma once

#include <glad/glad.h>
#include <glm/glm.hpp>

struct ImVec2;

namespace ScreenProjector
{
    bool WorldToScreen(
        const glm::vec3 &point,
        const glm::mat4 &viewMatrix,
        const glm::mat4 &projectionMatrix,
        int windowWidth,
        int windowHeight,
        ImVec2 &screenPos);
}
