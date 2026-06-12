#pragma once

#include <glm/glm.hpp>

struct ImDrawList;
typedef unsigned int ImU32;

namespace AxisOverlay
{
    void DrawAxisLabel(
        ImDrawList *drawList,
        const char *label,
        const glm::vec3 &worldPoint,
        const glm::mat4 &viewMatrix,
        const glm::mat4 &projectionMatrix,
        int windowWidth,
        int windowHeight,
        ImU32 textColor);
}
