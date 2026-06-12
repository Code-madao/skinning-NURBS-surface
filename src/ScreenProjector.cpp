#include "../include/ScreenProjector.h"

#include "imgui.h"

#include <glm/gtc/matrix_transform.hpp>

bool ScreenProjector::WorldToScreen(
    const glm::vec3 &point,
    const glm::mat4 &viewMatrix,
    const glm::mat4 &projectionMatrix,
    int windowWidth,
    int windowHeight,
    ImVec2 &screenPos)
{
    glm::vec3 projected = glm::project(
        point,
        viewMatrix,
        projectionMatrix,
        glm::vec4(0.0f, 0.0f, static_cast<float>(windowWidth), static_cast<float>(windowHeight)));

    if (projected.z < 0.0f || projected.z > 1.0f)
        return false;

    screenPos = ImVec2(projected.x, static_cast<float>(windowHeight) - projected.y);
    return true;
}
