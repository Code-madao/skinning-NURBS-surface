#pragma once

#include <glm/glm.hpp>

namespace GeometryDebugUtils
{
    void NormalizeVec3(glm::vec3 &input);
    glm::vec3 ComputeExtendedPoint(const glm::vec3 &point, const glm::vec3 &vec, float length);
}
