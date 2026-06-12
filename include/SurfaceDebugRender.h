#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "tinynurbs/tinynurbs.h"

namespace SurfaceDebugRender
{
    void BuildTangentLines(
        const tinynurbs::RationalSurface3d &surface,
        double u,
        double v,
        std::vector<std::vector<glm::vec3>> &uVertices,
        std::vector<std::vector<glm::vec3>> &uOffsetVertex,
        std::vector<std::vector<glm::vec3>> &vVertices,
        std::vector<std::vector<glm::vec3>> &vOffsetVertex);

    void BuildNormalLine(
        const tinynurbs::RationalSurface3d &surface,
        double u,
        double v,
        float length,
        std::vector<std::vector<glm::vec3>> &vertices,
        std::vector<std::vector<glm::vec3>> &offsetVertex);
}
