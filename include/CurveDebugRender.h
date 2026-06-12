#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "tinynurbs/tinynurbs.h"

namespace CurveDebugRender
{
    void BuildDerivativeLine(
        const tinynurbs::RationalCurve3d &curve,
        double u,
        float length,
        std::vector<std::vector<glm::vec3>> &vertices,
        std::vector<std::vector<glm::vec3>> &offsetVertex);
}
