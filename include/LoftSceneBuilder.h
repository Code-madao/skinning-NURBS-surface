#pragma once

#include <vector>

#include "TestCaseLoader.h"
#include "tinynurbs/tinynurbs.h"

struct LoftScene
{
    std::vector<tinynurbs::RationalCurve3d> curves;
    std::vector<tinynurbs::RationalSurface3d> surfaces;
};

namespace LoftSceneBuilder
{
    LoftScene Build(const TestCaseData &testCase);
    void PrintSurfaceInfo(const tinynurbs::RationalSurface3d &surface);
}
