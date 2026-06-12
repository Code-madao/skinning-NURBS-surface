#include "../include/LoftSceneBuilder.h"

#include <iostream>

#include "../include/mynurbs.h"

LoftScene LoftSceneBuilder::Build(const TestCaseData &testCase)
{
    LoftScene scene;
    scene.curves = testCase.curves;

    tinynurbs::RationalSurface3d loftSurface;
    myNurbs::createLoftSurface(scene.curves, loftSurface);

    assert(tinynurbs::surfaceIsValid(loftSurface));
    scene.surfaces.emplace_back(loftSurface);
    return scene;
}

void LoftSceneBuilder::PrintSurfaceInfo(const tinynurbs::RationalSurface3d &surface)
{
    using namespace std;
    cout << "degree u v:" << surface.degree_u << ' ' << surface.degree_v << endl;
    cout << "knots u:";
    for each (const auto &ele in surface.knots_u)
    {
        cout << ' ' << ele;
    }
    cout << endl;
    cout << "knots v:";
    for each (const auto &ele in surface.knots_v)
    {
        cout << ' ' << ele;
    }
    cout << endl;
    cout << "ctrpnts.rows cols: " << surface.control_points.rows() << ' ' << surface.control_points.cols() << endl;
}
