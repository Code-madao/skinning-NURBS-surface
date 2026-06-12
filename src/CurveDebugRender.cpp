#include "../include/CurveDebugRender.h"

#include "../include/GeometryDebugUtils.h"
#include "../include/mynurbs.h"

void CurveDebugRender::BuildDerivativeLine(
    const tinynurbs::RationalCurve3d &curve,
    double u,
    float length,
    std::vector<std::vector<glm::vec3>> &vertices,
    std::vector<std::vector<glm::vec3>> &offsetVertex)
{
    glm::vec3 point = myNurbs::myRationalCurvePoint(curve, u);
    std::vector<glm::vec<3, double>> derivatives = myNurbs::myRationalCurveDerivative(curve, 1, u);
    glm::vec3 direction = derivatives[1];
    GeometryDebugUtils::NormalizeVec3(direction);
    glm::vec3 end = GeometryDebugUtils::ComputeExtendedPoint(point, direction, length);

    vertices.resize(1);
    offsetVertex.resize(1);
    vertices[0].push_back(point);
    offsetVertex[0].push_back(end);
}
