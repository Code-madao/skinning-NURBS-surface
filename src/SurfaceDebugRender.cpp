#include "../include/SurfaceDebugRender.h"

#include "../include/GeometryDebugUtils.h"
#include "../include/mynurbs.h"

void SurfaceDebugRender::BuildTangentLines(
    const tinynurbs::RationalSurface3d &surface,
    double u,
    double v,
    std::vector<std::vector<glm::vec3>> &uVertices,
    std::vector<std::vector<glm::vec3>> &uOffsetVertex,
    std::vector<std::vector<glm::vec3>> &vVertices,
    std::vector<std::vector<glm::vec3>> &vOffsetVertex)
{
    glm::vec3 point = myNurbs::myRationalSurfPoint(surface, u, v);
    tinynurbs::array2<glm::vec<3, double>> derivatives = myNurbs::myRationalSurfDerivative(surface, 1, u, v);

    glm::vec3 uDirection = derivatives(1, 0);
    GeometryDebugUtils::NormalizeVec3(uDirection);
    glm::vec3 uEnd = GeometryDebugUtils::ComputeExtendedPoint(point, uDirection, 1.0f);

    uVertices.resize(1);
    uOffsetVertex.resize(1);
    uVertices[0].resize(1);
    uOffsetVertex[0].resize(1);
    uVertices[0][0] = point;
    uOffsetVertex[0][0] = uEnd;

    glm::vec3 vDirection = derivatives(0, 1);
    GeometryDebugUtils::NormalizeVec3(vDirection);
    glm::vec3 vEnd = GeometryDebugUtils::ComputeExtendedPoint(point, vDirection, 1.0f);

    vVertices.resize(1);
    vOffsetVertex.resize(1);
    vVertices[0].resize(1);
    vOffsetVertex[0].resize(1);
    vVertices[0][0] = point;
    vOffsetVertex[0][0] = vEnd;
}

void SurfaceDebugRender::BuildNormalLine(
    const tinynurbs::RationalSurface3d &surface,
    double u,
    double v,
    float length,
    std::vector<std::vector<glm::vec3>> &vertices,
    std::vector<std::vector<glm::vec3>> &offsetVertex)
{
    glm::vec3 point = myNurbs::myRationalSurfPoint(surface, u, v);
    glm::vec3 normal = myNurbs::myRationalSurfNormal(surface, u, v);
    GeometryDebugUtils::NormalizeVec3(normal);
    glm::vec3 end = GeometryDebugUtils::ComputeExtendedPoint(point, normal, length);

    vertices.resize(1);
    offsetVertex.resize(1);
    vertices[0].resize(1);
    offsetVertex[0].resize(1);
    vertices[0][0] = point;
    offsetVertex[0][0] = end;
}
