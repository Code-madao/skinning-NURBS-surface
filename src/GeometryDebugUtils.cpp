#include "../include/GeometryDebugUtils.h"

#include <cmath>

void GeometryDebugUtils::NormalizeVec3(glm::vec3 &input)
{
    double length = sqrt(pow(input.x, 2) + pow(input.y, 2) + pow(input.z, 2));
    input.x = input.x / length;
    input.y = input.y / length;
    input.z = input.z / length;
}

glm::vec3 GeometryDebugUtils::ComputeExtendedPoint(const glm::vec3 &point, const glm::vec3 &vec, float length)
{
    glm::vec3 result;
    result.x = point.x - vec.x * length;
    result.y = point.y - vec.y * length;
    result.z = point.z - vec.z * length;
    return result;
}
