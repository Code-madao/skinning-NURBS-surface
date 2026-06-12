#include "../include/ViewGizmo.h"

#include <algorithm>
#include <array>

#include "imgui.h"

namespace
{
    constexpr float kViewGizmoRadius = 34.0f;
    constexpr float kViewGizmoMargin = 18.0f;
    constexpr float kViewGizmoLabelFontSize = 14.0f;

    ImU32 makeDepthTintedColor(const ImVec4 &baseColor, float depth)
    {
        float clampedDepth = glm::clamp(depth, -1.0f, 1.0f);
        float brightness = 0.45f + 0.55f * ((clampedDepth + 1.0f) * 0.5f);
        int r = static_cast<int>(glm::clamp(baseColor.x * brightness, 0.0f, 1.0f) * 255.0f);
        int g = static_cast<int>(glm::clamp(baseColor.y * brightness, 0.0f, 1.0f) * 255.0f);
        int b = static_cast<int>(glm::clamp(baseColor.z * brightness, 0.0f, 1.0f) * 255.0f);
        return IM_COL32(r, g, b, 255);
    }
}

void ViewGizmo::Draw(ImDrawList *drawList, const glm::mat4 &viewMatrix, const glm::mat4 &modelMatrix, int windowWidth)
{
    struct GizmoAxis
    {
        const char *label;
        glm::vec3 direction;
        ImVec4 baseColor;
    };

    ImVec2 center(
        static_cast<float>(windowWidth) - kViewGizmoMargin - kViewGizmoRadius,
        kViewGizmoMargin + kViewGizmoRadius);

    drawList->AddCircleFilled(center, kViewGizmoRadius + 10.0f, IM_COL32(20, 24, 28, 190), 48);
    drawList->AddCircle(center, kViewGizmoRadius + 10.0f, IM_COL32(180, 180, 180, 180), 48, 1.0f);

    glm::mat3 orientation = glm::mat3(viewMatrix * modelMatrix);
    std::array<GizmoAxis, 3> axes = {
        GizmoAxis{"X", glm::normalize(orientation * glm::vec3(1.0f, 0.0f, 0.0f)), ImVec4(0.86f, 0.22f, 0.22f, 1.0f)},
        GizmoAxis{"Y", glm::normalize(orientation * glm::vec3(0.0f, 1.0f, 0.0f)), ImVec4(0.22f, 0.78f, 0.22f, 1.0f)},
        GizmoAxis{"Z", glm::normalize(orientation * glm::vec3(0.0f, 0.0f, 1.0f)), ImVec4(0.22f, 0.45f, 0.92f, 1.0f)}};

    std::sort(axes.begin(), axes.end(), [](const GizmoAxis &lhs, const GizmoAxis &rhs)
              { return lhs.direction.z < rhs.direction.z; });

    drawList->AddCircleFilled(center, 3.5f, IM_COL32(235, 235, 235, 255), 16);

    for (const GizmoAxis &axis : axes)
    {
        glm::vec2 screenDir(axis.direction.x, -axis.direction.y);
        ImVec2 end(
            center.x + screenDir.x * kViewGizmoRadius,
            center.y + screenDir.y * kViewGizmoRadius);

        ImU32 axisColor = makeDepthTintedColor(axis.baseColor, axis.direction.z);
        float endpointRadius = 4.5f + 1.5f * ((glm::clamp(axis.direction.z, -1.0f, 1.0f) + 1.0f) * 0.5f);

        drawList->AddLine(center, end, axisColor, 2.5f);
        drawList->AddCircleFilled(end, endpointRadius, axisColor, 16);

        ImVec2 textSize = ImGui::GetFont()->CalcTextSizeA(kViewGizmoLabelFontSize, 10000.0f, 0.0f, axis.label);
        ImVec2 labelPos(
            end.x + screenDir.x * (endpointRadius + 5.0f) - textSize.x * 0.5f,
            end.y + screenDir.y * (endpointRadius + 5.0f) - textSize.y * 0.5f);

        if (glm::length(screenDir) < 0.2f)
        {
            labelPos.x += 6.0f;
            labelPos.y -= 6.0f;
        }

        drawList->AddText(ImGui::GetFont(), kViewGizmoLabelFontSize, labelPos, IM_COL32(255, 255, 255, 240), axis.label);
    }
}
