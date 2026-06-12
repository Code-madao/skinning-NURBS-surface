#include "../include/AxisOverlay.h"

#include "imgui.h"

#include "../include/ScreenProjector.h"

void AxisOverlay::DrawAxisLabel(
    ImDrawList *drawList,
    const char *label,
    const glm::vec3 &worldPoint,
    const glm::mat4 &viewMatrix,
    const glm::mat4 &projectionMatrix,
    int windowWidth,
    int windowHeight,
    ImU32 textColor)
{
    ImVec2 screenPos;
    if (!ScreenProjector::WorldToScreen(worldPoint, viewMatrix, projectionMatrix, windowWidth, windowHeight, screenPos))
        return;

    ImVec2 textSize = ImGui::CalcTextSize(label);
    drawList->AddText(ImVec2(screenPos.x + 6.0f, screenPos.y - textSize.y - 6.0f), textColor, label);
}
