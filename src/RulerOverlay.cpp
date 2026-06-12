#include "../include/RulerOverlay.h"

#include <cmath>
#include <cstdio>

#include "imgui.h"

namespace
{
    constexpr float kRulerLabelFontSize = 16.0f;

    float computeNiceStep(float rawStep)
    {
        if (rawStep <= 0.0f)
            return 1.0f;

        float magnitude = pow(10.0f, floor(log10(rawStep)));
        float normalized = rawStep / magnitude;

        if (normalized <= 1.0f)
            return magnitude;
        if (normalized <= 2.0f)
            return 2.0f * magnitude;
        if (normalized <= 5.0f)
            return 5.0f * magnitude;
        return 10.0f * magnitude;
    }

    float roundUpToNiceStep(float rawStep)
    {
        if (rawStep <= 1.0f)
            return 1.0f;

        float step = computeNiceStep(rawStep);
        if (step < rawStep)
        {
            float magnitude = pow(10.0f, floor(log10(step)));
            float normalized = step / magnitude;
            if (normalized < 1.5f)
                step = 2.0f * magnitude;
            else if (normalized < 3.5f)
                step = 5.0f * magnitude;
            else
                step = 10.0f * magnitude;
        }
        return step;
    }

    void drawBorderRuler(
        ImDrawList *drawList,
        float visibleLength,
        float startCoordinate,
        float borderCoordinate,
        float availableLength,
        bool vertical,
        ImU32 color)
    {
        if (visibleLength <= 1e-4f || availableLength <= 1.0f)
            return;

        float unitStep = roundUpToNiceStep(visibleLength / 20.0f);
        if (unitStep <= 0.0f)
            return;
        float tickLength = 10.0f;
        float pixelsPerUnit = availableLength / visibleLength;

        for (float value = unitStep; value <= visibleLength + unitStep * 0.25f; value += unitStep)
        {
            float offset = value * pixelsPerUnit;
            float coord = vertical ? (startCoordinate - offset) : (startCoordinate + offset);

            ImVec2 point = vertical ? ImVec2(borderCoordinate, coord) : ImVec2(coord, borderCoordinate);
            ImVec2 tickEnd = vertical ? ImVec2(borderCoordinate + tickLength, coord) : ImVec2(coord, borderCoordinate - tickLength);
            drawList->AddLine(point, tickEnd, color, 1.0f);

            char label[32];
            snprintf(label, sizeof(label), "%.0f", value);
            ImVec2 textSize = ImGui::GetFont()->CalcTextSizeA(kRulerLabelFontSize, 10000.0f, 0.0f, label);
            ImVec2 textPos = vertical ? ImVec2(borderCoordinate + tickLength + 4.0f, coord - textSize.y * 0.5f)
                                      : ImVec2(coord - textSize.x * 0.5f, borderCoordinate - tickLength - textSize.y - 2.0f);
            drawList->AddText(ImGui::GetFont(), kRulerLabelFontSize, textPos, IM_COL32(255, 255, 255, 255), label);
        }
    }

    void drawRulerZeroLabel(ImDrawList *drawList, float leftBorderX, float bottomBorderY)
    {
        const char *label = "0";
        ImVec2 textSize = ImGui::GetFont()->CalcTextSizeA(kRulerLabelFontSize, 10000.0f, 0.0f, label);
        ImVec2 textPos(leftBorderX + 4.0f, bottomBorderY - textSize.y - 2.0f);
        drawList->AddText(ImGui::GetFont(), kRulerLabelFontSize, textPos, IM_COL32(255, 255, 255, 255), label);
    }
}

void RulerOverlay::Draw(ImDrawList *drawList, int windowWidth, int windowHeight, float cameraZoom)
{
    const float leftBorderX = 12.0f;
    const float bottomBorderY = static_cast<float>(windowHeight) - 12.0f;

    drawList->AddLine(ImVec2(leftBorderX, 0.0f), ImVec2(leftBorderX, static_cast<float>(windowHeight)), IM_COL32(110, 110, 110, 200), 1.0f);
    drawList->AddLine(ImVec2(0.0f, bottomBorderY), ImVec2(static_cast<float>(windowWidth), bottomBorderY), IM_COL32(110, 110, 110, 200), 1.0f);
    drawRulerZeroLabel(drawList, leftBorderX, bottomBorderY);

    float aspect = static_cast<float>(windowWidth) / static_cast<float>(windowHeight > 0 ? windowHeight : 1);
    float visibleHeight = cameraZoom;
    float visibleWidth = cameraZoom * aspect;

    drawBorderRuler(
        drawList,
        visibleWidth,
        leftBorderX,
        bottomBorderY,
        static_cast<float>(windowWidth) - leftBorderX,
        false,
        IM_COL32(170, 70, 70, 255));

    drawBorderRuler(
        drawList,
        visibleHeight,
        bottomBorderY,
        leftBorderX,
        bottomBorderY,
        true,
        IM_COL32(70, 70, 170, 255));
}
