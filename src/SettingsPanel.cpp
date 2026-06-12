#include "../include/SettingsPanel.h"

#include "imgui.h"

void SettingsPanel::Draw(RunUiState &state)
{
    ImGui::Begin("Settings");
    ImGui::Checkbox("show shade surface", &state.showShadeSurface);
    ImGui::Checkbox("show curves", &state.showCurves);
    ImGui::Checkbox("show grid", &state.showGrid);
    ImGui::Checkbox("show control mesh", &state.showControlMesh);
    ImGui::Checkbox("show curve control points line", &state.showCurveControlPoints);
    ImGui::Checkbox("show curves derivate", &state.showCurvesDerivate);
    
    ImGui::Checkbox("show tangent", &state.showTangent);
    ImGui::Checkbox("show normal", &state.showNormal);
    if (state.showTangent || state.showNormal)
    {
        ImGui::SliderFloat("u##1", &state.tangentU, 0.0f, 1.0f);
        ImGui::SliderFloat("v##1", &state.tangentV, 0.0f, 1.0f);
    }
    ImGui::Checkbox("show object axis", &state.showObjectAxis);
    ImGui::Checkbox("show ruler", &state.showRuler);
    ImGui::Checkbox("show view gizmo", &state.showViewGizmo);
    if (state.showCurvesDerivate)
    {
        ImGui::SliderFloat("u##2", &state.curvesDerivateU, 0.0f, 1.0f);
    }
    ImGui::End();
}
