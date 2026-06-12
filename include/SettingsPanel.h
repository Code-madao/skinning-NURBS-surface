#pragma once

struct RunUiState
{
    bool showGrid = false;
    bool showTangent = false;
    bool showNormal = false;
    float tangentU = 0.5f;
    float tangentV = 0.5f;
    bool showControlMesh = false;
    bool showShadeSurface = true;
    bool showCurves = false;
    bool showCurvesDerivate = false;
    float curvesDerivateU = 0.5f;
    bool showCurveControlPoints = false;
    bool showKnotInsertion = false;
    bool showKnotRemoval = false;
    bool showDegreeElevation = false;
    bool showDegreeReduction = false;
    bool showObjectAxis = true;
    bool showRuler = true;
    bool showViewGizmo = true;
};

namespace SettingsPanel
{
    void Draw(RunUiState &state);
}
