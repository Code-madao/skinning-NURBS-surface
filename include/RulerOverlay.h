#pragma once

struct ImDrawList;

namespace RulerOverlay
{
    void Draw(ImDrawList *drawList, int windowWidth, int windowHeight, float cameraZoom);
}
