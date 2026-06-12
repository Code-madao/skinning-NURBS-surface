#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#define MAX_ZOOM 100.0f
#define MIN_ZOOM 1.0f

enum CameraMovement
{
    UP,
    DOWN,
    LEFT,
    RIGHT,
    FITTING
};

class Camera
{
public:
    glm::vec3 position;
    glm::vec3 front;
    glm::vec3 up;
    glm::vec3 right;
    
    float movementSpeed;
    float zoom;
    float zoomSpeed;

    // unused
    float yaw;
    float pitch;

    Camera();
    void lookAt(const glm::vec3& target);
    void setPosition(const glm::vec3& _pos);
    glm::mat4 getViewMatrix(void);
    void processKeyboard(CameraMovement direction, float deltaTime);
    void processMouseScroll(float yOffset);
    void relocate(const glm::vec3& target, const glm::vec3& _pos);
};

#endif // CAMERA_H
