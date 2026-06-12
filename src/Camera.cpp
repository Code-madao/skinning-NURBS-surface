#include "../include/Camera.h"
#include "Camera.h"

// default constructor
Camera::Camera() : position(0.0f, 0.0f, 30.0f),
                   front(0.0f, 0.0f, -0.5f),
                   up(1.0f, 0.0f, 0.0f),
                   right(0.0f, -1.0f, 0.0f),
                   movementSpeed(15.0f),
                   zoom(45.0f),
                   zoomSpeed(2.0f)
{
}

void Camera::lookAt(const glm::vec3 &target)
{
    front = glm::normalize(target - position);

    glm::vec3 worldUp(0.0f, 0.0f, 1.0f);
    if (glm::abs(glm::dot(front, worldUp)) > 0.999f)
        worldUp = glm::vec3(0.0f, 1.0f, 0.0f);

    right = glm::normalize(glm::cross(front, worldUp));
    up = glm::normalize(glm::cross(right, front));

    pitch = glm::degrees(asin(front.y));
    yaw = glm::degrees(atan2(front.z, front.x));
}

void Camera::setPosition(const glm::vec3 &_pos)
{
    position = _pos;
}

glm::mat4 Camera::getViewMatrix(void)
{
    return glm::lookAt(position, position + front, up);
}

void Camera::processKeyboard(CameraMovement direction, float deltaTime)
{
    float velocity = movementSpeed * deltaTime;
    switch (direction)
    {
    case UP:
        position += up * velocity;
        break;
    case DOWN:
        position -= up * velocity;
        break;
    case LEFT:
        position -= right * velocity;
        break;
    case RIGHT:
        position += right * velocity;
        break;
    case FITTING:
        position += right * velocity;

        break;
    }
}

void Camera::processMouseScroll(float yOffset)
{
    zoom -= yOffset * zoomSpeed;

    // clamp between max and min zoom values
    if (zoom < MIN_ZOOM)
        zoom = MIN_ZOOM;
    if (zoom > MAX_ZOOM)
        zoom = MAX_ZOOM;
}

void Camera::relocate(const glm::vec3 &target, const glm::vec3 &_pos)
{
    setPosition(_pos);
    lookAt(target);
}
