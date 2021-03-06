//
// Created by iak on 14.01.17.
//

#include "Camera.hpp"

const GLfloat YAW        = -90.0f;
const GLfloat PITCH      =   0.0f;
const GLfloat SPEED      =   3.0f;
const GLfloat SENSITIVTY =   0.25f;
const GLfloat ZOOM       =  45.0f;

Camera::Camera(
        glm::vec3 position,
        glm::vec3 up,
        GLfloat yaw,
        GLfloat pitch):

        Position(position),
        Front(glm::vec3(0.0f, 0.0f, -1.0f)),
        WorldUp(up),
        Yaw(yaw),
        Pitch(pitch),
        MovementSpeed(SPEED),
        MouseSensitivity(SENSITIVTY),
        Zoom(ZOOM)
{
    updateCameraVectors();
}

Camera::Camera(
        GLfloat posX, GLfloat posY, GLfloat posZ,
        GLfloat upX, GLfloat upY, GLfloat upZ,
        GLfloat yaw, GLfloat pitch):

        Position(glm::vec3(posX, posY, posZ)),
        Front(glm::vec3(0.0f, 0.0f, -1.0f)),
        WorldUp(glm::vec3(upX, upY, upZ)),
        Yaw(yaw),
        Pitch(pitch),
        MovementSpeed(SPEED),
        MouseSensitivity(SENSITIVTY),
        Zoom(ZOOM)
{
    updateCameraVectors();
}

glm::mat4 Camera::GetViewMatrix() {
    return glm::lookAt(Position, Position + Front, Up);
}

void Camera::ProcessKeyboard(Camera_Movement direction, GLfloat deltaTime) {
    GLfloat velocity = MovementSpeed * deltaTime;
    if (direction == FORWARD)
        Position += Front * velocity;
    if (direction == BACKWARD)
        Position -= Front * velocity;
    if (direction == LEFT)
        Position -= Right * velocity;
    if (direction == RIGHT)
        Position += Right * velocity;
}

void Camera::ProcessMouseMovement(GLfloat xoffset, GLfloat yoffset, GLboolean constrainPitch) {
    xoffset *= MouseSensitivity;
    yoffset *= MouseSensitivity;

    Yaw   += xoffset;
    Pitch += yoffset;

    // Make sure that when pitch is out of bounds, screen doesn't get flipped
    if (constrainPitch)
    {
        if (Pitch > 89.0f)
            Pitch = 89.0f;
        if (Pitch < -89.0f)
            Pitch = -89.0f;
    }

    // Update Front, Right and Up Vectors using the updated Eular angles
    updateCameraVectors();
}

void Camera::ProcessMouseScroll(GLfloat yoffset) {
    if (Zoom >= 1.0f && Zoom <= 45.0f)
        Zoom -= yoffset;
    if (Zoom <= 1.0f)
        Zoom = 1.0f;
    if (Zoom >= 45.0f)
        Zoom = 45.0f;
}

void Camera::updateCameraVectors() {
    // Calculate the new Front vector
    Front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
    Front.y = sin(glm::radians(Pitch));
    Front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
    // Also re-calculate the Right and Up vector
    Right = glm::normalize(glm::cross(Front, WorldUp));  // Normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
    Up    = glm::normalize(glm::cross(Right, Front));
}
