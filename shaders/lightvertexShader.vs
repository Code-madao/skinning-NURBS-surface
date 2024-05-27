#version 450 core

layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 normal;

out vec3 fragPos;
out vec3 Normal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    Normal =  mat3(transpose(inverse(model))) * normal;      
    gl_Position = projection * view * model * vec4(pos, 1.0);
    fragPos = vec3(model * vec4(pos, 1.0));
}
