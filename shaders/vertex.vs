#version 410 core

layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 tCoord;

out vec2 texCoord;

uniform mat4 model;

void main() {
    gl_Position = model * vec4(pos, 1.0);
    texCoord = vec4(tCoord, 1).xy;
}