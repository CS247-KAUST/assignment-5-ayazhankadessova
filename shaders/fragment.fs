#version 410 core

in vec2 texCoord;

uniform sampler2D txtr;
uniform vec4 vertexColor;

uniform int   colormapMode; // -1 = solid overlay (use vertexColor), 0 = grayscale, 1 = rainbow, 2 = cool-warm
uniform float blendFactor;  // 0..1 blend between grayscale and colormap

out vec4 fragColor;

// Rainbow: blue -> cyan -> green -> yellow -> red
vec3 rainbow(float v) {
    v = clamp(v, 0.0, 1.0);
    vec3 c;
    if (v < 0.25) {
        c = mix(vec3(0.0, 0.0, 1.0), vec3(0.0, 1.0, 1.0), v / 0.25);
    } else if (v < 0.5) {
        c = mix(vec3(0.0, 1.0, 1.0), vec3(0.0, 1.0, 0.0), (v - 0.25) / 0.25);
    } else if (v < 0.75) {
        c = mix(vec3(0.0, 1.0, 0.0), vec3(1.0, 1.0, 0.0), (v - 0.5) / 0.25);
    } else {
        c = mix(vec3(1.0, 1.0, 0.0), vec3(1.0, 0.0, 0.0), (v - 0.75) / 0.25);
    }
    return c;
}

// Cool-warm: blue -> white -> red
vec3 coolwarm(float v) {
    v = clamp(v, 0.0, 1.0);
    if (v < 0.5) {
        return mix(vec3(0.23, 0.30, 0.75), vec3(0.95, 0.95, 0.95), v / 0.5);
    } else {
        return mix(vec3(0.95, 0.95, 0.95), vec3(0.71, 0.02, 0.15), (v - 0.5) / 0.5);
    }
}

void main() {
    // Overlay path: solid color, ignore texture
    if (colormapMode < 0) {
        fragColor = vertexColor;
        return;
    }

    float s = texture(txtr, texCoord).r;
    vec3 gray = vec3(s);
    vec3 mapped = gray;
    if (colormapMode == 1)      mapped = rainbow(s);
    else if (colormapMode == 2) mapped = coolwarm(s);

    vec3 finalColor = mix(gray, mapped, blendFactor);
    fragColor = vec4(finalColor, 1.0) + vertexColor;
}
