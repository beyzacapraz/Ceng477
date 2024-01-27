#version 330 core

out vec4 fragColor;

uniform vec3 black;  // Define the color for black squares
uniform vec3 white;  // Define the color for white squares
uniform float offset; // Control the starting point of the pattern
uniform float scale;  // Control the size of the squares
uniform float scaleZ;  // Control the size of the squares

in vec3 pos;

void main() {
    // Calculate checkerboard pattern
    bool x = int((pos.x + offset) * scale) % 2 == 0;
    bool y = int((pos.y + offset) * scale) % 2 == 0;
    bool z = int((pos.z + offset)  * scale + scaleZ) % 2 == 0;
    bool xorXY = x != y;

    // Assign color based on the checkerboard pattern
    if (xorXY != z) {
        fragColor = vec4(black, 1.0);
    } else {
        fragColor = vec4(white, 1.0);
    }

}