# Ceng477
Hardware and software components of graphics systems. Output and filled-data primitives.  2D and 3D geometric transformations. Two-dimensional viewing. Three-dimensional viewing: Viewing pipeline, viewing parameters, projections, viewing transformations, clipping. Visible surface detection. Introduction to illumination models and surface rendering.

# Homework 1
Ray tracing is a fundamental rendering algorithm. It is commonly used for animations and architectural simulations, in which the quality of the created images is more important than the time
it takes to create them. In this assignment, I am going to implement a basic ray tracer that
simulates the propagation of light in real world.  

An example output image in PNG format:

![example image](https://github.com/beyzacapraz/Ceng477/blob/main/HW1/inputs/Car.png?raw=true)

# Homework 2

In this assignment, I implemented Modeling Transformation, Viewing Transformation,
and Rasterization stages of the Forward Rendering Pipeline. Given a set of triangles
and the position and the color attributes of all of their vertices, my aim was to render the scene as a
two-dimensional (2D) image. I first transformed all of the vertices to the viewport and then used line
drawing and triangle rasterization algorithms to display the triangles in wireframe or solid modes.
I implemented Liang-Barsky  clipping algorithm for wireframe mode only. 
I also implemented the backface culling (for both modes) for correct
rendering of the scene. Implementations are in C++ language.  

An example output image for demonstrating line rasterization: 
![empty box](https://github.com/beyzacapraz/Ceng477/blob/main/HW2/inputs_outputs/culling_enabled_outputs/empty_box/empty_box_1.ppm.png?raw=true)

Another output example:

![horse and mug](https://github.com/beyzacapraz/Ceng477/blob/main/HW2/inputs_outputs/culling_enabled_outputs/horse_and_mug/horse_and_mug_1.ppm.png?raw=true)

# Homework 3

In this assignment, I implemented an OpenGL game called Bunny Run. In this game, the
user controls a hopping bunny in the horizontal direction to capture the yellow checkpoints while
avoiding the red checkpoints in an everlasting journey on a wide road to obtain a high score. The
bunny continuously hops forward at an increasing speed, which makes it hard to keep the bunny
alive as time passes in a round. Hitting one of the red checkpoints ends the round by making the
bunny faint and lie down. A and D buttons should make the bunny move to its left and right, respectively.
R button should restart the round by re-initializing the values of Bunny Run.  

**Additional Libraries:** GLM, GLEW, GLFW, stb image.h, and FreeType libraries.  

An example image of Bunny Run:



