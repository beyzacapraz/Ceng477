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

