# Shallow water equations solver

This is an example of a shallow water 2D equations solver in C++.

Numerical method approach is that of a common fluid motion solver:
a staggered grid for 2D-space discretization is used for water cylinder height and velocity vector unknowns.

Euler's method for time advancement is implemeted with two steps:
  1. Advect the unknowns with respect to the current velocity field
  2. Update the unknowns accordingly to the acceleration terms (right-hand side of the equation)
  
Solution can be visualized using VisualEngine class which is implemented with OpenGL&nbsp;3.3 core profile. Following libraries are used:
  * **GLEW** -- OpenGL extensions managing
  * **GLFW** -- windows and callbacks handler
  * **SOIL** -- texture loading
  * **GLM** -- linear algebra
  ---
  Based on <http://liris.cnrs.fr/alexandre.meyer/teaching/m2pro_charanim/papers_pdf/coursenotes_shallowWater.pdf>
