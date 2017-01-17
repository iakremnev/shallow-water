cmake_minimum_required(VERSION 3.6)
project(shallow_water)

set(CMAKE_CXX_STANDARD 11)

set(MY_TARGET shallow_water)

set(SOURCE_FILES main.cpp solver.cpp render.cpp Shader.cpp Camera.cpp)
add_executable(${MY_TARGET} ${SOURCE_FILES})

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/cmake-build-debug")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v")

find_package(glfw3 3.2 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(X11 REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

target_include_directories(${MY_TARGET} PUBLIC ${OPENGL_INCLUDE_DIR})
target_link_libraries(${MY_TARGET} glfw3 ${OPENGL_gl_LIBRARY} GLEW X11 Xxf86vm Xrandr Xi dl Xinerama Xcursor pthread SOIL)