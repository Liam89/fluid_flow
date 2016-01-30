CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

# DEAL_II_INITIALIZE_CACHED_VARIABLES() must be called before setting the project name
include(dealii)

# Must be set before qt/vtk
PROJECT(FluidFlow)

include(qvtk)

enable_testing()
include(gmock)

add_subdirectory(gui)
add_subdirectory(simulation)

SET(TARGET "fluid_flow")
ADD_EXECUTABLE(${TARGET} main.cpp)
target_link_libraries(${TARGET} gui simulation)
