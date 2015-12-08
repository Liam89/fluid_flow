CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(DEAL_II_DIR "/home/liam/dev/deallii/installed")
FIND_PACKAGE(deal.II 8.2 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Must be called after DEAL_II_INITIALIZE_CACHED_VARIABLES() but before qt/vtk
PROJECT(FluidFlow)

SET(VTK_DIR "/home/liam/dev/VTK/VTK-6.3.0-build")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET(QT_QMAKE_EXECUTABLE "/home/liam/dev/qt/qt5.5.1/5.5/gcc_64/bin/qmake")
set(CMAKE_AUTOMOC ON)
find_package(Qt5Widgets REQUIRED QUIET)

# qt generates source files, so include the following dirs to enable out-of-source builds
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)


SET(TARGET "simulation")
ADD_EXECUTABLE(${TARGET} main.cc)
DEAL_II_SETUP_TARGET(${TARGET})
target_link_libraries(${TARGET} ${VTK_LIBRARIES})
