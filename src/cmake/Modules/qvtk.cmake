SET(VTK_DIR "/home/liam/dev/VTK/VTK-6.3.0-build" CACHE PATH "qvtk directory")
SET(QT_QMAKE_EXECUTABLE "/home/liam/dev/qt/qt5.5.1/5.5/gcc_64/bin/qmake" CACHE PATH "qmake exe path")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(CMAKE_AUTOMOC ON)
find_package(Qt5Widgets REQUIRED QUIET)
# disable QTs 'emit' 'signals' and 'slots' macros because they clash with deal.ii and boost
ADD_DEFINITIONS(-DQT_NO_KEYWORDS)


INCLUDE_DIRECTORIES(SYSTEM
  # qt generates source files, so include the following dirs to enable out-of-source builds
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}

  ${VTK_INCLUDE_DIRS}
)
