file(GLOB CPP_FILES *.cpp)

add_library(mesh ${CPP_FILES})
DEAL_II_SETUP_TARGET(mesh)
