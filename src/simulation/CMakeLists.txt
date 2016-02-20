file(GLOB CPP_FILES *.cpp)

add_library(simulation ${CPP_FILES})
DEAL_II_SETUP_TARGET(simulation)

add_subdirectory(test)
