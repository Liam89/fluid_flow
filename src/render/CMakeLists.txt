file(GLOB CPP_FILES *.cpp)

add_library(render ${CPP_FILES})
DEAL_II_SETUP_TARGET(render)

#SET(TEST_AFTER_BUILD OFF)
add_subdirectory(test)
