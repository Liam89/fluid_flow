file(GLOB CPP_FILES *.cpp)

add_library(simulation ${CPP_FILES})
DEAL_II_SETUP_TARGET(simulation)
target_link_libraries(simulation mesh)

#SET(TEST_AFTER_BUILD OFF)
add_subdirectory(test)
