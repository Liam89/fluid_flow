SET(GMOCK_DIR "/home/liam/dev/gmock/gmock-1.7.0" CACHE PATH "gmock sources directory")

add_subdirectory(${GMOCK_DIR} ${CMAKE_BINARY_DIR}/gmock)
include_directories(SYSTEM ${GMOCK_DIR}/include ${GMOCK_DIR}/gtest/include)
