SET(GMOCK_DIR "../../gmock/gmock-1.7.0" CACHE PATH "gmock sources directory")
SET(TEST_AFTER_BUILD ON CACHE BOOL "Whether to automatically run tests after building")

add_subdirectory(${GMOCK_DIR} ${CMAKE_BINARY_DIR}/gmock)
include_directories(SYSTEM ${GMOCK_DIR}/include ${GMOCK_DIR}/gtest/include)

# Automatically run the tests after building
function(add_auto_test target)
  add_test(${target} ${target})

  IF(TEST_AFTER_BUILD)
    add_custom_command(
      TARGET ${target}
      POST_BUILD
      COMMAND ${target}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Running ${target}" VERBATIM
    )
  ENDIF()
endfunction()
