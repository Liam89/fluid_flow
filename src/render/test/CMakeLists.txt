add_executable(render_test renderer.cpp)

target_link_libraries(render_test gtest_main render)

add_auto_test(render_test)
