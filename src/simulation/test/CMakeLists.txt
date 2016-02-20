add_executable(sim_test simulation.cpp)

target_link_libraries(sim_test gtest_main simulation)

add_auto_test(sim_test)
