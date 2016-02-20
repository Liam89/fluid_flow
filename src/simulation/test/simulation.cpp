#include <iostream>

#include "gtest/gtest.h"

#include "simulation/poisson.h"
#include <deal.II/numerics/vector_tools.h>

#include <vector>

// Basic test of the poisson simulation with boundary of phi = x^2
TEST(Poisson, Saddle)
{
    struct result { dealii::Point<2> point; double value; };
    // boundary should be exact
    std::vector<result> expected_boundary = {
        //  x = -1         x = 0            x = 1
        { {-1,-1}, 1 }, { {0,-1}, 0 }, { {1,-1}, 1}, // y = -1
        { {-1, 0}, 1 },                { {1, 0}, 1}, // y =  0
        { {-1, 1}, 1 }, { {0, 1}, 0 }, { {1, 1}, 1}  // y =  1
    };
    // we know the result should look like an inverted saddle, so the center point should be ~0.5
    std::vector<result> expected_solution = { { {0,0}, 0.5 } };

    Simulation::Poisson poisson_problem;
    poisson_problem.run();

    for(auto expected_result : expected_boundary) {
        auto actual_value = poisson_problem.get_point_value(expected_result.point);
        EXPECT_DOUBLE_EQ(actual_value, expected_result.value);
    }

    for(auto expected_result : expected_solution) {
        auto actual_value = poisson_problem.get_point_value(expected_result.point);
        EXPECT_NEAR(actual_value, expected_result.value, 0.2);
    }
}
