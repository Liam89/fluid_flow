#include "simulation_base.h"
#include "mesh/grid.h"

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

namespace Simulation
{
    template <int dim>
    SimulationBase<dim>::SimulationBase()
        : grid{new Mesh::Grid<dim>{}},
          // brace style initialization broken with references in gcc <4.9
          // we rarely use grid directly and mostly use the dof_handler, so keep a reference to it for convenience
          dof_handler(grid->get_dof_handler())
    {}


    // explicit destructor so we can use unique_ptr with an incomplete type
    template <int dim>
    SimulationBase<dim>::~SimulationBase(){}


    // Get the value of the solution at the point
    template <int dim>
    double SimulationBase<dim>::get_point_value(const dealii::Point<dim> point) const
    {
        dealii::Vector<double> value{1};
        dealii::VectorTools::point_value(dof_handler,solution,point,value);
        return value[0];
    }

    // Convert the solution to a string in the vtu format
    template <int dim>
    std::string SimulationBase<dim>::get_vtu_solution()
    {
        std::stringstream ss;
        dealii::DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution, "solution");
        data_out.build_patches();
        data_out.write_vtu(ss);
        return ss.str();
    }

    template <int dim>
    const std::string SimulationBase<dim>::run()
    {
        setup_system();
        assemble_system();
        solve();
        return get_vtu_solution();
    }

    template class SimulationBase<2>;   // only care about 2d for now
}
