#include "grid.h"

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

namespace Mesh
{
    template <int dim>
    Grid<dim>::Grid()
        // dealii's dof_handler and triangulation are inextricably linked, so dof_handler is owned by the grid instead of the simulation
        : dof_handler{triangulation}, // associate degrees of freedom handler with triangulation
          boundary{new BoundaryValues<dim>{}}
    {
        setup_grid();
    }

    // Create the triangulation
    // For now a square mesh with 16 * 16 cells
    template <int dim>
    void Grid<dim>::setup_grid()
    {
        dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
        triangulation.refine_global(4); // 2^4 * 2^4 cells
    }

    // Convert the grid to a string in the vtu format
    template <int dim>
    const std::string Grid<dim>::get_vtu_grid() const
    {
        std::stringstream ss;
        dealii::GridOut grid_out;
        grid_out.write_vtu(triangulation, ss);
        return ss.str();
    }

    // Get the dof handler linked to the triangulation
    template <int dim>
    dealii::DoFHandler<dim>& Grid<dim>::get_dof_handler()
    {
        return dof_handler;
    }

    // Return dealii::Function describing the grid boundary
    template <int dim>
    dealii::Function<dim>& Grid<dim>::get_boundary()
    {
        return *boundary.get();
    }

    template class Grid<2>; // only care about 2d for now
}
