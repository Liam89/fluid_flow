#ifndef MESH_GRID_H
#define MESH_GRID_H

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include "boundary.h"

namespace Mesh
{
    template <int dim>
    class Grid
    {
    public:
        Grid();

        // Get the grid in vtu format so it can be understood by the Renderer (vtk)
        //  Note: the grid information is already contained in the simulation's get_vtu_solution
        const std::string get_vtu_grid() const;
        // Get the degrees of freedom handler associated with the triangulation
        dealii::DoFHandler<dim>& get_dof_handler();
        // Get the function decribing the grid's boundary conditions
        dealii::Function<dim>& get_boundary();
    private:
        void setup_grid();

        dealii::Triangulation<dim> triangulation;
        dealii::DoFHandler<dim> dof_handler;
        std::unique_ptr<dealii::Function<dim>> boundary;
    };
}

#endif // MESH_GRID_H
