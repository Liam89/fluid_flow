// Adapted from deall II tutorials: https://dealii.org/8.2.1/doxygen/deal.II/Tutorial.html
#include "poisson.h"
#include "deal.II/base/function_lib.h"
#include <iostream>
#include <fstream>
#include <sstream>

namespace Simulation
{
    Poisson::Poisson()
        : fe(1), // 2d bi-linear lagrange finite element (a square with degrees of freedom at the corners)
          dof_handler(triangulation) // associate degrees of freedom handler with triangulation (the mesh),
                                     //  triangulation doesn't need to be set up yet
    {}

    // square mesh with 16 * 16 cells
    void Poisson::setup_grid()
    {
        dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
        triangulation.refine_global(4); // 2^4 * 2^4 cells
    }

    const std::string Poisson::get_grid() const
    {
        std::stringstream ss;
        dealii::GridOut grid_out;
        grid_out.write_vtu(triangulation, ss);
        return ss.str();
    }

    std::string Poisson::get_solution()
    {
        std::stringstream ss;
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (solution, "solution");
        data_out.build_patches ();
        data_out.write_vtu(ss);
        return ss.str();
    }

    // Set up the dofs of the system and allocate memory for it
    void Poisson::setup_system()
    {
        // enumerate dofs using the mesh and finite element, the 2d bi-linear lagrange element means there
        //  will be 1 dof on each cell vertex of the mesh - i.e. 17*17 dof for a 16*16 grid.
        dof_handler.distribute_dofs(fe);

        // resize the system matrix, a sparse matrix is used since most elements are zero
        dealii::CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
        dealii::DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
        sparsity_pattern.copy_from(c_sparsity);
        system_matrix.reinit(sparsity_pattern);

        // resize the right hand side and solution vectors
        solution.reinit (dof_handler.n_dofs());
        system_rhs.reinit (dof_handler.n_dofs());
    }

    template <int dim>
    class BoundaryValues : public dealii::Function<dim>
    {
    public:
      BoundaryValues () : dealii::Function<dim>() {}
      virtual double value (const dealii::Point<dim>   &p,
                            const unsigned int  component = 0) const;
    };

    template <int dim>
    double BoundaryValues<dim>::value (const dealii::Point<dim> &p,
                                       const unsigned int /*component*/) const
    {
      return p(0)*p(0);
    }

    // Calculate the values of the system matrix and rhs
    void Poisson::assemble_system()
    {
        // quadrature to use when integrating over finite elements
        dealii::QGauss<2> quadrature_formula(2); // Gauss-Legendre quadrature, in 2d, with 2 points
        // fe_values will be used when calculating the shape function's value/gradient/quadrature weights at a point on each cell
        dealii::FEValues<2> fe_values(fe, quadrature_formula,
                                      dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);

        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points = quadrature_formula.size();

        // writing to the sparse global matrix is slow, so calculate the contribution of a cell on a local system then add it to
        //  the global system
        dealii::FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        dealii::Vector<double> cell_rhs(dofs_per_cell);
        std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

        for (auto cell : dof_handler.active_cell_iterators()) {
            // compute values, gradients, and weights of shape functions on cell
            fe_values.reinit(cell);
            // reset local system
            cell_matrix = 0;
            cell_rhs = 0;

            // Integrate over cell by looping over quadrature points and degrees of freedom, and assemble local system
            for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) {
                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                    for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                        cell_matrix(i,j) += (fe_values.shape_grad(i, q_index) *
                                            fe_values.shape_grad(j, q_index) *
                                            fe_values.JxW(q_index));
                    }
                }

                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                  cell_rhs(i) += (fe_values.shape_value(i, q_index) *
                                  0 *
                                  fe_values.JxW(q_index));
                }
            }

            // transfer cell contribution to global system
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    system_matrix.add(local_dof_indices[i],
                                      local_dof_indices[j],
                                      cell_matrix(i,j));
                }
            }
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
              system_rhs(local_dof_indices[i]) += cell_rhs(i);
            }

        }

        // apply boundary conditions
        std::map<dealii::types::global_dof_index,double> boundary_values;
        dealii::VectorTools::interpolate_boundary_values(dof_handler,
                                                        0,
                                                        BoundaryValues<2>(),
                                                        boundary_values);


        dealii::MatrixTools::apply_boundary_values(boundary_values,
                                                   system_matrix,
                                                   solution,
                                                   system_rhs);

    }




    void Poisson::solve()
    {
        dealii::SolverControl solver_control(1000, 1e-12);
        dealii::SolverCG<> solver(solver_control);
        solver.solve(system_matrix, solution, system_rhs, dealii::PreconditionIdentity());
    }

    const std::string Poisson::run()
    {
        setup_grid();
        setup_system();
        assemble_system();
        solve();
        return get_solution();
    }

}
