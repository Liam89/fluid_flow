#ifndef SIM_LAPLACE_H
#define SIM_LAPLACE_H

#include <string>

// grid setup
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

// setup of system of equations
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

// solving system of equations
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>

namespace Simulation
{
    class Simulation
    {
    public:
        Simulation();
        const std::string run();
        const std::string get_grid() const;

    private:
        void setup_grid();
        void setup_system();
        void assemble_system();
        void solve();
        std::string get_solution();

        dealii::Triangulation<2> triangulation;
        dealii::FE_Q<2> fe;
        dealii::DoFHandler<2> dof_handler;
        dealii::SparsityPattern sparsity_pattern;
        dealii::SparseMatrix<double> system_matrix;
        dealii::Vector<double> solution;
        dealii::Vector<double> system_rhs;
    };
}
#endif // SIM_LAPLACE_H
