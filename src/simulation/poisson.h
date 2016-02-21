#ifndef SIM_POISSON_H
#define SIM_POISSON_H

#include <string>

#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>

namespace Mesh
{
    template <int> class Grid;
}


namespace Simulation
{
    // Simulation solving the poisson equation (grad^2 phi = f).
    // Currently the grid is fixed to have 16*16 cells, the boundary is phi = x^2, and f = 0
    class Poisson
    {
    public:
        Poisson();
        ~Poisson();

        // Run an iteration of the simulation
        const std::string run();
        // Get the value of the solution at the point
        double get_point_value(const dealii::Point<2> point) const;

    private:
        void setup_system();
        void assemble_system();
        void solve();
        std::string get_vtu_solution();

        dealii::FE_Q<2> fe;
        std::unique_ptr<Mesh::Grid<2>> grid;
        dealii::DoFHandler<2>& dof_handler;
        dealii::Vector<double> solution;

        // sparsity pattern lifetime must be at least as long as the system_matrix that uses it
        dealii::SparsityPattern sparsity_pattern;
        dealii::SparseMatrix<double> system_matrix;
        dealii::Vector<double> system_rhs;
    };
}
#endif // SIM_POISSON_H
