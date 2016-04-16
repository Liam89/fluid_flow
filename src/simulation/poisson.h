#ifndef SIM_POISSON_H
#define SIM_POISSON_H

#include "simulation_base.h"

#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/sparse_matrix.h>

namespace Simulation
{
    // Simulation solving the poisson equation (grad^2 phi = f).
    // Currently the grid is fixed to have 16*16 cells, the boundary is phi = x^2, and f = 0
    class Poisson : public SimulationBase<2>
    {
    public:
        Poisson();
        ~Poisson();

    private:
        void setup_system();
        void assemble_system();
        void solve();

        dealii::FE_Q<2> fe;

        // sparsity pattern lifetime must be at least as long as the system_matrix that uses it
        dealii::SparsityPattern sparsity_pattern;
        dealii::SparseMatrix<double> system_matrix;
        dealii::Vector<double> system_rhs;
    };
}
#endif // SIM_POISSON_H
