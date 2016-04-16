#ifndef SIMBASE_H
#define SIMBASE_H

#include <string>

#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>

namespace Mesh
{
    template <int> class Grid;
}


namespace Simulation
{
    template <int dim>
    class SimulationBase
    {
    public:
        SimulationBase();
        ~SimulationBase();

        // Run an iteration of the simulation
        virtual const std::string run();
        // Get the value of the solution at the point
        double get_point_value(const dealii::Point<dim> point) const;

    protected:
        std::unique_ptr<Mesh::Grid<dim>> grid;
        // reference to DoFHandler associated with above Grid for convenience
        dealii::DoFHandler<dim>& dof_handler;
        dealii::Vector<double> solution;

    private:
        virtual void setup_system() = 0;
        virtual void assemble_system() = 0;
        virtual void solve() = 0;
        std::string get_vtu_solution();
    };
}
#endif // SIMBASE_H
