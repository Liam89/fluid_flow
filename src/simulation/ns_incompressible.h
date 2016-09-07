#ifndef NSINCOMPRESSIBLE_H
#define NSINCOMPRESSIBLE_H

#include "simulation_base.h"

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <cmath>
#include <iostream>

namespace Simulation
{
using namespace dealii;
    namespace RunTimeParameters
    {

    enum MethodFormulation
    {
      METHOD_STANDARD,
      METHOD_ROTATIONAL
    };

    class Data_Storage
    {
    public:
      Data_Storage();
      ~Data_Storage();
      MethodFormulation form;
      double initial_time,
             final_time,
             Reynolds;
      double dt;
      unsigned int n_global_refines,
               pressure_degree;
      unsigned int vel_max_iterations,
               vel_Krylov_size,
               vel_off_diagonals,
               vel_update_prec;
      double vel_eps,
             vel_diag_strength;
      bool verbose;
      unsigned int output_interval;
    protected:
      ParameterHandler prm;
    };

    }

    namespace EquationData
    {
        // As we have chosen a completely decoupled formulation, we will not take
        // advantage of deal.II's capabilities to handle vector valued
        // problems. We do, however, want to use an interface for the equation
        // data that is somehow dimension independent. To be able to do that, our
        // functions should be able to know on which spatial component we are
        // currently working, and we should be able to have a common interface to
        // do that. The following class is an attempt in that direction.
        template <int dim>
        class MultiComponentFunction: public Function<dim>
        {
        public:
          MultiComponentFunction (const double initial_time = 0.);
          void set_component (const unsigned int d);
        protected:
          unsigned int comp;
        };

        // With this class defined, we declare classes that describe the boundary
        // conditions for velocity and pressure:
        template <int dim>
        class Velocity : public MultiComponentFunction<dim>
        {
        public:
          Velocity (const double initial_time = 0.0);

          virtual double value (const Point<dim> &p,
                                const unsigned int component = 0) const;

          virtual void value_list (const std::vector< Point<dim> > &points,
                                   std::vector<double> &values,
                                   const unsigned int component = 0) const;
        };
    }

    using namespace dealii;
    // Navier Stokes incompressible flow
    template <int dim>
    class NSIncompressible : public SimulationBase<dim>
    {
    public:
        NSIncompressible();
        ~NSIncompressible();

        const std::string run();
      protected:
        RunTimeParameters::Data_Storage data;
        RunTimeParameters::MethodFormulation type;

        const unsigned int deg;
        const double       dt;
        const double       t_0, T, Re;

        EquationData::Velocity<dim>       vel_exact;
        std::map<types::global_dof_index, double>    boundary_values;
        std::vector<types::boundary_id> boundary_indicators;

        Triangulation<dim> triangulation;

        FE_Q<dim>          fe_velocity;
        FE_Q<dim>          fe_pressure;

        DoFHandler<dim>    dof_handler_velocity;
        DoFHandler<dim>    dof_handler_pressure;

        QGauss<dim>        quadrature_pressure;
        QGauss<dim>        quadrature_velocity;

        SparsityPattern    sparsity_pattern_velocity;
        SparsityPattern    sparsity_pattern_pressure;
        SparsityPattern    sparsity_pattern_pres_vel;

        SparseMatrix<double> vel_Laplace_plus_Mass;
        SparseMatrix<double> vel_it_matrix[dim];
        SparseMatrix<double> vel_Mass;
        SparseMatrix<double> vel_Laplace;
        SparseMatrix<double> vel_Advection;
        SparseMatrix<double> pres_Laplace;
        SparseMatrix<double> pres_Mass;
        SparseMatrix<double> pres_Diff[dim];
        SparseMatrix<double> pres_iterative;

        Vector<double> pres_n;
        Vector<double> pres_n_minus_1;
        Vector<double> phi_n;
        Vector<double> phi_n_minus_1;
        Vector<double> u_n[dim];
        Vector<double> u_n_minus_1[dim];
        Vector<double> u_star[dim];
        Vector<double> force[dim];
        Vector<double> v_tmp;
        Vector<double> pres_tmp;
        Vector<double> rot_u;

        SparseILU<double> prec_velocity[dim];
        SparseILU<double> prec_pres_Laplace;
        SparseDirectUMFPACK prec_mass;
        SparseDirectUMFPACK prec_vel_mass;

        DeclException2 (ExcInvalidTimeStep,
                        double, double,
                        << " The time step " << arg1 << " is out of range."
                        << std::endl
                        << " The permitted range is (0," << arg2 << "]");

        void create_triangulation_and_dofs (const unsigned int n_refines);

        void initialize();

        void interpolate_velocity ();

        void diffusion_step (const bool reinit_prec);

        void projection_step (const bool reinit_prec);

        void update_pressure (const bool reinit_prec);
    private:
        unsigned int vel_max_its;
        unsigned int vel_Krylov_size;
        unsigned int vel_off_diagonals;
        unsigned int vel_update_prec;
        double       vel_eps;
        double       vel_diag_strength;

        void initialize_velocity_matrices();

        void initialize_pressure_matrices();
        typedef std_cxx11::tuple< typename DoFHandler<dim>::active_cell_iterator,
                typename DoFHandler<dim>::active_cell_iterator
                > IteratorTuple;

        typedef SynchronousIterators<IteratorTuple> IteratorPair;

        void initialize_gradient_operator();

        struct InitGradPerTaskData
        {
          unsigned int              d;
          unsigned int              vel_dpc;
          unsigned int              pres_dpc;
          FullMatrix<double>        local_grad;
          std::vector<types::global_dof_index> vel_local_dof_indices;
          std::vector<types::global_dof_index> pres_local_dof_indices;

          InitGradPerTaskData (const unsigned int dd,
                               const unsigned int vdpc,
                               const unsigned int pdpc)
            :
            d(dd),
            vel_dpc (vdpc),
            pres_dpc (pdpc),
            local_grad (vdpc, pdpc),
            vel_local_dof_indices (vdpc),
            pres_local_dof_indices (pdpc)
          {}
        };

        struct InitGradScratchData
        {
          unsigned int  nqp;
          FEValues<dim> fe_val_vel;
          FEValues<dim> fe_val_pres;
          InitGradScratchData (const FE_Q<dim> &fe_v,
                               const FE_Q<dim> &fe_p,
                               const QGauss<dim> &quad,
                               const UpdateFlags flags_v,
                               const UpdateFlags flags_p)
            :
            nqp (quad.size()),
            fe_val_vel (fe_v, quad, flags_v),
            fe_val_pres (fe_p, quad, flags_p)
          {}
          InitGradScratchData (const InitGradScratchData &data)
            :
            nqp (data.nqp),
            fe_val_vel (data.fe_val_vel.get_fe(),
                        data.fe_val_vel.get_quadrature(),
                        data.fe_val_vel.get_update_flags()),
            fe_val_pres (data.fe_val_pres.get_fe(),
                         data.fe_val_pres.get_quadrature(),
                         data.fe_val_pres.get_update_flags())
          {}
        };

        void assemble_one_cell_of_gradient (const IteratorPair  &SI,
                                            InitGradScratchData &scratch,
                                            InitGradPerTaskData &data);

        void copy_gradient_local_to_global (const InitGradPerTaskData &data);

        // The same general layout also applies to the following classes and
        // functions implementing the assembly of the advection term:
        void assemble_advection_term();

        struct AdvectionPerTaskData
        {
          FullMatrix<double>        local_advection;
          std::vector<types::global_dof_index> local_dof_indices;
          AdvectionPerTaskData (const unsigned int dpc)
            :
            local_advection (dpc, dpc),
            local_dof_indices (dpc)
          {}
        };

        struct AdvectionScratchData
        {
          unsigned int                 nqp;
          unsigned int                 dpc;
          std::vector< Point<dim> >    u_star_local;
          std::vector< Tensor<1,dim> > grad_u_star;
          std::vector<double>          u_star_tmp;
          FEValues<dim>                fe_val;
          AdvectionScratchData (const FE_Q<dim> &fe,
                                const QGauss<dim> &quad,
                                const UpdateFlags flags)
            :
            nqp (quad.size()),
            dpc (fe.dofs_per_cell),
            u_star_local (nqp),
            grad_u_star (nqp),
            u_star_tmp (nqp),
            fe_val (fe, quad, flags)
          {}

          AdvectionScratchData (const AdvectionScratchData &data)
            :
            nqp (data.nqp),
            dpc (data.dpc),
            u_star_local (nqp),
            grad_u_star (nqp),
            u_star_tmp (nqp),
            fe_val (data.fe_val.get_fe(),
                    data.fe_val.get_quadrature(),
                    data.fe_val.get_update_flags())
          {}
        };

        void assemble_one_cell_of_advection (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                             AdvectionScratchData &scratch,
                                             AdvectionPerTaskData &data);

        void copy_advection_local_to_global (const AdvectionPerTaskData &data);

        // The final few functions implement the diffusion solve as well as
        // postprocessing the output, including computing the curl of the
        // velocity:
        void diffusion_component_solve (const unsigned int d);

        std::string output_results (const unsigned int step);

        void assemble_vorticity (const bool reinit_prec);
    };

}
#endif // NSINCOMPRESSIBLE_H
