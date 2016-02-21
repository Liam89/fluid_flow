#ifndef MESH_BOUNDARY_H
#define MESH_BOUNDARY_H

#include <deal.II/base/function.h>

namespace Mesh {
    template <int dim>
    class BoundaryValues : public dealii::Function<dim>
    {
    public:
        BoundaryValues();
        virtual ~BoundaryValues();
        virtual double value(const dealii::Point<dim>   &p,
                             const unsigned int  component = 0) const;
    };
}

#endif // MESH_BOUNDARY_H
