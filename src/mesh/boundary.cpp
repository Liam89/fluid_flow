#include "boundary.h"

namespace Mesh
{
    template <int dim>
    BoundaryValues<dim>::BoundaryValues() : dealii::Function<dim>() {}

    template <int dim>
    BoundaryValues<dim>::~BoundaryValues(){}

    template <int dim>
    double BoundaryValues<dim>::value(const dealii::Point<dim> &p, const unsigned int component) const
    {
        return p(0)*p(0);
    }

    template class BoundaryValues<2>; // only care about 2d for now
}

