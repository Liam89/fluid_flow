#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>

void create_grid()
{
	dealii::Triangulation<2> triangulation;
	dealii::GridGenerator::hyper_cube(triangulation);
	triangulation.refine_global(4);

	std::ofstream out("grid-1.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(triangulation, out);
	std::cout << "Grid written to grid-1.eps";
}

int main()
{
	create_grid();
}

