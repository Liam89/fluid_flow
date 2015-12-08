#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractEdges.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <iostream>
#include <fstream>
#include <sstream>

std::string create_grid()
{
    dealii::Triangulation<2> triangulation;
	dealii::GridGenerator::hyper_cube(triangulation);
	triangulation.refine_global(4);

	std::stringstream ss;
	dealii::GridOut grid_out;
	grid_out.write_vtu(triangulation, ss);
    std::ofstream out ("grid.xml");
    grid_out.write_eps(triangulation, out);
	std::cout << "Grid written to stringstream";
	return ss.str();
}

void render(std::string vtu)
{
	  //read all the data from the file
	  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
	    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	  reader->ReadFromInputStringOn();
	  reader->SetInputString(vtu);
	  reader->Update();

      vtkSmartPointer<vtkExtractEdges> extractEdges =
          vtkSmartPointer<vtkExtractEdges>::New();
      extractEdges->SetInputConnection(reader->GetOutputPort());
      extractEdges->Update();

	  //Create a mapper and actor
	  vtkSmartPointer<vtkDataSetMapper> mapper =
	    vtkSmartPointer<vtkDataSetMapper>::New();
      mapper->SetInputConnection(extractEdges->GetOutputPort());

	  vtkSmartPointer<vtkActor> actor =
	    vtkSmartPointer<vtkActor>::New();
	  actor->SetMapper(mapper);

	  //Create a renderer, render window, and interactor
	  vtkSmartPointer<vtkRenderer> renderer =
	    vtkSmartPointer<vtkRenderer>::New();
	  vtkSmartPointer<vtkRenderWindow> renderWindow =
	    vtkSmartPointer<vtkRenderWindow>::New();
      renderWindow->AddRenderer(renderer);
	  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();
	  renderWindowInteractor->SetRenderWindow(renderWindow);

	  //Add the actor to the scene
	  renderer->AddActor(actor);
	  renderer->SetBackground(.3, .6, .3); // Background color green

	  //Render and interact
	  renderWindow->Render();
	  renderWindowInteractor->Start();
}

int main()
{
    std::string ss = create_grid();
	render(ss);
}

